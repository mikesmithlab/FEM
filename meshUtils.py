import numpy as np

from scipy.spatial import Delaunay
from mayavi import mlab
import h5py
from Generic.filedialogs import load_filename, save_filename
import displayUtils as display



def generate_hex_mesh(width,length):
    #N must be an even number
    N = width
    x = np.linspace(0, width, N + 1 )
    y = np.linspace(0, length, N/2 + 1)
    
    xvals, yvals = np.meshgrid(x,y)
    a = (x[1]-x[0])/2.0
    b = (y[1]-y[0])/2.0
    
    xpts = np.append(xvals,xvals.copy()+a).reshape(1,-1)
    ypts = np.append(yvals,yvals.copy()+b).reshape(1,-1)
    node_pts = np.stack((xpts[0],ypts[0]),axis=-1)
    elements = define_element_connectivity(node_pts).simplices.copy()
    values = default_values(node_pts)
    return (node_pts,elements,values)

def default_values(node_pts):
    values = np.zeros(np.shape(node_pts)[0])
    return values

def generate_simple_square_mesh(width,length):
    x = np.linspace(0,width, width + 1)
    y = np.linspace(0, length, length + 1)
    xvals,yvals = np.meshgrid(x,y)
    
    
    node_pts = np.stack((xvals[0],yvals[0]),axis=-1)
    elements = define_element_connectivity(node_pts).simplices.copy()
    values = np.zeros(np.shape(node_pts)[0])
    return (node_pts, elements, values)

def save_mesh_to_file(node_pts,elements,filename=None, vals=None):
    #http://christopherlovell.co.uk/blog/2016/04/27/h5py-intro.html
    if filename is None:
        filename = save_filename(directory='/media/ppzmis/data/FEM/', file_filter='*.h5')

    hf = h5py.File(filename, 'w')
    hf.create_dataset('node_pts', data=str(node_pts))
    hf.create_dataset('elements', data=str(elements))
    if vals != None:
        hf.create_dataset('vals', data=vals)
    
    hf.close()    
    return filename

def save_mesh_constraints(known_displacements,known_forces,filename=None):
    if filename is None:
        filename = save_filename(directory='/media/ppzmis/data/FEM/', file_filter='*.h5')
    op_filename = filename[:-3] + 'constraints.h5'
    hf = h5py.File(op_filename, 'w')
    hf.create_dataset('displacements', data=known_displacements)
    hf.create_dataset('forces', data=known_forces)
    hf.close()
    return op_filename

def load_mesh_constraints(op_filename=None,directory='/home/ppzmis/Documents/PythonScripts/FEM'):
    if op_filename is None:
        op_filename = load_filename(directory=directory, file_filter='*.h5')
    hf = h5py.File(op_filename, 'r')
    known_displacements = hf.get('displacements')
    known_forces = hf.get('forces')
    return known_displacements, known_forces


def load_mesh_from_file(filename=None,directory='/home/ppzmis/Documents/PythonScripts/FEM'):
    if filename is None:
        filename=load_filename(directory=directory,file_filter='*.h5')
    hf = h5py.File(filename, 'r')
    
    node_pts = np.array(hf.get('node_pts'))
    elements = np.array(hf.get('elements'))
    if 'vals' in hf.keys():
        values = np.array(hf.get('vals'))
    else:
        values = np.ones(np.shape(elements)[0])
    
    return (node_pts,elements,values)



def define_element_connectivity(points):
    #find triangle elements nodal indices
    connectivity = Delaunay(points)#,incremental=True)
    return connectivity


def show_mesh(points, triangles, values, display_type = 'surface'):
    '''
    This function visualises the mesh
    displayType can be 'surface' or 'wireframe'
    '''
    t=np.ones(np.shape(points)[0])
    colors = np.random.uniform(0,1,size=np.shape(triangles)[0])
    
    mlab.figure(size = (1024,768),bgcolor = (0,0,0), fgcolor = (0.5, 0.5, 0.5))
    mesh = mlab.triangular_mesh(points[:,0],points[:,1],t, triangles, representation=display_type,opacity=0,color=colors)#displayType)
    
    mesh.mlab_source.dataset.cell_data.scalars = colors
    mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data'
    mesh.mlab_source.update()
    
    mesh2 = mlab.pipeline.set_active_attribute(mesh,cell_scalars='Cell data')
    mlab.pipeline.surface(mesh2)
    mesh2.scene.z_plus_view()
    mlab.show()
    
def length_tri_side(node1,node2):
    
    length_val = ((node1[0]-node2[0])**2+(node1[1]-node2[1])**2)**0.5
    return length_val
    
def refine_mesh_element(node_pts,elements,values,refine_element_id):
    '''
    We create a node on the longest triangular side we then create 4 new triangles in place of the 2 original ones.
    If the element is at the edge we only create 2.
        
    The value of the new node is the average of the two ends of the line which is divided.
    '''
    triangle_node_indices = elements[refine_element_id].copy()

    #Which side is longest?
    l = np.array([length_tri_side(node_pts[triangle_node_indices[0]],node_pts[triangle_node_indices[1]]),length_tri_side(node_pts[triangle_node_indices[1]],node_pts[triangle_node_indices[2]]),length_tri_side(node_pts[triangle_node_indices[2]],node_pts[triangle_node_indices[0]])])
    longest_side_index=np.argmax(l)
    if longest_side_index < 2:
        longest_side_b = longest_side_index + 1
    else:
        longest_side_b = 0
    
    #create new node coordinates
    new_node = [0.5*(node_pts[triangle_node_indices[longest_side_index]]-node_pts[triangle_node_indices[longest_side_b]])+node_pts[triangle_node_indices[longest_side_b]]]
    new_node_index = np.shape(node_pts)[0]
    
    #Find other triangle which shares this edge. This triangle shares 2 nodes with the element under refinement
    common_a = triangle_node_indices[longest_side_index]
    common_b = triangle_node_indices[longest_side_b]
    common_nodes = np.array([common_a,common_b])
    a = np.where(elements == common_a)[0]
    b = np.where(elements == common_b)[0]

    shared_edge_triangles = np.intersect1d(a,b)
    
    if np.shape(shared_edge_triangles)[0]==1:
        #This if statement is triggered for triangles at a boundary
        other_node_a = np.setdiff1d(elements[shared_edge_triangles[0]],common_nodes)[0]
        indices_triangle1 = np.array([common_a, new_node_index,other_node_a])
        indices_triangle2 = np.array([ other_node_a,common_b, new_node_index])
        elements[refine_element_id] = indices_triangle1
        indices_triangles = [indices_triangle2]
        #Give all the new elements the same value as the original.
        values = np.append(values,[values[refine_element_id]])
     
            

    else:
        #Index of other triangle that is being bisected
        other_element_id = np.setdiff1d(shared_edge_triangles,np.array([refine_element_id]))[0]
        
        #This is the normal condition
        other_node_a = np.setdiff1d(elements[shared_edge_triangles[0]],common_nodes)[0]
        other_node_b = np.setdiff1d(elements[shared_edge_triangles[1]],common_nodes)[0]
    
        indices_triangle1 = np.array([common_a,new_node_index, other_node_a ])
        indices_triangle2 = np.array([common_a, other_node_b,  new_node_index])
        indices_triangle3 = np.array([new_node_index,common_b,other_node_a,])
        indices_triangle4 = np.array([common_b,new_node_index, other_node_b])

        elements[refine_element_id] = indices_triangle1
        elements[other_element_id] = indices_triangle2
        indices_triangles = np.stack((indices_triangle3,indices_triangle4),axis=0)
        
        #Give all the new elements the same value as the original.
        values = np.append(values,[values[refine_element_id],values[other_element_id]])
    
    elements = np.concatenate((elements,indices_triangles),axis=0)
        
    #Add new node to the end of the list
    nodePts = np.concatenate((node_pts,new_node))
   
    
    return(nodePts,elements,values)
    


def find_tri_with_node(elements,node_id):
    '''
    Finds all elements containing nodeId
    '''
    indices = np.where(elements == node_id )
    return indices[0]
    
def min_length_side(nodepts,triangle_node_indices):
    '''
    Calculates the shortest side of the triangle and returns the node indices at either end and the mid point coords
    '''
    pts = nodepts[triangle_node_indices]

    length_vals = np.array([((pts[0,0]-pts[1,0])**2 + (pts[0,1]-pts[1,1])**2)**0.5,((pts[1,0]-pts[2,0])**2 + (pts[1,1]-pts[2,1])**2)**0.5,((pts[2,0]-pts[0,0])**2 + (pts[2,1]-pts[0,1])**2)**0.5])
    min_length_index = np.argmin(length_vals)

    if min_length_index < 2:
        node_indices = (triangle_node_indices[min_length_index],triangle_node_indices[min_length_index + 1])
    else:
        node_indices = (triangle_node_indices[min_length_index],triangle_node_indices[0])
    
    midpoint_coords = nodepts[node_indices[0]] + (nodepts[node_indices[1]] - nodepts[node_indices[0]])/2
    return (node_indices, midpoint_coords)
    
    
def coarsen_mesh_element(nodepts,elements,values,refine_element_id):
    '''
    coarsenMesh removes Elements that have been selected to be coarsened.
    It does this by picking a side from the Element = refineElementId and replacing the two terminating nodes
    with a node halfway along the side.
    
    '''
    #Find indices of nodes at corners of triangles
    triangle_node_indices = elements[refine_element_id].copy()
    node_indices_replace, new_node_coords = min_length_side(nodepts,triangle_node_indices)
    
    #Need to find triangles with the nodeIndicesReplace. If triangles contain one node 
    #nodeIndicesReplace then swap node for new node. If triangles contain both delete triangle and its associated value
    common_node_tri_indices1=find_tri_with_node(elements,node_indices_replace[0])
    common_node_tri_indices2=find_tri_with_node(elements,node_indices_replace[1])

    tri_contains_both_node = np.intersect1d(common_node_tri_indices1,common_node_tri_indices2)
    tri_contains_one_node = np.setdiff1d(np.append(common_node_tri_indices1,common_node_tri_indices2),tri_contains_both_node)
    
    #newNodeCoords in the nodepts at position of lowest nodeindex    
    if node_indices_replace[0] < node_indices_replace[1]:
        nodepts[node_indices_replace[0]] = new_node_coords
        #Make both node index values the same as the lower value.
        elements[np.where(elements == node_indices_replace[1])] = node_indices_replace[0]
        #Have to shift indices contained in elements to correspond to the fact we have one less node
        elements[np.where(elements > node_indices_replace[1])] = elements[np.where(elements > node_indices_replace[1])] - 1
        nodepts = np.delete(nodepts,node_indices_replace[1],0)
        
    else:
        nodepts[node_indices_replace[1]] = new_node_coords
        #Make both node index values the same as the lower value.
        elements[np.where(elements == node_indices_replace[0])] = node_indices_replace[1]
        elements[np.where(elements > node_indices_replace[0])] = elements[np.where(elements > node_indices_replace[0])] - 1
        nodepts = np.delete(nodepts,node_indices_replace[0],0)
   
    #Need to remove triangles which have duplicated nodes and the corresponding value in the value array
    elements=np.delete(elements,tri_contains_both_node,0)
    values=np.delete(values,tri_contains_both_node,0)
    return (nodepts, elements, values)
    

if __name__ == '__main__':
    #Create the intial mesh node pts
    nodepts,elements,values=generate_hex_mesh(20,30)
    save_mesh_to_file(nodepts,elements)
    #Define the element triangles that link the nodes

    nodepts,elements,values = load_mesh_from_file('/home/ppzmis/Documents/PythonScripts/FEM/test.h5')

    display.tri_mesh_plot(nodepts,elements,values)


    nodepts,elements,values = refine_mesh_element(nodepts,elements,values,6)

    for i in range(5):
        index = np.random.randint(0,np.shape(elements)[0])
        nodepts,elements,values = refine_mesh_element(nodepts,elements,values,index)
    
    for i in range(5):
        index = np.random.randint(0,np.shape(elements)[0])
        nodepts,elements,values=coarsen_mesh_element(nodepts,elements,values,index)
  


    
    #display.matplotlib_helper_plot(nodepts,elements)
    
    
