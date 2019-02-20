import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from mayavi import mlab

from meshpy.triangle import MeshInfo, build, write_gnuplot_mesh


def createInitialTriGridPts(width,length):
    #N must be an even number
    N = width
    x = np.linspace(0, width, N + 1 )
    y = np.linspace(0, length, N/2 + 1)
    
    xvals, yvals = np.meshgrid(x,y)
    a = (x[2]-x[1])/2.0
    b = (y[2]-y[1])/2.0
    
    xpts = np.append(xvals,xvals.copy()+a).reshape(1,-1)
    ypts = np.append(yvals,yvals.copy()+b).reshape(1,-1)
    nodePts = np.stack((xpts[0],ypts[0]),axis=-1)
    
    return nodePts


def defineElementConnectivity(points):
    #find triangle elements nodal indices
    connectivity = Delaunay(points,incremental=True)
    return connectivity


def showMesh(points,triangles,vals = None, displayType = 'surface'):
    '''
    This function visualises the mesh
    displayType can be 'surface' or 'wireframe'
    '''
    
    z=np.zeros(np.shape(points)[0])
    t=np.random.rand(np.shape(points)[0])
    
    mlab.figure(size = (1024,768),bgcolor = (0,0,0), fgcolor = (0.5, 0.5, 0.5))
    tri = mlab.triangular_mesh(points[:,0],points[:,1],z , triangles, scalars=t,representation=displayType)
    tri.scene.z_plus_view()
    mlab.show()
    


def refineMesh(nodepts,elements,refineElementId):
    '''
    We create 4 elements in place of 1 element by subdividing each vertex
    of triangle in 2 and creating a node at each point
    These additional nodes have numbers which are larger than any current values
    
    The central element of the 4 new elements retains the original element number
    and the node values are updated. The new element numbers start at the end of the current list.
    '''
    
    triangle_nodes = nodepts[elements.simplices[refineElementId]]
    numbElements = np.shape(elements.simplices)[0]
    
    #We generate the coordinates of the new nodes
    newnodea = np.array([triangle_nodes[0][0]-0.5*(triangle_nodes[0][0]-triangle_nodes[1][0]),triangle_nodes[0][1]-0.5*(triangle_nodes[0][1]-triangle_nodes[1][1])])
    newnodeb = np.array([triangle_nodes[1][0]-0.5*(triangle_nodes[1][0]-triangle_nodes[2][0]),triangle_nodes[1][1]-0.5*(triangle_nodes[1][1]-triangle_nodes[2][1])])
    newnodec = np.array([triangle_nodes[2][0]-0.5*(triangle_nodes[2][0]-triangle_nodes[0][0]),triangle_nodes[2][1]-0.5*(triangle_nodes[2][1]-triangle_nodes[0][1])])
    #Add those new coordinates to the node list
    newnodes = np.stack((newnodea,newnodeb,newnodec))
    nodepts = np.concatenate((nodepts,newnodes))

    #Add the triangular element to the list of simplices
    elements.add_points(newnodes)
    
    return(nodepts,elements)
    
    
    
def cleanMesh(nodepts,elements,scalarVals,threshold,condition='stress'):
    '''
    cleanMesh is run every once in a while to coarsen the mesh where it no longer needs
    to be fine by removing small triangles that do not match the condition
    
    scalarVals should be a numpy array with length = length of elements.simplices
    '''
    #Find those elements which need coarsening
    coarsen_Ids = np.where(scalarVals < threshold)
    
    #find one of the nodes in the element
    nodes_remove_ID = np.unique(elements.simplices[coarsenIds][0])
    
    #redraw the delaunay triangles
    
    
    
    
    
    
    
    
    
    return findDelaunayTriangles(nodepts)
    
    




    
if __name__ == '__main__':
    #Create the intial mesh node pts
    nodepts=createInitialTriGridPts(4,4)
    #Define the element triangles that link the nodes
    elements = defineElementConnectivity(nodepts)
    
    mesh_info = MeshInfo()
    print(nodepts)
    mesh_info.set_points(nodepts)
    print(elements.simplices)
    mesh_info.set_facets(elements.simplices)
    mesh = build(mesh_info)
    #mesh.triangle.
    write_gnuplot_mesh("test.vtk", mesh_info.points, facets=False)    
    
    #triangle=showMesh(nodepts,elements.simplices)
    print(elements)
    print(elements.simplices)
    print(elements.points)
    
    
    nodepts,elements = refineMesh(nodepts,elements,7)
    nodepts,elements = refineMesh(nodepts,elements,20)
    #nodepts,elements = refineMesh(nodepts,elements,100)
    #nodepts,elements = refineMesh(nodepts,elements,50)
    #nodepts,elements = refineMesh(nodepts,elements,2)
    #nodepts,elements = refineMesh(nodepts,elements,10)
    
    triangle=showMesh(nodepts,elements.simplices,displayType='wireframe')
    neighbours=elements.neighbors[7]
    
    
    nodepts,elements,scalarVals = cleanMesh(nodepts,elements,scalarVals,threshold,condition='stress')
    
    print(nodepts[elements.simplices[7]])
    print(neighbours)
    print(nodepts[elements.simplices[neighbours]])
    
    
    
    