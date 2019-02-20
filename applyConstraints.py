import femUtils as fem
import elementUtils as el
import meshUtils as mesh
import numpy as np

#nodepts,elements,values=mesh.generate_hex_mesh(20,30)
#filename = mesh.save_mesh_to_file(nodepts,elements)


#Define the element triangles that link the nodes
filename = '/home/ppzmis/Documents/PythonScripts/FEM/test.h5'

node_coords,connectivity,values = mesh.load_mesh_from_file(filename)

#Define connectivity of nodes
numb_nodes = np.shape(node_coords)[0]
node_dof = 2
numb_elements = np.shape(connectivity)[0]
default_vals = np.ones(np.shape(connectivity)[0])

E=[210e6]*numb_elements
t=[0.025]*numb_elements


#Assemble global stiffness matrix
K_global=np.zeros((numb_nodes*node_dof,numb_nodes*node_dof),dtype=float)
for i in range(numb_elements):
    K_global = el.linear_triangle_assemble(K_global,el.linear_triangle_element_stiffness(E[i],t[i],node_coords[connectivity[i]]),connectivity[i][0],connectivity[i][1],connectivity[i][2])


#Known Displacements
known_disp_x, known_disp_index_x = fem.known_node_val(node_coords, connectivity, 0,shape='rectangle')
known_disp_y, known_disp_index_y = fem.known_node_val(node_coords, connectivity, 0,shape='rectangle')

#Known Forces
known_force_x, known_force_index_x = fem.known_node_val(node_coords, connectivity, 0)
known_force_y, known_force_index_y = fem.known_node_val(node_coords, connectivity, 0)

displacements = [known_disp_x,known_disp_index_x,known_disp_y,known_disp_index_y]
#displacements = {'known_disp_x':known_disp_x, 'known_disp_index_x':known_disp_index_x, 'known_disp_y':known_disp_y, 'known_disp_index_y':known_disp_index_y}
forces = [known_force_x,known_force_index_x,known_force_y,known_force_index_y]
#forces = {'known_force_x':known_force_x, 'known_force_index_x':known_force_index_x, 'known_force_y':known_force_y, 'known_force_index_y':known_force_index_y}

mesh.save_mesh_constraints(displacements,forces,filename=filename)


#Initialise the Displacement Vector
U= fem.initialiseU(K_global,known_disp_index_x, known_disp_x,known_disp_index_y, known_disp_y)
#print(U)




F = fem.initialiseF(K_global,known_force_index_x, known_force_x,known_force_index_y, known_force_y)

#Calculate the partition matrix and perform gaussian elimination
u,f=fem.partition_matrix(K_global,U,F)
#print('partitionedK')
#print(u)


#Recombine to give nodal displacements and forces
U_global = fem.global_displacement(u,U)
F_global = fem.global_force(K_global,U_global)

#print(K_global)
#print(U_global)
#print(F_global)

#Calculate the element forces
u_element_store = np.zeros((numb_elements,1))
f_element = np.zeros((numb_elements,1))
stress_element = np.zeros((numb_elements,1))
for i in range(numb_elements):
    u = el.linear_triangle_element_displacements(U_global, connectivity, i)
    sigma = el.linear_triangle_element_stresses(E[i],node_coords[connectivity[i]],u)
    sigmaP = el.linear_triangle_element_pstresses(sigma)
    print('u_element')
    print(u)
    print('sigma_element')
    print(sigma)
    print('sigma principal')
    print(sigmaP)
    
    
    

print('\n ------------------------------------------------------------------------\n')
print( __file__)

#Print out the final answers
print('\nGlobal Stiffness Matrix')
print(K_global)
print('\nComplete Nodal Displacements')
print(U_global)
print('\nComplete list of Nodal Forces')
print(F_global)



