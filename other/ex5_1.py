from femUtils import *
from elementUtils import *




#Define connectivity of nodes
numbNodes = 3
node_dof = 2
connectivity = [[0,1],[0,2],[1,2]]
numbElements = np.shape(connectivity)[0]

E=[210e6]*numbElements
A=[1e-4]*numbElements
L = [4, PlaneTrussElementLength(0,0,2,3),PlaneTrussElementLength(4,0,2,3)]
theta = [0,56.3099,123.6901]



#Define spring elements
#Simpler to inject straight into assembly of stiffness matrix


#Assemble global stiffness matrix
K_global=np.zeros((numbNodes*node_dof,numbNodes*node_dof),dtype=float)
for i in range(numbElements):
    K_global = PlaneTrussAssemble(K_global,PlaneTrussElementStiffness(E[i],A[i],L[i],theta[i]),connectivity[i][0],connectivity[i][1])
    


#Known Displacements
knownDispIndex = np.array([0,1,3])
knownDisplacements = np.array([0,0,0])
#Initialise the Displacement Vector
U = initialiseU(K_global,knownDispIndex, knownDisplacements)


#Known Forces
knownForceIndex = np.array([2,4,5])
knownForces = np.array([0,5,-10])
F = initialiseF(K_global,knownForceIndex, knownForces)


#Calculate the partition matrix and perform gaussian elimination
u,f=partitionMatrix(K_global,U,F)


#Recombine to give nodal displacements and forces
U_global = globalDisplacement(u,U)
F_global = globalForce(K_global,U_global)


#Calculate the element forces
u_element_store = np.zeros((numbElements,1))
f_element = np.zeros((numbElements,1))
stress_element = np.zeros((numbElements,1))
for i in range(numbElements):
    k = PlaneTrussElementStiffness(E[i],A[i],L[i],theta[i])
    u_element = PlaneTrussElementDisplacements(U_global,connectivity,i)  
    u_element_store[i] = u_element[0]
    f_element[i] = PlaneTrussElementForce(E[i],A[i],L[i],theta[i],u_element)
    stress_element[i] = PlaneTrussElementStress(E[i],L[i],theta[i],u_element)
    

print('\n ------------------------------------------------------------------------\n')
print( __file__)

#Print out the final answers
print('\nGlobal Stiffness Matrix')
print(K_global)
print('\nComplete Nodal Displacements')
print(U_global)
print('\nComplete list of Nodal Forces')
print(F_global)

#Element displacements
print('\nElement Displacements')
print(u_element_store)
#Element forces
print('\nElement Forces')
print(f_element)
#Element stresses
print('\nElement Stresses')
print(stress_element)

