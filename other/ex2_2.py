from femUtils import *

#Define equilibrium positions of nodes in mesh
x = np.array([0,1,2])

#Define connectivity of nodes
numbNodes = 5
connectivity = [[0,2],[2,3],[2,4],[2,4],[4,3],[3,1]]
numbElements = np.shape(connectivity)[0]
kval=120
k_element = np.reshape([kval]*numbElements,(-1,1))




#Define spring elements
k = springElement(kval)


#Assemble global stiffness matrix
K_global=np.zeros((numbNodes,numbNodes),dtype=float)
for i in range(numbElements):
    K_global = springAssemble(K_global,k,connectivity[i][0],connectivity[i][1])


#Known Displacements
knownDispIndex = np.array([0,1])
knownDisplacements = np.array([0,0])


#Initialise the Displacement Vector
U = initialiseU(K_global,knownDispIndex, knownDisplacements)


#Known Forces
knownForceIndex = np.array([2,3,4])
knownForces = np.array([0,0,20])

F = initialiseF(K_global,knownForceIndex, knownForces)


#Calculate the partition matrix and perform gaussian elimination
u,f=partitionMatrix(K_global,U,F)


#Recombine to give nodal displacements and forces
U_global = globalDisplacement(u,U)
F_global = globalForce(K_global,U_global)


#Calculate the element forces
f_element = np.zeros((numbElements,1))
for i in range(numbElements):
	u_element = U_global[connectivity[i][0]] - U_global[connectivity[i][1]]
	f_element[i]=springElementForces(k_element[i],u_element)



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
print(u_element)
#Element forces
print('\nElement Forces')
print(f_element)


