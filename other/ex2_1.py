from femUtils import *

#Define equilibrium positions of nodes in mesh
x = np.array([0,1,2])

#Define connectivity of nodes
numbNodes = 3

#Define spring elements
k1 = springElement(100)
k2 = springElement(200)


#Assemble global stiffness matrix
K_global=np.zeros((numbNodes,numbNodes),dtype=float)
K_global = springAssemble(K_global,k1,0,1)
K_global = springAssemble(K_global,k2,1,2)



#Known Forces
knownForceIndex = np.array([1,2])
knownForces = np.array([0,15])

F = initialiseF(K_global,knownForceIndex, knownForces)



#Known Displacements
knownDispIndex = np.array([0])
knownDisplacements = np.array([0])


#Initialise the
U = initialiseU(K_global,knownDispIndex, knownDisplacements)





#Calculate the partition matrix and perform gaussian elimination
u,f=partitionMatrix(K_global,U,F)


#Recombine to give nodal displacements and forces
U_global = globalDisplacement(u,U)
F_global = globalForce(K_global,U_global)


#Print out the final answers
print('Global Stiffness Matrix')
print(K_global)
print('Complete Nodal Displacements')
print(U_global)
print('Complete list of Nodal Forces')
print(F_global)



