from femUtils import *
from elementUtils import *

#Close all open figures
plt.close('all')

#Define connectivity of nodes
numbNodes = 5
node_dof = 2
connectivity = [[0,1,4],[4,1,2],[3,4,2],[0,4,3]]
numbElements = np.shape(connectivity)[0]
nodeCoords = np.array([[0,0],[0.5,0],[0.5,0.25],[0,0.25],[0.25,0.125]])

E=[210e6]*numbElements
t=[0.025]*numbElements


#Assemble global stiffness matrix
K_global=np.zeros((numbNodes*node_dof,numbNodes*node_dof),dtype=float)

for i in range(numbElements):
    K_global = LinearTriangleAssemble(K_global,LinearTriangleElementStiffness(E[i],t[i],nodeCoords[connectivity[i]]),connectivity[i][0],connectivity[i][1],connectivity[i][2])

    
#Modify for inclined support
#T = np.identity(numbNodes*node_dof)
#T=PlaneTrussInclinedSupport(T,3,45)
#K_global = np.matmul(T,K_global)
#K_global = np.matmul(K_global,np.transpose(T))


#Known Displacements
knownDispIndex = np.array([0,1,6,7])
knownDisplacements = np.array([0.0,0.0,0.0,0.0])
#Initialise the Displacement Vector
U = initialiseU(K_global,knownDispIndex, knownDisplacements)


#Known Forces
knownForceIndex = np.array([2,3,4,5])
knownForces = np.array([9.375,0.0,9.375,0.0])
F = initialiseF(K_global,knownForceIndex, knownForces)


#Perform boundary conditions check
bad,errorMessage=checkBoundaryConditions(numbNodes,knownDispIndex,knownDisplacements,knownForceIndex,knownForces)
if bad:
    print(errorMessage)

#Calculate the partition matrix and perform gaussian elimination
u,f=partitionMatrix(K_global,U,F)


#Recombine to give nodal displacements and forces
U_global = globalDisplacement(u,U)
F_global = globalForce(K_global,U_global)

#print(K_global)
#print(U_global)
#print(F_global)

#Calculate the element forces
u_element_store = np.zeros((numbElements,1))
f_element = np.zeros((numbElements,1))
stress_element = np.zeros((numbElements,1))
for i in range(numbElements):
    u = LinearTriangleElementDisplacements(U_global, connectivity, i)
    sigma = LinearTriangleElementStresses(E[i],nodeCoords[connectivity[i]],u)
    sigmaP = LinearTriangleElementPStresses(sigma)
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



