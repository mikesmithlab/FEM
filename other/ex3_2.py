from femUtils import *
from elementUtils import *



#Define connectivity of nodes
numbNodes = 6
connectivity = [[0,1],[1,2],[2,3],[3,4],[4,5]]
numbElements = np.shape(connectivity)[0]

E=[210e6]*numbElements
A=[0.002+0.01*(0.6*a+0.3)/3 for a in range(numbElements)]
L = [0.6]*numbElements



#Define spring elements
#Simpler to inject straight into assembly of stiffness matrix


#Assemble global stiffness matrix
K_global=np.zeros((numbNodes,numbNodes),dtype=float)
for i in range(numbElements):
    K_global = linearBarAssemble(K_global,linearBarElement(E[i],A[i],L[i]),connectivity[i][0],connectivity[i][1])



#Known Displacements
knownDispIndex = np.array([5])
knownDisplacements = np.array([0])
#Initialise the Displacement Vector
U = initialiseU(K_global,knownDispIndex, knownDisplacements)


#Known Forces
knownForceIndex = np.array([0,1,2,3,4])
knownForces = np.array([-18,0,0,0,0])
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
    k=linearBarElement(E[i],A[i],L[i])
    u_element = linearBarElementDisplacements(U_global,connectivity,i)   
    u_element_store[i] = u_element[0]
    f_element[i] = linearBarElementForces(k,u_element)
    stress_element[i] = linearBarElementStresses(k,u_element,A[i])
    

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

