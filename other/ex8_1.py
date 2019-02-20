from femUtils import *
from elementUtils import *

#Close all open figures
plt.close('all')

#Define connectivity of nodes
numbNodes = 4
node_dof = 3
connectivity = [[0,1],[1,2],[2,3]]
numbElements = np.shape(connectivity)[0]

E=[210e6]*numbElements
A=[2e-2]*numbElements
I=[5e-5]*numbElements
L = [PlaneFrameElementLength(0,0,0,3),PlaneFrameElementLength(0,3,4,3),PlaneFrameElementLength(4,3,4,0)]
theta = [90,0,270]



#Assemble global stiffness matrix
K_global=np.zeros((numbNodes*node_dof,numbNodes*node_dof),dtype=float)

for i in range(numbElements):
    K_global = PlaneFrameAssemble(K_global,PlaneFrameElementStiffness(E[i],A[i],I[i],L[i],theta[i]),connectivity[i][0],connectivity[i][1])
   
#Modify for inclined support
#T = np.identity(numbNodes*node_dof)
#T=PlaneTrussInclinedSupport(T,3,45)
#K_global = np.matmul(T,K_global)
#K_global = np.matmul(K_global,np.transpose(T))


#Known Displacements
knownDispIndex = np.array([0,1,2,9,10,11])
knownDisplacements = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
#Initialise the Displacement Vector
U = initialiseU(K_global,knownDispIndex, knownDisplacements)


#Known Forces
knownForceIndex = np.array([3,4,5,6,7,8])
knownForces = np.array([-20.0,0.0,0.0,0.0,0.0,12.0])
F = initialiseF(K_global,knownForceIndex, knownForces)


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
    k = PlaneFrameElementStiffness(E[i],A[i],I[i],L[i],theta[i])
    u_element = PlaneFrameElementDisplacements(U_global,connectivity,i)  
    u_element_store[i] = u_element[0]
    f_element = PlaneFrameElementForces(E[i],A[i],I[i],L[i],theta[i],u_element)
    PlaneFrameElementPlots(f_element,L[i],i)    
   
    #stress_element[i] = PlaneFrameElementStress(E[i],L[i],theta[i],u_element)

    
    

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

