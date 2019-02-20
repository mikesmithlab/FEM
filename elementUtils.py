import numpy as np
import matplotlib.pyplot as plt

#Global constants
poissons = 0.3


'''
-------------------------------------------------------------------------------------------
The Spring Element - A simple spring connecting 2 nodes with displacements only in x
-------------------------------------------------------------------------------------------

'''

def springElement(k):
	'''
	This fn returns the element stiffness matrix for a spring with stiffness k
	'''
	return np.array([[k, -k],[-k, k]])


def springAssemble(K,k,i,j):
	'''
	This function assembles the element stiffness matrix k of the spring nodes i and j into global matrix K
	'''
	K[i,i] = K[i,i]+k[0,0]
	K[i,j] = K[i,j]+k[0,1]
	K[j,i] = K[j,i]+k[1,0]
	K[j,j] = K[j,j]+k[1,1]
	return K
	
def springElementForces(k,u):
	'''
	returns element nodal force vector given element stiffness amtrix k and element nodal displacements
	'''
	return k*u


'''
--------------------------------------------------------------------------------------------------
The Linear Bar Element - 1D Same as spring but stiffness given in terms of E, A, L
--------------------------------------------------------------------------------------------------
'''

def linearBarElement(E,A,L):
	'''
	This fn returns the element stiffness matrix for a spring with stiffness kalala
	'''
	n = E*A/L
	return np.array([[n, -n],[-n, n]])

def linearBarAssemble(K,k,i,j):
	'''
	This function assembles the element stiffness matrix k of the spring nodes i and j into global matrix K
	'''
	K[i,i] = K[i,i]+k[0,0]
	K[i,j] = K[i,j]+k[0,1]
	K[j,i] = K[j,i]+k[1,0]
	K[j,j] = K[j,j]+k[1,1]
	return K

def linearBarElementDisplacements(U_global, connectivity, i):
    displacement = [[U_global[connectivity[i][0]][0]],[U_global[connectivity[i][1]][0]]]
    return displacement
	
def linearBarElementForces(k,u):
	'''
	returns element nodal force vector given element stiffness amtrix k and element nodal displacements
	negative forces indicate tension.
	'''
	element_forces = np.matmul(k,u)
	return element_forces[0]

def linearBarElementStresses(k,u,A):
    '''
	A negative returned stress is a tensile stress
	'''
	
    sz = np.shape(u)
    Area = np.ones(sz)*A
    u_over_A = np.divide(u,Area)
    element_stresses = np.matmul(k,u_over_A)
    
    return element_stresses[0]

'''
--------------------------------------------------------------------------------------------------
The Plane Truss Element - 2D
--------------------------------------------------------------------------------------------------
'''
def PlaneTrussElementLength(x1,y1,x2,y2):
    '''
    Returns the length of the Truss element.
    '''
    length = np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) 
    return length

def PlaneTrussElementStiffness(E,A,L,theta):
    '''     
    Returns element stiffness matrix, Size of matrix is 4x4
    '''
    x=theta*np.pi/180.0
    C = np.cos(x)
    S = np.sin(x)
    k = (E*A/L)*np.array([[C*C, C*S, -C*C, -C*S],[C*S, S*S, -C*S, -S*S],[-C*C, -C*S, C*C, C*S],[-C*S, -S*S, C*S, S*S]])
    return k

def PlaneTrussAssemble(K,k,i,j):
    '''
    Assembles the element stiffness matrices k into global stiffness matrix K
    '''
    K[2*i,2*i] = K[2*i,2*i] + k[0,0]
    K[2*i,2*i+1] = K[2*i,2*i+1]     + k[0,1]
    K[2*i,2*j] = K[2*i,2*j]+ k[0,2]
    K[2*i,2*j+1] = K[2*i,2*j+1]     + k[0,3]
    K[2*i+1,2*i] = K[2*i+1,2*i]     + k[1,0]
    K[2*i+1,2*i+1] = K[2*i+1,2*i+1]         + k[1,1]
    K[2*i+1,2*j] = K[2*i+1,2*j]     + k[1,2]
    K[2*i+1,2*j+1] = K[2*i+1,2*j+1]         + k[1,3]
    K[2*j,2*i] = K[2*j,2*i] + k[2,0]
    K[2*j,2*i+1] = K[2*j,2*i+1]     + k[2,1]
    K[2*j,2*j] = K[2*j,2*j] + k[2,2]
    K[2*j,2*j+1] = K[2*j,2*j+1]     + k[2,3]
    K[2*j+1,2*i] = K[2*j+1,2*i]     + k[3,0]
    K[2*j+1,2*i+1] = K[2*j+1,2*i+1]         + k[3,1]
    K[2*j+1,2*j] = K[2*j+1,2*j]     + k[3,2]
    K[2*j+1,2*j+1] = K[2*j+1,2*j+1]         + k[3,3]
    return K
    
def PlaneTrussElementForce(E,A,L,theta,u):
    '''
    returns element force
    '''
    x = theta* np.pi/180.0
    C=np.cos(x)
    S=np.sin(x)
    return (E*A/L)*np.matmul(np.array([-C, -S, C, S]),u)
  

def PlaneTrussElementStress(E,L,theta,u):
    '''
    returns element stress
    '''
    x = theta* np.pi/180.0
    C=np.cos(x)
    S=np.sin(x)
    return (E/L)*np.matmul(np.array([-C, -S, C, S]),u)
    

def PlaneTrussElementDisplacements(U_global, connectivity, i):
    dof = 2
    displacement = [[U_global[dof*connectivity[i][0]][0]],[U_global[dof*connectivity[i][0]+1][0]],[U_global[dof*connectivity[i][1]][0]],[U_global[dof*connectivity[i][1]+1][0]]]
    return displacement
    
def PlaneTrussInclinedSupport(T,i,alpha):
    '''
    Calcs transformation matrix of inclined support at node i with angle inclination alpha in degrees
    '''
    x = alpha *np.pi / 180.0
    T[2*i,2*i] = np.cos(x)
    T[2*i,2*i+1]   = np.sin(x)
    T[2*i+1,2*i]   = -np.sin(x)
    T[2*i+1,2*i+1]     = np.cos(x)
    return T

'''
--------------------------------------------------------------------------------------------------
The Plane Frame Element - 2D
--------------------------------------------------------------------------------------------------
'''

def PlaneFrameElementLength(x1,y1,x2,y2):
    '''
    Returns the length of the Truss element.
    '''
    length = np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) 
    return length

def PlaneFrameElementStiffness(E,A,I,L,theta):
    '''     
    Returns element stiffness matrix, Size of matrix is 6x6
    '''
    x=theta*np.pi/180.0
    C = np.cos(x)
    S = np.sin(x)
    w1=A*C*C + 12*I*S*S/(L*L)
    w2=A*S*S + 12*I*C*C/(L*L)
    w3=(A-12*I/(L*L))*C*S
    w4=6*I*S/L
    w5=6*I*C/L
    k = (E/L)*np.array([[w1, w3, -w4, -w1, -w3, -w4],[w3, w2, w5, -w3, -w2, w5],[-w4, w5, 4*I, w4,-w5, 2*I],[-w1, -w3, w4, w1, w3, w4],[-w3, -w2, -w5, w3, w2, -w5],[-w4, w5, 2*I, w4, -w5, 4*I]])
    return k

def PlaneFrameAssemble(K,k,i,j):
    '''
    Assembles the element stiffness matrices k into global stiffness matrix K
    '''
    K[3*i,3*i] = K[3*i,3*i] + k[0,0]
    K[3*i,3*i+1] = K[3*i,3*i+1]     + k[0,1]
    K[3*i,3*i+2] = K[3*i,3*i+2]+ k[0,2]
    K[3*i,3*j] = K[3*i,3*j]     + k[0,3]
    K[3*i,3*j+1] = K[3*i,3*j+1]     + k[0,4]
    K[3*i,3*j+2] = K[3*i,3*j+2]         + k[0,5]

    K[3*i+1,3*i] = K[3*i+1,3*i] + k[1,0]
    K[3*i+1,3*i+1] = K[3*i+1,3*i+1]     + k[1,1]
    K[3*i+1,3*i+2] = K[3*i+1,3*i+2]+ k[1,2]
    K[3*i+1,3*j] = K[3*i+1,3*j]     + k[1,3]
    K[3*i+1,3*j+1] = K[3*i+1,3*j+1]     + k[1,4]
    K[3*i+1,3*j+2] = K[3*i+1,3*j+2]         + k[1,5]
    
    K[3*i+2,3*i] = K[3*i+2,3*i] + k[2,0]
    K[3*i+2,3*i+1] = K[3*i+2,3*i+1]     + k[2,1]
    K[3*i+2,3*i+2] = K[3*i+2,3*i+2]+ k[2,2]
    K[3*i+2,3*j] = K[3*i+2,3*j]     + k[2,3]
    K[3*i+2,3*j+1] = K[3*i+2,3*j+1]     + k[2,4]
    K[3*i+2,3*j+2] = K[3*i+2,3*j+2]         + k[2,5]
    
    K[3*j,3*i] = K[3*j,3*i] + k[3,0]
    K[3*j,3*i+1] = K[3*j,3*i+1]     + k[3,1]
    K[3*j,3*i+2] = K[3*j,3*i+2]+ k[3,2]
    K[3*j,3*j] = K[3*j,3*j]     + k[3,3]
    K[3*j,3*j+1] = K[3*j,3*j+1]     + k[3,4]
    K[3*j,3*j+2] = K[3*j,3*j+2]         + k[3,5]

    K[3*j+1,3*i] = K[3*j+1,3*i] + k[4,0]
    K[3*j+1,3*i+1] = K[3*j+1,3*i+1]  + k[4,1]
    K[3*j+1,3*i+2] = K[3*j+1,3*i+2]+ k[4,2]
    K[3*j+1,3*j] = K[3*j+1,3*j]     + k[4,3]
    K[3*j+1,3*j+1] = K[3*j+1,3*j+1]     + k[4,4]
    K[3*j+1,3*j+2] = K[3*j+1,3*j+2]         + k[4,5]
    
    K[3*j+2,3*i] = K[3*j+2,3*i] + k[5,0]
    K[3*j+2,3*i+1] = K[3*j+2,3*i+1]     + k[5,1]
    K[3*j+2,3*i+2] = K[3*j+2,3*i+2]+ k[5,2]
    K[3*j+2,3*j] = K[3*j+2,3*j]     + k[5,3]
    K[3*j+2,3*j+1] = K[3*j+2,3*j+1]     + k[5,4]
    K[3*j+2,3*j+2] = K[3*j+2,3*j+2]         + k[5,5]
    
    return K

def PlaneFrameElementDisplacements(U_global, connectivity, i):
    dof = 3
    displacements = [[U_global[dof*connectivity[i][0]][0]],[U_global[dof*connectivity[i][0]+1][0]],[U_global[dof*connectivity[i][0]+2][0]],[U_global[dof*connectivity[i][1]][0]],[U_global[dof*connectivity[i][1]+1][0]],[U_global[dof*connectivity[i][1]+2][0]]]
    return displacements
    
def PlaneFrameElementForces(E,A,I,L,theta,u):
    '''
    returns element force
    '''
    x = theta* np.pi/180.0
    C=np.cos(x)
    S=np.sin(x)
    w1=(E*A/L)
    w2=12*E*I/(L*L*L)
    w3=6*E*I/(L*L)
    w4=4*E*I/L
    w5=2*E*I/L
    kprime = np.array([[w1, 0, 0, -w1,0,0],[0,w2,w3,0,-w2,w3],[0,w3,w4,0,-w3,w5],[-w1,0,0,w1,0,0],[0,-w2,-w3,0,w2,-w3],[0,w3,w5,0,-w3,w4]])
    T=np.array([[C,S,0,0,0,0],[-S,C,0,0,0,0],[0,0,1,0,0,0],[0,0,0,C,S,0],[0,0,0,-S,C,0],[0,0,0,0,0,1]])
    a = np.matmul(kprime,T)
    return np.matmul(a,u)

def PlaneFrameElementPlots(f,L,i):
    '''
    plots the axial force diagram
    '''
    
    x = [0,L]
    
    zaxial=[-f[0],f[3]]
    zshear=[f[1],-f[4]]
    zbending = [-f[2],f[5]]
    
    fig, axarr = plt.subplots(3, 1)
    axarr[0,].plot(x, zaxial,'r-')
    axarr[0,].set_xlabel('L')
    axarr[0,].set_ylabel('Axial Force')
    axarr[1,].plot(x, zshear,'b-')
    axarr[1,].set_xlabel('L')
    axarr[1,].set_ylabel('Shear Force')
    axarr[2,].plot(x, zbending,'g-')
    axarr[2,].set_xlabel('L')
    axarr[2,].set_ylabel('Bending Moment')
    fig.canvas.set_window_title('Force Diagrams for Element ' + str(i))
    fig.tight_layout()
    plt.show()
    return np.shape(axarr)
    
def PlaneFrameInclinedSupport(T,i,alpha):
    '''
    Calcs transformation matrix of inclined support at node i with angle inclination alpha in degrees
    '''
    x = alpha *np.pi / 180.0
    T[3*i,3*i]     = np.cos(x)
    T[3*i,3*i+1]   = np.sin(x)
    T[3*i,3*i+2]   = 0
    T[3*i+1,3*i]   = -np.sin(x)
    T[3*i+1,3*i+1] = np.cos(x)
    T[3*i+1,3*i+2] = 0
    T[3*i+2,3*i]   = 0
    T[3*i+2,3*i+1] = 0
    T[3*i+2,3*i+2] = 1
    
    return T


'''
--------------------------------------------------------------------------------------------------
The Linear Triangular Element 2D
--------------------------------------------------------------------------------------------------
'''


def linear_triangle_element_area(xi,yi,xj,yj,xm,ym):
    '''
    returns the area of the triangular element whose coords are listed counter clockwise and should give a positive area
    '''
    area = (xi*(yj-ym) + xj*(ym-yi)+xm*(yi-yj))/2
    return area

def linear_triangle_element_stiffness(E,t,coords,NU=poissons,p=1):
    '''     
    Returns element stiffness matrix, Size of matrix is 6x6
    E is youngs modulus
    t is thickness
    NU is poissons ratio
    p=1 --> plane stress
    p=2 --> plane strain
    '''
    xi=coords[0,0]
    yi=coords[0,1]
    xj=coords[1,0]
    yj=coords[1,1]
    xm=coords[2,0]
    ym=coords[2,1]
    
    
    A = linear_triangle_element_area(xi,yi,xj,yj,xm,ym)
    betai = yj-ym
    betaj = ym-yi
    betam = yi-yj
    gammai = xm-xj
    gammaj = xi-xm
    gammam = xj-xi
    B = (1/(2*A))*np.array([[betai, 0, betaj, 0, betam, 0],[0, gammai, 0, gammaj, 0, gammam],[gammai, betai, gammaj, betaj, gammam, betam]])
    BTranspose = np.transpose(B)   
    if p == 1:
        D = (E/(1-NU*NU))*np.array([[1, NU, 0],[NU, 1, 0],[0, 0, (1-NU)/2]])
    else:
        #assumes p = 2
        D = (E/(1+NU)/(1-2*NU))*np.array([[1-NU, NU, 0],[NU, 1-NU, 0],[0, 0, (1-2*NU)/2]])
    BTransposeD = np.matmul(BTranspose,D)
    k = t*A*np.matmul(BTransposeD,B)
    return k

def linear_triangle_assemble(K,k,i,j,m):
    '''
    Assembles the element stiffness matrices k into global stiffness matrix K
    '''
    K[2*i,2*i]     = K[2*i,2*i]    + k[0,0]
    K[2*i,2*i+1]   = K[2*i,2*i+1]  + k[0,1]
    K[2*i,2*j]     = K[2*i,2*j]    + k[0,2]
    K[2*i,2*j+1]   = K[2*i,2*j+1]  + k[0,3]
    K[2*i,2*m]     = K[2*i,2*m]    + k[0,4]
    K[2*i,2*m+1]   = K[2*i,2*m+1]  + k[0,5]

    K[2*i+1,2*i]   = K[2*i+1,2*i]  + k[1,0]
    K[2*i+1,2*i+1] = K[2*i+1,2*i+1]+ k[1,1]
    K[2*i+1,2*j]   = K[2*i+1,2*j]  + k[1,2]
    K[2*i+1,2*j+1] = K[2*i+1,2*j+1]+ k[1,3]
    K[2*i+1,2*m]   = K[2*i+1,2*m]  + k[1,4]
    K[2*i+1,2*m+1] = K[2*i+1,2*m+1]+ k[1,5]
    
    K[2*j,2*i]     = K[2*j,2*i]    + k[2,0]
    K[2*j,2*i+1]   = K[2*j,2*i+1]  + k[2,1]
    K[2*j,2*j]     = K[2*j,2*j]    + k[2,2]
    K[2*j,2*j+1]   = K[2*j,2*j+1]  + k[2,3]
    K[2*j,2*m]     = K[2*j,2*m]    + k[2,4]
    K[2*j,2*m+1]   = K[2*j,2*m+1]  + k[2,5]
    
    K[2*j+1,2*i]   = K[2*j+1,2*i]  + k[3,0]
    K[2*j+1,2*i+1] = K[2*j+1,2*i+1]+ k[3,1]
    K[2*j+1,2*j]   = K[2*j+1,2*j]  + k[3,2]
    K[2*j+1,2*j+1] = K[2*j+1,2*j+1]+ k[3,3]
    K[2*j+1,2*m]   = K[2*j+1,2*m]  + k[3,4]
    K[2*j+1,2*m+1] = K[2*j+1,2*m+1]+ k[3,5]

    K[2*m,2*i]     = K[2*m,2*i]    + k[4,0]
    K[2*m,2*i+1]   = K[2*m,2*i+1]  + k[4,1]
    K[2*m,2*j]     = K[2*m,2*j]    + k[4,2]
    K[2*m,2*j+1]   = K[2*m,2*j+1]  + k[4,3]
    K[2*m,2*m]     = K[2*m,2*m]    + k[4,4]
    K[2*m,2*m+1]   = K[2*m,2*m+1]  + k[4,5]
    
    K[2*m+1,2*i]   = K[2*m+1,2*i]  + k[5,0]
    K[2*m+1,2*i+1] = K[2*m+1,2*i+1]+ k[5,1]
    K[2*m+1,2*j]   = K[2*m+1,2*j]  + k[5,2]
    K[2*m+1,2*j+1] = K[2*m+1,2*j+1]+ k[5,3]
    K[2*m+1,2*m]   = K[2*m+1,2*m]  + k[5,4]
    K[2*m+1,2*m+1] = K[2*m+1,2*m+1]+ k[5,5]
    
    return K

def linear_triangle_element_stresses(E,coords,u,NU=poissons,p=1):
    '''     
    Returns element stiffness matrix, Size of matrix is 6x6
    E is youngs modulus
    u is the element displacement vector
    NU is poissons ratio
    p=1 --> plane stress
    p=2 --> plane strain
    '''
    
    xi=coords[0,0]
    yi=coords[0,1]
    xj=coords[1,0]
    yj=coords[1,1]
    xm=coords[2,0]
    ym=coords[2,1]
    
    A = linear_triangle_element_area(xi,yi,xj,yj,xm,ym)
    betai = yj-ym
    betaj = ym-yi
    betam = yi-yj
    gammai = xm-xj
    gammaj = xi-xm
    gammam = xj-xi
    B = (1/(2*A))*np.array([[betai, 0, betaj, 0, betam, 0],[0, gammai, 0, gammaj, 0, gammam],[gammai, betai, gammaj, betaj, gammam, betam]])
    
    if p == 1:
        D = (E/(1-NU*NU))*np.array([[1, NU, 0],[NU, 1, 0],[0, 0, (1-NU)/2]])
    else:
        #assumes p = 2
        D = (E/(1+NU)/(1-2*NU))*np.array([[1-NU, NU, 0],[NU, 1-NU, 0],[0, 0, (1-2*NU)/2]])
    DB = np.matmul(D,B)
    stress = np.matmul(DB,u)
    return stress
        
def linear_triangle_element_pstresses(sigma):
    '''
    returns the element principal stresses and their angle given the stress vector
    '''
    R = (sigma[0] + sigma[1])/2
    Q = ((sigma[0] - sigma[1])/2)**2 + sigma[2]*sigma[2]
    M = 2*sigma[2]/(sigma[0]-sigma[1])
    s1 = R + np.sqrt(Q)
    s2 = R - np.sqrt(Q)
    theta = (np.arctan(M)/2)*180/np.pi
    return (s1,s2,theta)
    
    
    
def linear_triangle_element_displacements(U_global, connectivity, i):
    '''
    U_global is the final displacements of all nodes
    i is the element number
    connectivity tells you of which nodes an element consists
    
    '''
    dof = 2
    displacements = [[U_global[dof*connectivity[i][0]][0]],[U_global[dof*connectivity[i][0]+1][0]],[U_global[dof*connectivity[i][1]][0]],[U_global[dof*connectivity[i][1]+1][0]],[U_global[dof*connectivity[i][2]][0]],[U_global[dof*connectivity[i][2]+1][0]]]
    return displacements
    

def linear_triangle_element_strain(U_global, connectivity, i):
    '''
    U_global is the final displacements of all nodes
    i is the element number
    connectivity tells you of which nodes an element consists
    
    '''
    pass