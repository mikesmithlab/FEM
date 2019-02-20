import numpy as np
import scipy
import scipy.linalg
import displayUtils as disp

'''
k is a local spring constant stiffness
K is the global stiffness matrix
u is the nodal displacement vectors
i and j are the indices of nodes

'''
def known_node_val(nodepts,connectivity,value,known_index=None,shape='ellipse'):
    indices = np.where(disp.select_pts(nodepts, connectivity,shape=shape))
    print(np.shape(indices))
    if known_index is None:
        known_index = indices[0]
    else:
        known_index = np.append(known_index,indices[0])
    known = value*np.ones(np.shape(indices[0]),dtype=float)
    return known, known_index


def initialiseU(K,knownIndexX, knownValuesX,knownIndexY,knownValuesY):
    '''
    K is the global stiffness matrix
    knownIndex = numpy array of indices where displacement is known
    knownValue = numpy array of values corresponding
    '''
    
    U = np.ones((np.shape(K)[0],1),dtype=float)*np.nan
    for index,value in zip(knownIndexX,knownValuesX):
        U[int(2*index),0] = value
    for index,value in zip(knownIndexY,knownValuesY):
        U[int(2*index + 1),0]=value
    return U
    	
def initialiseF(K,knownIndexX, knownValuesX,knownIndexY,knownValuesY):
    '''
    K is the global stiffness matrix
    knownIndex = numpy array of indices where force is known
    knownValue = numpy array of values corresponding
    '''
    lenXvals = np.shape(K)[0]
    F = np.ones((np.shape(K)[0],1),dtype=float)*np.nan
    for index,value in zip(knownIndexX,knownValuesX):
        F[int(2*index),0] = value
    for index,value in zip(knownIndexY,knownValuesY):
        F[int(2*index + 1),0]=value
    return F
	
	
def partitionMatrix(K,U,F):
    ''' 
    K is the global stiffness matrix
    U is the nodal displacement's vector some known and some unknown
    F is the nodal force's vector some known and some unknown.
    
    
    '''
    nan_disp_index = np.where(np.isnan(U[:,0]))[0]
    non_zero_ornan_index = np.nonzero(U[:,0])[0]
    non_zero_disp_index = np.setdiff1d(non_zero_ornan_index,nan_disp_index)
    
    #If a value is unknown it will have the value np.nan. These are the rows and columns we want in our partitioned matrix
    partitionedK = K[np.reshape(nan_disp_index, (-1, 1)),[nan_disp_index]]
	
	#For each nonzero displacement want to multiply rows of stiffness matrix which are unknowns (ie nans in U) with column = row of non-zero displacement in U and subtract from f vector
    f=np.copy(F)
    f[nan_disp_index] = f[nan_disp_index] - np.matmul(K[np.reshape(nan_disp_index, (-1, 1)),[non_zero_disp_index]] , U[non_zero_disp_index])
    u = gaussElim(partitionedK, f[nan_disp_index])	
        
    return (u,f[nan_disp_index])



def gaussElim(a,b):
    '''
    Solves ax=b where a is square matrix mxm and b is vector of length m.
    returns solution x
    
    '''
    x = scipy.linalg.solve(a,b)
    return x
    
	
def globalDisplacement(u,U):
	'''
	When everything has been calculated we combine the originally known displacements
	with the ones we have just calculated to give the complete list of nodal displacements
	'''
	nans = np.where(np.isnan(U[:,0]))
	U[nans[0]] = u
	return U

def globalForce(K_global,U_global):
    '''
    We calculate the complete nodal forces vector using the complete
    list of nodal displacements and the global stiffness matrix.
    '''
    
    F_global = np.dot(K_global,U_global)
    return F_global
	
	
def checkBoundaryConditions(numbNodes,knownDispIndex,knownDisplacements,knownForceIndex,knownForces):
    #Perform some basic checks
    error = False
    #Are node indices same length as value arrays
    if np.shape(knownDispIndex)[0] != np.shape(knownDisplacements)[0]:
        return (True, 'Displacement indices different length values given')
    
    if np.shape(knownForceIndex)[0] != np.shape(knownForces)[0]:
        return (True, 'Force indices different length values given')
    
    nodenumbers = np.arange(numbNodes)
    unspecifiedBCs = np.setdiff1d(nodenumbers,np.append(knownDispIndex,knownForceIndex))
    if np.shape(unspecifiedBCs)[0] != 0:
        return (True, 'One node has not been assigned an initial force or displacement')
    
    return (False,'')
    
if __name__ ==  '__main__':
    pass
    
    
    