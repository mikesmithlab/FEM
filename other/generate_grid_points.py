import numpy as np
import matplotlib.pyplot as plt

Psi_t_arr = []

for rand_size in np.arange(0,0.6,0.005):
    
    x = []
    y = []
    a = [1,0]
    b = [np.cos(np.pi/3),np.sin(np.pi/3)]
    for n in np.arange(1,30,1):
        for m in np.arange(1,30,1):
    
            x.append((n*a[0] + m*b[0]) + np.random.uniform(-rand_size,rand_size))
            y.append((n*a[1] + m*b[1]) + np.random.uniform(-rand_size,rand_size) - 0.7)
            
    G = ((np.pi)*2)/np.array(b)
    Psi_t = []
    
    for o in np.arange(0,len(x),1):
        
        r = np.array([x[o],y[o]])
        
        Psi_t.append(np.exp(1j*np.dot(G,r)))
        
    Psi_t_arr.append(np.mean(Psi_t))
    
plt.plot(Psi_t_arr)
        