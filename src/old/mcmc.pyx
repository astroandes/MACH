STUFF = "Hi"

import numpy as np
cimport numpy as np

def chi2(np.ndarray p,model,np.ndarray x_obs,np.ndarray y_obs,np.ndarray sigma):
    cdef float chi_sq = np.sum(((model(x_obs,p) - y_obs)/sigma)**2)
    return chi_sq

def walker(dist,np.ndarray state_0,np.ndarray step_size,n_iterations,reestrictions = None):

    cdef int dims = len(state_0)
    cdef np.ndarray walk = np.empty((dims,n_iterations))
    cdef np.ndarray prob = np.empty(n_iterations)
    cdef np.ndarray t_init
    cdef np.ndarray t_prime
    cdef float P_init
    cdef float P_prime
    cdef float a
    cdef float b
    walk[:,0] = state_0
    prob[0] = dist(state_0)

    if reestrictions == None:
        reestrictions = [lambda x: True for i in range(dims)]
        
    for i in range(1,n_iterations):
        
        t_init = walk[:,i-1]
        t_prime = np.empty(dims)
        
        for n in range(dims):
            t_prime[n] = np.random.normal(t_init[n],step_size[n])
            while reestrictions[n](t_prime[n]) == False:
                t_prime[n] = np.random.normal(t_init[n],step_size[n])

        P_init = dist(t_init)
        P_prime = dist(t_prime)

        a = -(P_prime-P_init)

        if a > 0:
            walk[:,i] = t_prime
            prob[i] = P_prime
        else:
            b = np.log(np.random.random())
            if a > b:
                walk[:,i] = t_prime
                prob[i] = P_prime
            else:
                walk[:,i] = t_init
                prob[i] = P_init
                
    return walk,prob