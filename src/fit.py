import numpy as np, sys, os, nfw, emcee
from scipy.optimize import curve_fit
from emcee.utils import sample_ball

# Does a maximum-likelihood estimate using the standard emcee sampler for a model with two parameters.
# Requires:
#     data: free variable data
#     obs: observed data
#     mod: expected data from the model
#     n_iterations: number of steps for the random walk
# Returns:
#     best parameters, the random walk for each parameter and the chi squared for each step

def emcee_sampler(data,obs,mod,n_iterations):

    sys.stdout.write('\rRunning Emcee Sampling Algorithm... ')
    sys.stdout.flush()

    dim = 1
    walkers = 2

    guess = curve_fit(mod,data,obs,maxfev=n_iterations)[0]
    p0 = sample_ball(guess, np.ones(dim), walkers)
    
    def loglike(p):        
	return -np.sum(((mod(data,p[0],p[1]))-obs)**2)/dim

    sampler = emcee.EnsembleSampler(walkers, dim, loglike)
    sampler.run_mcmc(p0, n_iterations)
    
    a_walk = sampler.flatchain
    chisq = sampler.flatlnprobability

    max_pos = np.argmax(chisq)
    best_a = a_walk[max_pos]

    sys.stdout.write('Done\n')
    
    return best_a,a_walk,-chisq

# Gets the acceptance interval for a value of a parameter in a random walk
# Requires:
#     walk: array with the random walk
#     parameter: value of the parameter
#     opt: need to put 'log' if the random walk is in logarithmic scale
# Returns:
#     minimum and maximum boundary for the interval

def error_bars(walk,parameter,opt):

    n_iterations = len(walk)
    n_bins = int(n_iterations/100)

    if (opt=='log'):
        freq,bins = np.histogram(np.exp(walk), bins=n_bins)
        index = np.argmin(np.abs(bins-np.exp(parameter)))

    else:
        freq,bins = np.histogram(walk, bins=n_bins)
        index = np.argmin(np.abs(bins-parameter))

    max_c = [index+int(0.341*sum(freq)/n_bins),len(bins)-1]
    min_c = [index-int(0.341*sum(freq)/n_bins),0]

    M = max_c[np.argmin(max_c)]
    m = min_c[np.argmax(min_c)]
    
    if (opt=='log'):

        max = np.log(bins[M])
        min = np.log(bins[m])

    else:
        max = bins[M]
        min = bins[m]
    
    return max, min
