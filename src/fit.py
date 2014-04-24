import numpy as np, sys, os, nfw, emcee
from scipy.optimize import curve_fit
from emcee.utils import sample_ball
# Returns the chi squared estimate given a prediction and the real data
# Requires:
#     obs: observed data
#     mod: expected data from the model
# Returns:
#     The chi squared estimate

def chi2(obs,mod):
    return np.sum((obs-mod)**2)

# Does a maximum-likelihood estimate using the Metropolis-Hastings algorithm for a model with one parameters.
# Requires:
#     data: free variable data
#     obs: observed data
#     mod: expected data from the model
#     n_iterations: number of steps for the random walk
# Returns:
#     best parameter, the random walk for the parameter and the chi squared for each step

def metropolis_one(data,obs,mod,n_iterations):

    sys.stdout.write('\rRunning Metropolis-Hastings Algorithm... ')
    sys.stdout.flush()
    
    guess = [np.log(10)]
    
    a_walk = np.empty((n_iterations+1))
    chisq = np.empty((n_iterations+1))

    a_walk[0] = guess[0]
    chisq[0] = chi2(obs,mod(data,guess[0]))

    step_a = np.abs(np.log(1.01))    

    for i in range(n_iterations):
        
        a_prime = np.random.normal(a_walk[i],step_a)

        while a_prime < 0 :            
            a_prime = np.random.normal(a_walk[i],step_a)
        
        chi2_init = chi2(obs,mod(data,a_walk[i]))
        chi2_prime = chi2(obs,mod(data,a_prime))
        ratio = (chi2_init - chi2_prime)

        if (ratio > 0.0):
            a_walk[i+1] = a_prime
            chisq[i+1] = chi2_prime
        else:
            beta = np.random.random()
            if (ratio > np.log(beta)):
                a_walk[i+1] = a_prime
                chisq[i+1] = chi2_prime
            else:
                a_walk[i+1] = a_walk[i]
                chisq[i+1] = chisq[i]

    n = np.argmin(chisq)

    best_a = a_walk[n]

    sys.stdout.write('Done\n')

    return best_a,a_walk,chisq

# Does a maximum-likelihood estimate using the Metropolis-Hastings algorithm for a model with two parameters.
# Requires:
#     data: free variable data
#     obs: observed data
#     mod: expected data from the model
#     n_iterations: number of steps for the random walk
# Returns:
#     best parameters, the random walk for each parameter and the chi squared for each step

def metropolis(data,obs,mod,n_iterations, maxi,mini):

    sys.stdout.write('\rRunning Metropolis-Hastings Algorithm... ')
    sys.stdout.flush()
    
    guess = curve_fit(mod,data,obs,maxfev=n_iterations)[0]

    if guess[1] > maxi:
        guess[1] = maxi
    if guess[1] < mini:
        guess[1] = mini 

    a_walk = np.empty((n_iterations+1))
    b_walk = np.empty((n_iterations+1))
    chisq = np.empty((n_iterations+1))

    a_walk[0] = guess[0]
    b_walk[0] = guess[1]
    chisq[0] = chi2(obs,mod(data,guess[0],guess[1]))
    
    for i in range(n_iterations):

	step_a = 1.0
        step_b = 0.1
        a_prime = np.random.normal(a_walk[i],step_a)
        b_prime = np.random.normal(b_walk[i],step_b) 
        
        while(b_prime > maxi or b_prime < mini):
            b_prime = np.random.normal(b_walk[i],step_b)
        chi2_init = chi2(obs,mod(data,a_walk[i],b_walk[i]))
        chi2_prime = chi2(obs,mod(data,a_prime,b_prime))

        ratio = (chi2_init - chi2_prime)

        if (ratio > 0.0):
            a_walk[i+1] = a_prime
            b_walk[i+1] = b_prime
            chisq[i+1] = chi2_prime

        else:
            beta = np.random.random()
            if (ratio > np.log(beta)):
                a_walk[i+1] = a_prime
                b_walk[i+1] = b_prime
                chisq[i+1] = chi2_prime
            else:
                a_walk[i+1] = a_walk[i]
                b_walk[i+1] = b_walk[i]
                chisq[i+1] = chisq[i]

    n = np.argmin(chisq)

    best_a = a_walk[n]
    best_b = b_walk[n]

    sys.stdout.write('Done\n')

    return best_a,best_b,a_walk,b_walk,chisq

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
