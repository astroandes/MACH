import numpy as np, emcee, sys, os, nfw
from scipy.optimize import curve_fit
from emcee.utils import sample_ball

# Returns the chi squared estimate given a prediction and the real data
# Requires:
#     obs: observed data
#     mod: expected data from the model
# Returns:
#     The chi squared estimate

def chi2(obs,mod):
    return np.sum(-0.5*(obs-mod)**2)

# Does a maximum-likelihood estimate using the Metropolis-Hastings algorithm for a model with two parameters.
# Requires:
#     data: free variable data
#     obs: observed data
#     mod: expected data from the model
#     n_iterations: number of steps for the random walk
# Returns:
#     best parameters, the random walk for each parameter and the chi squared for each step

def metropolis(data,obs,mod,n_iterations):

    sys.stdout.write('\rRunning Metropolis-Hastings Algorithm... ')
    sys.stdout.flush()

    guess = curve_fit(mod,data,obs,maxfev=n_iterations)[0]

    a_walk = np.empty((0))
    b_walk = np.empty((0))
    chisq = np.empty((0))

    a_walk = np.append(a_walk,guess[0])
    b_walk = np.append(b_walk,guess[1])
    chisq = np.append(chisq,chi2(obs,mod(data,guess[0],guess[1])))
        
    for i in range(n_iterations):
        a_prime = np.random.normal(a_walk[i],0.001)
        b_prime = np.random.normal(b_walk[i],0.001) 
        
        chi2_init = chi2(obs,mod(data,a_walk[i],b_walk[i]))
        chi2_prime = chi2(obs,mod(data,a_prime,b_prime))

        alpha = np.exp(chi2_prime-chi2_init)
        ratio = chi2_init/chi2_prime
        if (ratio >= 1.0):
            a_walk = np.append(a_walk,a_prime)
            b_walk = np.append(b_walk,b_prime)
            chisq = np.append(chisq,chi2_prime)

        else:
            beta = np.random.random()
            if (alpha >= beta):
                a_walk = np.append(a_walk,a_prime)
                b_walk = np.append(b_walk,b_prime)
                chisq = np.append(chisq,chi2_prime)
            else:
                a_walk = np.append(a_walk,a_walk[i])
                b_walk = np.append(b_walk,b_walk[i])
                chisq = np.append(chisq,chi2_init)

    max_pos = np.argmax(chisq)
    best_a = a_walk[max_pos]
    best_b = b_walk[max_pos]

    sys.stdout.write('Done\n')

    return best_a,best_b,a_walk,b_walk,chisq

# Compiles the Metropolis-Hastings algorithm code in C for the loglogmass function in the nfw module.
def compile_c_metropolis():

    sys.stdout.write('\rCompiling Metropolis-Hastings Algorithm in C... ')
    sys.stdout.flush()
    os.system('cc metropolis.c -lm -o metropolis.out')
    sys.stdout.write('Done\n')

# Does a maximum-likelihood estimate using the Metropolis-Hastings algorithm in C for the loglogmass function in the nfw module.
# Requires:
#     data: free variable data
#     obs: observed data
#     mod: expected data from the model
#     n_iterations: number of steps for the random walk
#     Having excetuted the compile_c_metropolis function
# Returns:
#     best parameters, the random walk for each parameter and the chi squared for each step
def c_metropolis(data,obs,n_iterations):

    sys.stdout.write('\rRunning Metropolis-Hastings Algorithm in C... ')
    sys.stdout.flush()

    guess = curve_fit(nfw.loglogmass,data,obs,maxfev=n_iterations)[0]
    
    open('profile.dat', "w").write('\n'.join('%lf %lf' % (data[i],obs[i]) for i in range(len(data))))

    os.system('./metropolis.out profile.dat '+str(guess[0])+' '+str(guess[1])+' '+str(n_iterations))
    a_walk = np.loadtxt('a_walk.dat')
    b_walk = np.loadtxt('b_walk.dat')
    chisq = np.loadtxt('chi2.dat')

    first_max = np.argmax(chisq)
    first_min = np.argmin(chisq)
    first_a = a_walk[first_max]
    first_b = b_walk[first_max]
    last_a = a_walk[first_min]
    last_b = b_walk[first_min]
    p0 = np.array([first_a,first_b])
    p1 = np.array([last_a,last_b])
    distance = np.sqrt(np.sum((p0-p1)**2))

##    positions = sample_ball(np.array([first_a,first_b]), distance*np.ones(2), 5)
    
##    for i in range(5):
##        position = positions[i,:]
        
##        os.system('./metropolis.out profile.dat '+str(position[0])+' '+str(position[1])+' '+str(n_iterations))

##        a_walk = np.concatenate((a_walk,np.loadtxt('a_walk.dat')))
##        b_walk = np.concatenate((b_walk,np.loadtxt('b_walk.dat')))
##        chisq = np.concatenate((chisq,np.loadtxt('chi2.dat')))
 
    new_max = np.argmax(chisq)
    best_a = a_walk[new_max]
    best_b = b_walk[new_max]

    os.system('rm *.dat')

    sys.stdout.write('Done\n')

    return best_a,best_b,a_walk,b_walk,chisq

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

    dim = 2
    walkers = 4
    guess = curve_fit(mod,data,obs,maxfev=n_iterations)[0]
    print guess
    p0 = sample_ball(guess, 100*np.ones(dim), walkers)
    
    
    def loglike(p):        
	return -np.sum(((mod(data,p[0],p[1]))-obs)**2)/dim

    sampler = emcee.EnsembleSampler(walkers, dim, loglike)
    sampler.run_mcmc(p0, n_iterations)
    
    a_walk = sampler.flatchain[:,0]
    b_walk = sampler.flatchain[:,1]
    chisq = sampler.flatlnprobability
    for l in chisq:
        print l

    max_pos = np.argmax(chisq)
    best_a = a_walk[max_pos]
    best_b = b_walk[max_pos]

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
