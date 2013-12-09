import numpy as np

# Profiles in normal scale
# Parameters:
#     r: numpy array with the radial distances from the center of the halo
#     rho0: mean density
#     rs: scale radius
# Returns:
#     a numpy array with the mass/density profile

def mass(r,rho0,rs):
    return 4*np.pi*(rs**3)*rho0*(np.log(1+(r/rs))-(1+(rs/r))**(-1))

def density(r,rho0,rs):
    return rho0*((r/rs)*(1+(r/rs))**2)**(-1)

# Profiles in logarithmic scale
# Parameters:
#     logr: numpy array with the logarithm of the radial distances from the center of the halo
#     rho0: mean density
#     rs: scale radius
# Returns:
#     a numpy array with the logarithm of the mass/density profile

def logmass(logr,rho0,rs):
    return np.log(4*np.pi)+3*np.log(rho0)+3*np.log(rs)+np.log(np.log(1+np.exp(logr)/rs)-(1+rs*np.exp(-logr))**(-1))

def logdensity(logr,rho0,rs):
    return np.log(rho0)-np.log(rs)+logr+2*np.log(1+np.exp(logr)/rs)

# Profiles in logarithmic scale with logarithmic parameters
# Parameters:
#     logr: numpy array with the logarithm of the radial distances from the center of the halo
#     logrho0: logarithm of the mean density
#     logrs: logarithm of the scale radius
# Returns:
#     a numpy array with the logarithm of the mass/density profile

def loglogmass(logr,logrho0,logrs):
    return np.log(4*np.pi)+3*logrs+logrho0+np.log(np.log(1+np.exp(logr-logrs))-(1+np.exp(logrs-logr))**(-1))

def loglogdensity(logr,logrho0,logrs):
    return logrho0-logrs+logr+2*np.log(1+np.exp(logr-logrs))
