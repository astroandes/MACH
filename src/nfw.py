import numpy as np

# Profiles in normal scale
# Parameters:
#     r: numpy array with the radial distances from the center of the halo
#     rho0: mean density
#     rs: scale radius
# Returns:
#     a numpy array with the mass/density profile

def mass(r,rho0,rs):
    return 4*np.pi*(rs**3)*(np.log(1+(r/rs))+(1+(r/rs))**(-1))

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
    return np.log10(4*np.pi)+3*np.log10(rs)+np.log10(np.log(1+10.0**(logr-np.log10(rs)))+(1+10.0**(logr-np.log10(rs)))**(-1))

def logdensity(logr,rho0,rs):
    return np.log10(rho0)-np.log10(rs)+logr+2*np.log10(1+10.0**(logr-np.log10(rs)))

# Profiles in logarithmic scale with logarithmic parameters
# Parameters:
#     logr: numpy array with the logarithm of the radial distances from the center of the halo
#     logrho0: logarithm of the mean density
#     logrs: logarithm of the scale radius
# Returns:
#     a numpy array with the logarithm of the mass/density profile

def loglogmass(logr,logrho0,logrs):
    return np.log10(4*np.pi)+3*logrs+np.log10(np.log(1+10.0**(logr-logrs))+(1+10.0**(logr-logrs))**(-1))

def loglogdensity(logr,logrho0,logrs):
    return logrho0-logrs+logr+2*np.log10(1+10.0**(logr-logrs))
