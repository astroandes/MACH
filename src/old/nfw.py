import numpy as np

# Profiles in normal scale
# Parameters:
#     r: numpy array with the radial distances from the center of the halo
#     rho0: mean density
#     rs: scale radius
# Returns:
#     a numpy array with the mass/density profile

def mass(r,a,b):
    return 4*np.pi*a*(b**3)*(np.log(1+r/b)-(1+(b/r))**(-1))

def density(r,a,b):
    return a*((r/b)*(1+(r/b))**2)**(-1)

def mass_norm(r_norm,c):
    return (np.log(1.0+r_norm*c)-r_norm*c/(1.0+r_norm*c))/(np.log(1.0+c)-c/(1.0+c))

def density_norm(r_norm,c):
    return ((1+c)**2/(r_norm*(1+r_norm*c)**2))

def velocity_norm(r_norm,c):
    return np.sqrt(mass_norm(r_norm,c)/r_norm)

# Profiles in logarithmic scale
# Parameters:
#     logr: numpy array with the logarithm of the radial distances from the center of the halo
#     a: mean density
#     b: scale radius
# Returns:
#     a numpy array with the logarithm of the mass/density profile

def logmass(logr,a,b):
    return np.log(4*np.pi)+ np.log(a)+3*np.log(b)+np.log(np.log(1+np.exp(logr)/b)-(1+b/np.exp(logr))**(-1))

def logdensity(logr,a,b):
    return np.log(a)+np.log(b)-logr-2*np.log(1+np.exp(logr)/b)

def logdensity_norm(logr_norm,c):
    return (2.0*(np.log(1+c)-np.log(1+c*np.exp(logr_norm)))-np.log(logr_norm))

# Profiles in logarithmic scale with logarithmic parameters
# Parameters:
#     logr: numpy array with the logarithm of the radial distances from the center of the halo
#     loga: logarithm of the mean density
#     logb: logarithm of the scale radius
# Returns:
#     a numpy array with the logarithm of the mass/density profile

def loglogmass(logr,loga,logb):
    return np.log(4*np.pi)+ loga+3*logb+np.log(np.log(1+np.exp(logr-logb))-(1+np.exp(logb-logr))**(-1))

def loglogdensity(logr,loga,logb):
    return loga+logb-logr-2*np.log(1+np.exp(logr-logb))

def loglogmass_norm(logr_norm,logc):
    return np.log((np.log(1.0+np.exp(logr_norm+logc))-np.exp(logr_norm+logc)/(1.0+np.exp(logr_norm+logc)))/(np.log(1.0+np.exp(logc))-np.exp(logc)/(1.0+np.exp(logc))))

def loglogdensity_norm(logr_norm,logc):
    return 2.0*(np.log(1.0+np.exp(logc))-np.log(1.0+np.exp(logr_norm+logc)))-logr_norm