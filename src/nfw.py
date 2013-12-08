# This code has the density and mass functions of the NFW profile
import numpy as np

# density profile in normal scale
def density(r,rho0,rs):
    return rho0*((r/rs)*(1+(r/rs))**2)**(-1)

# density profile in logarithmic scale
def logdensity(logr,rho0,rs):
    return np.log10(rho0)-np.log10(rs)+logr+2*np.log10(1+10.0**(logr-np.log10(rs)))

# density profile in logarithmic scale with logarithmic parameters
def loglogdensity(logr,logrho0,logrs):
    return logrho0-logrs+logr+2*np.log10(1+10.0**(logr-logrs))

# mass profile in normal scale
def mass(r,rho0,rs):
    return 4*np.pi*(rs**3)*(np.log(1+(r/rs))+(1+(r/rs))**(-1))

# mass profile in logarithmic scale
def logmass(logr,rho0,rs):
    return np.log10(4*np.pi)+3*np.log10(rs)+np.log10(np.log(1+10.0**(logr-np.log10(rs)))+(1+10.0**(logr-np.log10(rs)))**(-1))

# mass profile in logarithmic scale with logarithmic parameters
def logmass(logr,rho0,rs):
    return np.log10(4*np.pi)+3*logrs+np.log10(np.log(1+10.0**(logr-logrs))+(1+10.0**(logr-logrs))**(-1))
