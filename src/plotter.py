import numpy as np, pylab,sys, nfw
from matplotlib.mlab import griddata

angle = np.arange(0,2*np.pi,0.01)

# Plots a halo
# Requires:
#     x,y,z: arrays with the positions for each particle
#     x_center,y_center,z_center: position of the halo center
# Returns:
#     Three files with graphics in the planes x-y, y-z and z-x of the halo

def halo(x,y,z,x_center,y_center,z_center,r_bdmv,r_bdmw):
    
    sys.stdout.write('\rPlotting halo... ')
    sys.stdout.flush()
    pylab.plot(x , y, 'k.')
    pylab.plot(x_center , y_center, 'ko')
    pylab.plot(r_bdmv*np.cos(angle)+x_center,r_bdmv*np.sin(angle)+y_center,'--b',label='BDMV')
    pylab.plot(r_bdmw*np.cos(angle)+x_center,r_bdmw*np.sin(angle)+y_center,'--r',label='BDMW')
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xlabel('x (Mpc/H)')
    pylab.ylabel('y (Mpc/H)')
    pylab.title('Halo x-y')
    pylab.savefig('halo_xy.png',format='png',dpi=300)
    pylab.close()

    pylab.plot(y , z, 'k.')
    pylab.plot(y_center , z_center, 'ko')
    pylab.plot(r_bdmv*np.cos(angle)+y_center,r_bdmv*np.sin(angle)+z_center,'--b',label='BDMV')
    pylab.plot(r_bdmw*np.cos(angle)+y_center,r_bdmw*np.sin(angle)+z_center,'--r',label='BDMW')
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xlabel('y (Mpc/H)')
    pylab.ylabel('z (Mpc/H)')
    pylab.title('Halo y-z')
    pylab.savefig('halo_yz.png',format='png',dpi=300)
    pylab.close()

    pylab.plot(z , x, 'k.')
    pylab.plot(z_center , x_center, 'ko')
    pylab.plot(r_bdmv*np.cos(angle)+z_center,r_bdmv*np.sin(angle)+x_center,'--b',label='BDMV')
    pylab.plot(r_bdmw*np.cos(angle)+z_center,r_bdmw*np.sin(angle)+x_center,'--r',label='BDMW')
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xlabel('z (Mpc/H)')
    pylab.ylabel('x (Mpc/H)')
    pylab.title('Halo z-x')
    pylab.savefig('halo_zx.png',format='png',dpi=300)
    pylab.close()
    sys.stdout.write('Done\n')

# Plots the NFW mass profile 
# Requires:
#     radius: array with the radial distances
#     mass: array with the mass values
#     parameters: an array with the parameters of the NFW profile (mean density and scale radius)
# Returns:
#     A file with agraphic of the mass profile

def mass(radius, mass, parameters):
    
    pylab.plot(radius , mass,'.r',label="Real Mass")
    pylab.plot(radius , nfw.mass(radius,parameters[0],parameters[1]),'k',label="NFW profile")
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xlabel('Radius (Mpc/H)')
    pylab.ylabel('Mass (10^11 Solar Masses)')
    pylab.title('Mass')
    pylab.savefig('mass.png',format='png',dpi=300)
    pylab.close()

# Plots the NFW mass profile in logarithmic scale 
# Requires:
#     logR: array with the logarithm of the radial distances
#     logM: array with the logarithm of the mass values
#     parameters: an array with the parameters of the NFW profile (mean density and scale radius)
# Returns:
#     A file with a graphic of the mass profile in logarithmic scale
    
def logmass(radius, mass, parameters,lims_a,lims_b):
    
    pylab.plot(radius, mass,'.r',label="Real Mass")
    pylab.plot(radius, nfw.mass(radius,parameters[0],parameters[1]),'k',label="NFW profile")

    pylab.plot(radius , nfw.mass(radius,lims_a[0],lims_b[0]),'--r',label="++")
    pylab.plot(radius , nfw.mass(radius,lims_a[1],lims_b[1]),'--b',label="--")

    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Radius (Mpc/H)')
    pylab.ylabel('Mass (10^11 Solar Masses)')
    pylab.title('Mass (log-log)')
    pylab.savefig('log_mass.png',format='png',dpi=300)
    pylab.close()

def mass_norm(radius,mass,c,c_max,c_min,name):
    
    pylab.plot(radius, mass,'.r',label="Real Norm Mass")
    pylab.plot(radius, nfw.mass_norm(radius,c),'k',label="NFW profile")
    pylab.plot(radius , nfw.mass_norm(radius,c_max),'--r',label="Max param")
    pylab.plot(radius , nfw.mass_norm(radius,c_min),'--b',label="Min param")
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xlabel('Radius (Normalized)')
    pylab.ylabel('Mass (Normalized)')
    pylab.savefig('mass_norm_'+name+'.png',format='png',dpi=300)
    pylab.close()

# Plots the NFW density profile in logarithmic scale 
# Requires:
#     r_density: array with the radial distances
#     density: array with the density values
#     logR: array with the logarithm of the radial distances
#     parameters: an array with the parameters of the NFW profile (mean density and scale radius)
# Returns:
#     A file with a graphic of the density profile in logarithmic scale
    
def logdensity(r_density, density, parameters):
    pylab.plot(r_density , density,'.r',label="Real Density")
    pylab.plot(r_density , nfw.density(r_density,parameters[0],parameters[1]),'k',label="NFW profile")
    pylab.legend(loc=4, borderaxespad=0.5)
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Radius (Mpc/H)')
    pylab.ylabel('Density (10^11 Solar Masses/(Mpc/H)^3)')
    pylab.title('Density (log-log)')
    pylab.savefig('log_density.png',format='png',dpi=300)
    pylab.close()

# Makes a contour plot of the chi squared 
# Requires:
#     a_walk: array with the random walk of the first parameter
#     b_walk: array with the random walk of the second parameter
#     chi2: array with the chi squared estimate for each step in the random walk
# Returns:
#     A file with a contour plot of the chi squared in function of both parameters
    
def rainbow_chi2(a_walk,b_walk,chi2):
    
    N = 2000j
    extent = (a_walk[np.argmin(a_walk)],a_walk[np.argmax(a_walk)],b_walk[np.argmin(b_walk)],b_walk[np.argmax(b_walk)])

    my_xs,my_ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

    my_resampled = griddata(a_walk, b_walk, -chi2, my_xs, my_ys)
    pylab.imshow(my_resampled.T,extent=extent,origin='lower',interpolation='bicubic',cmap='spectral',aspect='auto')
    pylab.title(r'$\ln(\cal{L})$')
    pylab.xlabel(r'$\log(\rho_{0})$')
    pylab.ylabel(r'$\log(R_s)$')
    pylab.colorbar()
    pylab.savefig('chi2.png',format='png',dpi=300)
    pylab.close()

# Makes a contour plot of the likelihood 
# Requires:
#     a_walk: array with the random walk of the first parameter
#     b_walk: array with the random walk of the second parameter
#     chi2: array with the chi squared estimate for each step in the random walk
# Returns:
#     A file with a contour plot of the likelihood in function of both parameters
    
def rainbow_likelihood(a_walk,b_walk,chi2):
    
    N = 2000j
    extent = (a_walk[np.argmin(a_walk)],a_walk[np.argmax(a_walk)],b_walk[np.argmin(b_walk)],b_walk[np.argmax(b_walk)])

    my_xs,my_ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

    my_resampled = griddata(a_walk, b_walk, np.exp(chi2), my_xs, my_ys)
    pylab.imshow(my_resampled.T,extent=extent,origin='lower',interpolation='bicubic',cmap='spectral',aspect='auto')
    pylab.plot(a_walk[0],b_walk[0],'or')
    pylab.plot(a_walk[-1],b_walk[-1],'ob')
    pylab.title(r'$\cal{L}$')
    pylab.xlabel(r'$\ln(R_s)$')
    pylab.ylabel(r'$\ln(\rho_{0})$')
    pylab.colorbar()
    pylab.savefig('likelihood.png',format='png',dpi=300)
    pylab.close()

# Makes a contour plot of the likelihood in logarithmic scale
# Requires:
#     a_walk: array with the random walk of the first parameter
#     b_walk: array with the random walk of the second parameter
#     chi2: array with the chi squared estimate for each step in the random walk
# Returns:
#     A file with a contour plot of the likelihood in logarithmic scale in function of both parameters
    
def rainbow_loglikelihood(a_walk,b_walk,chi2):
    
    N = 2000j
    extent = (a_walk[np.argmin(a_walk)],a_walk[np.argmax(a_walk)],b_walk[np.argmin(b_walk)],b_walk[np.argmax(b_walk)])

    my_xs,my_ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

    my_resampled = griddata(a_walk, b_walk, chi2, my_xs, my_ys)
    pylab.imshow(my_resampled.T,extent=extent,origin='lower',interpolation='bicubic',cmap='spectral',aspect='auto')
    pylab.title(r'$\ln(\cal{L})$')
    pylab.xlabel(r'$\rho_{0}$')
    pylab.ylabel(r'$\R_s$')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.colorbar()
    pylab.savefig('loglikelihood.png',format='png',dpi=300)
    pylab.close()

def random_walk(a_walk,b_walk,n_iterations):
    
    pylab.plot(a_walk,b_walk,'-k')
    pylab.plot(a_walk,b_walk,',r')
    pylab.plot(a_walk[0],b_walk[0],'ob')
    for i in range(1,5):
        pylab.plot(a_walk[n_iterations*i-1],b_walk[n_iterations*i-1],'oy')
        pylab.plot(a_walk[n_iterations*i],b_walk[n_iterations*i],'oc')
    pylab.xlabel(r'$\ln(R_s)$')
    pylab.ylabel(r'$\ln(\rho_{0})$')    
    pylab.savefig('randomwalk.png',format='png',dpi=300)
    pylab.close()
