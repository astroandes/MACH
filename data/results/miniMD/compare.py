import numpy as np, pylab,sys

mcmc = np.loadtxt(sys.argv[1],delimiter=',')
vel = np.loadtxt(sys.argv[2],delimiter=',')
dens = np.loadtxt(sys.argv[3],delimiter=',')

id_mcmc = mcmc[:,0]
id_dens = dens[:,0]

c_mcmc = mcmc[:,4]
c_vel = vel[:,4]
c_dens = dens[:,4]

c_mcmc_dens = np.empty(0)
c_vel_dens = np.empty(0)

for id in id_dens:
    index = np.where(id_mcmc == id)
    c_mcmc_dens = np.append(c_mcmc_dens,c_mcmc[index])
    c_vel_dens = np.append(c_vel_dens,c_vel[index])

pylab.scatter(c_vel,c_mcmc)
pylab.plot(c_vel,c_vel,'-r')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel("$Velocity\ Method$")
pylab.ylabel("$Mass\ Method$")
pylab.savefig('velocity-mass.png',dpi=200)
pylab.close()

pylab.scatter(c_dens,c_mcmc_dens)
pylab.plot(c_dens,c_dens,'-r')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel("$Density\ Method$")
pylab.ylabel("$Mass\ Method$")
pylab.savefig('density-mass.png',dpi=200)
pylab.close()

pylab.scatter(c_dens,c_vel_dens)
pylab.plot(c_dens,c_dens,'-r')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel("$Density\ Method$")
pylab.ylabel("$Velocity\ Method$")
pylab.savefig('density-velocity.png',dpi=200)
pylab.close()
