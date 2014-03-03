import numpy as np, pylab, sys

orig = np.loadtxt(sys.argv[1],delimiter=',')
mark = np.loadtxt(sys.argv[2]) 

mcmc = mark[mark[:,0].argsort()]

orig_rho0 = orig[:,1]
mcmc_rho0 = mcmc[:,4]
mcmc_rho0_min = mcmc[:,5]
mcmc_rho0_max = mcmc[:,6]
rho0_sigma = np.array([mcmc_rho0_max-mcmc_rho0,mcmc_rho0-mcmc_rho0_min])

orig_rs = orig[:,2]
mcmc_rs = mcmc[:,7]
mcmc_rs_min = mcmc[:,8]
mcmc_rs_max = mcmc[:,9]
rs_sigma = np.array([mcmc_rs_max-mcmc_rs,mcmc_rs-mcmc_rs_min])

pylab.errorbar(orig_rho0,mcmc_rho0,yerr=rho0_sigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.plot(orig_rs,orig_rs,'--k')
pylab.title('Density (Test)')
pylab.xlabel('Original')
pylab.ylabel('Montecarlo-NFW')
pylab.savefig('test_rho0.png',dpi=150)
pylab.close()

pylab.errorbar(orig_rs,mcmc_rs,yerr=rs_sigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.plot(orig_rs,orig_rs,'--k')
pylab.title('Scale Radius (Test)')
pylab.xlabel('Original')
pylab.ylabel('Montecarlo-NFW')
pylab.savefig('test_rs.png',dpi=150)
pylab.close()
