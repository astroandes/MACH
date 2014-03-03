import numpy as np, pylab, sys

orig = np.loadtxt(sys.argv[1],delimiter=',')
mark = np.loadtxt(sys.argv[2]) 

mcmc = mark[mark[:,0].argsort()]

orig_rs = orig[:,2]
mcmc_rs = mcmc[:,7]
mcmc_rs_min = mcmc[:,8]
mcmc_rs_max = mcmc[:,9]
rs_sigma = np.array([mcmc_rs_max-mcmc_rs,mcmc_rs-mcmc_rs_min])

pylab.errorbar(orig_rs,mcmc_rs,yerr=rs_sigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.plot(orig_rs,orig_rs,'--k')
pylab.title('Scale Radius (Test)')
pylab.xlabel('Original')
pylab.ylabel('Montecarlo-NFW')
pylab.savefig('test.png',dpi=150)
pylab.close()
