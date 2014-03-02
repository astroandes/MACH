import numpy as np, pylab, sys

orig = np.loadtxt(sys.argv[1],delimiter=',')
mark = np.loadtxt(sys.argv[2]) 

mcmc = mark[mark[:,0].argsort()]

pylab.plot(orig[:,1],mcmc[:,4])
pylab.show()

pylab.plot(orig[:,2],mcmc[:,7])
pylab.show()
