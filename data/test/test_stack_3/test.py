import numpy as np, pylab

data = np.loadtxt('gets_center.csv',delimiter=',')
orig = np.loadtxt('parameters.dat',delimiter=',')
data = data[data[:,0].argsort()]
print data
c_bdmv = data[:,4]
c_bdmw = data[:,7]

r_bdmv = data[:,10]
r_bdmw = data[:,11]

rs_bdmv = r_bdmv/c_bdmv
rs_bdmw = r_bdmw/c_bdmw

rs = orig[:,1]

pylab.plot(rs,rs,'-k')
pylab.plot(rs,rs_bdmv,'ob')
pylab.plot(rs,rs_bdmw,'or')
pylab.show()
