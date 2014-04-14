import numpy as np, pylab

data = np.loadtxt('gets_center.csv',delimiter=',')

orig = np.loadtxt('parameters.dat',delimiter=',')
data = data[data[:,0].argsort()]

c_bdmv = data[:,4]
c_bdmw = data[:,7]

r_bdmv = data[:,10]
r_bdmw = data[:,11]

rs_bdmv = r_bdmv/c_bdmv
rs_bdmw = r_bdmw/c_bdmw

rs = orig[:,1]

pylab.title('$r_s$')
pylab.plot(rs,rs,'-k')
pylab.plot(rs,rs_bdmv,'ob')
pylab.xlabel('Real')
pylab.ylabel('BDMV')
pylab.savefig('gets_bdmv.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs-rs,'-k')
pylab.plot(rs,rs_bdmv-rs,'ob')
pylab.xlabel('Real')
pylab.ylabel('BDMV-Real')
pylab.savefig('gets_bdmv-real.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,np.abs(rs_bdmv-rs)/rs,'ob')
pylab.xlabel('Real')
pylab.ylabel('|BDMV-Real|/Real')
pylab.savefig('gets_bdmv_relative.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,np.abs(rs_bdmw-rs)/rs,'or')
pylab.xlabel('Real')
pylab.ylabel('|BDMW-Real|/Real')
pylab.savefig('gets_bdmw_relative.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs,'-k')
pylab.plot(rs,rs_bdmw,'or')
pylab.xlabel('Real')
pylab.ylabel('BDMW')
pylab.savefig('gets_bdmw.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs-rs,'-k')
pylab.plot(rs,rs_bdmw-rs,'or')
pylab.xlabel('Real')
pylab.ylabel('BDMW-Real')
pylab.savefig('gets_bdmw-real.png',dpi=300)
pylab.close()

data = np.loadtxt('auto_center.csv',delimiter=',')

orig = np.loadtxt('parameters.dat',delimiter=',')
data = data[data[:,0].argsort()]

c_bdmv = data[:,4]
c_bdmw = data[:,7]

r_bdmv = data[:,10]
r_bdmw = data[:,11]

rs_bdmv = r_bdmv/c_bdmv
rs_bdmw = r_bdmw/c_bdmw

rs = orig[:,1]

pylab.title('$r_s$')
pylab.plot(rs,rs,'-k')
pylab.plot(rs,rs_bdmv,'ob')
pylab.xlabel('Real')
pylab.ylabel('BDMV')
pylab.savefig('auto_bdmv.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs-rs,'-k')
pylab.plot(rs,rs_bdmv-rs,'ob')
pylab.xlabel('Real')
pylab.ylabel('BDMV-Real')
pylab.savefig('auto_bdmv-real.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs,'-k')
pylab.plot(rs,rs_bdmw,'or')
pylab.xlabel('Real')
pylab.ylabel('BDMW')
pylab.savefig('auto_bdmw.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,rs-rs,'-k')
pylab.plot(rs,rs_bdmw-rs,'or')
pylab.xlabel('Real')
pylab.ylabel('BDMW-Real')
pylab.savefig('auto_bdmw-real.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,np.abs(rs_bdmv-rs)/rs,'ob')
pylab.xlabel('Real')
pylab.ylabel('|BDMV-Real|/Real')
pylab.savefig('auto_bdmv_relative.png',dpi=300)
pylab.close()

pylab.title('$r_s$')
pylab.plot(rs,np.abs(rs_bdmw-rs)/rs,'or')
pylab.xlabel('Real')
pylab.ylabel('|BDMW-Real|/Real')
pylab.savefig('auto_bdmw_relative.png',dpi=300)
pylab.close()
