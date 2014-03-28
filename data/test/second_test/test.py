import numpy as np, pylab

real = np.loadtxt('parameters.dat',delimiter=',')
auto = np.loadtxt('auto_center.csv',delimiter=',')
gets = np.loadtxt('gets_center.csv',delimiter=',')

auto = auto[auto[:,0].argsort()]
gets =gets[gets[:,0].argsort()]

real_a = real[:,0]
auto_a = auto[:,4]
gets_a = gets[:,4]

real_b = real[:,1]
auto_b = auto[:,7]
gets_b = gets[:,7]

pylab.plot(real_a,auto_a,'or')
pylab.plot(real_a,real_a,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Real')
pylab.ylabel('Obtained')
pylab.title('$\\rho_0$ with center at 0')
pylab.savefig('rho0_auto.png',dpi=300)
pylab.close()

pylab.plot(real_b,auto_b,'or')
pylab.plot(real_b,real_b,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Real')
pylab.ylabel('Obtained')
pylab.title('$R_s$ with center at 0')
pylab.savefig('Rs_auto.png',dpi=300)
pylab.close()

pylab.plot(real_a,gets_a,'or')
pylab.plot(real_a,real_a,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Real')
pylab.ylabel('Obtained')
pylab.title('$\\rho_0$')
pylab.savefig('rho0_gets.png',dpi=300)
pylab.close()

pylab.plot(real_b,gets_b,'or')
pylab.plot(real_b,real_b,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Real')
pylab.ylabel('Obtained')
pylab.title('$R_s$')
pylab.savefig('Rs_gets.png',dpi=300)
pylab.close()

