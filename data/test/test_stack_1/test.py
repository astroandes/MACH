import numpy as np, pylab

real = np.loadtxt('parameters.dat',delimiter=',')
auto = np.loadtxt('auto_center.csv',delimiter=',')
gets = np.loadtxt('gets_center.csv',delimiter=',')

auto = auto[auto[:,0].argsort()]
gets =gets[gets[:,0].argsort()]

real_a = real[:,0]
auto_a,auto_amin,auto_amax = auto[:,4],auto[:,5],auto[:,6]
gets_a,gets_amin,gets_amax = gets[:,4],gets[:,5],gets[:,6]
auto_asigma = np.array([auto_amax-auto_a,auto_a-auto_amin])
gets_asigma = np.array([gets_amax-gets_a,gets_a-gets_amin])

real_b = real[:,1]
auto_b,auto_bmin,auto_bmax = auto[:,7],auto[:,8],auto[:,9]
gets_b,gets_bmin,gets_bmax = gets[:,7],gets[:,8],gets[:,9]
auto_bsigma = np.array([auto_bmax-auto_b,auto_b-auto_bmin])
gets_bsigma = np.array([gets_bmax-gets_b,gets_b-gets_bmin])


pylab.plot(real_a,real_a,'-k')
pylab.errorbar(real_a,auto_a,yerr=auto_asigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$Obtained$')
pylab.title('$\\rho_0$ with center at 0')
pylab.savefig('rho0_auto.png',dpi=300)
pylab.close()

pylab.plot(real_b,real_b,'-k')
pylab.errorbar(real_b,auto_b,yerr=auto_bsigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$Obtained$')
pylab.title('$R_s$ with center at 0')
pylab.savefig('Rs_auto.png',dpi=300)
pylab.close()

pylab.plot(real_a,real_a,'-k')
pylab.errorbar(real_a,gets_a,yerr=gets_asigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$Obtained$')
pylab.title('$\\rho_0$')
pylab.savefig('rho0_gets.png',dpi=300)
pylab.close()

pylab.plot(real_b,real_b,'-k')
pylab.errorbar(real_b,gets_b,yerr=gets_bsigma,fmt='.r',elinewidth=1,capsize=2,capthick=1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$Obtained$')
pylab.title('$R_s$')
pylab.savefig('Rs_gets.png',dpi=300)
pylab.close()

pylab.plot(real_a,np.abs(auto_a-real_a)/real_a,'or')
pylab.xscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$|Obtained-Real|/Real$')
pylab.title('$\\rho_0$ with center at 0')
pylab.savefig('rho0_auto_relative.png',dpi=300)
pylab.close()

pylab.plot(real_b,np.abs(auto_b-real_b)/real_b,'or')
pylab.xscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$|Obtained-Real|/Real$')
pylab.title('$R_s$ with center at 0')
pylab.savefig('Rs_auto_relative.png',dpi=300)
pylab.close()

pylab.plot(real_a,np.abs(gets_a-real_a)/real_a,'or')
pylab.xscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$|Obtained-Real|/Real$')
pylab.title('$\\rho_0$')
pylab.savefig('rho0_gets_relative.png',dpi=300)
pylab.close()

pylab.plot(real_b,np.abs(gets_b-real_b)/real_b,'or')
pylab.xscale('log')
pylab.xlabel('$Real$')
pylab.ylabel('$|Obtained-Real|/Real$')
pylab.title('$R_s$')
pylab.savefig('Rs_gets_relative.png',dpi=300)
pylab.close()


