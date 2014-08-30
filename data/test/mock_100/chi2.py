import numpy as np, pylab

vlct = np.loadtxt('mock_100_velocity/table.csv',delimiter=',')
mass = np.loadtxt('mock_100_mass/table.csv',delimiter=',')
dens = np.loadtxt('mock_100_density/table.csv',delimiter=',')

orig = vlct[:,0]

def chi2(x_obs,x_org):
    return np.sum((x_obs-x_org)**2)

pylab.title('$Velocity$')
pylab.plot(orig,orig,label='$real$')
pylab.plot(orig,vlct[:,1],label='$n=20$')
pylab.plot(orig,vlct[:,2],label='$n=200$')
pylab.plot(orig,vlct[:,3],label='$n=2000$')
pylab.plot(orig,vlct[:,4],label='$n=20000$')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('vlct.png',dpi=300)
pylab.close()


chi_vlct = [chi2(orig,vlct[:,1]),chi2(orig,vlct[:,2]),chi2(orig,vlct[:,3]),chi2(orig,vlct[:,4])]

orig = dens[:,0]

pylab.title('$Density$')
pylab.plot(orig,orig,label='$real$')
pylab.plot(orig,dens[:,1],label='$n=20$')
pylab.plot(orig,dens[:,2],label='$n=200$')
pylab.plot(orig,dens[:,3],label='$n=2000$')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('dens.png',dpi=300)
pylab.close()

chi_dens =[chi2(orig,dens[:,1]),chi2(orig,dens[:,2]),chi2(orig,dens[:,3])]

orig = mass[:,0]

pylab.title('$Mass$')
pylab.plot(orig,orig,label='$real$')
pylab.plot(orig,mass[:,1],label='$n=20$')
pylab.plot(orig,mass[:,2],label='$n=200$')
pylab.plot(orig,mass[:,3],label='$n=2000$')
pylab.plot(orig,mass[:,4],label='$n=20000$')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('mass.png',dpi=300)
pylab.close()

chi_mass = [chi2(orig,mass[:,1]),chi2(orig,mass[:,2]),chi2(orig,mass[:,3]),chi2(orig,mass[:,4])]

pylab.plot([20,200,2000,20000],chi_mass,'or',label='$Mass$')
pylab.plot([20,200,2000,20000],chi_mass,'-r')
pylab.plot([20,200,2000,20000],chi_vlct,'sg',label='$Velocty$')
pylab.plot([20,200,2000,20000],chi_vlct,'-g')
pylab.plot([200,2000,20000],chi_dens,'^b',label='$Density$')
pylab.plot([200,2000,20000],chi_dens,'-b')
pylab.xscale('log')
pylab.yscale('log')
pylab.legend()
pylab.xlabel('$Number\ of\ Points$',fontsize=20)
pylab.ylabel('$\chi^{2}$',fontsize=20)
pylab.savefig('chi2.png',dpi=200)
pylab.close()
