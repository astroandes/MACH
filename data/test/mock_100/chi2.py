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
pylab.xlabel=('$Original$')
pylab.xlabel=('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('vlct.png',dpi=300)
pylab.close()

orig = dens[:,0]

pylab.title('$Density$')
pylab.plot(orig,orig,label='$real$')
pylab.plot(orig,dens[:,1],label='$n=20$')
pylab.plot(orig,dens[:,2],label='$n=200$')
pylab.plot(orig,dens[:,3],label='$n=2000$')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel=('$Original$')
pylab.xlabel=('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('dens.png',dpi=300)
pylab.close()

orig = mass[:,0]

pylab.title('$Mass$')
pylab.plot(orig,orig,label='$real$')
pylab.plot(orig,mass[:,1],label='$n=20$')
pylab.plot(orig,mass[:,2],label='$n=200$')
pylab.plot(orig,mass[:,3],label='$n=2000$')
pylab.plot(orig,mass[:,4],label='$n=20000$')
pylab.xlabel=('$Original$')
pylab.xlabel=('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('mass.png',dpi=300)
pylab.close()
