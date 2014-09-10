import numpy as np, pylab

vlct = np.loadtxt('velocity/table.csv',delimiter=',')
mass = np.loadtxt('mass/table.csv',delimiter=',')
dens = np.loadtxt('density/table.csv',delimiter=',')

orig = vlct[:,0]

def chi2(x_obs,x_org):
    return np.sum(np.abs((x_obs-x_org)/x_org))/len(x_obs)

def best_fit(original,results):
    ans = np.empty(len(original))
    for i in range(len(original)):
        concentrations = results[i,:]
        index = np.argmin(np.abs(original[i]-concentrations))
        ans[i] = concentrations[index]
        print original[i],concentrations,index
    return ans
    
pylab.title('$Velocity$')
pylab.plot(orig,orig,'--k',label='$real$')
#pylab.plot(orig,vlct[:,1],'.b',label='$n=20$')
#pylab.plot(orig,vlct[:,2],'.r',label='$n=200$')
#pylab.plot(orig,vlct[:,3],'.g',label='$n=2000$')
#pylab.plot(orig,vlct[:,4],'.c',label='$n=20000$')
pylab.plot(orig,best_fit(orig,vlct[:,1:5]),'-r',label='$best\ fit$')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('vlct.png',dpi=300)
pylab.close()


chi_vlct = np.array([chi2(orig,vlct[:,1]),chi2(orig,vlct[:,2]),chi2(orig,vlct[:,3]),chi2(orig,vlct[:,4])])

orig = mass[:,0]

pylab.title('$Mass$')
pylab.plot(orig,orig,'--k',label='$real$')
#pylab.plot(orig,mass[:,1],'.',label='$n=20$',)
#pylab.plot(orig,mass[:,2],'.',label='$n=200$')
#pylab.plot(orig,mass[:,3],'.',label='$n=2000$',)
#pylab.plot(orig,mass[:,4],'.',label='$n=20000$')
pylab.plot(orig,best_fit(orig,mass[:,1:5]),'-r',label='$best\ fit$')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('mass.png',dpi=300)
pylab.close()

chi_mass = np.array([chi2(orig,mass[:,1]),chi2(orig,mass[:,2]),chi2(orig,mass[:,3]),chi2(orig,mass[:,4])])

orig = dens[:,0]

pylab.title('$Density$')
pylab.plot(orig,orig,label='$real$')
#pylab.plot(orig,dens[:,1],'.',label='$n=20$')
#pylab.plot(orig,dens[:,2],'.',label='$n=200$')
#pylab.plot(orig,dens[:,3],'.',label='$n=2000$')
pylab.plot(orig,best_fit(orig,dens[:,1:5]),'-r',label='$best\ fit$')
pylab.xlabel('$Original$')
pylab.ylabel('$Obtained$')
pylab.legend(loc=4, borderaxespad=0.5)
pylab.savefig('dens.png',dpi=300)
pylab.close()

chi_dens = np.array([chi2(orig,dens[:,1]),chi2(orig,dens[:,2]),chi2(orig,dens[:,3]),chi2(orig,dens[:,4])])

pylab.plot([20,200,2000,20000],chi_mass,'-k',lw=2)
pylab.plot([20,200,2000,20000],chi_mass,'ok',label='$Mass$',ms=7)
pylab.plot([20,200,2000,20000],chi_vlct,'-k',lw=2)
pylab.plot([20,200,2000,20000],chi_vlct,'sk',label='$Velocity$',ms=7)
pylab.plot([20,200,2000,20000],chi_dens,'-k',lw=2)
pylab.plot([20,200,2000,20000],chi_dens,'^k',label='$Density$',ms=7)
pylab.xscale('log')
pylab.yscale('log')
pylab.legend(loc=3)
pylab.xlabel('$Number\ of\ Points$',fontsize=20)
pylab.ylabel('$\\xi$',fontsize=20)
pylab.savefig('error.png',dpi=200)
pylab.savefig('error.pdf',dpi=200)
pylab.close()

c = orig
c_dens = dens[:,3]
c_mass = np.empty(0)
c_vlct = np.empty(0)

for id in dens[:,0]:
    index = np.where(dens[:,0] == id)
    c_mass = np.append(c_mass,mass[index,4])
    c_vlct = np.append(c_vlct,vlct[index,4])
