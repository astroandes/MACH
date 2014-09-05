import numpy as np, pylab,sys

def get_stats(log_mass,conc):

    bins = np.arange(0,np.amax(log_mass)+0.5,0.5)
    bins_index = np.digitize(log_mass,bins)
    median = np.empty(len(bins))
    quar1 = np.empty(len(bins))
    quar2 = np.empty(len(bins))

    for i in range(len(bins)):
        tempo = []
        for j in range(len(bins_index)):
            if (i == bins_index[j]):
                tempo.append(conc[j])
        tempo = np.array(tempo)
        if (tempo.size):
            median[i] = np.median(tempo)
            quar1[i],quar2[i] = np.percentile(tempo,[25,75])
        else:
            median[i],quar1[i],quar2[i]=0,0,0
    return bins,median,quar1,quar2

mcmc = np.loadtxt(sys.argv[1],delimiter=',')
vel = np.loadtxt(sys.argv[2],delimiter=',')
dens = np.loadtxt(sys.argv[3],delimiter=',')

id_mcmc = mcmc[:,0]
id_dens = dens[:,0]

c_mcmc = mcmc[:,4]
c_vel = vel[:,4]
c_dens = dens[:,4]

m_mcmc = mcmc[:,8]
m_vel = vel[:,7]
m_dens = dens[:,8]

c_mcmc_dens = np.empty(0)
c_vel_dens = np.empty(0)

'''
pylab.scatter(m_mcmc,c_mcmc,c='b')
pylab.scatter(m_vel,c_vel,c='y')
pylab.scatter(m_dens,c_dens,c='r')
pylab.xscale('log')
pylab.yscale('log')
pylab.show()
'''

bins_mcmc,median_mcmc,quar1_mcmc,quar2_mcmc = get_stats(np.log10(m_mcmc),c_mcmc)
bins_vel,median_vel,quar1_vel,quar2_vel = get_stats(np.log10(m_vel),c_vel)
bins_dens,median_dens,quar1_dens,quar2_dens = get_stats(np.log10(m_dens),c_dens)

pylab.plot(10**(bins_mcmc),median_mcmc,'-g',lw=2,label='Mass')
pylab.fill_between(10**(bins_mcmc),quar1_mcmc,quar2_mcmc,color='g',alpha=0.4)
pylab.plot(10**(bins_vel),median_vel,'-b',lw=2,label='Velocity')
pylab.fill_between(10**(bins_vel),quar1_vel,quar2_vel,color='b',alpha=0.4)
pylab.plot(10**(bins_dens),median_dens,'-r',lw=2,label='Density')
pylab.fill_between(10**(bins_dens),quar1_dens,quar2_dens,color='r',alpha=0.4)
pylab.legend(loc=1, borderaxespad=0.5)
pylab.xscale('log')
pylab.xlim([1E2,1E5])
pylab.xlabel('$Number\ of\ Particles$')
pylab.ylabel('$Concentration$')
pylab.savefig('concentration.pdf',dpi=300)
pylab.close()


for id in id_dens:
    index = np.where(id_mcmc == id)
    c_mcmc_dens = np.append(c_mcmc_dens,c_mcmc[index])
    c_vel_dens = np.append(c_vel_dens,c_vel[index])

pylab.scatter(c_vel,c_mcmc,c='k',alpha=0.1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(xmin=np.amin(c_vel))
pylab.ylim(ymin=np.amin(c_mcmc))
pylab.xlabel("$Velocity\ Method$")
pylab.ylabel("$Mass\ Method$")
pylab.savefig('velocity-mass.png',dpi=200)
pylab.close()

pylab.scatter(c_dens,c_mcmc_dens,c='k',alpha=0.1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(xmin=np.amin(c_dens))
pylab.ylim(ymin=np.amin(c_mcmc_dens))
pylab.xlabel("$Density\ Method$")
pylab.ylabel("$Mass\ Method$")
pylab.savefig('density-mass.png',dpi=200)
pylab.close()

pylab.scatter(c_dens,c_vel_dens,c='k',alpha=0.1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim(xmin=np.amin(c_dens))
pylab.ylim(ymin=np.amin(c_vel_dens))
pylab.xlabel("$Density\ Method$")
pylab.ylabel("$Velocity\ Method$")
pylab.savefig('density-velocity.png',dpi=200)
pylab.close()
