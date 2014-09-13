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

index = [int(x) for x in np.argwhere(dens[:,9] > 100)]

id_mcmc = mcmc[:,0][index]
id_dens = dens[:,0][index]

c_mcmc = mcmc[:,4][index]
c_vel = vel[:,4][index]
c_dens = dens[:,4][index]

m_mcmc = mcmc[:,8][index]
m_vel = vel[:,7][index]
m_dens = dens[:,8][index]


bins_mcmc,median_mcmc,quar1_mcmc,quar2_mcmc = get_stats(np.log10(m_mcmc),c_mcmc)
bins_vel,median_vel,quar1_vel,quar2_vel = get_stats(np.log10(m_vel),c_vel)
bins_dens,median_dens,quar1_dens,quar2_dens = get_stats(np.log10(m_dens),c_dens)

pylab.plot(10**(bins_dens),median_dens,'-g',lw=2,label='Density')
pylab.fill_between(10**(bins_dens),quar1_dens,quar2_dens,color='g',alpha=0.4)

pylab.plot(10**(bins_mcmc),median_mcmc,'-r',lw=2,label='Mass')
pylab.fill_between(10**(bins_mcmc),quar1_mcmc,quar2_mcmc,color='r',alpha=0.4)

pylab.plot(10**(bins_vel),median_vel,'-b',lw=2,label='Velocity')
pylab.fill_between(10**(bins_vel),quar1_vel,quar2_vel,color='b',alpha=0.4)

pylab.legend(loc=1, borderaxespad=0.5)
pylab.xscale('log')
pylab.xlim([1E2,1E5])
pylab.ylim([0,10])
pylab.xlabel('$\mathrm{Number\ of\ Particles}$',fontsize=20)
pylab.ylabel('$\mathrm{Concentration}$',fontsize=20)
pylab.savefig('concentration.pdf',dpi=300)
pylab.close()


fig, (ax1,ax2,ax3) = pylab.subplots(1,3)

ax1.scatter(c_vel,c_mcmc,c='k',alpha=0.1)
ax1.plot(c_vel,c_vel,'-r',alpha=0.4,lw=2)
ax1.set_xlim(xmin=np.amin(c_vel),xmax=np.amax(c_vel))
ax1.set_ylim(ymin=np.amin(c_vel),ymax=np.amax(c_vel))
ax1.set_xlabel("$\mathrm{Velocity\ Method}$",fontsize=15)
ax1.set_ylabel("$\mathrm{Mass\ Method}$",fontsize=15)
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.scatter(c_mcmc,c_dens,c='k',alpha=0.1)
ax2.plot(c_mcmc,c_mcmc,'-r',alpha=0.4,lw=2)
ax2.set_xlim(xmin=np.amin(c_dens),xmax=np.amax(c_mcmc))
ax2.set_ylim(ymin=np.amin(c_dens),ymax=np.amax(c_mcmc))
ax2.set_xlabel("$\mathrm{Mass\ Method}$",fontsize=15)
ax2.set_ylabel("$\mathrm{Density\ Method}$",fontsize=15)
ax2.set_xscale('log')
ax2.set_yscale('log')

ax3.scatter(c_dens,c_vel,c='k',alpha=0.1)
ax3.plot(c_vel,c_vel,'-r',alpha=0.4,lw=2)
ax3.set_xlim(xmin=np.amin(c_vel),xmax=np.amax(c_vel))
ax3.set_ylim(ymin=np.amin(c_vel),ymax=np.amax(c_vel))
ax3.set_xlabel("$\mathrm{Density\ Method}$",fontsize=15)
ax3.set_ylabel("$\mathrm{Velocity\ Method}$",fontsize=15)
ax3.set_xscale('log')
ax3.set_yscale('log')

fig.subplots_adjust(wspace=.3)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
fig.set_size_inches(15,6)
pylab.savefig('mass-density-velocity.pdf',dpi=200,bbox_inches='tight')
pylab.savefig('mass-density-velocity.png',dpi=200,bbox_inches='tight')
pylab.close()
