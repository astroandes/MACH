import numpy as np, pylab,sys

dm = 8720999991.75 #https://raw.githubusercontent.com/forero/TrainingHalos/master/data/halo_ID/MDmini_all_halo-ID.csv

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


fig, ax1 = pylab.subplots()

ax1.plot(10**(bins_dens)*dm,median_dens,'-g',lw=2,label='$\mathrm{Density\ Method}$')
ax1.fill_between(10**(bins_dens)*dm,quar1_dens,quar2_dens,color='g',alpha=0.4)

ax1.plot(10**(bins_mcmc)*dm,median_mcmc,'-r',lw=2,label='$\mathrm{Our\ Method}$')
ax1.fill_between(10**(bins_mcmc)*dm,quar1_mcmc,quar2_mcmc,color='r',alpha=0.4)

ax1.plot(10**(bins_vel)*dm,median_vel,'-b',lw=2,label='$\mathrm{Velocity\ Method}$')
ax1.fill_between(10**(bins_vel)*dm,quar1_vel,quar2_vel,color='b',alpha=0.4)

ax1.legend(loc=1, borderaxespad=0.5)
ax1.set_xscale('log')
ax1.set_xlim([1E2*dm,1E5*dm])
ax1.set_ylim([0,10])
ax1.set_xlabel('$\mathrm{Halo\ Mass\ (M_{\odot})}$',fontsize=20)
ax1.set_ylabel('$\mathrm{Concentration}$',fontsize=20)

ax2 = ax1.twiny()
ax2.plot(10**(bins_dens),median_dens,lw=0)
ax2.set_xscale('log')
ax2.set_xlim([1E2,1E5])
ax2.set_ylim([0,10])
ax2.set_xlabel('$\mathrm{Number\ of\ particles}$',fontsize=20)
pylab.tight_layout()
pylab.savefig('concentration.pdf',dpi=300)
pylab.close()


fig, ax1 = pylab.subplots()
ax1.scatter(c_mcmc,c_vel,c='k',alpha=0.1)
ax1.plot(c_mcmc,c_mcmc,'-r',alpha=0.4,lw=2)
ax1.set_xlim(xmin=np.amin(c_vel),xmax=np.amax(c_mcmc))
ax1.set_ylim(ymin=np.amin(c_vel),ymax=np.amax(c_mcmc))
ax1.set_xlabel("$\mathrm{Concentration\ (Our\ Method)}$",fontsize=15)
ax1.set_ylabel("$\mathrm{Concentration\ (Velocity\ Method)}$",fontsize=15)
ax1.set_xscale('log')
ax1.set_yscale('log')
pylab.savefig('mass-velocity.pdf',dpi=200,bbox_inches='tight')
pylab.savefig('mass-velocity.png',dpi=200,bbox_inches='tight')
pylab.close()

fig, ax2 = pylab.subplots()
ax2.scatter(c_mcmc,c_dens,c='k',alpha=0.1)
ax2.plot(c_mcmc,c_mcmc,'-r',alpha=0.4,lw=2)
ax2.set_xlim(xmin=np.amin(c_dens),xmax=np.amax(c_mcmc))
ax2.set_ylim(ymin=np.amin(c_dens),ymax=np.amax(c_mcmc))
ax2.set_xlabel("$\mathrm{Concentration\ (Our\ Method)}$",fontsize=25)
ax2.set_ylabel("$\mathrm{Concentration\ (Density\ Method)}$",fontsize=25)
ax2.set_xscale('log')
ax2.set_yscale('log')
pylab.savefig('mass-density.pdf',dpi=200,bbox_inches='tight')
pylab.savefig('mass-density.png',dpi=200,bbox_inches='tight')
pylab.close()

fig, ax3 = pylab.subplots()
ax3.scatter(c_dens,c_vel,c='k',alpha=0.1)
ax3.plot(c_vel,c_vel,'-r',alpha=0.4,lw=2)
ax3.set_xlim(xmin=np.amin(c_vel),xmax=np.amax(c_vel))
ax3.set_ylim(ymin=np.amin(c_vel),ymax=np.amax(c_vel))
ax3.set_xlabel("$\mathrm{Concentration\ (Density\ Method)}$",fontsize=25)
ax3.set_ylabel("$\mathrm{Concentration\ (Velocity\ Method)}$",fontsize=25)
ax3.set_xscale('log')
ax3.set_yscale('log')
pylab.savefig('density-velocity.pdf',dpi=200,bbox_inches='tight')
pylab.savefig('density-velocity.png',dpi=200,bbox_inches='tight')
pylab.close()
