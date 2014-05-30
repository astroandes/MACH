import numpy as np, pylab, sys, matplotlib
from matplotlib.mlab import griddata
from scipy.optimize import curve_fit

def distance(x1,y1,z1,x2,y2,z2):
    return np.sqrt(((x1-x2)**2.0)+((y1-y2)**2.0)+((z1-z2)**2.0))

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
        median[i] = np.median(tempo)
        if (tempo.size):
            quar1[i],quar2[i] = np.percentile(tempo,[25,75])            
        else:
            quar1[i],quar2[i]=0,0
    return bins,median,quar1,quar2

dm = 8.721e9

mcmc = np.loadtxt(sys.argv[1],delimiter=',')
bdm = np.loadtxt(sys.argv[2],delimiter=',')

mcmc_id = mcmc[:,0]
mcmc_x = mcmc[:,1]
mcmc_y = mcmc[:,2]
mcmc_z = mcmc[:,3]
mcmc_c = mcmc[:,4]
mcmc_c_max = mcmc[:,5]
mcmc_c_min = mcmc[:,6]
mcmc_r = mcmc[:,7]
mcmc_m = mcmc[:,8]*dm

bdm_id = bdm[:,0]
bdm_x = bdm[:,4]
bdm_y = bdm[:,5]
bdm_z = bdm[:,6]
bdm_m = bdm[:,11]
bdm_r = bdm[:,13]
bdm_c = bdm[:,16]

used_id = np.empty(0)
used_x = np.empty(0)
used_y = np.empty(0)
used_z = np.empty(0)
used_c = np.empty(0)
used_c_max = np.empty(0)
used_c_min = np.empty(0)
used_r = np.empty(0)
used_m = np.empty(0)

for i in range(len(bdm_id)):
	r = distance(bdm_x[i],bdm_y[i],bdm_z[i],mcmc_x,mcmc_y,mcmc_z)
	j = np.argmin(r)
	used_id = np.append(used_id,mcmc_id[j])
	used_x = np.append(used_x,mcmc_x[j])
	used_y = np.append(used_y,mcmc_y[j])
	used_z = np.append(used_z,mcmc_z[j])
	used_c = np.append(used_c,mcmc_c[j])
	used_c_max = np.append(used_c_max,mcmc_c_max[j])
	used_c_min = np.append(used_c_min,mcmc_c_min[j])
	used_r = np.append(used_r,mcmc_r[j])
	used_m = np.append(used_m,mcmc_m[j])

pylab.scatter(bdm_c,used_c)
pylab.plot(bdm_c,bdm_c,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$BDM$')
pylab.ylabel('$MCMC$')
pylab.title('$Concentration$')
pylab.savefig('concentration.pdf',dpi=300)
pylab.close()

log_m_mcmc = np.log10(used_m)
log_m_bdm = np.log10(bdm_m)

bins_mcmc,median_c_mcmc,quartile1_c_mcmc,quartile2_c_mcmc = get_stats(log_m_mcmc,used_c)
bins_bdm,median_c_bdm,quartile1_c_bdm,quartile2_c_bdm = get_stats(log_m_bdm,bdm_c)

pylab.plot(10**(bins_mcmc),median_c_mcmc,'-r',lw=2,label='MCMC')
pylab.fill_between(10**(bins_mcmc),quartile1_c_mcmc,quartile2_c_mcmc,color='r',alpha=0.3)
pylab.plot(10**(bins_bdm),median_c_bdm,'-b',lw=2,label='BDM')
pylab.fill_between(10**(bins_bdm),quartile1_c_bdm,quartile2_c_bdm,color='b',alpha=0.3)
pylab.legend(loc=1, borderaxespad=0.5)
pylab.xscale('log')
pylab.xlim([1E11,1E15])
#pylab.ylim([0,15])
pylab.xlabel('$Mass\ (M_{sun}/h )$')
pylab.ylabel('$Concentration$')
pylab.savefig('concentration2.pdf',dpi=300)
pylab.close()
