import numpy as np, pylab, sys, matplotlib, plotter
from matplotlib.mlab import griddata
from scipy.optimize import curve_fit

def distance(x1,y1,z1,x2,y2,z2):
    return np.sqrt(((x1-x2)**2.0)+((y1-y2)**2.0)+((z1-z2)**2.0))

def get_stats(log_mass,conc):

    bins = np.arange(0,np.amax(log_mass)+0.5,0.5)
    bins_index = np.digitize(log_mass,bins)
    median = np.empty(len(bins))
    quartile = np.empty(len(bins))

    for i in range(len(bins)):
        tempo = []
        for j in range(len(bins_index)):
            if (i == bins_index[j]):
                tempo.append(conc[j])
        tempo = np.array(tempo)
        median[i] = np.median(tempo)
        quartile[i] = 0
    return bins,median,quartile

dm = 8.721e9

mcmc = np.loadtxt(sys.argv[1],delimiter=',')
bdmv = np.loadtxt(sys.argv[2],delimiter=',')
bdmw = np.loadtxt(sys.argv[3],delimiter=',')

mcmc_id = mcmc[:,0]
mcmc_x = mcmc[:,1]
mcmc_y = mcmc[:,2]
mcmc_z = mcmc[:,3]
mcmc_c_bdmv = mcmc[:,4]
mcmc_c_bdmv_max = mcmc[:,5]
mcmc_c_bdmv_min = mcmc[:,6]
mcmc_c_bdmw = mcmc[:,7]
mcmc_c_bdmw_max = mcmc[:,8]
mcmc_c_bdmw_min = mcmc[:,9]
mcmc_r_bdmv = mcmc[:,10]
mcmc_r_bdmw = mcmc[:,11]
mcmc_n_bdmv = mcmc[:,12]*dm
mcmc_n_bdmw = mcmc[:,13]*dm

bdmv_id = bdmv[:,0]
bdmv_x = bdmv[:,1]
bdmv_y = bdmv[:,2]
bdmv_z = bdmv[:,3]
bdmv_m = bdmv[:,6]
bdmv_r = bdmv[:,7]
bdmv_c = bdmv[:,8]

bdmw_id = bdmw[:,0]
bdmw_x = bdmw[:,1]
bdmw_y = bdmw[:,2]
bdmw_z = bdmw[:,3]
bdmw_m = bdmw[:,6]
bdmw_r = bdmw[:,7]
bdmw_c = bdmw[:,8]

used_id_bdmv = np.empty(0)
used_x_bdmv = np.empty(0)
used_y_bdmv = np.empty(0)
used_z_bdmv = np.empty(0)
used_id_bdmw = np.empty(0)
used_x_bdmw = np.empty(0)
used_y_bdmw = np.empty(0)
used_z_bdmw = np.empty(0)
used_c_bdmv = np.empty(0)
used_c_bdmv_max = np.empty(0)
used_c_bdmv_min = np.empty(0)
used_c_bdmw = np.empty(0)
used_c_bdmw_max = np.empty(0)
used_c_bdmw_min = np.empty(0)
used_r_bdmv = np.empty(0)
used_r_bdmw = np.empty(0)
used_n_bdmv = np.empty(0)
used_n_bdmw = np.empty(0)

for i in range(len(bdmv_id)):
	r = distance(bdmv_x[i],bdmv_y[i],bdmv_z[i],mcmc_x,mcmc_y,mcmc_z)
	j = np.argmin(r)
	used_id_bdmv = np.append(used_id_bdmv,mcmc_id[j])
	used_x_bdmv = np.append(used_x_bdmv,mcmc_x[j])
	used_y_bdmv = np.append(used_y_bdmv,mcmc_y[j])
	used_z_bdmv = np.append(used_z_bdmv,mcmc_z[j])
	used_c_bdmv = np.append(used_c_bdmv,mcmc_c_bdmv[j])
	used_c_bdmv_max = np.append(used_c_bdmv_max,mcmc_c_bdmv_max[j])
	used_c_bdmv_min = np.append(used_c_bdmv_min,mcmc_c_bdmv_min[j])
	used_r_bdmv = np.append(used_r_bdmv,mcmc_r_bdmv[j])
	used_n_bdmv = np.append(used_n_bdmv,mcmc_n_bdmv[j])

for i in range(len(bdmw_id)):
	r = distance(bdmw_x[i],bdmw_y[i],bdmw_z[i],mcmc_x,mcmc_y,mcmc_z)
	j = np.argmin(r)
	used_id_bdmw = np.append(used_id_bdmw,mcmc_id[j])
	used_x_bdmw = np.append(used_x_bdmw,mcmc_x[j])
	used_y_bdmw = np.append(used_y_bdmw,mcmc_y[j])
	used_z_bdmw = np.append(used_z_bdmw,mcmc_z[j])
	used_c_bdmw = np.append(used_c_bdmw,mcmc_c_bdmw[j])
	used_c_bdmw_max = np.append(used_c_bdmw_max,mcmc_c_bdmw_max[j])
	used_c_bdmw_min = np.append(used_c_bdmw_min,mcmc_c_bdmw_min[j])
	used_r_bdmw = np.append(used_r_bdmw,mcmc_r_bdmw[j])
	used_n_bdmw = np.append(used_n_bdmw,mcmc_n_bdmw[j])

pylab.plot(bdmv_c,used_c_bdmv,'.b')
pylab.plot(bdmv_c,bdmv_c,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$BDMV$')
pylab.ylabel('$MCMC$')
pylab.title('$Concentration$')
pylab.savefig('bdmv.png',dpi=300)
pylab.close()

pylab.plot(bdmw_c,used_c_bdmw,'.b')
pylab.plot(bdmw_c,bdmw_c,'-k')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$BDMW$')
pylab.ylabel('$MCMC$')
pylab.title('$Concentration$')
pylab.savefig('bdmw.png',dpi=300)
pylab.close()

log_m_bdmv_mcmc = np.log10(used_n_bdmv)
log_m_bdmw_mcmc = np.log10(used_n_bdmw)
log_m_bdmv_orig = np.log10(bdmv_m)
log_m_bdmw_orig = np.log10(bdmw_m)

bins_bdmv_mcmc,median_c_bdmv_mcmc,quartile_c_bdmv_mcmc  = get_stats(log_m_bdmv_mcmc,used_c_bdmv)
bins_bdmw_mcmc,median_c_bdmw_mcmc,quartile_c_bdmw_mcmc = get_stats(log_m_bdmw_mcmc,used_c_bdmw)
bins_bdmv_orig,median_c_bdmv_orig,quartile_c_bdmv_orig = get_stats(log_m_bdmv_orig,bdmv_c)
bins_bdmw_orig,median_c_bdmw_orig,quartile_c_bdmw_orig = get_stats(log_m_bdmw_orig,bdmw_c)

pylab.plot(10**(bins_bdmv_mcmc),median_c_bdmv_mcmc,'-r',lw=2,label='MCMC-V')
pylab.plot(10**(bins_bdmw_mcmc),median_c_bdmw_mcmc,'-b',lw=2,label='MCMC-W')
pylab.plot(10**(bins_bdmv_orig),median_c_bdmv_orig,'--r',lw=2,label='BDMV')
pylab.plot(10**(bins_bdmw_orig),median_c_bdmw_orig,'--b',lw=2,label='BDMW')
pylab.legend(loc=1, borderaxespad=0.5)
pylab.xscale('log')
pylab.xlim([1E11,1E15])
pylab.xlabel('$Mass\ (M_{sun}/h )$')
pylab.ylabel('$Concentration$')
pylab.savefig('concentration.png',dpi=300)
pylab.close()

pylab.plot(used_n_bdmv,used_c_bdmv,'.b')
pylab.plot(10**(bins_bdmv_mcmc),median_c_bdmv_mcmc,'-r',lw=2)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Mass\ (M_{sun}/h )$')
pylab.ylabel('$Concentration$')
pylab.title('$BDMV$')
pylab.savefig('n_vs_c_bdmv.png',dpi=300)
pylab.close()

pylab.plot(used_n_bdmw,used_c_bdmw,'.b')
pylab.plot(10**(bins_bdmw_mcmc),median_c_bdmw_mcmc,'-r',lw=2)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('$Mass\ (M_{sun}/h )$')
pylab.ylabel('$Concentration$')
pylab.title('$BDMW$')
pylab.savefig('n_vs_c_bdmw.png',dpi=300)
pylab.close()
