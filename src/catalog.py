import numpy as np, pylab, sys, matplotlib, plotter
from matplotlib.mlab import griddata
from scipy.optimize import curve_fit

def distance(x1,y1,z1,x2,y2,z2):
    return np.sqrt(((x1-x2)**2.0)+((y1-y2)**2.0)+((z1-z2)**2.0))

mcmc = np.loadtxt(sys.argv[1],delimiter=',')
bdmv = np.loadtxt(sys.argv[2],delimiter=',')

mcmc_id = mcmc[:,0]
mcmc_x = mcmc[:,1]
mcmc_y = mcmc[:,2]
mcmc_z = mcmc[:,3]
mcmc_rho0 = mcmc[:,4]
mcmc_rho0_max = mcmc[:,5]
mcmc_rho0_min = mcmc[:,6]
mcmc_rs = mcmc[:,7]
mcmc_rs_max = mcmc[:,8]
mcmc_rs_min = mcmc[:,9]
mcmc_rvir = mcmc[:,10]
mcmc_c = mcmc[:,11]

bdmv_id = bdmv[:,0]
bdmv_x = bdmv[:,1]
bdmv_y = bdmv[:,2]
bdmv_z = bdmv[:,3]
bdmv_rvir = bdmv[:,7]
bdmv_c = bdmv[:,8]
bdmv_rs = bdmv_rvir/bdmv_c

used_id = np.empty((0))
used_x = np.empty((0))
used_y = np.empty((0))
used_z = np.empty((0))
used_rho0 = np.empty((0))
used_rho0_max = np.empty((0))
used_rho0_min = np.empty((0))
used_rs = np.empty((0))
used_rs_max = np.empty((0))
used_rs_min = np.empty((0))
used_rvir = np.empty((0))
used_c = np.empty((0))

FA = 0
NV = 0

for i in range(len(bdmv_x)):
    d = distance(bdmv_x[i],bdmv_y[i],bdmv_z[i],mcmc_x,mcmc_y,mcmc_z)
    j = np.argmin(d)
    d_min = d[j]

    if (d_min <= (mcmc_rvir[j] + bdmv_rvir[i]) and mcmc_rvir[j] > 0):
        used_id = np.append(used_id , mcmc_id[j])
        used_x = np.append(used_x , mcmc_x[j])
        used_y = np.append(used_y , mcmc_y[j])
        used_z = np.append(used_z , mcmc_z[j])
        used_rho0 = np.append(used_rho0 , mcmc_rho0[j])
        used_rho0_max = np.append(used_rho0_max , mcmc_rho0_max[j])
        used_rho0_min = np.append(used_rho0_min , mcmc_rho0_min[j])
        used_rs = np.append(used_rs , mcmc_rs[j])
        used_rs_max = np.append(used_rs_max , mcmc_rs_max[j])
        used_rs_min = np.append(used_rs_min , mcmc_rs_min[j])
        used_rvir = np.append(used_rvir , mcmc_rvir[j])
        used_c = np.append(used_c , mcmc_c[j])
    else:
        if(mcmc_rvir[j] > 0):
            FA = FA + 1
        else:
            NV = NV + 1
        used_id = np.append(used_id , 0)
        used_x = np.append(used_x , 0)
        used_y = np.append(used_y , 0)
        used_z = np.append(used_z , 0)
        used_rho0 = np.append(used_rho0 , 0)
        used_rho0_max = np.append(used_rho0_max , 0)
        used_rho0_min = np.append(used_rho0_min , 0)
        used_rs = np.append(used_rs , 0)
        used_rs_max = np.append(used_rs_max , 0)
        used_rs_min = np.append(used_rs_min , 0)
        used_rvir = np.append(used_rvir , 0)
        used_c = np.append(used_c , 0)

good_bdmv_x = np.empty((0))
good_bdmv_y = np.empty((0))
good_bdmv_z = np.empty((0))
good_bdmv_rs = np.empty((0))

good_used_x = np.empty((0))
good_used_y = np.empty((0))
good_used_z = np.empty((0))
good_used_rs = np.empty((0))
good_used_rs_max = np.empty((0))
good_used_rs_min = np.empty((0))

for i in range(len(bdmv_rs)):
    if(used_rs_min[i] <= bdmv_rs[i]  <= used_rs_max[i]):

        good_bdmv_x = np.append(good_bdmv_x,bdmv_x[i])
        good_bdmv_y = np.append(good_bdmv_y,bdmv_y[i])
        good_bdmv_z = np.append(good_bdmv_z,bdmv_z[i])
        good_bdmv_rs = np.append(good_bdmv_rs,bdmv_rs[i])

        good_used_x = np.append(good_used_x,used_x[i])
        good_used_y = np.append(good_used_y,used_y[i])
        good_used_z = np.append(good_used_z,used_z[i])
        good_used_rs = np.append(good_used_rs,used_rs[i])
        good_used_rs_max = np.append(good_used_rs_max,used_rs_max[i])
        good_used_rs_min = np.append(good_used_rs_min,used_rs_min[i])

sigma = np.array([good_used_rs_max-good_used_rs,good_used_rs-good_used_rs_min])

pylab.plot(bdmv_rs, used_rs,'.r')
pylab.plot(bdmv_rs, bdmv_rs,'-k')
pylab.plot(good_bdmv_rs,good_used_rs,'.b')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim([0,bdmv_rs[np.argmax(bdmv_rs)]])
pylab.ylim([0,bdmv_rs[np.argmax(bdmv_rs)]])
pylab.title('Scale Radius')
pylab.xlabel('BDMV')
pylab.ylabel('MCMC')
pylab.savefig('rs.png',format='png',dpi=600)
pylab.close()

pylab.plot(bdmv_rs, used_rs,'.r')
pylab.plot(bdmv_rs, bdmv_rs,'-k')
pylab.errorbar(good_bdmv_rs,good_used_rs,yerr=sigma,fmt='.b',elinewidth=1,capsize=2,capthick=1)
pylab.xscale('log')
pylab.yscale('log')
pylab.xlim([0,bdmv_rs[np.argmax(bdmv_rs)]])
pylab.ylim([0,bdmv_rs[np.argmax(bdmv_rs)]])
pylab.title('Scale Radius')
pylab.xlabel('BDMV')
pylab.ylabel('MCMC')
pylab.savefig('rs_errorbars.png',format='png',dpi=600)
pylab.close()

print 'Good Haloes: '+str(len(good_used_x))
print 'Bad Haloes: '+str(len(bdmv_x)-len(good_used_x)-FA)
print 'Forever Alone Haloes: '+str(FA)
print 'Non-Virialized Haloes: '+str(NV)
print 'Total Number of Haloes: '+str(len(bdmv_x))
