import numpy as np, emcee
import sys, os
from time import clock

def nfw_log_mass(log_r,log_c):
    return np.log((np.log(1.0+np.exp(log_c+log_r))-np.exp(log_c+log_r)/(1.0+np.exp(log_c+log_r)))/(np.log(1.0+np.exp(log_c))-np.exp(log_c)/(1.0+np.exp(log_c))))

def log_likelihood(log_c,log_r,log_m):
    if log_c >= 0:
        s = np.sum(log_m-nfw_log_mass(log_r,log_c))
        return -.5*s*s
    return -np.inf

def center(data):
    n = len(data)

    dr = data.reshape(n,1,3)-data
    dr = np.sum(dr*dr,axis=2)
    dr = np.sqrt(dr)

    U = -1.0/np.sum(dr,axis=1)
    n_min = np.argmin(U)

    return data[n_min]

def main():
    results = []
    outfile = open('output.csv','w')
    try:
        root_path = sys.argv[1]
        directory_list = [root_path+'/'+x for x in os.listdir(root_path)]
        directory_list.sort()

        guess = 2
        n_dimensions = 1
        n_walkers = 2

        for filename in directory_list:
            t_init = clock()
            data = np.loadtxt(filename,delimiter=',')
            #x_center,y_center,z_center = center(data)
            x_center,y_center,z_center = 0,0,0
            data = np.transpose(data)
            x,y,z = data
            del data

            x -= x_center
            y -= y_center
            z -= z_center

            r = np.sqrt(x*x+y*y+z*z)
            r = np.sort(r)
            r = np.delete(r,0)
            r = r/r[-1]
            del x,y,z

            n_particles = len(r)

            m = np.arange(n_particles)+1.0
            m /= n_particles

            log_r,log_m = np.log(r),np.log(m)
            del m,r

            position = np.linspace(np.log(guess),np.log(2*guess),n_walkers).reshape((n_walkers,1))
            
            sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, log_likelihood, args=(log_r, log_m),threads=8)
            sampler.run_mcmc(position, 500)
            chain = np.exp(sampler.flatchain)
            
            c_low,c_mid,c_high = np.percentile(chain,[16,50,84])
            del sampler, chain
            
            interval = clock()-t_init
            name = filename.split('/')[-1]
            print(name,c_low,c_mid,c_high,n_particles,interval)

            results.append([name,str(c_low),str(c_mid),str(c_high),str(n_particles),str(interval)])
            del c_low,c_mid,c_high

    finally:
        print('writing results')
        results = [','.join(x) for x in results]
        outfile.write('\n'.join(results))
        outfile.close()

if __name__=='__main__':
    main()
