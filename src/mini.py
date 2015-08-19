import numpy as np, matplotlib.pyplot as plt, emcee
import sys, os

def nfw_log_mass(log_r,c):
    return np.log((np.log(1.0+c*np.exp(log_r))-c*np.exp(log_r)/(1.0+c*np.exp(log_r)))/(np.log(1.0+c)-c/(1.0+c)))

def log_likelihood(c,log_r,log_m):
    if c >= 1:
        s = np.sum(log_m-nfw_log_mass(log_r,c))
        return -.5*s*s
    else:
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

            position = np.linspace(guess,2*guess,n_walkers)
            position = [[x] for x in position]

            sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, log_likelihood, args=(log_r, log_m),threads=4)
            sampler.run_mcmc(position, 500)
            chain = sampler.flatchain
            c_low,c_mid,c_high = np.percentile(chain,[16,50,84])
            del sampler, chain

            print(filename,str(c_low),str(c_mid),str(c_high))
            results.append([filename,str(c_low),str(c_mid),str(c_high)])
            del c_low,c_mid,c_high

    finally:
        print('writing results')
        results = [','.join(x) for x in results]
        outfile.write('\n'.join(results))
        outfile.close()

if __name__=='__main__':
    main()
