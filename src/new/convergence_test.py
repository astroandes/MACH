import argparse # For parsing arguments

parser = argparse.ArgumentParser(__file__, description="A simple MCMC adjuster for concentration in halos")
parser.add_argument("--file", "-f", help="Path to the directory with the halos", type=str)
parser.add_argument("--outputdir","-o",help="Folder output", type=str)
args = parser.parse_args()


import numpy as np, emcee # For numerical computation and MCMC
import os # For handling directories


# Navarro-Frenk-White mass profile with logarithmic scale on mass, radius and
# concentration. It assumes the radius is normalized respect to the virial
# radius.
def nfw_log_mass(log_r,log_c):
    return np.log((np.log(1.0+np.exp(log_c+log_r))-np.exp(log_c+log_r)/(1.0+np.exp(log_c+log_r)))/(np.log(1.0+np.exp(log_c))-np.exp(log_c)/(1.0+np.exp(log_c))))

# Loglikelihood using the chi squared. If the concentration is less than one it
# returns -inf.
def log_likelihood(log_c,log_r,log_m):
    if log_c >= 0:
        s = np.sum(log_m-nfw_log_mass(log_r,log_c))
        return -.5*s*s
    return -np.inf

def fit_data(x, y, z):
    guess = 2
    n_dimensions = 1
    n_walkers = 2
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
        
    sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, log_likelihood, args=(log_r, log_m), threads=4)
    sampler.run_mcmc(position, 500)

    chain = np.exp(sampler.flatchain).flatten()
    chisq = -sampler.flatlnprobability.flatten()
    c_low,c_mid,c_high = np.percentile(chain,[16,50,84])
    return c_low, c_mid, c_high

def main():
    results = []
    filename = args.file
    
    halo_name =  filename.split('/')[-1]
    output_filename = args.outputdir+'/downsample_'+ halo_name
    fout = open(output_filename, 'w')
    print(output_filename)




    data = np.loadtxt(filename,delimiter=',')
    data = np.transpose(data)
    x,y,z = data
    del data


    n_points = len(x)
    print('{} points in total'.format(n_points))


    n_factors = 16
    downsample_factors = np.logspace(1.0, 3.0, n_factors)
    id_list = np.arange(n_points)

    n_selection = 5 

    for i in range(n_factors):
        for j in range(n_selection):
            factor = downsample_factors[i]
            n_select = int(n_points/factor)
            sample_id = np.random.choice(id_list, n_select)

            c_low, c_mid, c_high = fit_data(x[sample_id], y[sample_id], z[sample_id])
            print('{} {} {} {} {}'.format(n_select, factor, c_low, c_mid, c_high))
            fout.write('{} {} {} {} {}\n'.format(n_select, factor, c_low, c_mid, c_high))
        
    fout.close()

if __name__=='__main__':
    main()
