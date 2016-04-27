import argparse # For parsing arguments

parser = argparse.ArgumentParser(__file__, description="A simple script to downsample halo particles")
parser.add_argument("--file", "-f", help="Path to the directory with the halos", type=str)
parser.add_argument("--outputdir","-o",help="Folder output", type=str)
args = parser.parse_args()

import numpy as np # For numerical computation and MCMC
import os # For handling directories
from scipy.optimize import fsolve

# Defines the relation between the concentration and the greatest velocity


def f(c,v_max):
    return v_max - np.sqrt(0.216*c/(np.log(1+c)-c/(1+c)))

def fit_data_vel(x,y,z):
    r = np.sqrt(x*x+y*y+z*z)
    r = np.sort(r)
    r = np.delete(r,0)
    r = r/r[-1]
    del x,y,z
    
    n_particles = len(r)
    
    m = np.arange(n_particles)+1.0
    m /= n_particles
    
    v = np.sqrt(m/r)
    v_max = np.max(v)
    c = fsolve(f,10.0,args=(v_max))[0]
    return c

def main():
    results = []
    filename = args.file
    
    halo_name =  filename.split('/')[-1]

    data = np.loadtxt(filename,delimiter=',')
    data = np.transpose(data)
    x,y,z = data
    del data


    n_points = len(x)
    print('{} points in total'.format(n_points))


    n_factors = 7
    downsample_factors = np.logspace(1.0, 3.0, n_factors)[::-1]


    print(downsample_factors)
    id_list = np.arange(n_points)

    n_iterations = 100
    i_sample = 0
    for i in range(n_factors):
        for j in range(n_iterations):
            factor = downsample_factors[i]
            n_select = int(n_points/factor)
            sample_id = np.random.choice(id_list, n_select)

            output_filename = args.outputdir+'/downsample_{}_'.format(i_sample)+ halo_name    
            print(output_filename)
            data = x[sample_id], y[sample_id], z[sample_id]
            data = np.transpose(data)
            np.savetxt(output_filename, data, fmt='%f')            
            i_sample += 1
            

if __name__=='__main__':
    main()
