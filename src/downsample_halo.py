import argparse # For parsing arguments

parser = argparse.ArgumentParser(__file__, description="A simple script to downsample halo particles")
parser.add_argument("--file", "-f", help="Path to the directory with the halos", type=str)
parser.add_argument("--outputdir","-o",help="Folder output", type=str)
parser.add_argument("--csv", "-c", help="Is input in csv format?")
args = parser.parse_args()

import numpy as np # For numerical computation and MCMC
import os # For handling directories
from scipy.optimize import fsolve


# Defines the relation between the concentration and the greatest velocity

def main():
    results = []
    filename = args.file
    
    halo_name =  filename.split('/')[-1]

    if args.csv:
        data = np.loadtxt(filename,delimiter=',')
    else:
        data = np.loadtxt(filename)

    data = np.transpose(data)
    x,y,z = data
    del data


    n_points = len(x)
    print('{} points in total'.format(n_points))
    max_sample_factor = np.log10(n_points/200.0)
    min_sample_factor = np.log10(2.0)
    print('(sampling) max {} min {}'.format(10**max_sample_factor, 10**min_sample_factor),)

    n_factors = 10
    downsample_factors = np.logspace(0.3, 2.0, n_factors)[::-1]
    print('{}'.format(downsample_factors))

    print(downsample_factors)
    id_list = np.arange(n_points)
    for i in range(5):
        id_list = np.random.permutation(id_list)

    n_iterations = 10000
    i_sample = 0
    for i in range(n_iterations):
        log_factor = (max_sample_factor - min_sample_factor) * np.random.random() + min_sample_factor 
        factor = 10**log_factor
        n_select = int(n_points/factor)
        sample_id = np.random.choice(id_list, n_select)

        output_filename = args.outputdir+'/sample_{:05d}.dat'.format(i_sample)
        print(output_filename, factor, n_select)
        data = x[sample_id], y[sample_id], z[sample_id]
        data = np.transpose(data)
        np.savetxt(output_filename, data, fmt='%f')            
        i_sample += 1
            

if __name__=='__main__':
    main()
