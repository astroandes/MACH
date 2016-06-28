"""
Simple script to select halos to downsample
"""
import numpy as np

data = np.loadtxt("../data/results/Bolshoi/density/results.dat")
mass = data[:,4]*1E10
xoff = data[:,5]
ID = data[:,0]

ii = (mass > 4E13) & (xoff<0.1)
for i in ID[ii]:
    haloid = "%06d"%(int(i))
    output_dir = "/srv/cosmusr/scratch/jforero/Bolshoi/downsample/{}/".format(haloid)
    input_dir = "/srv/cosmusr/scratch/jforero/Bolshoi/particle_subset/code/CenteredHalos/largest/"
    input_file = "{}/halo_{}.dat".format(input_dir, haloid)
    command_a = "mkdir -p {}".format(output_dir)
    command_b = "python downsample_halo.py -f {} -o {}".format(input_file, output_dir)
    print(command_a)
    print(command_b)
