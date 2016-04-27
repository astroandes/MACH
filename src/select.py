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
    print(int(i))
