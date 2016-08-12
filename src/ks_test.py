import numpy as np
from scipy import stats

def make_ks_test(halo_id=0):
    n_samples = 700
    fileout = open('ks_test_p_values_{:06d}.dat'.format(halo_id), 'w')
    fileout.write("#id_i id_j n_part_i n_part_j ks_test_value pvalue\n")
    print(fileout)

    filename = '/scratch/jforero/Bolshoi/allhalos/largest/halo_{:06}.dat'.format(halo_id)
    data = np.loadtxt(filename, delimiter=',')
    r = np.sqrt(data[:,0]**2 + data[:,1]**2+ data[:,2]**2)
    for j in range(n_samples):
        print halo_id, j
        filename_j = '//scratch/jforero/Bolshoi/jcdownsample/halos/{:06d}/sample_{:03d}.dat'.format(halo_id,j)
        data_j = np.loadtxt(filename_j)
        r_j = np.sqrt(data_j[:,0]**2 + data_j[:,1]**2+ data_j[:,2]**2)

        ks_value, p_value = stats.ks_2samp(r, r_j)
        fileout.write("{} {} {} {}\n".format(j, len(r_j), ks_value, p_value))

    fileout.close()



halo_id_list = [3438, 4497, 14976, 19876, 26975, 30491, 36320, 55806, 59095, 59875, 60662, 69962, 73222, 78002]

for halo_id in halo_id_list:
    make_ks_test(halo_id=halo_id)
