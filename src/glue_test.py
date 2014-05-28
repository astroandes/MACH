import sys, os, numpy as np

directory = str(sys.argv[1])
total_list = os.listdir(directory)
filenames = sorted([directory+'/'+element for element in total_list if 'result' in element])


export = open(directory+'/results.csv', 'w')
head='Id,x_center,y_center,z_center,conc,c_max,c_min'
for fname in filenames:
    with open(fname) as infile:
        for line in infile:
            export.write(line)
export.close()

data = np.loadtxt(directory+'/results.csv',delimiter=',')
data = data[data[:,0].argsort()]

np.savetxt(directory+'/results.csv',data,fmt='%d,%f,%f,%f,%f,%f,%f',header=head)
