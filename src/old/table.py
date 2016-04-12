import numpy as np

data = np.loadtxt('results.csv',delimiter=',')
orig = np.loadtxt('parameters.csv',delimiter=',')

conc = orig[:,0]

n = int(len(conc)/4)
print n
conc = conc[0:n]
c_20 = data[0:n,4]
c_200 = data[n:2*n,4]
c_2000 = data[2*n:3*n,4]
c_20000 = data[3*n:4*n,4]
print len(conc),len(c_20),len(c_200),len(c_2000),len(c_20000)
export = np.dstack((conc,c_20,c_200,c_2000,c_20000))[0]
head = 'c,c_20,c_200,c_2000,c_20000'
np.savetxt('table.csv',export,fmt='%f,%f,%f,%f,%f',header=head)