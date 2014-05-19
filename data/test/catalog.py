import numpy as np, pylab

expected = np.loadtxt('parameters.csv',delimiter = ',')
obtained = np.loadtxt('results.csv',delimiter = ',')

c_exp = expected[:,0]
r_vir = expected[:,1]
m_vir = expected[:,2]
c_obt = obtained[:,4]
c_max = obtained[:,5]
c_min = obtained[:,6]
count = 0

for i in range(len(c_obt)):
    if c_exp[i] > c_max[i] or c_exp[i] < c_min[i]:
        count += 1

print count

pylab.plot(c_exp,c_max,'.r')
pylab.plot(c_exp,c_min,'.g')
pylab.plot(c_exp,c_obt,'.b')
pylab.plot(np.sort(c_exp),np.sort(c_exp),'--k',lw=2)
pylab.xscale('log')
pylab.yscale('log')
pylab.title('Concentration')
pylab.xlabel('Expected')
pylab.ylabel('Obtained')
pylab.savefig('concentration.png',dpi=300)
pylab.close()

pylab.plot(m_vir,np.abs(c_exp - c_obt),'.k')
pylab.xlabel('Number of Particles')
pylab.ylabel('Relative Error')
pylab.savefig('relative_error.png',dpi=300)
pylab.close()


pylab.plot(m_vir,100.0*np.abs(c_exp - c_obt)/c_exp,'.k')
pylab.xlabel('Number of Particles')
pylab.ylabel('Porcentual Error')
pylab.savefig('porcentual_error.png',dpi=300)
pylab.close()
