import numpy as np, pylab,sys

data = np.loadtxt(sys.argv[1],delimiter=',')
v = data[:,-3:]
print v
r = np.empty(len(v))

for i in range(len(r)):
    r[i] = np.sum(v[i]*np.ones(3))/(np.sqrt(3*np.sum(v[i]*v[i])))

r -= (3)**(-.5)
r /= 1-(3)**(-.5)

pylab.plot(data[:,8]/data[:,9],r,'.k')
pylab.plot([0,1],[1,1],'--',label='$\mathrm{Bun\~uelo}$')
pylab.plot([0,1],[0.93689,0.93689],'--',label='$\mathrm{Empanada}$')
pylab.plot([0,1],[2/(3*(np.sqrt(3)-1)),2/(3*(np.sqrt(3)-1))],'--',label='$\mathrm{Arepa}$')
pylab.plot([0,1],[(np.sqrt(2)-1)/(np.sqrt(3)-1),(np.sqrt(2)-1)/(np.sqrt(3)-1)],'--',label='$\mathrm{Palito\ de\ queso}$')
pylab.ylabel('$\mathrm{Coeficiente\ de\ bu\~nuelidad}$',fontsize=16)
pylab.xlabel('$n_{vir}/n_{tot}$',fontsize=16)
pylab.xlim([0,1])
pylab.ylim([0,1.1])
pylab.legend(loc=4)
pylab.savefig('plot.png',dpi=200)
pylab.show()

pylab.hist(r,bins=100)
pylab.plot([1,1],[0,len(r)],'--',label='$\mathrm{Bu\~nuelo}$')
pylab.plot([0.93689,0.93689],[0,len(r)],'--',label='$\mathrm{Empanada}$')
pylab.plot([2/(3*(np.sqrt(3)-1)),2/(3*(np.sqrt(3)-1))],[0,len(r)],'--',label='$\mathrm{Arepa}$')
pylab.plot([(np.sqrt(2)-1)/(np.sqrt(3)-1),(np.sqrt(2)-1)/(np.sqrt(3)-1)],[0,len(r)],'--',label='$\mathrm{Palito\ de\ queso}$')
pylab.xlim([0,1.1])
pylab.ylim([0,len(r)])
pylab.xlabel('$\mathrm{Coeficiente\ de\ bu\~nuelidad}$',fontsize=16)
pylab.legend(loc=2)
pylab.savefig('hist.png',dpi=200)
pylab.show()
