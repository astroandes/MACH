import pylab, numpy as np, nfw, matplotlib

r = np.arange(0,1,0.001)
C = np.linspace(1,20,6)

for c in C:
    pylab.plot(r,nfw.mass_norm(r,c),color='black',lw= np.log(c+1)/2,label="$c="+str(c)+"$")
pylab.legend(loc=0, scatterpoints=1, prop={'size':20})
filename = 'nfw_normalized'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$r_{norm}$',fontsize=20)
pylab.ylabel('$m_{norm}$',fontsize=20)
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()
