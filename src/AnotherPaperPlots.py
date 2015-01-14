import pylab, numpy as np, nfw, matplotlib

r = np.arange(0,1,0.001)
C = np.linspace(1,20,4)

for c in C:
    c = int(c)
    pylab.plot(r,nfw.mass_norm(r,c),color='black',lw = .2*c+1,label="$c="+str(c)+"$")
pylab.legend(loc=0, scatterpoints=1, prop={'size':20})
filename = 'nfw_normalized'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$x$',fontsize=20)
pylab.ylabel('$m(<x)$',fontsize=20)
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()


for c in C:
    c = int(c)
    pylab.plot(r,np.sqrt(nfw.mass_norm(r,c)/r),color='black',lw = .2*c+1,label="$c="+str(c)+"$")
pylab.legend(loc=0, scatterpoints=1, prop={'size':20})
filename = 'vel_normalized'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$x$',fontsize=20)
pylab.ylabel('$v(<x)$',fontsize=20)
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()
