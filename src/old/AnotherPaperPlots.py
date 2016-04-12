import pylab, numpy as np, nfw, matplotlib

r = np.linspace(0,1,1000)
C = [1,2,6,12,20]

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
    pylab.plot(r,nfw.velocity_norm(r,c),color='black',lw = .2*c+1,label="$c="+str(c)+"$")
pylab.legend(loc=0, scatterpoints=1, prop={'size':20})
filename = 'vel_normalized'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$x$',fontsize=20)
pylab.ylabel('$v(<x)$',fontsize=20)
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()

for c in C:
    c = int(c)
    dv_dx = (nfw.velocity_norm(np.delete(r,0),c)-nfw.velocity_norm(np.delete(r,-1),c))/(r[1]-r[0])
    pylab.plot(np.delete(r,0),dv_dx,color='black',lw = .2*c+1,label="$c="+str(c)+"$")
pylab.legend(loc=0, scatterpoints=1, prop={'size':20})
filename = 'dv_dx'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$x$',fontsize=20)
pylab.ylabel('$dv/dx$',fontsize=20)
pylab.yscale('log')
pylab.ylim(ymin=1E-2)
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()

def f(x,c):
    return ((2*c*x+1)*x*c/(c*x+1)**2 - np.log(c*x+1))
c = np.linspace(1,20,1000)
X,Y = pylab.meshgrid(r,c)

pylab.imshow(f(X,Y),cmap='spectral',origin='lower',extent=[1,20,0,1],aspect=20)
filename = 'zeros'
matplotlib.rcParams.update({'font.size': 18})
pylab.xlabel('$c$',fontsize=20)
pylab.ylabel('$x$',fontsize=20)
pylab.title('$2Ax^{2}v (dv/dx)$',fontsize=20)
pylab.colorbar()
pylab.savefig(filename+'.eps',transparent=True)
pylab.savefig(filename+'.pdf',transparent=True)
pylab.savefig(filename+'.png',transparent=True)
pylab.close()
