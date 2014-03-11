import numpy as np, pylab, math, sys, timeit, os, random, datetime, nfw, fit, plotter
from scipy.optimize import fsolve

# Gets the mean-density, scale radius and virial radius for several haloes
# It must be executed with the following command line:
#
#      python main.py directory #x #y #x noplot*
#
# Where:
#      directory: is the path of the directory with the files with the positions of the haloes
#      #x, #y, #z: are the column position of each coordinate in the files (for Multidark is 2 3 4)
#      noplot: is an optional parameter, if it is added the code will not make any graphics
#
# Have fun! 
#
# Done by: Christian Poveda (cn.poveda542@uniandes.edu.co)

totaltime = 0.0
now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")

plt = 1
try:
    if (sys.argv[5]=='noplot'):
        print 'No Plotting mode enabled'
        plt = 0
except:
    plt = 1

os.system('mkdir results_'+now)

sys.stdout.write('\rCompiling the code used to calculate the center of each halo... ')
sys.stdout.flush()
os.system('cc potential.c -lm -o  potential.out')
sys.stdout.write('Done\n')

fit.compile_c_metropolis()

count = 1
for filename in os.listdir('./'+str(sys.argv[1])):

    start = timeit.default_timer()
    print ''
    sys.stdout.write('\rLoading data from "'+str(filename)+'" (file '+str(count)+' of '+str(len(os.listdir('./'+str(sys.argv[1]))))+') ... ')
    count = count + 1
    sys.stdout.flush()

    path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
    data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=16)

    x = data[:,int(sys.argv[2])]
    y = data[:,int(sys.argv[3])]
    z = data[:,int(sys.argv[4])]
    sys.stdout.write('Done\n')
    n_points = len(x)
    
    results_folder ='./results_'+now+'/'+str(filename)

    mass_element = 1.7

    open('positions.dat', "w").write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(n_points)))

    sys.stdout.write('\rRunning c script for potential... ')
    sys.stdout.flush()
    os.system('./potential.out positions.dat')
    sys.stdout.write('Done\n')

    sys.stdout.write('\rCalculating the halo center... ')
    sys.stdout.flush()
    potential = np.loadtxt(open('potential.dat', 'r'))
    maximum = np.argmax(potential)

    x_center,y_center,z_center=x[maximum],y[maximum],z[maximum]
   
    os.system('rm potential.dat positions.dat')
    sys.stdout.write('Done\n')

    x_new,y_new,z_new = x-x_center,y-y_center,z-z_center

    r = np.sort(np.sqrt((x_new)**2 + (y_new)**2 + (z_new)**2), kind='quicksort')
    m = mass_element*np.arange(1,n_points+1,1)

    avg_dens = m/((4.0/3.0)*np.pi*(radius**3))
    r_vir = radius[np.argmin(np.abs(avg_dens-200*(mass_element*(2048.0/1000.0)**3)))]
    
    norm_r = r/r_vir
    norm_m = np.arange(1,n_points+1,1)

    n_iterations = 20000
    logc,logc_walk,chi2 = fit.emcee_sampler(np.log(norm_r),np.log(norm_m),nfw.loglogmass,n_iterations)
    conc = np.exp(logc)
    c_walk = np.exp(c_walk)

    if (plt == 1):

        os.system('mkdir '+results_folder)
        os.chdir(results_folder)
        sys.stdout.write('\rPlotting results... ')
        
        sys.stdout.write('Done\n')
        os.chdir('../')

    else:
        os.chdir('./results_'+now+'/')

    export = open('results_'+now+'.dat', 'a')
    c_max, c_min = np.exp(fit.error_bars(logc_walk,logc,'log'))
    line = [[int(filename.split('_')[1]),x_center,y_center,z_center,conc,c_max,c_min,r_vir]]
    np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf'])

    stop = timeit.default_timer()
    time = stop - start

    print 'Done in ' + str(int(time/60.0))+' min. '+ str(time%60.0)+' sec.'
    os.chdir('../')
    totaltime= totaltime+int(time)
    
print ''
print 'Total time: ' + str(int(totaltime/60.0))+' min. '+ str(int(totaltime%60.0))+' sec. '
os.system('rm *.pyc *.out *.dat')
