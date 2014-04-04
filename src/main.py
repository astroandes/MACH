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

export = open('./results_'+now+'/results_'+now+'.csv', 'w')

export.write('#Id,x_center,y_center,z_center,c_bdmv,c_bdmv_max,c_bdmv_min,c_bdmw,c_bdmw_max,c_bdmw_min,r_vir_bdmv,r_vir_bdmw,c\n')
export.close()

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

    r_values = np.sort(np.sqrt((x-x_center)**2 + (y-y_center)**2 + (z-z_center)**2), kind='quicksort')

    mass = mass_element*np.arange(1,n_points,1)
    radius = np.delete(r_values,0)

    avg_density = mass/((4.0/3.0)*np.pi*(radius**3))
    rho_back = mass_element*(2048.0/1000.0)**3

    bdmv = 360.0
    bdmw = 740.0

    bdmv_index = np.argmin(np.abs(avg_density-bdmv*rho_back))
    bdmw_index = np.argmin(np.abs(avg_density-bdmw*rho_back))

    if np.argmin(np.abs(avg_density-bdmv*rho_back)) != len(avg_density)-1:
        r_bdmv = radius[bdmv_index]
    else:
        r_bdmv = radius[-1]
        bdmv_index = -1

    if np.argmin(np.abs(avg_density-bdmw*rho_back)) != len(avg_density)-1:
        r_bdmw = radius[bdmw_index]
    else:
        r_bdmw = radius[-1]
        bdmw_index = -1
    
    bdmv_mass = np.split(mass,bdmv_index)[0]
    bdmv_radius = np.split(radius,bdmv_index)[0]

    bdmv_mass = bdmv_mass/bdmv_mass[-1]
    bdmv_radius = bdmv_radius/bdmv_radius[-1]
    
    bdmw_mass = np.split(mass,bdmw_index)[0]
    bdmw_radius = np.split(radius,bdmw_index)[0]

    bdmw_mass = bdmw_mass/bdmw_mass[-1]
    bdmw_radius = bdmw_radius/bdmw_radius[-1]

    pylab.scatter(bdmv_radius,bdmv_mass,c='b')
    pylab.scatter(bdmw_radius,bdmw_mass,c='r')
    pylab.show()

    n_iterations = 60000
    log_bdmv,bdmv_walk,bdmv_chi2 = fit.metropolis_one(np.log(bdmv_radius),np.log(bdmv_mass),nfw.loglogmass_norm,n_iterations)
    log_bdmw,bdmw_walk,bdmw_chi2 = fit.metropolis_one(np.log(bdmw_radius),np.log(bdmw_mass),nfw.loglogmass_norm,n_iterations)
    
    c_bdmv = np.exp(log_bdmv)
    c_bdmw = np.exp(log_bdmw)

    bdmv_max, bdmv_min = np.exp(fit.error_bars(bdmv_walk,c_bdmv,'log'))
    bdmw_max, bdmw_min = np.exp(fit.error_bars(bdmw_walk,c_bdmw,'log'))
    
    if (plt == 1):

        os.system('mkdir '+results_folder)
        os.chdir(results_folder)
        plotter.halo(x,y,z,x_center,y_center,z_center)
        sys.stdout.write('\rPlotting results... ')
        plotter.mass_norm(bdmv_radius,bdmv_mass,c_bdmv,bdmv_max,bdmv_min,'BDMV')
        plotter.mass_norm(bdmw_radius,bdmw_mass,c_bdmw,bdmw_max,bdmw_min,'BDMW')
        sys.stdout.flush()
        sys.stdout.write('Done\n')
        os.chdir('../')

    else:
        os.chdir('./results_'+now+'/')

    export = open('results_'+now+'.csv', 'a')
    line = [[int(filename.split('_')[1]),x_center,y_center,z_center,c_bdmv,bdmv_max,bdmv_min,c_bdmw,bdmw_max,bdmw_min,r_bdmv,r_bdmw]]
    np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf'],delimiter=',')
    
    stop = timeit.default_timer()
    time = stop - start

    print 'Done in ' + str(int(time/60.0))+' min. '+ str(time%60.0)+' sec.'
    os.chdir('../')
    totaltime= totaltime+int(time)
    
print ''
print 'Total time: ' + str(int(totaltime/60.0))+' min. '+ str(int(totaltime%60.0))+' sec. '
os.system('rm *.pyc *.out *.dat')
