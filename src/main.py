import numpy as np, pylab, math, sys, timeit, os, random, datetime, nfw, mcmc,fit , plotter, multiprocessing, time
from scipy.optimize import fsolve

# Gets the concentration and virial radius for several haloes
# It must be executed with the following command line:
#
#      python main.py directory #x #y #x skip processes noplot*
#
# Where:
#      directory: is the path of the directory with the files with the positions of the haloes
#      #x, #y, #z: are the column position of each coordinate in the files (for Multidark is 2 3 4)
#      skip: number of rows to skip (For Multidark is 16) 
#      processes: number of child processes
#      noplot: is an optional parameter, if it is added the code will not make any graphics
#
# Have fun! 
#
# Done by: Christian Poveda (cn.poveda542@uniandes.edu.co)
#
# Thanks to Diva Martinez (dm.martinez831@uniandes.edu.co) for the multiprocessing idea

processes = int(sys.argv[6])
now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
mass_element = 1.0   
plt = 1
total_list = os.listdir('./'+str(sys.argv[1]))
len_list = len(total_list)
lists = [total_list[i*len_list // processes: (i+1)*len_list //processes] for i in range(processes)]
jobs = []

try:
    if (sys.argv[7]=='noplot'):
        print 'No Plotting mode enabled'
        plt = 0
except:
    plt = 1
 
os.system('mkdir results_'+now)
sys.stdout.write('\rCompiling the code used to calculate the center of each halo... ')
sys.stdout.flush()
os.system('cc potential.c -lm -o  potential.out')
sys.stdout.write('Done\n')

def run(directories,process_number):

    process_number+=1
    count = 0
    filename_export = './results_'+now+'/results_'+str(process_number)+'.csv'
    os.system('touch '+filename_export)

    for filename in directories:

        count = count + 1
        print '\rWorking with file '+str(count)+' of '+str(len(directories))+' in process '+str(process_number)

    # Creates the "data" array with the information from the file in "path"
        path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
        data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=int(sys.argv[5]))

    # Gets the cartesian coordinates for each particle in the halo and the number of particles
        x = data[:,int(sys.argv[2])]
        y = data[:,int(sys.argv[3])]
        z = data[:,int(sys.argv[4])]
        n_points = len(x)

    # Exports a file with the cartesian coordinates of each particle
        file_id = int(filename.split('_')[1])
        positions_name ='positions_'+str(file_id)+'.dat'
        potential_name ='potential_'+str(file_id)+'.dat'
        results_folder ='./results_'+now+'/'+str(file_id)
        open(positions_name, "w").write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(n_points)))

    # Runs the executable "potential.out" that will get the potential energy for each particle
        os.system('./potential.out '+positions_name+' '+potential_name)

    # Finds the particle with the lowest potential and puts it as the new origin
        potential = np.loadtxt(open(potential_name, 'r'))
        maximum = np.argmax(potential)
        x_center,y_center,z_center = x[maximum],y[maximum],z[maximum]

    # Gets the radial distance from the new origin to each particle and sorts it
        r_values = np.sort(np.sqrt((x-x_center)**2 + (y-y_center)**2 + (z-z_center)**2), kind='quicksort')
    
    # Removes the files used by "potential.out"
        os.system('rm '+potential_name+' '+positions_name)

    # Gets the mass for each radius and removes the particle in the origin
        mass = mass_element*np.arange(1,n_points,1)
        radius = np.delete(r_values,0)

    # Gets the average density in function of the radius
        avg_density = mass/((4.0/3.0)*np.pi*(radius**3))
        rho_back = mass_element*(2048.0/1000.0)**3

    # Gets the virial radius
        crit = 740.0
        vir_index = np.argmin(np.abs(avg_density-crit*rho_back))

        if np.argmin(np.abs(avg_density-crit*rho_back)) > 1:
            r_vir = radius[vir_index]
        else:
            r_vir = radius[-1]
            vir_index = len(avg_density)-1

    # Removes the particles that have a greater radius than the virial radius
        mass = np.resize(mass,vir_index+1)
        radius = np.resize(radius,vir_index+1)

    # Normalizes the mass and radius
        mass = mass/mass[-1]
        radius = radius/radius[-1]

    # Does the Metropolis algorithm in order to find the concentration
        
        n_iterations = 50000

        step = np.array([np.log(1.03)])
        guess = np.array([np.log(10)])
        reest = [lambda x: x >= 0]

        chi_2 = lambda p : mcmc.chi2(p,nfw.loglogmass_norm,np.log(radius),np.log(mass),np.ones(len(radius)))    
        walk,chi2 = mcmc.walker(chi_2,guess,step,n_iterations,reest)
        walk = walk[0]
        log_c = walk[np.argmin(chi2)]
        c = np.exp(log_c)
        c_max, c_min = np.exp(fit.error_bars(walk,log_c,'log'))

    # Generates plots
        if (plt == 1):

            os.system('mkdir '+results_folder)
            os.chdir(results_folder)

            plotter.halo(x,y,z,x_center,y_center,z_center,r_vir)
            plotter.mass_norm(radius,mass,c,c_max,c_min)
            pylab.scatter(np.exp(walk),chi2)
            pylab.xlabel('$c$')
            pylab.ylabel('$\chi ^2$')
            pylab.savefig('chi2.pdf',dpi=300)
            pylab.close()

            pylab.scatter(np.exp(walk),np.exp(-chi2/2.0),c='r',label='BDMW')
            pylab.xlabel('$c$')
            pylab.ylabel('$\cal{L}$')
            pylab.savefig('likelihood.pdf',dpi=300)
            pylab.close()

            pylab.hist(np.exp(walk),1000,normed=True)
            pylab.xlabel('c')
            pylab.ylabel('P(c)')
            pylab.savefig('histogram.pdf',dpi=300)
            pylab.close()

            os.chdir('../../')
            
    # Writes the results for the halo
        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c,c_max,c_min,r_vir,len(mass),n_points]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%d','%d'],delimiter=',')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()

