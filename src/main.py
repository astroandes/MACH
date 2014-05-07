import numpy as np, pylab, math, sys, timeit, os, random, datetime, nfw, fit, plotter
from scipy.optimize import fsolve
import multiprocessing

# Gets the concentration and virial radius for several haloes
# It must be executed with the following command line:
#
#      python main.py directory #x #y #x skip noplot*
#
# Where:
#      directory: is the path of the directory with the files with the positions of the haloes
#      #x, #y, #z: are the column position of each coordinate in the files (for Multidark is 2 3 4)
#      skip: number of rows to skip (For Multidark is 16) 
#      noplot: is an optional parameter, if it is added the code will not make any graphics
#
# Have fun! 
#
# Done by: Christian Poveda (cn.poveda542@uniandes.edu.co)
#
# Thanks to Diva Martinez (dm.martinez831@uniandes.edu.co) for the multiprocessing idea

processes = 4
now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
mass_element = 1.0   
plt = 1
total_list = os.listdir('./'+str(sys.argv[1]))
len_list = len(total_list)
lists = [total_list[i*len_list // processes: (i+1)*len_list //processes] for i in range(processes)]
jobs = []

try:
    if (sys.argv[6]=='noplot'):
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
        sys.stdout.write('\rWorking with file '+str(count)+' of '+str(len(directories))+' in process '+str(process_number)+'... ')
        sys.stdout.flush()

        path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
        data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=int(sys.argv[5]))

        x = data[:,int(sys.argv[2])]
        y = data[:,int(sys.argv[3])]
        z = data[:,int(sys.argv[4])]

        n_points = len(x)
        results_folder ='./results_'+now+'/'+str(filename)
        
        file_id = int(filename.split('_')[1])
        positions_name ='positions_'+str(file_id)+'.dat'
        potential_name ='potential_'+str(file_id)+'.dat'
        open(positions_name, "w").write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(n_points)))

        os.system('./potential.out '+positions_name+' '+potential_name)

        potential = np.loadtxt(open(potential_name, 'r'))
        maximum = np.argmax(potential)

        x_center,y_center,z_center = x[maximum],y[maximum],z[maximum]
        
        os.system('rm '+potential_name+' '+positions_name)
        
        r_values = np.sort(np.sqrt((x-x_center)**2 + (y-y_center)**2 + (z-z_center)**2), kind='quicksort')
        
        mass = mass_element*np.arange(1,n_points,1)
        radius = np.delete(r_values,0)

        avg_density = mass/((4.0/3.0)*np.pi*(radius**3))
        rho_back = mass_element*(2048.0/1000.0)**3

        bdmv = 360.0
        bdmw = 740.0

        bdmv_index = np.argmin(np.abs(avg_density-bdmv*rho_back))
        bdmw_index = np.argmin(np.abs(avg_density-bdmw*rho_back))

        if np.argmin(np.abs(avg_density-bdmv*rho_back)) > 1:
            r_bdmv = radius[bdmv_index]
        else:
            r_bdmv = radius[-1]
            bdmv_index = len(avg_density)-1

        if np.argmin(np.abs(avg_density-bdmw*rho_back)) > 1:
            r_bdmw = radius[bdmw_index]
        else:
            r_bdmw = radius[-1]
            bdmw_index = len(avg_density)-1

        bdmv_mass = np.resize(mass,bdmv_index)
        bdmv_radius = np.resize(radius,bdmv_index)
        
        bdmv_mass = bdmv_mass/bdmv_mass[-1]
        bdmv_radius = bdmv_radius/bdmv_radius[-1]
        
        bdmw_mass = np.resize(mass,bdmw_index)
        bdmw_radius = np.resize(radius,bdmw_index)

        bdmw_mass = bdmw_mass/bdmw_mass[-1]
        bdmw_radius = bdmw_radius/bdmw_radius[-1]

        n_iterations = 25000
        log_bdmv,bdmv_walk,bdmv_chi2 = fit.metropolis_one(np.log(bdmv_radius),np.log(bdmv_mass),nfw.loglogmass_norm,n_iterations)
        log_bdmw,bdmw_walk,bdmw_chi2 = fit.metropolis_one(np.log(bdmw_radius),np.log(bdmw_mass),nfw.loglogmass_norm,n_iterations)

        c_bdmv = np.exp(log_bdmv)
        c_bdmw = np.exp(log_bdmw)
    
        bdmv_max, bdmv_min = np.exp(fit.error_bars(bdmv_walk,log_bdmv,'log'))
        bdmw_max, bdmw_min = np.exp(fit.error_bars(bdmw_walk,log_bdmw,'log'))
    
        if (plt == 1):

            os.system('mkdir '+results_folder)
            os.chdir(results_folder)

            plotter.halo(x,y,z,x_center,y_center,z_center,r_bdmv,r_bdmw)
            plotter.mass_norm(bdmv_radius,bdmv_mass,c_bdmv,bdmv_max,bdmv_min,'BDMV')
            plotter.mass_norm(bdmw_radius,bdmw_mass,c_bdmw,bdmw_max,bdmw_min,'BDMW')

            pylab.scatter(np.exp(bdmv_walk),bdmv_chi2,label='BDMV')
            pylab.scatter(np.exp(bdmw_walk),bdmw_chi2,c='r',label='BDMW')
            pylab.legend(loc=4, borderaxespad=0.5)
            pylab.xlabel('$c$')
            pylab.ylabel('$\chi ^2$')
            pylab.savefig('chi2.png',dpi=300)
            pylab.close()

            os.chdir('../../')

        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c_bdmv,bdmv_max,bdmv_min,c_bdmw,bdmw_max,bdmw_min,r_bdmv,r_bdmw,len(bdmv_mass),len(bdmw_mass),n_points]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%d','%d','%d'],delimiter=',')

        sys.stdout.write('Done\n')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()


for p in jobs: 
    while(p.is_alive()):
        pass

export = open('./results_'+now+'/results_'+now+'.csv', 'w')

filenames = ['./results_'+now+'/results_'+str(i)+'.csv' for i in range(processes)]
with export as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)    


sys.stdout.write('Done\n')

def run(directories,process_number):
    
    count = 0
    filename_export = './results_'+now+'/results_'+str(process_number)+'.csv'
    os.system('touch '+filename_export)

    for filename in directories:

        count = count + 1
        sys.stdout.write('\rWorking with file '+str(count)+' of '+str(len(directories))+' in process '+str(process_number)+'... ')
        sys.stdout.flush()

        path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
        data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=int(sys.argv[5]))

        x = data[:,int(sys.argv[2])]
        y = data[:,int(sys.argv[3])]
        z = data[:,int(sys.argv[4])]

        n_points = len(x)
        results_folder ='./results_'+now+'/'+str(filename)
        
        file_id = int(filename.split('_')[1])
        positions_name ='positions_'+str(file_id)+'.dat'
        potential_name ='potential_'+str(file_id)+'.dat'
        open(positions_name, "w").write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(n_points)))

        os.system('./potential.out '+positions_name+' '+potential_name)

        potential = np.loadtxt(open(potential_name, 'r'))
        maximum = np.argmax(potential)

        x_center,y_center,z_center = x[maximum],y[maximum],z[maximum]
        
        os.system('rm '+potential_name+' '+positions_name)
        
        r_values = np.sort(np.sqrt((x-x_center)**2 + (y-y_center)**2 + (z-z_center)**2), kind='quicksort')
        
        mass = mass_element*np.arange(1,n_points,1)
        radius = np.delete(r_values,0)

        avg_density = mass/((4.0/3.0)*np.pi*(radius**3))
        rho_back = mass_element*(2048.0/1000.0)**3

        bdmv = 360.0
        bdmw = 740.0

        bdmv_index = np.argmin(np.abs(avg_density-bdmv*rho_back))
        bdmw_index = np.argmin(np.abs(avg_density-bdmw*rho_back))

        if np.argmin(np.abs(avg_density-bdmv*rho_back)) > 1:
            r_bdmv = radius[bdmv_index]
        else:
            r_bdmv = radius[-1]
            bdmv_index = len(avg_density)-1

        if np.argmin(np.abs(avg_density-bdmw*rho_back)) > 1:
            r_bdmw = radius[bdmw_index]
        else:
            r_bdmw = radius[-1]
            bdmw_index = len(avg_density)-1

        bdmv_mass = np.resize(mass,bdmv_index)
        bdmv_radius = np.resize(radius,bdmv_index)
        
        bdmv_mass = bdmv_mass/bdmv_mass[-1]
        bdmv_radius = bdmv_radius/bdmv_radius[-1]
        
        bdmw_mass = np.resize(mass,bdmw_index)
        bdmw_radius = np.resize(radius,bdmw_index)

        bdmw_mass = bdmw_mass/bdmw_mass[-1]
        bdmw_radius = bdmw_radius/bdmw_radius[-1]

        n_iterations = 25000
        log_bdmv,bdmv_walk,bdmv_chi2 = fit.metropolis_one(np.log(bdmv_radius),np.log(bdmv_mass),nfw.loglogmass_norm,n_iterations)
        log_bdmw,bdmw_walk,bdmw_chi2 = fit.metropolis_one(np.log(bdmw_radius),np.log(bdmw_mass),nfw.loglogmass_norm,n_iterations)

        c_bdmv = np.exp(log_bdmv)
        c_bdmw = np.exp(log_bdmw)
    
        bdmv_max, bdmv_min = np.exp(fit.error_bars(bdmv_walk,log_bdmv,'log'))
        bdmw_max, bdmw_min = np.exp(fit.error_bars(bdmw_walk,log_bdmw,'log'))
    
        if (plt == 1):

            os.system('mkdir '+results_folder)
            os.chdir(results_folder)

            plotter.halo(x,y,z,x_center,y_center,z_center,r_bdmv,r_bdmw)
            plotter.mass_norm(bdmv_radius,bdmv_mass,c_bdmv,bdmv_max,bdmv_min,'BDMV')
            plotter.mass_norm(bdmw_radius,bdmw_mass,c_bdmw,bdmw_max,bdmw_min,'BDMW')

            pylab.scatter(np.exp(bdmv_walk),bdmv_chi2,label='BDMV')
            pylab.scatter(np.exp(bdmw_walk),bdmw_chi2,c='r',label='BDMW')
            pylab.legend(loc=4, borderaxespad=0.5)
            pylab.xlabel('$c$')
            pylab.ylabel('$\chi ^2$')
            pylab.savefig('chi2.png',dpi=300)
            pylab.close()

            os.chdir('../../')

        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c_bdmv,bdmv_max,bdmv_min,c_bdmw,bdmw_max,bdmw_min,r_bdmv,r_bdmw,len(bdmv_mass),len(bdmw_mass),n_points]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%d','%d','%d'],delimiter=',')

        sys.stdout.write('Done\n')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()


for p in jobs: 
    while(p.is_alive()):
        pass

export = open('./results_'+now+'/results_'+now+'.csv', 'w')

filenames = ['./results_'+now+'/results_'+str(i)+'.csv' for i in range(processes)]
with export as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)    

sys.stdout.write('Done\n')

def run(directories,process_number):
    
    count = 0
    filename_export = './results_'+now+'/results_'+str(process_number)+'.csv'
    os.system('touch '+filename_export)

    for filename in directories:

        count = count + 1
        sys.stdout.write('\rWorking with file '+str(count)+' of '+str(len(directories))+' in process '+str(process_number)+'... ')
        sys.stdout.flush()

        path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
        data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=int(sys.argv[5]))

        x = data[:,int(sys.argv[2])]
        y = data[:,int(sys.argv[3])]
        z = data[:,int(sys.argv[4])]

        n_points = len(x)
        results_folder ='./results_'+now+'/'+str(filename)
        
        file_id = int(filename.split('_')[1])
        positions_name ='positions_'+str(file_id)+'.dat'
        potential_name ='potential_'+str(file_id)+'.dat'
        open(positions_name, "w").write('\n'.join('%lf,%lf,%lf' % (x[i],y[i],z[i]) for i in range(n_points)))

        os.system('./potential.out '+positions_name+' '+potential_name)

        potential = np.loadtxt(open(potential_name, 'r'))
        maximum = np.argmax(potential)

        x_center,y_center,z_center = x[maximum],y[maximum],z[maximum]
        
        os.system('rm '+potential_name+' '+positions_name)
        
        r_values = np.sort(np.sqrt((x-x_center)**2 + (y-y_center)**2 + (z-z_center)**2), kind='quicksort')
        
        mass = mass_element*np.arange(1,n_points,1)
        radius = np.delete(r_values,0)

        avg_density = mass/((4.0/3.0)*np.pi*(radius**3))
        rho_back = mass_element*(2048.0/1000.0)**3

        bdmv = 360.0
        bdmw = 740.0

        bdmv_index = np.argmin(np.abs(avg_density-bdmv*rho_back))
        bdmw_index = np.argmin(np.abs(avg_density-bdmw*rho_back))

        if np.argmin(np.abs(avg_density-bdmv*rho_back)) > 1:
            r_bdmv = radius[bdmv_index]
        else:
            r_bdmv = radius[-1]
            bdmv_index = len(avg_density)-1

        if np.argmin(np.abs(avg_density-bdmw*rho_back)) > 1:
            r_bdmw = radius[bdmw_index]
        else:
            r_bdmw = radius[-1]
            bdmw_index = len(avg_density)-1

        bdmv_mass = np.resize(mass,bdmv_index)
        bdmv_radius = np.resize(radius,bdmv_index)
        
        bdmv_mass = bdmv_mass/bdmv_mass[-1]
        bdmv_radius = bdmv_radius/bdmv_radius[-1]
        
        bdmw_mass = np.resize(mass,bdmw_index)
        bdmw_radius = np.resize(radius,bdmw_index)

        bdmw_mass = bdmw_mass/bdmw_mass[-1]
        bdmw_radius = bdmw_radius/bdmw_radius[-1]

        n_iterations = 25000
        log_bdmv,bdmv_walk,bdmv_chi2 = fit.metropolis_one(np.log(bdmv_radius),np.log(bdmv_mass),nfw.loglogmass_norm,n_iterations)
        log_bdmw,bdmw_walk,bdmw_chi2 = fit.metropolis_one(np.log(bdmw_radius),np.log(bdmw_mass),nfw.loglogmass_norm,n_iterations)

        c_bdmv = np.exp(log_bdmv)
        c_bdmw = np.exp(log_bdmw)
    
        bdmv_max, bdmv_min = np.exp(fit.error_bars(bdmv_walk,log_bdmv,'log'))
        bdmw_max, bdmw_min = np.exp(fit.error_bars(bdmw_walk,log_bdmw,'log'))
    
        if (plt == 1):

            os.system('mkdir '+results_folder)
            os.chdir(results_folder)

            plotter.halo(x,y,z,x_center,y_center,z_center,r_bdmv,r_bdmw)
            plotter.mass_norm(bdmv_radius,bdmv_mass,c_bdmv,bdmv_max,bdmv_min,'BDMV')
            plotter.mass_norm(bdmw_radius,bdmw_mass,c_bdmw,bdmw_max,bdmw_min,'BDMW')

            pylab.scatter(np.exp(bdmv_walk),bdmv_chi2,label='BDMV')
            pylab.scatter(np.exp(bdmw_walk),bdmw_chi2,c='r',label='BDMW')
            pylab.legend(loc=4, borderaxespad=0.5)
            pylab.xlabel('$c$')
            pylab.ylabel('$\chi ^2$')
            pylab.savefig('chi2.png',dpi=300)
            pylab.close()

            os.chdir('../../')

        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c_bdmv,bdmv_max,bdmv_min,c_bdmw,bdmw_max,bdmw_min,r_bdmv,r_bdmw,len(bdmv_mass),len(bdmw_mass),n_points]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%lf','%d','%d','%d'],delimiter=',')

        sys.stdout.write('Done\n')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()


for p in jobs: 
    while(p.is_alive()):
        pass

export = open('./results_'+now+'/results_'+now+'.csv', 'w')

filenames = ['./results_'+now+'/results_'+str(i)+'.csv' for i in range(processes)]
with export as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)    
os.system('rm  potential.out')
sys.stdout.write('Fresh data from the oven!\n')

