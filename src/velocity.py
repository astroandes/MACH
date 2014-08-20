import numpy as np, pylab, math, sys, timeit, os, random, datetime, nfw, mcmc,fit , plotter, multiprocessing, time
from scipy.optimize import fsolve


G = 4.49E-15 # Mpc^3/(solar masses Gyr^2)
processes = int(sys.argv[6])
now = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")
mass_element = 	8.721E9    
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

        path = os.path.expanduser('./'+str(sys.argv[1])+'/'+filename)
        data = np.loadtxt(open(path, 'r'), delimiter=",",skiprows=int(sys.argv[5]))

        x = data[:,int(sys.argv[2])]
        y = data[:,int(sys.argv[3])]
        z = data[:,int(sys.argv[4])]

        n_points = len(x)
        
        file_id = int(filename.split('_')[1])
        positions_name ='positions_'+str(file_id)+'.dat'
        potential_name ='potential_'+str(file_id)+'.dat'
        results_folder ='./results_'+now+'/'+str(file_id)
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

        crit = 740.0

        vir_index = np.argmin(np.abs(avg_density-crit*rho_back))

        if np.argmin(np.abs(avg_density-crit*rho_back)) > 1:
            r_vir = radius[vir_index]
        else:
            r_vir = radius[-1]
            vir_index = len(avg_density)-1

        mass = np.resize(mass,vir_index+1)
        radius = np.resize(radius,vir_index+1)
	
	velocity = np.sqrt(G*mass/radius)
	v_max = np.max(velocity)
	v_vir = velocity[-1]
	def f(c):
		return (v_max/v_vir) - 0.216*c/(np.log(1+c)-c/(1+c)) 
	
	c = fsolve(f,10.0)
        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c,f(c),r_vir,len(mass),n_points]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf','%d','%d'],delimiter=',')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()

