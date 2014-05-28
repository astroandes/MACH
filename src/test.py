import numpy as np, pylab, math, sys, timeit, os, random, datetime, nfw, fit, mcmc, plotter, multiprocessing, time
from scipy.optimize import fsolve

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
        
        mass = mass/mass[-1]
        radius = radius/radius[-1]

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
    
        export = open(filename_export, 'a')
        line = [[file_id,x_center,y_center,z_center,c,c_max,c_min]]
        np.savetxt(export,line,fmt=['%d','%lf','%lf','%lf','%lf','%lf','%lf'],delimiter=',')

for i in range(processes):
    p = multiprocessing.Process(target=run, args=(lists[i],i))
    jobs.append(p)
    p.start()

