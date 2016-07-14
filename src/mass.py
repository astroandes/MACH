import argparse # For parsing arguments

parser = argparse.ArgumentParser(__file__, description="A simple MCMC adjuster for concentration in halos")
parser.add_argument("--dir", "-d", help="Path to the directory with the halos", type=str)
parser.add_argument("--verbose","-v",help="Turns on verbose mode", action="store_true")
parser.add_argument("--time", "-t", help="It measures the time of execution", action="store_true")
parser.add_argument("--plot","-p",help="Exports useful plots for each halo", action="store_true")
parser.add_argument("--output","-o",help="File output", type=str)
parser.add_argument("--threads","-T",help="Number of threads (8 by default)", type=int)
parser.add_argument("--radius", "-r", help="Desired value for the virial radius. If empty the code takes the largest radius.", type=float)
parser.add_argument("--csv", "-c", help="Is input in csv format?")
args = parser.parse_args()
if args.dir == None:
    print("Please specify the path of the directory with the option -d\nUse --help for more information")
    exit(0)

if args.output == None:
    print("Saving output on output.csv")

if args.threads == None:
    threads = 1;
else:
    threads = args.threads

print("Working with %d threads" % threads)

import numpy as np, emcee # For numerical computation and MCMC
import os # For handling directories
if args.time:
    from time import clock # For measuring time
if args.plot:
    import matplotlib.pyplot as plt # For plotting

# Navarro-Frenk-White mass profile with logarithmic scale on mass, radius and
# concentration. It assumes the radius is normalized respect to the virial
# radius.
def nfw_log_mass(log_r,log_c):
    return np.log((np.log(1.0+np.exp(log_c+log_r))-np.exp(log_c+log_r)/(1.0+np.exp(log_c+log_r)))/(np.log(1.0+np.exp(log_c))-np.exp(log_c)/(1.0+np.exp(log_c))))

# Loglikelihood using the chi squared. If the concentration is less than one it
# returns -inf.
def log_likelihood(log_c,log_r,log_m):
    n = 1.0*len(log_m)
    if log_c >= 0:
        s = np.sum(log_m-nfw_log_mass(log_r,log_c))
        sigma2 = (1.0 - np.exp(log_r))/(np.exp(log_r)) * (n **1.15)/4500.0
        return -0.5*s*s/sigma2
    return -np.inf

def main():
    results = []
    if args.output == None:
        outfile = open('output.csv','w')
    else:
        outfile = open(args.output,'w')
    try:



        root_path = args.dir
        directory_list = [root_path+'/'+x for x in os.listdir(root_path)]
        directory_list.sort()

        guess = 2
        n_dimensions = 1
        n_walkers = 2

        for filename in directory_list:
            if args.time:
                t_init = clock()
            if args.csv:
                data = np.loadtxt(filename,delimiter=',')
            else:
                data = np.loadtxt(filename)

            x_center,y_center,z_center = 0,0,0
            data = np.transpose(data)
            x,y,z = data
            del data

            x -= x_center
            y -= y_center
            z -= z_center

            r = np.sqrt(x*x+y*y+z*z)
            r = np.sort(r)
            r = np.delete(r,0)
            if args.radius == None:
                virial_radius = r[-1]
            else:
                virial_radius = args.radius
            r = r/virial_radius

            r = np.delete(r,-1)

            del x,y,z

            n_particles = len(r)

            m = np.arange(n_particles)+1.0
            m /= n_particles

            log_r,log_m = np.log(r),np.log(m)
            del m,r

            position = np.linspace(np.log(guess),np.log(2*guess),n_walkers).reshape((n_walkers,1))

            sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, log_likelihood, args=(log_r, log_m),threads=threads)
            sampler.run_mcmc(position, 500)

            chain = np.exp(sampler.flatchain).flatten()
            chisq = -sampler.flatlnprobability.flatten()
            c_low,c_mid,c_high = np.percentile(chain,[16,50,84])

            name = filename.split('/')[-1].split('.')[0]
            name = name.split('_')[-1]

            if args.plot:
                x_min = c_mid - 2*(c_mid - c_low)
                x_max = c_mid + 2*(c_high - c_mid)

                chisq = chisq[np.logical_and(chain > x_min, chain < x_max)]
                chain = chain[np.logical_and(chain > x_min, chain < x_max)]
                a,b,c = np.polyfit(chain,chisq,2)
                x = np.linspace(x_min,x_max,100)
                y = a*x**2+b*x+c
                plt.plot(chain, chisq, 'ok')
                plt.plot(x,y,'g')

                plt.axvspan(c_low,c_high, color='b',label='16-84',alpha=0.5)
                plt.axvline(x=c_mid, color='r',label='50')

                plt.legend()
                plt.xlabel('Concentration')
                plt.ylabel('Chi Squared')
                plt.savefig('%s.png'%name,dpi=200)
                plt.close()
                del x_min, x_max, a, b, c, x, y

            if args.time:
                interval = clock()-t_init
                results.append([name,str(c_low),str(c_mid),str(c_high),str(n_particles),str(virial_radius),str(interval)])
                del interval
            else:
                results.append([name,str(c_low),str(c_mid),str(c_high),str(n_particles),str(virial_radius)])

            if args.verbose:
                print(name,c_low,c_mid,c_high,n_particles,virial_radius)

            del name,sampler,chain,c_low,c_mid,c_high,n_particles,virial_radius

    finally:
        print('writing results')
        results = [','.join(x) for x in results]
        if args.time:
            results = ['haloID,c_low,c_mid,c_high,n_particles,virial_r,time'] + results
        else:
            results = ['haloID,c_low,c_mid,c_high,n_particles,virial_r'] + results
        outfile.write('\n'.join(results))
        outfile.close()

if __name__=='__main__':
    main()
