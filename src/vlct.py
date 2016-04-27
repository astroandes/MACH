import argparse # For parsing arguments

parser = argparse.ArgumentParser(__file__, description="Velocity adjuster for concentration in halos")
parser.add_argument("--dir", "-d", help="Path to the directory with the halos", type=str)
parser.add_argument("--verbose","-v",help="Turns on verbose mode", action="store_true")
parser.add_argument("--time", "-t", help="It measures the time of execution", action="store_true")
parser.add_argument("--output","-o",help="File output", type=str)
args = parser.parse_args()
if args.dir == None:
    print("Please specify the path of the directory with the option -d\nUse --help for more information")
    exit(0)

if args.output == None:
    print("Saving output on output.csv")

import numpy as np # For numerical computation
import os # For handling directories
from scipy.optimize import fsolve

if args.time:
    from time import clock # For measuring time

# Defines the relation between the concentration and the greatest velocity
def f(c,v_max):
    return v_max - np.sqrt(0.216*c/(np.log(1+c)-c/(1+c)))


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
            data = np.loadtxt(filename,delimiter=',')
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
            r = r/r[-1]
            del x,y,z

            n_particles = len(r)

            m = np.arange(n_particles)+1.0
            m /= n_particles

            v = np.sqrt(m/r)
            v_max = np.max(v)
            c = fsolve(f,10.0,args=(v_max))[0]
            name = filename.split('/')[-1].split('.')[0]

            if args.time:
                interval = clock()-t_init
                results.append([name,str(c),str(n_particles),str(interval)])
                del interval
            else:
                results.append([name,str(c),str(n_particles)])

            if args.verbose:
                print(name,c,n_particles)

            del name,c,n_particles

    finally:
        print('writing results')
        results = [','.join(x) for x in results]
        if args.time:
            results = ['haloID,c,n_particles,time'] + results
        else:
            results = ['haloID,c,n_particles'] + results
        outfile.write('\n'.join(results))
        outfile.close()

if __name__=='__main__':
    main()
