import random
import os
import sys
import time
import numpy
import math
import subprocess as sp
from scipy.linalg import hilbert
from matplotlib.lines import Line2D

import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "font.family": "serif", # use serif/main font for text elements
    "font.size": 12,        # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
#         r"\usepackage[utf8x]{inputenc}",
#         r"\usepackage[T1]{fontenc}",
#         r"\usepackage{cmbright}",
         r"\usepackage{times}",
         ]
}
mpl.rcParams.update(pgf_with_pdflatex)
import matplotlib.pyplot as plt

markers = Line2D.filled_markers

num_run = 4
results = {}

#nvalues = [50, 75, 100, 125] #, 1000, 1500, 2000]
nvalues = [500, 1000, 1500, 2000, 2500, 3000]

#solvers = [[1, 4], [7, 4], [7, 8], [9, 4], [9, 8]]
solvers = [[1, 4], [7, 4], [7, 8], [7, 16], [7, 32], [9, 4], [9, 8], [9, 16], [9, 32]]

def generate_type4_matrix(m, n):
     x = numpy.zeros((m,n));
     for i in range(m):
         for j in range(n):
             x[i][j]=random.uniform(-1,1)
     numpy.savetxt('test_type4.out', x) 

def generate_type3_matrix(m,n):
    if m!=n:
        print "error : m!=n"
    else:
        c = 0.2
        s = math.sqrt(1 - math.pow(c,2))
        d = numpy.zeros((m,m));
        for i in range(m):
            d[i][i]=math.pow(s,i-1)
        r = numpy.zeros((m,m))
        numpy.fill_diagonal(r,1)
        for i in range(m):
            for j in range(m):
                if j>i:
                    r[i][j]=(-c)
        t = numpy.dot(d,r)
        numpy.savetxt('test_type3.out',t)
    
def generate_type1_matrix(m,n):    
    if m!=n:
        print "Hilbert matrices are squared matrix"
    else:
        h = numpy.zeros((m,m))        
        h = hilbert(m)
        numpy.savetxt('test_type1.out',h)

def run_and_average(n, s, k):
    #generate_type4_matrix(3, 3)
    cmd = "./jacobi -m %d -n %d -s %d -k %d" % (n, n, s, k)
    print cmd
    timing = 0.0;
    updates = 0.0;
    solver = str(s)+"."+str(k)

    for i in range(num_run):
        output = sp.check_output(cmd, shell=True).split('\n')
        for line in output:
            dat = line.split(',')
            if (dat[0] == '@@@@'):
               updates += float(dat[3])
               timing += float(dat[4])
    print "%d\t%d\t%f\t%f" % (n, k, updates/num_run, timing/num_run)
    results[solver]['updates'].append(updates/num_run)
    results[solver]['timing'].append(timing/num_run)

                

def plot(to_plot, yname, ylabel, xlabel):
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel(xlabel) 
    plt.ylabel(ylabel) 
    plt.tight_layout(1.25)
    plt.grid(b=True, which='both', color='0.15')
    plots = []
    for [s, config] in to_plot:
        line, = plt.plot(nvalues, results[s][yname], label = config, linewidth=3, marker=markers[len(plots)-1], markersize=5)
        plots.append(line)
    plt.legend(plots, loc=2, labelspacing=0.05, numpoints=1, borderpad=0, frameon=False, handlelength=1)
    plt.savefig(yname+".pdf")

update_plot = []
time_plot = []

for [s, k] in solvers:
    solver = str(s)+"."+str(k)
    results[solver] = {}
    results[solver]['updates'] = []
    results[solver]['timing'] = []
    for n in nvalues:
        run_and_average(n, s, k)
    update_plot.append([solver, "s"+str(s)+".k"+str(k)])
    time_plot.append([solver, "s"+str(s)+".k"+str(k)])

plot(update_plot, 'updates', 'Number of updates', 'Matrix size')
plot(time_plot, 'timing', 'Time taken', 'Matrix size')

