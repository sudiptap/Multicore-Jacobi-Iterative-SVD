import random
import os
import sys
import time
import numpy
import subprocess as sp

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


num_run = 5
results = {}

nvalues = [100, 200, 300]
solvers = [1, 7]


def run_and_average(n, s):
    cmd = "./jacobi -m %d -n %d -s %d" % (n, n, s)
    print cmd
    timing = 0.0;
    updates = 0.0;

    for i in range(num_run):
        output = sp.check_output(cmd, shell=True).split('\n')
        for line in output:
            dat = line.split(',')
            if (dat[0] == '@@@@'):
               updates += float(dat[3])
               timing += float(dat[4])
    print "%d\t%f\t%f" % (n, updates/num_run, timing/num_run)
    results[s]['updates'].append(updates/num_run)
    results[s]['timing'].append(timing/num_run)

                

def plot(to_plot, yname, ylabel, xlabel):
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.xlabel(xlabel) 
    plt.ylabel(ylabel) 
    plt.tight_layout(1.25)
    plt.grid(b=True, which='both', color='0.15')
    plots = []
    for [s, config] in to_plot:
        print len(nvalues)
        print len(results[s][yname])
        line, = plt.plot(nvalues, results[s][yname], label = config, linewidth=4, marker='o')
        plots.append(line)
    #plt.legend(handles=plots, loc=1)
    plt.savefig(yname+".pdf")

update_plot = []
time_plot = []

for s in solvers:
    results[s] = {}
    results[s]['updates'] = []
    results[s]['timing'] = []
    for n in nvalues:
        run_and_average(n, s)
    update_plot.append([s, "solver "+str(s)])
    time_plot.append([s, "solver "+str(s)])

plot(update_plot, 'updates', 'Number of updates', 'Matrix size')
plot(time_plot, 'timing', 'Time taken', 'Matrix size')

