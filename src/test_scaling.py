import random
import os
import sys
import time
import numpy
import math
import subprocess as sp
#from scipy.linalg import hilbert
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

num_run = 2
results = []
script_dir = os.path.dirname(os.path.abspath(__file__))
top_dir = os.path.abspath(script_dir + "/..")
data_dir = top_dir + "/data"

#nvalues = [30, 50, 100]
nvalues = [500, 1000, 1500, 2000, 2500, 3000]

solvers = [[9, 1], [9, 2], [9, 4], [9, 8], [9, 16]]

def run_and_average(n, s, k, p):
    base_cmd = "./jacobi -m %d -n %d -s %d -k %d -t %d" % (n, n, s, k, p)
    print base_cmd
    timing = 0.0;
    updates = 0.0;
    solver = str(s)+"."+str(k)

    for i in range(num_run):
        mat_name = "sm%dx%d" % (n, n)
        cmd = "%s -f %s/%s/%s-%d.txt" % (base_cmd, data_dir, mat_name, mat_name, i+1) 
        print cmd
        output = sp.check_output(cmd, shell=True).split('\n')
        for line in output:
            dat = line.split(',')
            if (dat[0] == '@@@@'):
               print line
               updates += float(dat[3])
               timing += float(dat[4])
    print "%d\t%d\t%d\t%d\t%f\t%f" % (s, n, k, p, updates/num_run, timing/num_run)
    results.append([s, n, k, p, updates/num_run, timing/num_run])
    with open("results_scaling.py", "wt") as f:
        print  >>f, "results =", results
        f.close()

for [s, p] in solvers:
    for n in nvalues:
        run_and_average(n, s, 4, p)


