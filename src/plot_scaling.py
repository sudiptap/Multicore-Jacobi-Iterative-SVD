import random
import os
import sys
import time
import numpy
import math
import subprocess as sp
#from scipy.linalg import hilbert
from matplotlib.lines import Line2D
import itertools as it
import pandas as pd
import numpy as np

import matplotlib as mpl
mpl.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "font.family": "serif", # use serif/main font for text elements
    "font.size": 24,        # use serif/main font for text elements
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

from results_scaling import *
df = pd.DataFrame(results, columns=['solver', 'n', 'tau', 'threads', 'update', 'time'])

num_run = 2
script_dir = os.path.dirname(os.path.abspath(__file__))
top_dir = os.path.abspath(script_dir + "/..")
data_dir = top_dir + "/data"

nvalues = [500, 1000, 1500, 2000, 2500, 3000]
nvalues = [30, 50, 100]
solvers = [[9, 1], [9, 2], [9, 4], [9, 8], [9, 16]]

solver_names = {
        1: "Cyclic",
        7: "SeqJTS",
        9: "ParJTS",
}


def get_col(cols, values):
  select = [True for x in df.index]
  for (c,v) in it.izip(cols, values):
    #print c, v
    select = select & (df[c] == v)
    #print select
  print df[select]
  return df[select]


def plot(to_plot, yname, ylabel, xlabel, ycol, xcol):
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    num_plots = 0
    for [n, config] in to_plot:
        newdf = get_col(['n'], [n])
        thread1_time = newdf[ycol][newdf[xcol] == 1].iloc[0]
        newdf[ycol] = thread1_time / newdf[ycol]
        print newdf
        newax = newdf.plot(x=xcol, y=ycol, label = config, linewidth=3, marker=markers[num_plots], markersize=5)
        num_plots += 1
    lines, labels = newax.get_legend_handles_labels() 
    #print "*************", labels, "***************"
    plt.legend(lines, labels, loc=1, labelspacing=0.07, numpoints=1, borderpad=0, frameon=False, handlelength=1)
    plt.xlabel(xlabel) 
    plt.ylabel(ylabel) 
    fig.tight_layout()
    plt.subplots_adjust(left=0.01, right=0.95, top=0.95, bottom=0.01)
    plt.grid(b=True, which='both', color='0.15')
    plt.savefig(yname+".pdf", bbox_inches='tight')

update_plot = []
time_plot = []

print results

print df

#ax1 = df2.plot(x='n', y='time', label = "haha", linewidth=3, marker=markers[len(plots)-1], markersize=5)
#lines, lab = ax1.get_legend_handles_labels()
#print lines
#print labels
#plt.legend(plots, labels, loc=2, labelspacing=0.05, numpoints=1, borderpad=0, frameon=False, handlelength=1)
#plt.savefig("test.pdf")


#read_data("mylog.txt")

for n in nvalues:
    config = "$n$=%d" % (n)
    time_plot.append([n, config])

plot(time_plot, 'scaling', 'Scaling factor', 'Number of threads ($p$)', 'time', 'threads')

