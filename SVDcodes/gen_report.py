import random
import math
import os
import sys
import subprocess as sp
from test_config import *

sweeps = {}

def read_log(logfile):
    log = open(logfile, "rt")
    for line in log:
        dat = line.strip().split(",");
        mat = "%sx%s" % (dat[1], dat[2])
        solver = dat[4]
        if mat in sweeps.keys():
            if solver in sweeps[mat].keys():
                sweeps[mat][solver].append(int(dat[5]))
            else:
                sweeps[mat][solver] = [int(dat[5])]
        else:
            sweeps[mat] = {}
            sweeps[mat][solver] = [int(dat[5])]
    log.close()


def gen_report(m,n,prefix,solvers,f):
  if (prefix == "m"):
    sol_codes = [100+x for x in solvers]
  else: 
    sol_codes = [200+x for x in solvers]
  mat = "%dx%d" % (m, n)
  print >>f, mat,
  for code in sol_codes:
      l = sweeps[mat][str(code)]
      print >>f, "& %.1f" % (float(sum(l))/len(l)),
  print >>f, "\\\\"

def copy_file(fname):
    sp.call('cp '+fname+' '+paper_dir+"/", shell=True)


read_log("mylog.txt")

with open("twosided.tex", "wt") as f:
    for n in sym_mat_sizes:
        gen_report(n,n,"sm",solvers,f)

with open("onesided.tex", "wt") as f:
    for (m, n) in mat_sizes:
      gen_report(m,n,"m", solvers, f);

copy_file("twosided.tex")
copy_file("onesided.tex")
