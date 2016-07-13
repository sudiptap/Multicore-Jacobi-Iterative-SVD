import random
import math
import os
import sys
import subprocess as sp
from test_config import *

sweeps = {}
solvers = []

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
        s = int(solver) % 100
        if (not s in solvers):
            solvers.append(s)
    log.close()

def gen_header(solvers, f):
  print >>f, "\\begin{tabular}{r%s}" % ('r' * len(solvers))
  print >>f, "\\toprule"
  print >>f, "Matrix",
  for s in solvers:
      print >>f, "& %s" % (solver_names[s]),
  print >>f, "\\\\"
  print >>f, "\\midrule"

def gen_footer(solvers, f):
  print >>f, "\\bottomrule"
  print >>f, "\\end{tabular}"

def gen_report_old(m,n,prefix,solvers,f):
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

def gen_report(prefix, solvers, f):
  base = 2 if (prefix == "sm") else 1
  for mat in sorted(sweeps.keys(), cmp=lambda x,y: 1 if (int(x.split('x')[0]) > int(y.split('x')[0])) else -1):
      valid_solvers = filter(lambda x: int(x)/100 == base, sweeps[mat].keys())
      print mat, len(valid_solvers)
      if (len(valid_solvers) < 1):
          continue
      print >>f, mat,
      for solver in solvers:
          code = str(base*100 + solver)
          if code in sweeps[mat].keys():
              l = sweeps[mat][code]
              print >>f, "& %.1f" % (float(sum(l))/len(l)),
          else:
              print >>f, "& ",
      print >>f, "\\\\"

def copy_file(fname):
    sp.call('cp '+fname+' '+paper_dir+"/", shell=True)


#read_log("mylog.gold.txt")
read_log("mylog.txt")

print solvers
print sweeps

with open("twosided.tex", "wt") as f:
    gen_header(solvers, f)
#   for n in sym_mat_sizes:
#        gen_report(n,n,"sm",solvers,f)
    gen_report("sm", solvers,f)
    gen_footer(solvers, f)

with open("onesided.tex", "wt") as f:
    gen_header(solvers, f)
#    for (m, n) in mat_sizes:
#      gen_report(m,n,"m", solvers, f);
    gen_report("m", solvers,f)
    gen_footer(solvers, f)

copy_file("twosided.tex")
copy_file("onesided.tex")
