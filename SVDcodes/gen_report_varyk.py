import random
import math
import os
import sys
import subprocess as sp
from test_config import *

sweeps = {}
solvers = []
kvalues = []

def read_log(logfile):
    log = open(logfile, "rt")
    for line in log:
        dat = line.strip().split(",");
        mat = "%sx%s" % (dat[1], dat[2])
        solver = dat[4]
        topk = int(dat[6])
        if mat in sweeps.keys():
            if solver in sweeps[mat].keys():
              if topk in sweeps[mat][solver].keys():
                sweeps[mat][solver][topk].append(int(dat[5]))
              else:
                sweeps[mat][solver][topk] = [int(dat[5])]
            else:
                sweeps[mat][solver] = {}
                sweeps[mat][solver][topk] = [int(dat[5])]
        else:
            sweeps[mat] = {}
            sweeps[mat][solver] = {}
            sweeps[mat][solver][topk] = [int(dat[5])]
        s = int(solver) % 100
        if (not s in solvers):
            solvers.append(s)
        if (not topk in kvalues):
            kvalues.append(topk)
    log.close()

def gen_header(solvers, kvalues, f):
  print >>f, "\\begin{tabular}{r%s}" % ('r' * len(solvers) * len(kvalues))
  print >>f, "\\toprule"
  print >>f, "Matrix",
  for s in solvers:
    print >>f, "& \multicolumn{%d}{c}{%s}" % (len(kvalues), solver_names[s]),
  print >>f, "\\\\"
  print >>f, "\\cmidrule(lr){%d-%d} \\cmidrule(lr){%d-%d}" % (2, 2+len(kvalues)-1, 2+len(kvalues), 2+2*len(kvalues)-1 )
  for s in solvers:
    for k in kvalues:
      print >>f, "& $\\tau$=%d" % (k),
  print >>f, "\\\\"
  print >>f, "\\midrule"

def gen_footer(solvers, f):
  print >>f, "\\bottomrule"
  print >>f, "\\end{tabular}"

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
              print sweeps[mat][code]
              for topk in kvalues:
                  print mat, code, topk
                  l = sweeps[mat][code][topk]
                  print >>f, "& %.1f" % (float(sum(l))/len(l)),
          else:
              print >>f, "& ",
      print >>f, "\\\\"

def copy_file(fname):
    sp.call('cp '+fname+' '+paper_dir+"/", shell=True)


read_log("mylog_varyk.txt")

print solvers
print sweeps

solvers.sort()
kvalues.sort()
print solvers
print kvalues

with open("varyktwo.tex", "wt") as f:
    gen_header(solvers, kvalues, f)
    gen_report("sm", solvers,f)
    gen_footer(solvers, f)


with open("varykone.tex", "wt") as f:
    gen_header(solvers, kvalues, f)
    gen_report("m", solvers,f)
    gen_footer(solvers, f)


copy_file("varykone.tex")
copy_file("varyktwo.tex")
