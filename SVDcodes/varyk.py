import random
import math
import os
import sys
import subprocess as sp
from test_config import *

kvalues = [1, 2, 4]
solvers = [3, 6, 9]

def run_cmd_tee(cmd, log):
    proc = sp.Popen(cmd, shell=True, stdout=sp.PIPE, bufsize=1)
    for line in iter(proc.stdout.readline, ''):
        sys.stdout.write(line)
        if (line[0:4] == '@@@@'):
          log.write(line)
          log.flush()
    if proc.wait() != 0:
      print "Error running code!! exiting"
      exit(-1)

def run_solvers(m,n,prefix,solvers,topk,log):
  if (prefix == "m"):
    lamda = 1.0 - 2.2019*math.pow(m, -0.3382)
    sol_codes = [100+x for x in solvers]
  else: 
    lamda = 1.0 - 2.9267*math.pow(m, -0.4284)
    sol_codes = [200+x for x in solvers]
  mat = prefix + str(m) + "x" + str(n)
  for i in range(num_instances):
    mat_file = data_dir + "/" + mat + "/" + mat + "-" + str(i+1) + ".txt"
    for s in sol_codes:
      cmd = "%s %s %d %e %e %d" % (cmd_exe, mat_file, s, tol, lamda, topk)
      print cmd
      run_cmd_tee(cmd, log)

log = open("mylog_varyk.txt", "at")
for n in sym_mat_sizes:
  for k in kvalues: 
    run_solvers(n, n, "sm", solvers, k, log);

for (m, n) in mat_sizes:
  for k in kvalues: 
    run_solvers(m, n,"m", solvers, k, log);

log.close();
