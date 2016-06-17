import random
import math
import os
import sys
import subprocess as sp

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(script_dir + "/../data")
cmd_exe = script_dir + "/svd"

num_instances = 2
tol = 1e-15
lamda = 0.0

sym_mat_sizes = [500, 750, 1000, 1250, 1500, 1750, 2000]
mat_sizes = [[500, 300], [750, 500], [1000, 700], [1250, 900], [1500, 1300], [1750, 1500], [2000,1800]]
solvers = [1,2,3,4,5,6,7,8,9,10]

sym_mat_sizes = [30, 50]
mat_sizes = [[30, 25], [50, 20]]
solvers = [1,2]

def run_cmd_tee(cmd, log):
    proc = sp.Popen(cmd, shell=True, stdout=sp.PIPE, bufsize=1)
    for line in iter(proc.stdout.readline, ''):
        sys.stdout.write(line)
        if (line[0:4] == '@@@@'):
          log.write(line)
    if proc.wait() != 0:
      print "Error running code!! exiting"
      exit(-1)

def run_solvers(m,n,prefix,solvers,log):
  if (prefix == "m"):
    lamda = 1.0 - 2.2019*math.pow(m, -0.3382)
    sol_codes = [100+x for x in solvers]
  else: 
    lamda = 1.0 - 2.2019*math.pow(m, -0.3382)
    sol_codes = [200+x for x in solvers]
  mat = prefix + str(m) + "x" + str(n)
  for i in range(num_instances):
    mat_file = data_dir + "/" + mat + "/" + mat + "-" + str(i+1) + ".txt"
    for s in sol_codes:
      cmd = "%s %s %d %e %e" % (cmd_exe, mat_file, s, tol, lamda)
      print cmd
      run_cmd_tee(cmd, log)

log = open("mylog.txt", "at")
for n in sym_mat_sizes:
  run_solvers(n,n,"sm", solvers, log);

for (m, n) in mat_sizes:
  run_solvers(m,n,"m", solvers, log);

log.close();
