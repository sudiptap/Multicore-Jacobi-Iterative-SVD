import random
import math
import subprocess as sp

def read_vector(fname):
  v = []
  with open(fname, "r") as f:
    for l in f:
        v.append(float(l.rstrip()))
  return v



m = 300 
n = 300
fname = "mymatrix.txt"
t = [0 for i in range(n)]
with open(fname, "w") as f:
  print >>f, m
  for i in range(m):
    for j in range(n):
      t[j] = "%.12f" % (random.uniform(1, 10))
    print >>f, '\t'.join(t)
  f.close()

#sp.call("matlab -nodesktop -nosplash -nojvm < 1.bat", shell=True)
#sp.call("./jacobi -m %d -n %d -s 1 -f mymatrix.txt" % (m, n), shell=True)
#
#S1 = read_vector("matlab.dat")
#S2 = read_vector("JacobiGSL.dat")
#
#diff_count = 0
#for i in range(len(S1)):
#  if (abs(S1[i]-S2[i]) > 1e-3):
#    print S1[i], S2[i], S1[i]-S2[i]
#    diff_count += 1
#print diff_count
#
#print S1
#print S2
#
lamda = 1.0 - 2.2019*math.pow(m, -0.3382)
tol = 1e-15
print lamda
cmd = "./a.out mymatrix.txt 77 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
cmd = "./a.out mymatrix.txt 7 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
