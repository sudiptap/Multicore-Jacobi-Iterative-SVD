import random
import math
import subprocess as sp

def read_vector(fname):
  v = []
  with open(fname, "r") as f:
    for l in f:
        v.append(float(l.rstrip()))
  return v

def gen_matrix(m, n):
    A = [[0 for y in range(n)] for x in range(m)]
    for i in range(m):
        for j in range(n):
            A[i][j] = random.uniform(1, 10)
    return A


def gen_symm_matrix(n):
    A = [[0 for y in range(n)] for x in range(n)]
    for i in range(n):
        for j in range(i, n):
            A[i][j] = A[j][i] = random.uniform(1, 10)
    return A

def print_matrix(A, fname):
    m = len(A)
    n = len(A[0])
    with open(fname, "w") as f:
      print >>f, m, n
      for i in range(m):
        print >>f, '\t'.join(["%.12e" % (x) for x in A[i]])
      f.close()


print_matrix(gen_matrix(10, 5), "mymatrix.txt")
print_matrix(gen_symm_matrix(10), "mymatrixsymm.txt")

exit(0)

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
cmd = "./svd mymatrix.txt 5 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
cmd = "./svd mymatrix.txt 7 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
cmd = "./svd mymatrix.txt 77 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
cmd = "./svd mymatrix.txt 8 %e %e" % (tol, lamda)
print cmd
sp.call(cmd, shell=True) 
