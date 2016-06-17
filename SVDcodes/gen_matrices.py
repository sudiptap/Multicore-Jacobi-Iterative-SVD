import random
import math
import os
import subprocess as sp

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.abspath(script_dir + "/../data")

print data_dir

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

sym_mat_sizes = [30, 50,
    500, 750, 1000, 1250, 1500, 1750, 2000]
mat_sizes = [[30, 25], [50, 20],
    [500, 300], [750, 500], [1000, 700], [1250, 900], [1500, 1300], [1750, 1500], [2000,1800]]
num_instances = 10

for n in sym_mat_sizes:
    mat = "sm" + str(n) + "x" + str(n)
    mat_dir = data_dir + "/" + mat
    if not os.path.exists(mat_dir):
        os.mkdir(mat_dir)
    for i in range(num_instances):
        mat_ins = mat + "-" + str(i+1) + ".txt"
        mat_file = mat_dir + "/" + mat_ins
        print mat_file
        print_matrix(gen_symm_matrix(n), mat_file)

for (m, n) in mat_sizes:
    mat = "m" + str(m) + "x" + str(n)
    mat_dir = data_dir + "/" + mat
    if not os.path.exists(mat_dir):
        os.mkdir(mat_dir)
    for i in range(num_instances):
        mat_ins = mat + "-" + str(i+1) + ".txt"
        mat_file = mat_dir + "/" + mat_ins
        print mat_file
        print_matrix(gen_matrix(m, n), mat_file)
