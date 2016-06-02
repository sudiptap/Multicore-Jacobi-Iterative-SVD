import random
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
  for i in range(m):
    for j in range(n):
      t[j] = "%.12f" % (random.random())
    print >>f, '\t'.join(t)
  f.close()

sp.call("matlab -nodesktop -nosplash -nojvm < 1.bat", shell=True)
sp.call("./jacobi -m %d -n %d -s 1 -f mymatrix.txt" % (m, n), shell=True)

S1 = read_vector("matlab.dat")
S2 = read_vector("JacobiGSL.dat")

diff_count = 0
for i in range(len(S1)):
  if (abs(S1[i]-S2[i]) > 1e-3):
    print S1[i], S2[i], S1[i]-S2[i]
    diff_count += 1
print diff_count

print S1
print S2

