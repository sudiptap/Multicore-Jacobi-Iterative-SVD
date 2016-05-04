set term png
set output "convergence.png"
plot \
     "parallel.dat" using 1:2 title "parallel" with linespoints, \
     "sequential.dat" using 1:2 title "sequential" with linespoints
