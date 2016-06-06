A=dlmread('mymatrix.txt');
[U,S,V]=svd(A);
dlmwrite('matlab.dat',diag(S))
