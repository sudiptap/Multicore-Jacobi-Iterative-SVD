import sys
import re
from itertools import izip
i=0;
for s,t,s1,u in izip(sys.stdin,sys.stdin,sys.stdin,sys.stdin):
    s1 = re.split(' = |,|\n',s)
    t1 = re.split(' = |,|\n',t)
    u1 = re.split(' = |,|\n',u)
    print s1[1], t1[1], u1[1]
