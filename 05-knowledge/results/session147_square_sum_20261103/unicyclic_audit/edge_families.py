#!/usr/bin/env python3
from math import isqrt

def square(m): return isqrt(m)**2==m

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

for k in range(1,100001):
    t=isqrt(2*k+1)
    consecutive=(t%2==1 and t>=3 and t*t==2*k+1)
    require(square(2*k+1)==consecutive, f'consecutive {k}')
    t=isqrt((k+1)//2)
    gap=(2*t*t-1==k and t>=1)
    require(square(2*k+2)==gap, f'gap-two {k}')
    if k%2==0: require(not square(2*k+2),f'even gap {k}')
# Directly search for common neighbors in bounded G_1000.
for k in range(1,999):
    for gap in (1,2):
        if k+gap<=1000 and square(2*k+gap):
            for w in range(1,1001):
                if w not in (k,k+gap):
                    require(not (square(k+w) and square(k+gap+w)),f'triangle {k}, {k+gap}, {w}')
try: require(False,'intentional negative control')
except RuntimeError: pass
else: raise RuntimeError('failure detection disabled')
print('families checked for k=1..100000; no common neighbors for qualifying pairs in G_1000; explicit failure detected')
