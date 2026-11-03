#!/usr/bin/env python3
"""All one-vertex orientations of R4: polynomial and independent literal DP."""
from collections import Counter

def need(p,msg):
 if not p:raise ValueError(msg)
def H(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1<<n):
  for v in range(n):
   c=D[mask][v];z=rows[v]&~mask
   while c and z:
    b=z&-z;j=b.bit_length()-1;D[mask|b][j]+=c;z-=b
 return sum(D[-1])
def add(rows,mask):
 n=len(rows)
 return tuple(r|((1<<n) if not(mask>>i&1) else 0) for i,r in enumerate(rows))+(mask,)

r=(996,993,994,999,975,912,800,576,128,256);n=len(r)
f0=H(add(r,0));fi=[H(add(r,1<<i)) for i in range(n)]
fij={(i,j):H(add(r,(1<<i)|(1<<j))) for i in range(n) for j in range(i+1,n)}
a=[x-f0 for x in fi]; b={(i,j):v-fi[i]-fi[j]+f0 for (i,j),v in fij.items()}
values=Counter()
for mask in range(1<<n):
 D=[i for i in range(n) if mask>>i&1]
 prediction=f0+sum(a[i] for i in D)+sum(b[i,j] for ix,i in enumerate(D) for j in D[ix+1:])
 direct=H(add(r,mask));need(prediction==direct,f'quadratic mismatch {mask}')
 values[direct]+=1
need(not values[875] and not values[877],'unexpected ternary branch')
print('PASS: 1024 orientations, polynomial = literal DP; no H=875 or 877')
print('values near targets:',sorted(v for v in values if 700<=v<=1000))
