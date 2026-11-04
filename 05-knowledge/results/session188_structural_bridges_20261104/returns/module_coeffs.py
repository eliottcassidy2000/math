"""Nonnegative exterior skeleton coefficients for arbitrary module replacements."""
from itertools import permutations
from collections import Counter
from math import comb
from skeleton import gaps,outer,W

def module_coeffs():
 C={z:Counter() for z in ('h','Jempty','E','S')}
 for p in permutations(outer):
  A,R=gaps(p)
  if R<=A:
   for k in range(len(A)+1):
    if 0<=k-len(R)<=len(A)-len(R):C['h'][k]+=comb(len(A)-len(R),k-len(R))
  for cut in range(6):
   x,y=p[:cut],p[cut:];ax,rx=gaps(x);ay,ry=gaps(y)
   for mode in ('Jempty','E','S'):
    a,b=set(ax),set(ay)
    if mode in ('Jempty','S'):
     if not x or x[-1] not in W:continue
     a.discard(len(x))
    if mode in ('Jempty','E'):
     if not y or y[0] not in W:continue
     b.discard(0)
    forced=0
    if mode=='E':
     if len(x) not in a:continue
     a.remove(len(x));forced=1
    if mode=='S':
     if 0 not in b:continue
     b.remove(0);forced=1
    if not(rx<=a and ry<=b):continue
    n=len(a)+len(b);r=len(rx)+len(ry)
    for k in range(n+1):
     if 0<=k-r<=n-r:C[mode][k+forced]+=comb(n-r,k-r)
 return {mode:{str(k):v for k,v in t.items() if v} for mode,t in C.items()}
if __name__=='__main__':print(module_coeffs())
