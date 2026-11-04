from itertools import permutations
from collections import Counter
from math import comb
from carrier import base,outer,v
W={1,5}

def gaps(seq):
 L=len(seq);allowed=set();need=set()
 for q in range(L+1):
  if (q==0 or base[seq[q-1]]>>v&1) and (q==L or base[v]>>seq[q]&1):allowed.add(q)
  if 0<q<L and not(base[seq[q-1]]>>seq[q]&1):need.add(q)
 return allowed,need

def coeffs():
 ans={z:Counter() for z in ['h','J0','end','start']}
 for pi in permutations(outer):
  A,R=gaps(pi)
  if R<=A:
   for s in range(len(R)+1):ans['h'][len(A)-s]+=(-1)**s*comb(len(R),s)
  for c in range(6):
   x,y=pi[:c],pi[c:];ax,rx=gaps(x);ay,ry=gaps(y)
   for mode in ['J0','end','start']:
    a,b=set(ax),set(ay)
    if mode in ['J0','start']:
     if not x or x[-1] not in W:continue
     a.discard(len(x))
    if mode in ['J0','end']:
     if not y or y[0] not in W:continue
     b.discard(0)
    if mode=='end':
     if len(x) not in a:continue
     a.discard(len(x))
    if mode=='start':
     if 0 not in b:continue
     b.discard(0)
    if not(rx<=a and ry<=b):continue
    k=len(a)+len(b);r=len(rx)+len(ry)
    for s in range(r+1):
     ab=(k+1-s,k-s) if mode=='end' else (k-s,k+1-s)
     key=k-s if mode=='J0' else ab
     ans[mode][key]+=(-1)**s*comb(r,s)
 return {mode:{str(k):v for k,v in C.items() if v} for mode,C in ans.items()}
if __name__=='__main__':print(coeffs())
