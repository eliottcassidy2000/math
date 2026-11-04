"""Finite occupied-gap coefficient certificate for the first-marked clone family."""
from itertools import permutations
from collections import Counter
from math import comb
from skeleton import gaps,outer,W

def family_coeffs():
 ans={z:Counter() for z in ['h','p','Fouter','Gouter','Ffirst','Gfirst','eouter','souter','efirst','sfirst']}
 def add(mode,k,r):
  for q in range(r+1):ans[mode][k-q]+=(-1)**q*comb(r,q)
 for pi in permutations(outer):
  for mode in ['h','eouter','souter','efirst','sfirst']:
   A,R=gaps(pi);A=set(A)
   if mode=='eouter':
    if pi[-1] not in W:continue
    A.discard(5)
   if mode=='souter':
    if pi[0] not in W:continue
    A.discard(0)
   if mode=='efirst':
    if 5 not in A:continue
    A.discard(5) # Every other clone is greater and cannot precede first clone in its final run.
   if mode=='sfirst':
    if 0 not in A:continue
    # The forced first clone's gap still accepts all greater clones.
   if R<=A:add(mode,len(A),len(R))
  for cut in range(6):
   x,y=pi[:cut],pi[cut:];ax,rx=gaps(x);ay,ry=gaps(y)
   for mode in ['p','Fouter','Gouter','Ffirst','Gfirst']:
    a,b=set(ax),set(ay)
    if mode=='Fouter':
     if not x or x[-1] not in W:continue
     a.discard(len(x))
    if mode=='Gouter':
     if not y or y[0] not in W:continue
     b.discard(0)
    if mode=='Ffirst':
     if len(x) not in a:continue
     a.discard(len(x))
    if mode=='Gfirst':
     if 0 not in b:continue
    if rx<=a and ry<=b:add(mode,len(a)+len(b),len(rx)+len(ry))
 return {mode:{str(k):c for k,c in C.items() if c} for mode,C in ans.items()}
if __name__=='__main__':print(family_coeffs())
