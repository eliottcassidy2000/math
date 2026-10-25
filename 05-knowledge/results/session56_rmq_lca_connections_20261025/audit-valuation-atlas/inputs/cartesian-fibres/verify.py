#!/usr/bin/env python3
from itertools import permutations
from collections import Counter
from fractions import Fraction
from math import factorial,comb
from pathlib import Path
import json, argparse
P=Path(__file__).resolve().parent

def ck(b,m):
 if not b:raise RuntimeError(m)
def tree(a):
 if not a:return ()
 k=min(range(len(a)),key=lambda i:(a[i],i))
 return (tree(a[:k]),tree(a[k+1:]))
def size(t):return 0 if not t else 1+size(t[0])+size(t[1])
def hooks(t):return 1 if not t else size(t)*hooks(t[0])*hooks(t[1])
def signature(a):return tuple(min(range(i,j+1),key=lambda k:(a[k],k)) for i in range(len(a)) for j in range(i,len(a)))
def run():
 rows=[]
 for n in range(9):
  counts=Counter(tree(p) for p in permutations(range(n)))
  ck(len(counts)==comb(2*n,n)//(n+1),'Catalan')
  weights=Fraction()
  for t,c in counts.items():
   ck(c*hooks(t)==factorial(n),'hook fibres')
   weights+=Fraction(1,hooks(t))
  ck(weights==1,'normalization')
  rows.append({'n':n,'shapes':len(counts),'permutations':sum(counts.values()),'min_fibre':min(counts.values()),'max_fibre':max(counts.values()),'weighted_sum':str(weights)})
 a=(1,0,2);b=(2,0,1)
 ck(tree(a)==tree(b) and signature(a)==signature(b) and (a[0]<a[2])!=(b[0]<b[2]),'cross-branch hostile')
 # Every RMQ signature determines exactly one stable tree, including ties.
 from itertools import product
 sigs={};cases=0
 for n in range(1,7):
  for a in product(range(3),repeat=n):
   sig=signature(a);t=tree(a)
   if (n,sig) in sigs:ck(sigs[n,sig]==t,'RMQ determines tree')
   sigs[n,sig]=t;cases+=1
 capacity=[]
 for n in [7,8,10,20,50,100]:
  C=comb(2*n,n)//(n+1);N=n*(n-1)//2;B=max(0,N-(C*factorial(n)).bit_length()+1)
  ck(C*factorial(n)*2**B>=2**N,'bound ceiling')
  if B:ck(C*factorial(n)*2**(B-1)<2**N,'bound sharp from cardinality')
  capacity.append({'n':n,'ordered_shapes':C,'necessary_fixed_sidecar_bits_from_orbit_bound':B})
 return {'fibre_counts':rows,'ternary_array_cases':cases,'cross_branch_arrays':[[1,0,2],[2,0,1]],'capacity_bounds':capacity,'caveat':'Encoding capacity only; no counting-time lower bound.'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();s=json.dumps(run(),indent=2,sort_keys=True)+'\n'
 if a.write:(P/'certificate.json').write_text(s)
 else:ck((P/'certificate.json').read_text()==s,'certificate')
 print('PASS: all permutation fibres n<=8, tied RMQ signatures, exact capacity bounds')
if __name__=='__main__':main()
