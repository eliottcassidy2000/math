#!/usr/bin/env python3
"""Independent audit: exact checks survive -O; no source edits or imported controls."""
from pathlib import Path
import importlib.util, json, hashlib
from fractions import Fraction as F
from math import factorial, gcd, comb
from collections import Counter
from itertools import combinations, product
ROOT=Path(__file__).resolve().parents[1]
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
src=ROOT/'session22-enumeration/local_scale.py'
candidate=module('candidate',src)
baseline=module('baseline',ROOT/'session22-enumeration/baseline.py')
checks=Counter()
def need(x,label):
 checks[label]+=1
 if not x: raise RuntimeError(label)
def parts(n,low=1):
 if n==0: yield ()
 for j in range(low,n+1,2):
  for tail in parts(n-j,j): yield (j,)+tail
def direct(N,core,nofixed):
 out=[F(0) for _ in range(N+1)]
 for n in range(N+1):
  for p in parts(n,3 if nofixed else 1):
   mult=Counter(p); z=1
   for j,k in mult.items(): z*=j**k*factorial(k)
   exponent=sum((j-1)//2 for j in p)+sum(gcd(a,b) for a,b in combinations(p,2))
   poly={n:F(2**exponent,z)}
   if core:
    for j,k in mult.items():
     new={}
     for d,w in poly.items():
      for e in range((N-d)//j+1):
       nd=d+j*e; new[nd]=new.get(nd,F(0))+w*(-1)**e*comb(k+e-1,e)
     poly=new
   for d,w in poly.items():out[d]+=w
 return out
# Literal permutation edge action. No gcd/exponent formula is used here.
def fixed_counts(p):
 n=sum(p); perm=[]; pos=0
 for length in p:
  perm += list(range(pos+1,pos+length))+[pos]; pos+=length
 edges=list(combinations(range(n),2)); seen=set(); orbits=[]
 for edge in edges:
  if edge in seen:continue
  orbit=[]; a,b=edge
  while tuple(sorted((a,b))) not in seen:
   seen.add(tuple(sorted((a,b)))); orbit.append((a,b)); a,b=perm[a],perm[b]
  if (a,b)!=edge: return (0,0)
  orbits.append(orbit)
 total=2**len(orbits); reduced=0
 for bits in product((0,1),repeat=len(orbits)):
  out=[0]*n
  for bit,orbit in zip(bits,orbits):
   for a,b in orbit:
    if bit:a,b=b,a
    out[a]|=1<<b
  has_pair=any(((out[a]^out[b]) & ~((1<<a)|(1<<b)))==0 for a,b in edges)
  reduced+=not has_pair
 return total,reduced
for N in range(21):
 for core,nofixed in product((False,True),repeat=2):
  got,stats=candidate.local_dp(N,core=core,no_fixed=nofixed)
  need(got==direct(N,core,nofixed),'independent partition and substitution all four modes')
  old,oldstats=baseline.packed_dp(N,core=core,no_fixed=nofixed)
  need(got==old,'incoming baseline outputs')
  need(stats['rows']==[{k:r[k] for k in ('cycle_length','states','transitions')} for r in oldstats['rows']], 'same state and transition counts')
for n in range(7):
 raw=F(0); core=F(0)
 for p in parts(n):
  z=1
  for j,k in Counter(p).items():z*=j**k*factorial(k)
  a,b=fixed_counts(p); raw+=F(a,z);core+=F(b,z)
 need(candidate.local_dp(n)[0][n]==raw,'literal edge-action Burnside raw n<=6')
 need(candidate.local_dp(n,core=True)[0][n]==core,'literal edge-action Burnside pair-free n<=6')
for d in range(41):
 for ell in range(1,32,2):
  for k in range(9):
   for excess in (0,1,3):
    q=d+ell*(k+excess)
    odd=lambda x:x//(x&-x)
    num=odd(factorial(q));den=odd(factorial(d))*ell**k*odd(factorial(k))
    need(num%den==0,'universal local transfer finite grid')
    need(k*(ell-1)//2+ell*k*(k-1)//2-(factorial(k)&-factorial(k)).bit_length()+1>=0,'nonnegative binary shift')
for bad in (-1,True,2.0,'2',None):
 try:candidate.local_dp(bad)
 except (ValueError,TypeError):need(True,'bad bound rejected')
 else:need(False,'bad bound rejected')
for mode in ('core','no_fixed'):
 for bad in (0,1,None,'yes'):
  try:candidate.local_dp(4,**{mode:bad})
  except (ValueError,TypeError):need(True,'bad mode rejected')
  else:need(False,'bad mode rejected')
need(candidate.local_dp(3,core=True,no_fixed=True)[0][3]==F(2,3),'observer remains rational mass')
print(json.dumps({'status':'PASS','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'checks':dict(checks),'total':sum(checks.values())},indent=2))
