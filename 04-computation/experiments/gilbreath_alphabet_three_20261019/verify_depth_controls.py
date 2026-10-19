#!/usr/bin/env python3
from itertools import product
from pathlib import Path
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def D(w):return [abs(x-y)for x,y in zip(w,w[1:])]
def triangle(w):
 while w:
  yield w
  w=D(w)
def depth(w):
 for r,row in enumerate(triangle(w)):
  if r and row[0]==1 and max(row[1:],default=0)<=2:return r
 return None
# Fixed maxgap, all-unit edge, arbitrarily slow sea-entry family.
for N in range(5,1001):
 seq=[2,3,5,9,15]+list(range(17,17+2*(N-5),2))
 need(len(seq)==N,'family length')
 rows=list(triangle(seq))
 need(all(r[0]==1 for r in rows[1:]),'family unsafe edge')
 need(max(D(seq))==6 and depth(seq)==N-3,'family depth')
# Safe later row does not certify earlier edge.
rows=list(triangle([2,3,7,13]))
need([r[0]for r in rows]==[2,1,3,1],'recovery hostile')
need(depth([2,3,7,13])==3,'recovery depth')
# Independent standard-library prime sieve and scalar CA replay.
M=400000;sieve=bytearray(b'\1')*(M+1);sieve[0:2]=b'\0\0'
for p in range(2,int(M**.5)+1):
 if sieve[p]:sieve[p*p:M+1:p]=b'\0'*(((M-p*p)//p)+1)
primes=[p for p in range(M+1)if sieve[p]]
records=[]
for N in [1000,2814,31569]:
 seed=primes[:N];dep=depth(seed);G=max(D(seed));records.append({'N':N,'pN':seed[-1],'depth':dep,'maxgap':G})
need(records==[{'N':1000,'pN':7919,'depth':35,'maxgap':34},{'N':2814,'pN':25561,'depth':51,'maxgap':52},{'N':31569,'pN':370609,'depth':69,'maxgap':112}],'prime replay')
print(json.dumps({'family_exact_checks':996,'family_formula':'k(N)=N-3, maxgap6, allpositive leftedges1','recovery_prime_seed':[2,3,7,13],'prime_scalar_replays':records},indent=2,sort_keys=True))
