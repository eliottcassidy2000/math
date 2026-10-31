#!/usr/bin/env python3
"""Independent exact prime-coordinate covariance tests. Standard library only."""
import json,math
from fractions import Fraction
P=(2,3,5,7,11,13,17,19,23,29)
Q=(31,37,41,43,47,53,59,61,67,71)
F=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def require(x,msg):
 if not x:raise RuntimeError(msg)
def vec(n):
 v=[]
 for p in P:
  e=0
  while n%p==0:n//=p;e+=1
  v.append(e)
 require(n==1,'outside active monoid');return v
def enc(n):return math.prod(p**e for p,e in zip(Q,vec(n)))
def step(n,fs):
 for i,(a,b) in enumerate(fs):
  if n%b==0:return i,n//b*a
 return None,n
G=tuple((enc(a),enc(b)) for a,b in F)
checks=0
for seed in range(1,501):
 # This finite test uses active-monoid seeds only. A global extension
 # permutes P union Q and may fix primes outside that union.
 n=seed
 for p in P:
  while n%p==0:n//=p
 if n!=1:continue
 n=seed;m=enc(n)
 for t in range(100):
  i,n=step(n,F);j,m=step(m,G)
  require(i==j and m==enc(n),'conjugacy');checks+=1
  if i is None:break
N=enc(55)
r=math.isqrt(8*N+1)
a,b=0,1
while b<N:a,b=b,a+b
require(r*r!=8*N+1 and b!=N,'figurate hostile')
# Reordering the same fractions changes the first step at 2: L versus N.
require(step(2,F)[1]==15 and step(2,tuple(reversed(F)))[1]==110,'priority hostile')
print(json.dumps({'status':'passed','paired_steps':checks,'prime_image_of_55':N,'image_is_triangular':False,'image_is_fibonacci':False,'priority_hostile_seed':2,'original_next':15,'reversed_next':110,'renamed_fractions':G},sort_keys=True,indent=2))
