#!/usr/bin/env python3
"""Exact controls for finite absolute-difference strips; no prime inference."""
from itertools import product
from collections import Counter
import json

def check(p,msg):
    global checks
    checks+=1
    if not p: raise RuntimeError(msg)

def step(x,u):return tuple(abs(a-b) for a,b in zip(x,x[1:]+(u,)))
def apex(x):
    while len(x)>1:x=tuple(abs(a-b) for a,b in zip(x,x[1:]))
    return x[0]
checks=0
rigidity=0
for M in range(2,6):
 for n in range(1,7):
  for x in product(range(M+1),repeat=n):
   z=apex(x)
   check(z!=M or all(a in(0,M) for a in x),('rigidity',x,z))
   rigidity+=1
reset=[]
for M in range(2,7):
 for F in range(1,8):
  if (M+1)**F>350000:continue
  worst=0;witness=None
  for x in product(range(M+1),repeat=F):
   y=x;t=0
   for _ in range(F):y=step(y,1)
   check(max(y)<M,('max_drop',x,y))
   y=x
   while max(y)>1:
    y=step(y,1);t+=1
    check(t<=(M-1)*F,('entry',x,t))
   if t>worst:worst=t;witness=x
  reset.append(dict(M=M,F=F,max_entry=worst,witness=witness,bound=(M-1)*F))
single_hit_checks=0
for M in range(2,6):
 for F in range(1,6):
  for x in product(range(M+1),repeat=F):
   for tail in product((0,1),repeat=F-1):
    y=x
    for u in (1,)+tail:y=step(y,u)
    check(max(y)<M,('single_hit',M,F,x,tail,y))
    single_hit_checks+=1
for F in (1,2,4,8,16,32):
 for M in range(2,20):
  x=(M,)+(0,)*(F-1);t=0
  while max(x)>1:x=step(x,1);t+=1
  check(t==(M-1)*F,('dyadic_sharpness',M,F,t))
mixing=[]
for F in range(1,8):
 for x in product((0,1),repeat=F):
  outputs=Counter()
  for us in product((0,1),repeat=F):
   y=x
   for u in us:y=step(y,u)
   outputs[y]+=1
  check(len(outputs)==2**F and set(outputs.values())=={1},('uniform',F,x))
 # Every fixed word preserves all distinctions among binary states.
 for us in product((0,1),repeat=min(F,4)):
  ys=set()
  for x in product((0,1),repeat=F):
   y=x
   for u in us:y=step(y,u)
   ys.add(y)
  check(len(ys)==2**F,('bijection',F,us))
 mixing.append(F)
print(json.dumps(dict(status='PASS',checks=checks,rigidity_states=rigidity,single_hit_checks=single_hit_checks,reset=reset,exact_mixing_widths=mixing),indent=2))
