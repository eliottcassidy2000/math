#!/usr/bin/env python3
"""Exact three-strand compiler for V at an internal root of a final block."""
from itertools import product
import json
import core as v
D=(18,8,3,5,14); A=v.attach(D,9)
def table(T,case=None,r=1):
 H,end=v.hp_subsets(T);RH,start=v.hp_subsets(v.reverse(T));N=len(T);out=[]
 for labels in product(range(3),repeat=N-1):
  masks=[0]*3
  for i,j in enumerate(labels,1):masks[j]|=1<<i
  if case is not None and not masks[case]>>r&1:continue
  options=[]
  for j,S in enumerate(masks):
   if not S:options.append([(0,1)]);continue
   total=H[S]
   if case==j:total=end(S,r)
   if j==0:
    good=sum(end(S,i) for i in range(1,N) if S>>i&1 and T[i]&1)
    if case==j:good=total if T[r]&1 else 0
    options.append([(1,good),(2,total-good)])
   elif j==1:
    # Q must start at an outneighbor of z; if it must end r use direct endpoint DP.
    good=0
    for i in range(1,N):
     if S>>i&1 and T[0]>>i&1:
      if case!=j:good+=start(S,i)
      else:
       ids=[x for x in range(N) if S>>x&1]
       dp={(1<<i,i):1}
       for size in range(1,len(ids)):
        for (mask,last),count in list(dp.items()):
         if mask.bit_count()!=size:continue
         for x in ids:
          if not mask>>x&1 and T[last]>>x&1:dp[(mask|1<<x,x)]=dp.get((mask|1<<x,x),0)+count
       good+=dp.get((S,r),0)
    options.append([(1,good),(2,total-good)])
   else:options.append([(1,total)])
  for z in product(*options):
   w=z[0][1]*z[1][1]*z[2][1]
   if w:out.append((tuple(x[0] for x in z),w))
 return out
def step(states,tab):
 out={}
 for (p,q,r),n in states.items():
  for (a,b,c),w in tab:
   t=(a or p,q or b,int(r or c));out[t]=out.get(t,0)+n*w
 return out
def finish(states,tab):
 return sum(n for (p,q,r),n in step(states,tab).items() if p!=2 and q!=2)
