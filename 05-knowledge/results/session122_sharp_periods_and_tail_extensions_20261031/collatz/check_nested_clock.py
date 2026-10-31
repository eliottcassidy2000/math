#!/usr/bin/env python3
"""Necessary nested-return grammar, never a word realizability certificate."""
from fractions import Fraction as F
import json

def need(x,m):
    if not x: raise ValueError(m)

def floor_sum(n,m,a,b):
    ans=0
    while True:
        qa,a=divmod(a,m);qb,b=divmod(b,m)
        ans+=qa*n*(n-1)//2+qb*n
        y=a*n+b
        if y<m:return ans
        n,b=divmod(y,m);m,a=a,m

def count(N,M,a,b,r,step=8):
    first=r if r else step
    n=max(0,(N-first)//step+1)
    return (floor_sum(n,M,a*step,a*first+b-1)
            -floor_sum(n,M,(a-1)*step,(a-1)*first+b))

tests=0
for M in range(2,41):
 for a in (0,1,M//2,M-1,M+3):
  for b in (-M-1,-1,0,1,M//2,M+1):
   for N in (0,1,M//2,M-1):
    for r in range(8):
     literal=sum(1 for A in range(1,N+1) if A%8==r and 1<=(a*A+b)%M<A)
     need(count(N,M,a,b,r)==literal,'floor sum')
     tests+=1
eta=F(2832031,781250);gamma=eta-1;D=8*5**10;T=2*5**9
need(D*(F(21,8)-gamma)==25,'small cap')
small=0
for A in range(1,26):
 for L in range(1,A):
  need(pow(3,L,5**10)!=pow(2,A,5**10),'small multiplier');small+=1
need(2*F(22,8)-1>3,'two cores excluded')
q=5**14;M=4*5**13;log2=2883827743
need(pow(3,log2,q)==2,'log2')
need(pow(3,M,q)==1 and pow(3,M//2,q)!=1 and pow(3,M//5,q)!=1,'order')
h=pow(3*pow(2,-1,q)%q,T,q)
rows=[]
for j,oldcap in ((23,48828467),(24,97656592)):
 for a in range((j+8-22)//2+1):
  jc=j+8-2*a
  cap_rational=D*(F(jc,8)-gamma)
  need(cap_rational.denominator==1,'integer cap')
  cap=min(int(cap_rational),oldcap-a*T)
  residue=(5*jc)%8
  b=(a*T*(log2-1))%M
  target=pow(h,-a,q)
  need(pow(3,b,q)==target,'core multiplier')
  total=count(cap,M,log2,b,residue)
  def select(rank):
   lo,hi=1,cap
   while lo<hi:
    mid=(lo+hi)//2
    if count(mid,M,log2,b,residue)>=rank:hi=mid
    else:lo=mid+1
   length=(log2*lo+b)%M
   need(lo%8==residue and 1<=length<lo,'selection')
   need(pow(3,length,q)==target*pow(2,lo,q)%q,'selected multiplier')
   need(count(lo,M,log2,b,residue)==rank and count(lo-1,M,log2,b,residue)==rank-1,'selected rank')
   return [lo,length]
  placements=[p for p in range(a+1) if p<5 and a-p<5]
  # Exact affine composition for every permitted prefix/suffix placement.
  for p in placements:
   hp=pow(h,p,q);hs=pow(h,a-p,q)
   m=hs*target*hp%q
   t=(hs*(target*(hp-1)+(target-1))+hs-1)%q
   need((m,t)==(1,0),'affine composition')
  rows.append(dict(j11=j,all_one10_blocks=a,j10=jc,
       core_clock_cap=cap,core_clock_mod8=residue,
       prefix_block_counts=placements,required_core_map=[target,(target-1)%q],
       pairs=total,first_pair=select(1) if total else None,
       last_pair=select(total) if total else None))
print(json.dumps(dict(status='PASS',scope='conditional necessary grammar; E14 remains open',
 inherited_E9_R2_rerun=False,floor_sum_controls=tests,small_clock_pairs=small,
 primitive10_quantized_H='11/4',patterns=rows,
 grammar_forms=sum(len(r['prefix_block_counts']) for r in rows),
 multiplier_pairs_without_placement=sum(r['pairs'] for r in rows)),indent=2,sort_keys=True))
