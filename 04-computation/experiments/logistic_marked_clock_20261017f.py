#!/usr/bin/env python3
"""Exact 2-adic logistic conjugacy tables with marked integer-clock controls.

This is a representation check, not a convergence or Haar-to-source theorem.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from fractions import Fraction


def need(x,msg):
 if not x:raise RuntimeError(msg)
def parity_table(k,logistic):
 x=np.arange(1<<k,dtype=np.int64);q=np.zeros_like(x)
 for j in range(k):
  q|=(x&1)<<j
  x=(x*(1-x)//2 if logistic else np.where(x&1,(3*x+1)//2,x//2))% (1<<(k-j-1))
 return q

def run(K):
 previous=None;rows=[]
 for k in range(1,K+2):
  n=1<<k;qt=parity_table(k,False);ql=parity_table(k,True)
  need(np.array_equal(np.sort(qt),np.arange(n)) and np.array_equal(np.sort(ql),np.arange(n)),'both itinerary maps are permutations')
  inv=np.zeros(n,dtype=np.int64);inv[ql]=np.arange(n);H=inv[qt]
  need(np.array_equal(ql[H],qt),'conjugacy preserves all k parity bits')
  if previous is not None:
   m=n//2;src=np.arange(m,dtype=np.int64)
   need(np.array_equal(H[:m]%m,previous),'conjugacy survives refinement')
   target=np.where(src&1,(3*src+1)//2,src//2)%m
   need(np.array_equal((H[:m]*(1-H[:m])//2)%m,previous[target]),'every integer shortcut clock step survives conjugacy')
  need(all((int(H[a])**2-3*int(H[a])+6)%n==0 for a in (1,2%n)) if n>2 else True,'ROOT cycle polynomial')
  need(H[0]==0 and H[-1]==n-1,'zero and negative fixed point')
  if k>=2:
   source=(-pow(3,-1,n))%n
   need(H[source]==1,'braid accumulation point is the logistic preimage of zero')
  odd=np.arange(1,n,2);B=(4*odd+1)%n
  need(np.array_equal(ql[H[B]],(4*ql[H[odd]]-3)%n),'braid inserts exactly two even clocks')
  rows.append(dict(bits=k,states=n,ROOT_images=[int(H[1]),int(H[2%n])],conjugacy_sha256=hashlib.sha256(H.tobytes()).hexdigest()))
  if k==K:table=H.copy()
  previous=H
 # Exact rooted integer paths are a finite positive control. Their endpoints,
 # all intermediate integers and actual shortcut clocks remain attached.
 paths=0;steps=0;longest=0;digest=hashlib.sha256();mod=1<<K
 for source in range(1,2049):
  x=source;clock=0
  while x!=1:
   y=x//2 if x%2==0 else (3*x+1)//2
   hx=int(previous[x%(2*mod)])
   need((hx*(1-hx)//2)%mod==int(table[y%mod]),'marked path clock check')
   digest.update(f'{source},{clock},{x},{y};'.encode());x=y;clock+=1
   need(clock<10000,'finite positive control failed to reach ROOT within declared cap')
  paths+=1;steps+=clock;longest=max(longest,clock)
 # An algebraic identity, independently evaluated with rational arithmetic.
 for x in map(Fraction,range(-20,21)):
  L=lambda z:z*(1-z)/2
  need(L(L(x))-x==(L(x)-x)*(x*x-3*x+6)/4,'logistic period-two polynomial identity')
 return dict(status='FINITE-EXACT conjugacy, refinement and clock checks; uniform current OPEN',max_source_bits=K,rows=rows,
  rooted_control_sources=paths,rooted_control_steps=steps,maximum_control_clock=longest,clock_sha256=digest.hexdigest(),
  ROOT_cycle_polynomial=[6,-3,1],ROOT_cycle_discriminant=-15,braid_itinerary='Q(Bn)=4 Q(n)-3 for odd n',braid_fixed_point='H(-1/3)=1')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--bits',type=int,default=16);p.add_argument('--output',required=True);a=p.parse_args()
 r=run(a.bits);Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(r,sort_keys=True))
