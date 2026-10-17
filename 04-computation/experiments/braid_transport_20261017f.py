#!/usr/bin/env python3
"""Source-preserving odd transport after the two exact current folds.

Retains every primitive source and its actual accelerated shortcut clock.
No eventual ROOT census is used. Opposite signed weights cancel only at the
same primitive braid fibre, whose coalescence identity is proved separately.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np


def need(x,m):
 if not x:raise RuntimeError(m)
def v2(n):
 z=n.copy();a=np.zeros(len(n),dtype=np.int64)
 while True:
  even=z%2==0
  if not np.any(even):return a,z
  z[even]//=2;a[even]+=1

def primitive(n):
 s=n.copy();j=np.zeros(len(n),dtype=np.int64)
 while True:
  b=s%8==5
  if not np.any(b):return s,j
  s[b]=(s[b]-1)//4;j[b]+=1

def run(N,stages):
 sources=np.arange(1,4*N+1,2,dtype=np.int64);sources=sources[sources%8!=5]
 weights=np.zeros(len(sources),dtype=np.int64);m=sources.copy()
 while np.any(m<=4*N):
  inside=m<=2*N;rim=(m>2*N)&(m<=4*N)
  weights+=np.where(inside,np.where(m%3==1,-3,np.where(m%3==2,3,0)),np.where(rim,np.where(m%3==2,2,-1),0))
  m=4*m+1
 x=sources.copy();clock=np.zeros(len(x),dtype=np.int64);rows=[];peak=int(x.max());initial=int(np.abs(weights).sum())
 for depth in range(max(stages)+1):
  seed,j=primitive(x)
  # Every member B^j(1) reaches1 in an even number 2+2j of shortcut
  # steps (j=0 is stopped). Retain that literal ROOT clock.
  hit=(seed==1)&(x!=1);clock[hit]+=2+2*j[hit];x[hit]=1;j[hit]=0
  if depth in stages:
   signed=weights*np.where(clock%2,-1,1);ends,inv=np.unique(seed,return_inverse=True)
   total=np.zeros(len(ends),dtype=np.int64);np.add.at(total,inv,signed)
   root=int(total[ends==1].sum());res=int(np.abs(total[ends!=1]).sum())
   need(res+abs(root)<=initial,'signed transport never increases total variation')
   rows.append(dict(depth=depth,root_current=root,residual_l1=res,total_l1=res+abs(root),nonzero_terminals=int(np.count_nonzero(total)),
     certified_sources=int(np.count_nonzero(x==1)),max_actual_clock=int(clock.max()),
     source_endpoint_sha256=hashlib.sha256(x.tobytes()).hexdigest(),source_clock_sha256=hashlib.sha256(clock.tobytes()).hexdigest()))
  if depth==max(stages):break
  active=x!=1;before=x[active].copy();need(np.all(before<=(2**63-2)//3),'checked int64 arithmetic')
  a,y=v2(3*before+1);x[active]=y;clock[active]+=a;peak=max(peak,int(x.max()))
  # Primitive coordinate transition and exact fibre clock, independently.
  s,oldj=primitive(before);aa=(s%4==3).astype(np.int64)+2*(s%8==1).astype(np.int64)
  need(np.all((aa==1)|(aa==2)),'primitive seeds have odd-step cost1 or2')
  need(np.array_equal(a,aa+2*oldj) and np.array_equal(y,(3*s+1)//(1<<aa)),'skew product retains the actual integer clock')
 return dict(N=N,sources=len(sources),initial_l1=initial,peak_odd_endpoint=peak,source_sha256=hashlib.sha256(sources.tobytes()).hexdigest(),rows=rows)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--N',nargs='+',type=int,default=[1024,16384,262144]);p.add_argument('--output',required=True);a=p.parse_args()
 r=dict(status='FINITE-EXACT marked signed transport; uniform contraction OPEN',runs=[])
 for N in a.N:
  row=run(N,[0,1,2,4,8,12,20,28,40,64]);r['runs'].append(row);print(json.dumps(row,sort_keys=True),flush=True)
  Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
