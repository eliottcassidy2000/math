#!/usr/bin/env python3
"""Exact interval-union tails and source-owned weighted Voronoi receipts."""
import argparse,hashlib,json
from array import array
from itertools import product
from math import isqrt
from pathlib import Path
import collatz_interval_atlas_20261016e as A
from collatz_clock_tail_20261017b import tail
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def union_tail(X,anchors,t):
 intervals=[]
 for a,B,stage in anchors:
  if B<t:
   r=(t-1-B)//4;l=max(X,a-r);h=min(2*X-1,a+r)
   if l<=h:intervals.append((l,h))
 covered=0;end=X-1
 for l,h in sorted(intervals):
  if h>end:covered+=h-max(l,end+1)+1;end=h
 return X-covered

def cells(X,anchors):
 # Ordered ties choose the smallest literal anchor, then its stage.
 anchors=sorted(anchors);inf=10**100;values=[inf]*X;owner=array('i',[-1])*X
 def better(i,v,a):
  if v<values[i] or v==values[i] and (owner[i]<0 or a<owner[i]):values[i]=v;owner[i]=a
 for a,(n,B,stage) in enumerate(anchors):better(n-X,B,a)
 for i in range(1,X):
  if owner[i-1]>=0:better(i,values[i-1]+4,owner[i-1])
 for i in range(X-2,-1,-1):
  if owner[i+1]>=0:better(i,values[i+1]+4,owner[i+1])
 need(values==A.cone_envelope(X,anchors),'owner pass agrees with independent scalar envelope')
 pieces=[];lo=0
 def powers(n):return n*(n+1)//2,n*(n+1)*(2*n+1)//6
 while lo<X:
  j=owner[lo];hi=lo+1
  while hi<X and owner[hi]==j:hi+=1
  a,B,stage=anchors[j];L=X+lo;R=X+hi
  need(L<=a<R,'every live cell contains its original anchor')
  l1,l2=powers(a-L);r1,r2=powers(R-1-a);q=(R-L)*B*B;d=l2+r2
  pieces.append([L,R,a,B,stage,q,l1+r1,d,q+8*B*(l1+r1)+16*d]);lo=hi
 need(all(a[2]<b[2] for a,b in zip(pieces,pieces[1:])),'cell owners strictly increase')
 E=sum(v*v for v in values);Q=sum(p[5] for p in pieces);R2=sum(p[7] for p in pieces)
 need(sum(p[8] for p in pieces)==E and Q<=E and 16*R2<=E,'exact graded cell identity and nonnegative constituent bounds')
 return values,dict(cells=pieces,energy=E,threshold_energy=Q,distance_energy=R2)

def controls():
 universes=0
 for X in range(1,6):
  for raw in product((-1,0,3,8),repeat=X):
   anchors=[[X+i,b,0] for i,b in enumerate(raw) if b>=0]
   if not anchors:continue
   v,c=cells(X,anchors);universes+=1
   need(v==[min(B+4*abs(n-a) for a,B,s in anchors) for n in range(X,2*X)],'exhaustive independent direct minimum')
   for t in range(max(v)+2):need(union_tail(X,anchors,t)==sum(z>=t for z in v),'exact strict integer interval endpoints')
 hostiles=[]
 for m in (1,2,3):
  X=1<<(6*m);s=1<<(4*m);v=[isqrt(X*X//(j+s)) for j in range(X)];q=tail(v)
  need(all(abs(a-b)<=1 for a,b in zip(v,v[1:])),'logarithmic hostile is one-Lipschitz')
  need(q['weak_numerator']<=X*X,'logarithmic hostile has critical weak constant at most one')
  hostiles.append(dict(m=m,X=X,s=s,energy=q['energy'],weak_numerator=q['weak_numerator']))
 return dict(exhaustive_small_anchor_profiles=universes,log_hostiles=hostiles)

def run(path):
 raw=Path(path).read_bytes();obj=json.loads(raw);out=[]
 for r in obj['runs']:
  X=r['X']
  for stage in r['stages']:
   for record in stage['selections']:
    anchors=record['anchors'];v,c=cells(X,anchors);q=record['cone_tail'];need(c['energy']==q['energy'],'native cone energy retained')
    heights=sorted({1,q['critical_height'],q['maximum']+1}|{1<<j for j in range(max(v).bit_length())})
    hist=q['histogram'];samples=[]
    for t in heights:
     u=union_tail(X,anchors,t);need(u==sum(count for value,count in hist if value>=t),'interval coverage matches independent saved histogram')
     samples.append([t,u,t*t*u])
    kd=max((score for t,u,score in samples if t>0 and t&(t-1)==0),default=0)
    need(kd<=q['weak_numerator']<=4*kd,'dyadic threshold tests control every height within factor four')
    out.append(dict(k=r['k'],stage=stage['stage'],bin_width=record['bin_width'],graded=c,tail_queries=samples,dyadic_weak_numerator=kd))
 return dict(source_sha256=hashlib.sha256(raw).hexdigest(),records=out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--receipt',nargs='*',default=[]);a=p.parse_args();control=controls();runs=[run(p) for p in a.receipt]
 print(json.dumps(dict(status='PROVED elementary extraction with FINITE-EXACT controls; uniform native arithmetic bound OPEN',controls=control,runs=runs,checks=checks),sort_keys=True,separators=(',',':')))
