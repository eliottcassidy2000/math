#!/usr/bin/env python3
"""Exact star-bipartition/independent-phase optimization, m<=5."""
import itertools as it,collections,json
from pathlib import Path

def optimize(m,q):
 es=list(it.combinations(range(m),2));E=len(es);idx={frozenset(e):k for k,e in enumerate(es)}
 ctrl=[v if q>>k&1 else u for k,(u,v) in enumerate(es)]
 stars=[[e for e,c in enumerate(ctrl) if c==v] for v in range(m)]
 free=[e for star in stars for e in star[1:]]
 terms=[]
 for p in it.permutations(range(m)):
  wedges=[]
  for j in range(1,m-1):
   a=idx[frozenset((p[j-1],p[j]))];b=idx[frozenset((p[j],p[j+1]))]
   if ctrl[a]==p[j]==ctrl[b]:wedges.append((a,b))
  boundary=None
  if m>1:
   e=idx[frozenset(p[-2:])]
   if ctrl[e]==p[-1]:boundary=(p[-1],e)
  terms.append((wedges,boundary))
 best=-1;bestbases=[]
 for mask in range(1<<len(free)):
  signs=[1]*E
  for k,e in enumerate(free):signs[e]=-1 if mask>>k&1 else 1
  A=0;B=[0]*m
  for wedges,boundary in terms:
   w=1
   for a,b in wedges:
    if signs[a]==signs[b]:w=0;break
    w*=2
   A+=w
   if boundary:
    v,e=boundary;B[v]+=signs[e]*w
  f=A+sum(abs(x) for x in B)
  if f>best:best=f;bestbases=[{'base':mask,'A':A,'B':B}]
  elif f==best:bestbases.append({'base':mask,'A':A,'B':B})
 return best,bestbases,len(free)

def main():
 result=[]
 for m in range(1,6):
  hist=collections.Counter();witnesses=[];best=0;partitions=0
  for q in range(1<<(m*(m-1)//2)):
   f,bases,d=optimize(m,q);hist[f]+=1;partitions+=1<<d
   if f>best:best=f;witnesses=[q]
   elif f==best:witnesses.append(q)
  result.append({'m':m,'max':best,'max_per_Q_histogram':dict(sorted(hist.items())),'bipartitions_examined':partitions,'maximizing_Q':witnesses})
 result.append({'m5_Q2_detail':optimize(5,2)})
 Path(__file__).with_name('phase_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
