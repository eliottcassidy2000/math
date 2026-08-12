#!/usr/bin/env python3
"""THM-3350 uniform high-forest tail for non-dense connected-low six-rays.

For each primitive high channel, this exact referee maximizes the rigorous
all-channel midpoint error over the 180 ordered endpoint-label contexts that
occur in the 649 reflected bodies (using the smallest ruler in each context).
Positive certified gains 1/105-error are assembled by maximum-forest
Kruskal.  The forest extends to a complete-graph spanning tree.  Comparison
with the exact universal singleton-debt envelope proves a body- and
labelling-uniform physical tail.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from math import gcd,lcm
from pathlib import Path
import os

ROOT=Path(os.environ.get('LRC_REPO',Path(__file__).resolve().parents[1])).resolve()
SELECTOR=ROOT/'04-computation/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.py'
EXPECTED_SELECTOR='cd8b08087f0f7e1e0c0c7d0be673629c7c2702c170c5c1e771e1d76df1d3cd1c'
EXPECTED_CONTEXT='ee75c24771753b676cbc296a283c5f4efd87053f7266e716cf3ac36d8a4cedaf'
EXPECTED_HIST=((1,101280),(2,159228),(3,36631),(4,7204),(5,1046),(6,175),(7,329),(8,7),(11,7))
EXPECTED_SEMANTIC='da1f0ec439307701c0773368f9dd2c343dcc274b6c18ae89f81ff62266829f8a'
HIGH_FLOOR=F(1,105)
DEBT_MAX=F(186636088362,11773143757375)
DENSE={(1,2,3,4,6,12),(2,3,4,6,8,12)}
UP=frozenset((F(4,3),F(3,2),F(2),F(5,2),F(3),F(4),F(5),F(6)))
RATIOS=tuple(sorted(UP|{1/x for x in UP}))

def req(x,d):
 if not x:raise RuntimeError(d)

def load(path):
 s=spec_from_file_location('tail_selector',path);req(s and s.loader,path)
 m=module_from_spec(s);s.loader.exec_module(m);return m

req(sha256(SELECTOR.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==EXPECTED_SELECTOR,('selector',SELECTOR))
SEL=load(SELECTOR)

def low(x,y):
 a,b=sorted((x,y));return b/a in UP

def connected(row):
 unseen=set(row[1:]);front={row[0]}
 while front:
  x=front.pop();hit={y for y in unseen if low(x,y)};unseen-=hit;front|=hit
 return not unseen

def parent(row):
 for y in reversed(row[1:]):
  c=tuple(x for x in row if x!=y)
  if connected(c):return c
 raise RuntimeError(row)

def primitive(row):
 d=lcm(*(x.denominator for x in row));a=tuple(int(x*d) for x in row);g=gcd(*a)
 return tuple(x//g for x in a)

def rows():
 layer=[(F(1),)]
 for size in range(2,7):
  out=[]
  for row in layer:
   have=set(row)
   for y in sorted({x*r for x in row for r in RATIOS if x*r>=1 and x*r not in have}):
    child=tuple(sorted((*row,y)))
    if parent(child)==row:out.append(child)
  req(len(out)==len(set(out)),('duplicates',size));layer=out
 ans=sorted(primitive(r) for r in layer);req(len(ans)==len(set(ans))==305909,len(ans));return ans

def high_edges(row):
 out=[]
 for i in range(6):
  for j in range(i+1,6):
   d=gcd(row[i],row[j]);P,Q=row[i]//d,row[j]//d
   if P+Q>=8:out.append((i,j,d,P,Q))
 return out

def gamma(k):return F(1,2) if k%2==0 else F(k*k+1,2*k*k)

contexts={}
for body,L in SEL.MS.body_universe():
 for e in body:
  for f in body:
   if e!=f:contexts[e,f]=min(contexts.get((e,f),10**100),L)
req(len(contexts)==180,len(contexts))
ch=sha256()
for k,v in sorted(contexts.items()):ch.update((repr((k,v))+'\n').encode())
req(ch.hexdigest()==EXPECTED_CONTEXT,ch.hexdigest())

@cache
def error(scale,d,P,Q):
 h=scale*d;vals=[]
 for (e,f),L in contexts.items():
  C=Q*e-P*f
  vals.append(gamma(P)*F(e*P,h*L*P-e)+gamma(Q)*F(f*Q,h*L*Q-f)
              +F(abs(C)*(abs(C)//L+1),2*h*h*L*P*Q))
 return max(vals)

def forest(row,scale):
 par=list(range(6));chosen=[];total=F()
 def root(x):
  while par[x]!=x:par[x]=par[par[x]];x=par[x]
  return x
 weighted=[]
 for i,j,d,P,Q in high_edges(row):
  E=error(scale,d,P,Q);weighted.append((HIGH_FLOOR-E,(i,j,d,P,Q),E))
 for gain,edge,E in sorted(weighted,reverse=True):
  if gain<=0:break
  a,b=root(edge[0]),root(edge[1])
  if a!=b:par[a]=b;total+=gain;chosen.append((gain,edge,E))
 return total,tuple(chosen)

def threshold(row):
 for scale in range(1,100):
  total,T=forest(row,scale)
  if total>DEBT_MAX/scale:return scale,T,total
 raise RuntimeError(('threshold',row))

def main():
 hist=Counter();sem=sha256();maxrows=[];atlas=rows();highhist=Counter()
 for row in atlas:
  H=high_edges(row);highhist[len(H)]+=1
  if row in DENSE:req(len(H)==1,(row,H));continue
  scale,T,total=threshold(row);hist[scale]+=1
  rec=(scale,row,T,total,total-DEBT_MAX/scale)
  if not maxrows or scale>maxrows[0][0]:maxrows=[rec]
  elif scale==maxrows[0][0]:maxrows.append(rec)
  sem.update((repr(rec)+'\n').encode())
 actual=tuple(sorted(hist.items()));req(actual==EXPECTED_HIST,actual)
 req(sem.hexdigest()==EXPECTED_SEMANTIC,sem.hexdigest());req(len(maxrows)==7 and maxrows[0][0]==11,maxrows)
 print('LRC14 CONNECTED-LOW UNIFORM HIGH-FOREST TAIL')
 print('status=PROVED analytic transport + VERIFIED-EXACT intrinsic/body-context atlases')
 print('selector_sha256',EXPECTED_SELECTOR)
 print('atlas',len(atlas),'non_dense',sum(hist.values()),'dense_exceptions',len(DENSE))
 print('high_edge_count_hist',tuple(sorted(highhist.items())))
 print('body_endpoint_contexts',len(contexts),'context_semantic_sha256',ch.hexdigest())
 print('threshold_hist',actual)
 print('shape_scale_heads',sum((s-1)*c for s,c in hist.items()))
 print('max_threshold',maxrows[0][0],'max_count',len(maxrows),'max_rows',tuple(maxrows))
 print('semantic_sha256',sem.hexdigest())
 print('conclusion=every non-dense connected-low ray closes on every body/labelling for all common scales s>=11')

if __name__=='__main__':main()
