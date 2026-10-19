#!/usr/bin/env python3
"""Exact finite reduction certificate for all one-sided exteriors of one core.
Use --table for the complete small profile/count table, otherwise JSON receipt.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
import hashlib,json,sys
checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
Q=(0,0,0,6,0,19,2,18,86,151);roots=(0,3);n=10;k=7
def inc(p,a,b):return a!=b and not(p[a]>>b&1 or p[b]>>a&1)
delta=tuple(k-sum(inc(Q,a,b) for b in range(n)) for a in range(n))
require(delta==(0,4,1,0,2,1,0,1,2,3),'core slack')
def closure(n,ed):
 p=[0]*n
 for a,b in ed:p[b]|=1<<a
 for _ in range(n):
  old=p.copy()
  for b in range(n):
   for a in range(n):
    if p[b]>>a&1:p[b]|=p[a]
  if p==old:break
 require(all(not(p[a]>>a&1) for a in range(n)),'acyclic model')
 return tuple(p)
def count(p,orders=()):
 p=list(p)
 for a,b in orders:p[b]|=1<<a
 full=(1<<len(p))-1
 @lru_cache(None)
 def E(s):
  if s==full:return 1
  return sum(E(s|1<<a) for a in range(len(p)) if not(s>>a&1) and not(p[a]&~s))
 return E(0)
def model(side,pr):
 h=max(pr);chain=tuple(range(n,n+h));ed=[(a,b) for a in range(n) for b in range(n) if Q[b]>>a&1]+list(zip(chain,chain[1:]))
 ed += [(c,q) for i,c in enumerate(chain) for q in range(n) if i<h-pr[q]] if side=='lower' else [(q,c) for i,c in enumerate(chain) for q in range(n) if i>=pr[q]]
 return closure(n+h,ed)
def explicit_profiles(side):
 if side=='lower':return sorted((0,a,b,0,c,0,0,d,0,0) for a in range(5) for b in range(2) for c in range(3) for d in range(min(1,a,c)+1))
 return sorted((0,0,0,0,c,e,0,f,g,h) for c in range(2) for e in range(c,2) for f in range(c,2) for g in range(c,3) for h in range(f,4))
def general_profiles(side):
 eligible=[q for q in range(n) if q not in roots and (not(any(Q[q]>>r&1 for r in roots)) if side=='lower' else not(any(Q[r]>>q&1 for r in roots)))]
 out=[]
 for vals in product(*(range(delta[q]+1) for q in eligible)):
  pr=[0]*n
  for q,v in zip(eligible,vals):pr[q]=v
  if sum(v>0 for v in pr)>k:continue
  if any((pr[a]<pr[b] if side=='lower' else pr[a]>pr[b]) for a in range(n) for b in range(n) if Q[b]>>a&1):continue
  out.append(tuple(pr))
 return sorted(out)
table=['side\tprofile\tvertices\textensions\tTT\tTF\tFT\tFF\tminimum_margin'];reports=[]
for side,expected_size,first,second,sharp,simple in [('lower',46,(8,9),(3,8),F(5031,38734),F(1,8)),('upper',48,(2,4),(1,4),F(95,796),F(1,9))]:
 profiles=explicit_profiles(side);require(len(profiles)==expected_size,'profile universe size');require(profiles==general_profiles(side),'explicit universe equals all allowable profiles')
 worst=F(1);witness=None;min_slack=None
 for pr in profiles:
  p=model(side,pr);E=count(p);cc=[count(p,(a,b)) for a,b in product((first,first[::-1]),(second,second[::-1]))]
  require(sum(cc)==E,'joint cells partition ambient extensions')
  require(all(sum(inc(p,a,b) for b in range(len(p)))<=k for a in range(len(p))),'all profile models seven-thin')
  require(tuple(v&1023 for v in p[:n])==Q,'same induced core')
  require(set(roots)|{a for a in range(len(p)) if any(inc(p,a,r) for r in roots)}==set(range(n)),'same full critical-root core')
  mm=[E-sum(cc[:2])-cc[j] if j<2 else E-sum(cc[2:])-cc[j] for j in range(4)]
  require(all(F(v,E)>=sharp for v in mm),'sharp normalized Gold margin')
  require(all(F(v,E)>simple for v in mm),'simple strict Gold margin')
  ratio=F(min(mm),E);slack=min(v*sharp.denominator-E*sharp.numerator for v in mm)
  if ratio<worst:worst=ratio;witness={'profile':pr,'ambient_predecessors':p,'extensions':E,'cells':cc,'margins':mm}
  min_slack=slack if min_slack is None else min(min_slack,slack)
  table.append('\t'.join([side,','.join(map(str,pr)),str(len(p)),str(E),*map(str,cc),str(min(mm))]))
 require(worst==sharp and min_slack==0,'claimed sharp model exists')
 reports.append({'side':side,'profile_count':len(profiles),'first_comparison':first,'second_comparison':second,'sharp_margin':str(sharp),'simple_strict_margin':str(simple),'minimum_integer_sharp_slack':min_slack,'sharp_witness':witness})
text='\n'.join(table)+'\n'
if '--table' in sys.argv:print(text,end='')
else:print(json.dumps({'status':'PASS','guards':checks,'core_predecessors':Q,'roots':roots,'degree_bound':k,'slack':delta,'cases':reports,'profile_table_sha256':hashlib.sha256(text.encode()).hexdigest(),'scope':'Exact full finite profile universes for one fixed core. The compression theorem extends these linear Gold margins to arbitrary-size active exteriors on one side with the opposite exterior comparable to every core vertex. No two-sided arbitrary-exterior claim, morphology exclusion, or general conjecture resolution.'},indent=2))
