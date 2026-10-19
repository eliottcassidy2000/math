#!/usr/bin/env python3
"""Exact certificate: no common first core pair, with adaptive Gold repair.
Only integer ideal DP and rational arithmetic; no numerical solver or dependencies.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations,product
import json
checks=0
def require(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
def closure(n,edges):
 p=[0]*n
 for a,b in edges:p[b]|=1<<a
 for _ in range(n):
  old=p.copy()
  for b in range(n):
   for a in range(n):
    if p[b]>>a&1:p[b]|=p[a]
  if p==old:break
 require(all(not(p[a]>>a&1) for a in range(n)),'acyclic')
 return tuple(p)
def inc(p,a,b):return a!=b and not(p[a]>>b&1 or p[b]>>a&1)
def count(p,orders=()):
 p=list(p)
 for a,b in orders:p[b]|=1<<a
 full=(1<<len(p))-1
 @lru_cache(None)
 def E(s):
  if s==full:return 1
  return sum(E(s|1<<x) for x in range(len(p)) if not(s>>x&1) and not(p[x]&~s))
 return E(0)
Q=(0,0,0,6,0,19,2,18,86,151);roots=(0,3)
profiles=[([0,1,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0]),([0,1,1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,1,1,2])]
expected_E=[29193,84707]
expected_N=[[5740,13424,25049,8066,19671,21775,27596,22912,17187,9068,26602,21015,22807,1557,15736,7920,8917,20750,19674,24774,6959,7461,19752,18789,16327,25258,25090,12589],[15784,26618,64855,32660,52558,67846,78742,53285,64166,51075,81015,67321,76072,12104,61462,27334,45243,66490,74485,62599,12781,25544,46552,61196,59524,80128,63098,57740]]
pairs=[(a,b) for a,b in combinations(range(10),2) if inc(Q,a,b)]
def atom(alpha,beta):
 h=max(alpha);j=max(beta);lo=tuple(range(10,10+h));up=tuple(range(10+h,10+h+j))
 ed=[(a,b) for a in range(10) for b in range(10) if Q[b]>>a&1]+list(zip(lo,lo[1:]))+list(zip(up,up[1:]))
 ed += [(c,q) for i,c in enumerate(lo) for q in range(10) if i<h-alpha[q]]
 ed += [(q,c) for i,c in enumerate(up) for q in range(10) if i>=beta[q]]+[(c,d) for c in lo for d in up]
 return closure(10+h+j,ed)
reports=[];balanced=[];cells=[]
A=(0,2);B=(5,8)
for i,(alpha,beta) in enumerate(profiles):
 p=atom(alpha,beta);E=count(p);N=[count(p,(pair,)) for pair in pairs]
 require(E==expected_E[i],'ambient extension total')
 require(N==expected_N[i],'all core comparison numerators')
 require(all(sum(inc(p,a,b) for b in range(len(p)))<=7 for a in range(len(p))),'seven-thin ambient model')
 require(tuple(x&1023 for x in p[:10])==Q,'core induced order unchanged')
 actual_core=set(roots)|{a for a in range(len(p)) if any(inc(p,a,r) for r in roots)}
 require(actual_core==set(range(10)),'same full root core')
 require(not(p[0]&~p[3]),'critical lower-set inclusion')
 require(all(not(p[b]>>3&1) or p[b]>>0&1 for b in range(len(p))),'critical upper-set inclusion')
 bb=[pair for pair,v in zip(pairs,N) if E<=3*v<=2*E];balanced.append(set(bb))
 # Cell order is A yes/B yes, A yes/B no, A no/B yes, A no/B no.
 cc=[count(p,(a,b)) for a,b in product((A,A[::-1]),(B,B[::-1]))]
 require(sum(cc)==E,'four joint cells partition all extensions');cells.append(cc)
 reports.append({'lower_profile':alpha,'upper_profile':beta,'ambient_predecessors':p,'extensions':E,'comparison_counts':[{'pair':pair,'numerator':v} for pair,v in zip(pairs,N)],'balanced_core_pairs':bb,'joint_cells_A02_B58':cc})
require(not(balanced[0]&balanced[1]),'no common balanced core pair')
require(cells==[[10042,3382,9710,6059],[16922,9696,29630,28459]],'prescribed joint cells')
def margins(cc,first):
 if first=='B':cc=[cc[0],cc[2],cc[1],cc[3]]
 return [1-sum(cc[:2])-cc[j] if j<2 else 1-sum(cc[2:])-cc[j] for j in range(4)]
def law(t):return [t*F(cells[0][j],expected_E[0])+(1-t)*F(cells[1][j],expected_E[1]) for j in range(4)]
endpoint_rows=[]
for first,interval in [('B',(F(0),F(1,2))),('A',(F(1,2),F(1)))]:
 for t in interval:
  mm=margins(law(t),first)
  require(all(v>F(1,25) for v in mm),'uniform strict Gold margin at segment endpoints')
  endpoint_rows.append({'first':first,'t':str(t),'normalized_margins':list(map(str,mm))})
require(min(F(v) for row in endpoint_rows for v in row['normalized_margins'])==F(113042071,2472851451),'prescribed global half-interval minimum')
# Affinity implies the same bound on each whole half interval.
a0=F(13424,29193);a1=F(26618,84707);b0=F(19752,29193);b1=F(46552,84707)
switch=(1-a1-b1)/(a0+b0-a1-b1)
require(switch==F(336799641,674187622),'more-balanced switch')
require(0<switch<1,'switch lies inside segment')
for first,interval in [('B',(F(0),switch)),('A',(switch,F(1)))]:
 for t in interval:require(all(v>F(1,25) for v in margins(law(t),first)),'more-balanced-first rule also has uniform strict margin')
# Exact simple four-cell hostile/edge controls for the general ordering lemma.
ordering_rows=0
for cc_int in product(range(7),repeat=4):
 if sum(cc_int)!=6:continue
 cc=[F(x,6) for x in cc_int];ha=max(sum(cc[:2]),sum(cc[2:]));hb=max(cc[0]+cc[2],cc[1]+cc[3]);ma=min(margins(cc,'A'));mb=min(margins(cc,'B'))
 require((ma>=0)==(ha+max(cc)<=1),'exact weak row scalar criterion')
 require((mb>=0)==(hb+max(cc)<=1),'exact weak column scalar criterion')
 if ha<=hb:
  require(not(mb>=0) or ma>=0,'more-balanced order preserves weak Gold')
  require(not(mb>0) or ma>0,'more-balanced order preserves strict Gold')
 if hb<=ha:
  require(not(ma>=0) or mb>=0,'dual more-balanced order preserves weak Gold')
  require(not(ma>0) or mb>0,'dual more-balanced order preserves strict Gold')
 ordering_rows+=1
require(min(margins([F(3,10),F(3,10),F(4,10),F(0)],'A'))>0 and F(6,10)+F(4,10)==1,'strict scalar converse fails')
print(json.dumps({'status':'PASS','guards':checks,'core_predecessors':Q,'critical_roots':roots,'models':reports,'adaptive_repair':{'A':A,'B':B,'rule':'Ask B first for t<=1/2; ask A first for t>=1/2. Both questions are asked, with the same second question on both branches.','uniform_normalized_Gold_margin_lower_bound':'1/25','half_interval_endpoint_controls':endpoint_rows,'more_balanced_rule_switch':str(switch),'simplex_ordering_controls':ordering_rows},'scope':'Two actual models exclude any common first comparison restricted to their common core. The entire convex segment of their projected laws has strict adaptive Gold. This does not cover all attachments, prove general seven-thin Gold, or assert every mixed law is realizable.'},indent=2))
