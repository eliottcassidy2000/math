#!/usr/bin/env python3
"""Standalone integer-only factorized feasibility replay; no scan imports."""
from itertools import permutations
from pathlib import Path
import json
from fractions import Fraction
ROOT=Path(__file__).resolve().parent
def require(condition, message):
 if not condition: raise ValueError(message)
P=(0,28,0,0,12,12,13,125,125)
ORDER=(2,3,0,4,5,6,1,7,8)
U={():1,(2,):1,(2,3,4,5):6}
V={():1,(8,):1,(5,6,7,8):6}
witness=json.loads((ROOT/'witness.json').read_text())
require(tuple(witness['core'])==P and tuple(witness['order'])==ORDER, 'Witness core or order changed')
require({tuple(t['word']):t['increment'] for t in witness['U']}==U and len(witness['U'])==len(U), 'Prefix witness changed')
require({tuple(t['word']):t['increment'] for t in witness['V']}==V and len(witness['V'])==len(V), 'Suffix witness changed')
require(witness['all_unlisted_increments']==0, 'Unlisted increments must be zero')
pairs=[(a,b) if ORDER.index(a)<ORDER.index(b) else (b,a) for a in range(9) for b in range(a+1,9) if not(P[a]>>b&1 or P[b]>>a&1)]
E=0; mass=0; before=[0]*len(pairs); gens={}; words=[]
for w in permutations(range(9)):
 mask=0
 for v in w:
  if P[v]&~mask:break
  mask|=1<<v
 else:
  E+=1; r={v:i for i,v in enumerate(w)}; h=min(r[0],r[1]); k=max(r[0],r[1])+1
  L=sum(U.get(w[:i],0) for i in range(h+1)); R=sum(V.get(w[j:],0) for j in range(k,10)); weight=L*R
  mass+=weight
  for t,(a,b) in enumerate(pairs):before[t]+=weight*int(r[a]<r[b])
  for i in range(h+1):
   for j in range(k,10):
    key=(w[:i],w[j:]);gens[key]=gens.get(key,0)+1
  words.append({'word':w,'L':L,'R':R,'weight':weight})
require(E==408 and len(gens)==133, 'Core extension or support count changed')
require(U[()]==V[()]==1 and all(isinstance(x,int) and x>=0 for x in [*U.values(),*V.values()]), 'Nonnegative integer increments with positive empty normalization required')
require(mass==948, 'Exact mass must be 948')
margins=[3*n-2*mass for n in before]
require(pairs==[(0,1),(2,0),(3,0),(0,4),(0,5),(5,1),(6,1),(1,7),(1,8),(2,3),(4,5),(4,6),(5,6),(7,8)], 'Forward pair list changed')
require(margins==[705,300,84,3,246,246,3,84,300,12,3,228,3,12], 'Exact margin list changed')
require(min(margins)>0, 'All forward probabilities must strictly exceed two thirds')
unsupported=((2,3,4,5),(5,6,7,8))
require(U[unsupported[0]]*V[unsupported[1]]==36, 'Unsupported product must remain 36')
require(gens.get(unsupported,0)==0, 'Overlapping words must have zero support')
require(min(Fraction(n,mass) for n in before)==Fraction(211,316), 'Minimum probability changed')
active=[{'prefix':u,'suffix':v,'coefficient':U.get(u,0)*V.get(v,0),'extensions':e} for (u,v),e in sorted(gens.items()) if U.get(u,0)*V.get(v,0)]
require(len(active)==8, 'Eight supported products required')
report={'unsupported_product':{'prefix':unsupported[0],'suffix':unsupported[1],'coefficient':36,'extensions':0},'status':'PASS','scope':'Exact nonnegative factorized boundary array; actual exterior realizability is not asserted','core':P,'order':ORDER,'extensions':E,'supported_generators':len(gens),'empty_normalizations':[U[()],V[()]],'prefix_increments':[{'word':w,'value':v} for w,v in U.items()],'suffix_increments':[{'word':w,'value':v} for w,v in V.items()],'all_unlisted_increments':0,'total_mass':mass,'minimum_scaled_margin':min(margins),'minimum_forward_probability':str(min(Fraction(n,mass) for n in before)),'pairs':[{'pair':p,'before_mass':n,'scaled_margin':m,'probability':str(Fraction(n,mass))} for p,n,m in zip(pairs,before,margins)],'active_supported_products':active}
(ROOT/'replay_report.json').write_text(json.dumps(report,indent=2)+'\n');(ROOT/'weighted_extensions.json').write_text(json.dumps(words,indent=2)+'\n');print(json.dumps(report,indent=2))
