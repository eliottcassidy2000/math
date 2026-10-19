#!/usr/bin/env python3
"""Exact standard-library replay of a fixed-core exterior-cap certificate."""
from pathlib import Path
from itertools import permutations
import json,hashlib
ROOT=Path(__file__).resolve().parent

def require(condition,message):
 if not condition: raise ValueError(message)

c=json.loads((ROOT/'fixed_certificate.json').read_text())
p=tuple(c['core']);o=tuple(c['order']);n=len(p)
require(p==(0,28,0,0,12,12,13,125,125),'unexpected core')
require(o==(2,3,0,4,5,6,1,7,8),'unexpected projected order')
for b in range(n):
 require(not(p[b]>>b&1),'reflexive relation')
 for a in range(n):
  if p[b]>>a&1: require(p[a]&~p[b]==0,'nontransitive relation')
deg=[sum(a!=b and not(p[a]>>b&1 or p[b]>>a&1) for b in range(n)) for a in range(n)]
require(deg==[5,5,2,2,3,4,3,2,2],'degree mismatch')
pairs=[];weights=[]
for t in c['comparisons']:
 a,b=t['pair'];w=t['weight']
 require(isinstance(w,int) and w>0,'invalid comparison weight')
 require(a!=b and not(p[a]>>b&1 or p[b]>>a&1),'comparable certificate pair')
 require(o.index(a)<o.index(b),'backward certificate pair')
 pairs.append((a,b));weights.append(w)
require(len(set(pairs))==len(pairs),'duplicate pair')
cols={};count=0
for word in permutations(range(n)):
 seen=0
 for q in word:
  if p[q]&~seen:break
  seen|=1<<q
 else:
  count+=1;r={q:i for i,q in enumerate(word)}
  first=min(r[0],r[1]);last=max(r[0],r[1])
  indicators=[int(r[a]<r[b]) for a,b in pairs]
  for i in range(first+1):
   for j in range(last+1,n+1):
    key=(word[:i],word[j:]);v=cols.setdefault(key,[0]*(len(pairs)+1));v[0]+=1
    for k,ind in enumerate(indicators):v[k+1]+=ind
require(count==408 and len(cols)==133,'enumeration mismatch')
res={key:sum(w*(3*v[k+1]-2*v[0]) for k,w in enumerate(weights)) for key,v in cols.items()}
seen_caps=set()
for t in c['caps']:
 side=t['side'];u=tuple(t['prefix']);v=tuple(t['suffix']);key=(u,v);w=t['weight'];bound=t['bound']
 require(side in ('A','B'),'unknown side');require(key in cols,'unsupported child')
 require(isinstance(w,int) and w>0,'invalid cap weight')
 require((side,key) not in seen_caps,'duplicate cap');seen_caps.add((side,key))
 if side=='A':
  require(bool(u),'empty prefix cap');q=u[-1];parent=(u[:-1],v)
 else:
  require(bool(v),'empty suffix cap');q=v[0];parent=(u,v[1:])
 require(parent in cols,'unsupported parent')
 require(bound==7-deg[q],'invalid exterior degree cap')
 # C row is lambda_child - bound*lambda_parent <= 0.
 # Verify comparison row minus a nonnegative combination of C rows <= 0.
 res[key]-=w;res[parent]+=w*bound
require(max(res.values())<=-1,'certificate residual is not strictly negative')
require((min(res.values()),max(res.values()))==(-540,-1),'unexpected residual range')
report={'status':'PASS','core':p,'order':o,'degrees':deg,'extensions':count,'supported_generators':len(cols),'comparison_terms':len(pairs),'cap_terms':len(c['caps']),'residual_range':[min(res.values()),max(res.values())],'certificate_sha256':hashlib.sha256((ROOT/'fixed_certificate.json').read_bytes()).hexdigest(),'scope':'Every actual maximum-incomparability-degree-seven exterior of this fixed full root core and projected strict majority order is excluded. No whole-morphology claim.'}
print(json.dumps(report,indent=2))
