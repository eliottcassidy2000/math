#!/usr/bin/env python3
"""Independent quadratic controls; no producer imports and no native orbit census."""
from itertools import product,combinations
from pathlib import Path
from hashlib import sha256
import argparse,json
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo-root",type=Path,required=True,help="Path to the math repository checkout")
REPO=parser.parse_args().repo_root.resolve()
if not (REPO/"05-knowledge/results/future_tail_blocks_20261021b/receipt.json").is_file():
 parser.error("--repo-root does not contain the incoming stopped-tail receipt")
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok: raise ValueError(msg)
def plus(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
def q(v): return v[0]**2-v[0]*v[1]+v[1]**2
def d(v): return q(plus(v[:2],v[2:]))-q(v[:2])
def k(a,b): return d(plus(a,b))-d(a)-d(b)
def energy(C): return sum(map(d,C))
def U(C,B):
 return energy(C)+sum(max(0,d(b)+max([0]+[abs(k(c,b)) for c in C])) for b in B)+sum(abs(k(a,b)) for a,b in combinations(B,2))
def assignments(C,B):
 choices=[None]+list(product(range(len(C)),(-1,1)))
 for route in product(choices,repeat=len(B)):
  out=list(C)
  for b,r in zip(B,route):
   if r is not None: out[r[0]]=plus(out[r[0]],b if r[1]==1 else neg(b))
  yield energy(out)
def signed(F):
 return [tuple(sum(s*v[i] for s,v in zip(S,F)) for i in (0,1)) for S in product((-1,1),repeat=len(F))]
def beta(a,b):return q(plus(a,b))-q(a)-q(b)
def Bgate(C,Fo,Fn):
 G=max(q(n)+beta(o,n) for o in signed(Fo) for n in signed(Fn))
 L=sum(max([0]+[abs(beta(c[2:],v)) for c in C]) for v in Fo)
 L+=sum(max([0]+[abs(beta(plus(c[:2],c[2:]),v)) for c in C]) for v in Fn)
 return energy(C)+G+L
C=[(1,0,0,1),(0,1,-1,0),(0,0,0,0)]
V=list(product((-1,0,1),repeat=4))
for a,b in product(V,repeat=2):
 upper=U(C,[a,b])
 need(max(assignments(C,[a,b]))<=upper,'optional signed placements')
 need(U(C,[plus(a,b)])<=upper,'merge')
 for sign in (-1,1):
  aa=a if sign==1 else neg(a)
  need(U([plus(C[0],aa),*C[1:]],[b])<=upper,'signed arrival')
 need(U(C,[neg(a),b])==upper,'signed transport')
 need(d(plus(a,b))==d(a)+d(b)+k(a,b),'polarization')
# Cross-comparator examples with exact full independent source-sign search.
Z=[(0,0,0,0)]
tri=[(1,0),(0,1),(-1,-1)]
tri_blocks=[(0,0,*v) for v in tri]
need(U(Z,tri_blocks)==6 and Bgate(Z,[],tri)==4,'zonotope sharper than abs-pair envelope')
need(U(Z,[])==0 and Bgate(Z,[],[(1,0),(-1,0)])==4,'coalesced zero sharper than individual source zonotope')
first=Bgate(Z,[],tri)
after=Bgate([(0,0,1,0)],[],tri[1:])
need((first,after)==(4,6),'raw zonotope gate can increase on arrival')
need(min(first,after)==4,'running minimum retains prior all-future validity')
mixed=(1,0,-1,0)
need(d(mixed)==-1 and U(Z,[mixed])==0,'full-sign coupled-block maximum misses absence')
need(k(mixed,(1,0,0,0))==-2,'mixed total zero retains baseline')
# Mature zero new packet removes source norm but never energy or R-prime.
for w,h in product(product(range(-2,3),repeat=2),repeat=2):
 for A,Znorm in [(8,3),(13,8)]:
  physical=q(plus(w,h))-q(w)
  Gamma=physical-(A-Znorm)
  need(Gamma-Znorm==physical-A,'compression leaves physical interaction invariant')
# Literal pending-pair receipt arithmetic, intentionally not source-route certification.
packet=REPO/'05-knowledge/results/future_tail_blocks_20261021b'
receipt=json.loads((packet/'receipt.json').read_text())
rows=[]
for case in receipt['cases']:
 ports=[tuple(row[1:5]) for row in case['arrived_ports']]
 blocks=[tuple(row[1:5]) for row in case['pending']]
 current=energy(ports)
 local=sum(max(0,d(b)+max([0]+[abs(k(c,b)) for c in ports])) for b in blocks)
 pair=sum(abs(k(a,b)) for a,b in combinations(blocks,2))
 need([current,local,pair,current+local+pair]==case['envelope'],'saved literal pair arithmetic')
 need(len({row[0] for row in case['pending']})==len(blocks),'literal current node uniqueness')
 need(case['early_cuts_pass']==all(r[1]<=case['budget'] for r in case['history']),'earlier cuts gate')
 rows.append({'N':case['N'],'D':case['D'],'pending_nodes':len(blocks),'envelope':current+local+pair})
result={'checks':checks,'native_receipt_arithmetic_only':rows,'synthetic_bound_comparison':{'triangle_U':6,'triangle_B':4,'B_before_arrival':first,'B_after_arrival':after},'scope':'Quadratic controls and saved arithmetic only; no independent native prefix replay.'}
print(json.dumps(result,sort_keys=True))
