#!/usr/bin/env python3
"""Converse-equivariant multiplicative F64 carrier, represented by exponents mod63."""
import verify_base as v, json, collections as co,itertools as it
from pathlib import Path
D=Path(__file__).resolve().parent
R=json.loads((D/'blockers.json').read_text());B=[r['blocker'] for r in R['bijection']]
def inv(b):return v.am([((-w)%7,(-u)%7) for i,(u,w) in enumerate(v.E) if b>>i&1])
def encode(b,baseline=10):
 _,typ=v.encode(b)
 remain=[e for i,e in enumerate(v.E) if not (b>>i&1)]
 out={u:{w for a,w in remain if a==u} for u in v.V}
 inc={u:{a for a,w in remain if w==u} for u in v.V}
 if typ in ('two-source','two-sink'):
  u,w=[x for x in v.V if not (inc if typ=='two-source' else out)[x]]
  if (u,w) not in v.EI:u,w=w,u
  k={1:0,2:1,4:2}[(w-u)%7]
  return (9*u+baseline*pow(16,k,63))%63 if typ=='two-source' else (9*w-baseline*pow(16,k,63))%63
 if typ=='star':
  c=next(x for x in v.V if not inc[x] and not out[x]);return 9*c
 di=out if typ=='common-out' else inc
 c=next(x for x in v.V if sum(di[y]=={x} for y in v.V)==3)
 return (9*c+(21 if typ=='common-out' else 42))%63
P={b:encode(b) for b in B}
v.require(set(P.values())==set(range(63)),'multiplicative bijection')
for b in B:
 v.require(P[inv(b)]==(-P[b])%63,'inversion law')
 for a,c in v.G:
  k={1:0,2:1,4:2}[a]
  v.require(P[v.ga(b,(a,c))]==(9*c+pow(16,k,63)*P[b])%63,'affine action')
left=set(B);orbits=[]
while left:
 b=min(left);O={v.ga(z,g) for z in (b,inv(b)) for g in v.G};left-=O
 orbits.append({'size':len(O),'types':sorted({v.encode(z)[1] for z in O}),'phase_orders':sorted({9//__import__('math').gcd(P[z],9) for z in O})})
fixed=[b for b in B if inv(b)==b]
v.require(len(fixed)==1,'converse fixed count')
v.require(sorted(o['size'] for o in orbits)==[7,14,42],'extended orbit count')
# Fixed points of the ordinary vertex-negation action on even subsets.
evenfixed=[s for s in range(1,128) if s.bit_count()%2==0 and v.vm([(-x)%7 for x in v.vertices(s)])==s]
v.require(len(evenfixed)==7,'linear hostile')
def rho(e):
 for c in v.V:
  if e%21==9*c%21:return v.FULL^(1<<c)
  if e%21==(9*c+7)%21:return v.vm({c}|v.IN[c])
  if e%21==(9*c-7)%21:return v.vm({c}|v.OUT[c])
 raise RuntimeError('quotient decoder')
for b in B:v.require(rho(P[b])==v.boundary(b),'commutative square')
baseline_hostile=next((a,b,encode(a,1),encode(b,1)) for a,b in it.combinations(B,2) if v.boundary(a)==v.boundary(b) and encode(a,1)%21!=encode(b,1)%21)
boundary_pairs={(v.boundary(b),P[b]%21) for b in B}
v.require(len(boundary_pairs)==21,'projective boundary graph')
for a in B:
 for b in B:v.require((v.boundary(a)==v.boundary(b))==(P[a]%21==P[b]%21),'fiber equivalence')
print(json.dumps({'status':'PASS','boundary_graph':sorted(boundary_pairs),'baseline_one_failure':baseline_hostile,'map':[{'blocker':b,'boundary':v.boundary(b),'exponent':P[b],'type':v.encode(b)[1]} for b in B]},indent=2,sort_keys=True))
