#!/usr/bin/env python3
"""Converse-equivariant multiplicative F64 carrier, represented by exponents mod63."""
import verify as v, json, collections as co,itertools as it
from pathlib import Path
D=Path(__file__).resolve().parent
R=json.loads((D/'result.json').read_text());B=[r['blocker'] for r in R['bijection']]
def inv(b):return v.am([((-w)%7,(-u)%7) for i,(u,w) in enumerate(v.E) if b>>i&1])
def encode(b):
 _,typ=v.encode(b)
 remain=[e for i,e in enumerate(v.E) if not (b>>i&1)]
 out={u:{w for a,w in remain if a==u} for u in v.V}
 inc={u:{a for a,w in remain if w==u} for u in v.V}
 if typ in ('two-source','two-sink'):
  u,w=[x for x in v.V if not (inc if typ=='two-source' else out)[x]]
  if (u,w) not in v.EI:u,w=w,u
  k={1:0,2:1,4:2}[(w-u)%7]
  return (9*u+pow(16,k,63))%63 if typ=='two-source' else (9*w-pow(16,k,63))%63
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
print(json.dumps({'status':'PASS','orbits':orbits,'fixed_blockers':fixed,'ordinary_even_fixed_count':len(evenfixed),'map':[{'blocker':b,'exponent':P[b],'type':v.encode(b)[1]} for b in B]},indent=2,sort_keys=True))
