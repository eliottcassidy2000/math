#!/usr/bin/env python3
"""Independent exact audit, with no imports from the author's scripts.
Enumerate Hasse polynomial roots in F_{p^2}, obtain j-values, enumerate
X0(2) u-values, and quotient by u -> 4096/u. This does NOT begin by
assuming supersingular j-values are rational. It includes p=3 explicitly.
"""
from math import comb, prod
from pathlib import Path
import json

BASE=Path(__file__).resolve().parent
SOURCE=BASE.parent/'monster'/'level_two_exact.json'
M={2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
B={2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}
source={r['p']:r for r in json.loads(SOURCE.read_text())['rows']}

def row(p):
 d=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
 Z=(0,0); O=(1,0)
 def plus(x,y): return ((x[0]+y[0])%p,(x[1]+y[1])%p)
 def minus(x,y): return ((x[0]-y[0])%p,(x[1]-y[1])%p)
 def times(x,y): return ((x[0]*y[0]+d*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
 def inverse(x):
  denom=pow((x[0]*x[0]-d*x[1]*x[1])%p,-1,p)
  return (x[0]*denom%p,-x[1]*denom%p)
 def cube(x): return times(times(x,x),x)
 def scalar(c): return (c%p,0)
 def evaluate(cs,x):
  z=Z
  for c in reversed(cs): z=plus(times(z,x),scalar(c))
  return z
 field=[(a,b) for a in range(p) for b in range(p)]
 h=[comb((p-1)//2,k)**2%p for k in range((p-1)//2+1)]
 L={x for x in field if evaluate(h,x)==Z}
 assert len(L)==(p-1)//2
 assert not (L&{Z,O})
 J=set()
 for x in L:
  q=plus(minus(times(x,x),x),O)
  den=times(times(x,x),times(minus(x,O),minus(x,O)))
  J.add(times(times(scalar(256),cube(q)),inverse(den)))
 U=set(); S=set(); fibers={}
 for u in field[1:]:
  j=times(cube(plus(u,scalar(256))),inverse(times(u,u)))
  if j in J:
   U.add(u)
   s=plus(u,times(scalar(4096),inverse(u)))
   S.add(s); fibers.setdefault(s,set()).add(u)
 # Completeness in u follows independently from its cubic j-map: generic
 # fibers have 3 points, j=0 has 1, and j=1728 has 2 (p>3).
 if p>3:
  assert len(U)==3*len(J)-2*(Z in J)-(scalar(1728) in J)
 else:
  assert p==3 and J=={Z} and U=={(2,0)} and S=={(1,0)}
 assert all(times(scalar(4096),inverse(u)) in U for u in U)
 assert all((a,-b%p) in S for a,b in S)
 assert all(len(us) in [1,2] for us in fibers.values())
 assert all(len(us)==1 if next(iter(us)) in {(64%p,0),(-64%p,0)} else len(us)==2 for us in fibers.values())
 edge_endpoints=[]
 for s,us in sorted(fibers.items()):
  if s[1]==0: continue
  ends={times(cube(plus(u,scalar(256))),inverse(times(u,u))) for u in us}
  edge_endpoints.append({'s':list(s),'j_endpoints':sorted([list(j) for j in ends])})
 if p in {29,41,59,71}:
  expected={29:[2,25],41:[3,32],59:[28,48],71:[41,66]}[p]
  assert len(edge_endpoints)==2
  assert all(x['j_endpoints']==[[j,0] for j in expected] for x in edge_endpoints)
  assert not set(expected)&{0,1728%p}
 pairs={tuple([1,-2*a%p,(a*a-d*b*b)%p]) for a,b in S if b}
 assert 2*len(pairs)==sum(bool(b) for a,b in S)
 assert all(pow((c[1]*c[1]-4*c[2])%p,(p-1)//2,p)==p-1 for c in pairs)
 if p in source:
  r=source[p]
  assert (len(J),len(U),len(S),sum(bool(b) for a,b in S))==(r['ss_j_count'],r['ss_X0_2_count'],r['ss_X0plus_2_count'],r['nonrational_X0plus_2_points'])
  coeff=r['X0plus_2_polynomial_coefficients_high_to_low']
  assert all(evaluate(coeff[::-1],s)==Z for s in S)
  assert len(coeff)==len(S)+1
 return {'p':p,'quadratic_field_nonsquare':d,'hasse_root_count':len(L),
  'j_count':len(J),'nonrational_j_count':sum(bool(b) for a,b in J),'u_count':len(U),'s_count':len(S),
  'fricke_fixed_u':sorted([list(u) for u in U if times(u,u)==scalar(4096)]),
  'nonrational_s_count':sum(bool(b) for a,b in S),'irreducible_quadratics_high_to_low':sorted(pairs),
  'rational_s':sorted(a for a,b in S if b==0),'nonrational_s_endpoints':edge_endpoints,'PASS':True}

rows=[row(p) for p in sorted((set(M)|{37,43})-{2})]
orderM=prod(p**e for p,e in M.items()); orderB=prod(p**e for p,e in B.items())
assert orderM==808017424794512875886459904961710757005754368000000000
assert orderB==4154781481226426191177580544000000
assert orderM%orderB==0
strata={'equal_positive':[p for p in M if M[p]==B.get(p,0)],'Monster_only':[p for p in M if p not in B],'increased_common':[p for p in B if M[p]>B[p]]}
assert strata=={'equal_positive':[17,19,23,31,47],'Monster_only':[29,41,59,71],'increased_common':[2,3,5,7,11,13]}
quotient=orderM//orderB
core=prod(p**(M[p]-B[p]) for p in B)
exceptional=prod(p for p in M if p not in B)
assert core==39046339332000 and exceptional==4980721 and quotient==core*exceptional
f={p:M[p]-B[p] for p in B if M[p]>B[p]}
mod_witness={m: [(p,q,f[p],f[q]) for p in f for q in f if p<q and (p-q)%m==0 and (f[p]-f[q])%m] for m in [3,9]}
assert all(mod_witness.values())
cycles=[]
for correction in [0,1]:
 g={p:(4 if p==2 and correction else e) for p,e in f.items()}
 q=lambda x:13 if x==4 else x
 x=13; path=[]
 while x not in path:
  path.append(x); x=q(g[x])
 assert x==13
 cycles.append(path)
assert cycles==[[13,2,5,3,7],[13,2]]
result={'method':'Independent exhaustive Fp2 Hasse roots, exact j/u/s maps, and set quotient; no imported author code',
 'orders':{'Monster':orderM,'Baby':orderB,'M_over_B':quotient,'M_over_2B':quotient//2,'core':core,'exceptional_product':exceptional},
 'valuation_strata':strata,'valuation_cycles':cycles,'modular_descent_witnesses':mod_witness,'rows':rows,'PASS':True}
print(json.dumps(result,indent=2))
