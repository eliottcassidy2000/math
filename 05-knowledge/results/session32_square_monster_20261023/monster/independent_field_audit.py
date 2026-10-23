#!/usr/bin/env python3
"""Independent small finite-field enumeration: no symbolic resultants.
For each Monster prime p>=5, detect rational supersingular j by point counts,
then enumerate every nonzero u in F_{p^2}, compute j=(u+256)^3/u^2 and
s=u+4096/u. Compare coarse-point counts to the polynomial audit JSON.
The premise that all supersingular j are rational is independently checked
in audit_level_two.py, and is also Ogg's cited characterization.
"""
import json
from pathlib import Path
from audit_level_two import ss_j_from_counts, MONSTER
BASE=Path(__file__).resolve().parent
rows=json.loads((BASE/'level_two_exact.json').read_text())['rows']
out=[]
for row in rows:
 p=row['p']
 if p not in MONSTER: continue
 d=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
 def add(x,y): return ((x[0]+y[0])%p,(x[1]+y[1])%p)
 def mul(x,y): return ((x[0]*y[0]+d*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
 def inv(x):
  den=pow((x[0]**2-d*x[1]**2)%p,-1,p)
  return (x[0]*den%p,-x[1]*den%p)
 ss=set(ss_j_from_counts(p))
 U=set();S=set();nonrational_edges={}
 for a in range(p):
  for b in range(p):
   if (a,b)==(0,0): continue
   u=(a,b);v=add(u,(256%p,0));ui=inv(u)
   j=mul(mul(mul(v,v),v),mul(ui,ui))
   if j[1]==0 and j[0] in ss:
    v=mul((4096%p,0),ui)
    sval=add(u,v)
    U.add(u);S.add(sval)
    if sval[1]!=0:
     w=add(v,(256%p,0));vi=inv(v)
     jright=mul(mul(mul(w,w),w),mul(vi,vi))
     assert jright[1]==0 and jright[0] in ss
     endpoints=sorted([j[0],jright[0]])
     if sval in nonrational_edges: assert nonrational_edges[sval]==endpoints
     nonrational_edges[sval]=endpoints
 nrat=sum(b!=0 for a,b in S)
 assert len(U)==row['ss_X0_2_count']
 assert len(S)==row['ss_X0plus_2_count']
 assert nrat==row['nonrational_X0plus_2_points']
 # Frobenius(a+b sqrt(d)) = a-b sqrt(d).
 assert all((a,(-b)%p) in S for a,b in S)
 out.append({'p':p,'extension_relation':f'z^2={d}', 'enumerated_nonzero_u':p*p-1,
   'supersingular_j_from_pointcounts':sorted(ss),'u_point_count':len(U),
   's_point_count':len(S),'s_nonrational_count':nrat,
   'nonrational_s_coordinates':sorted([list(x) for x in S if x[1]!=0]),
   'nonrational_unordered_isogenies':[{'s_coordinates':list(x),'unordered_j_endpoints':nonrational_edges[x]} for x in sorted(nonrational_edges)],'PASS':True})
print(json.dumps(out,indent=2))
