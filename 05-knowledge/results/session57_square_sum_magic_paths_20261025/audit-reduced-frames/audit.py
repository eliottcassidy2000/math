#!/usr/bin/env python3
"""Independent exact checks; no source writes, and no assertions."""
import itertools as it,json,hashlib,argparse
from pathlib import Path
import sympy as s

def need(condition,label):
 if not condition: raise RuntimeError(label)
def edges(p):return frozenset(tuple(sorted(e)) for e in zip(p,p[1:]))
def components(nodes,es):
 left=set(nodes);out=[]
 while left:
  seen={min(left)};stack=list(seen)
  while stack:
   u=stack.pop()
   for a,b in es:
    if u in (a,b):
     v=b if u==a else a
     if v not in seen:seen.add(v);stack.append(v)
  left-=seen;out.append(sorted(seen))
 return out
ap=argparse.ArgumentParser();ap.add_argument('--expected-covers',type=int,default=36);args=ap.parse_args()
p=Path(__file__).resolve().parent/'input'/'two-path-frame';c=json.loads((p/'certificate.json').read_text())
P=(0,7,1,5,2,4,3,6);Q=(0,7,2,5,1,6,3,4);U=sorted(edges(P)|edges(Q))
E=s.Matrix([[int(v in e) for v in range(8)] for e in U]);A=s.Matrix(c['A']);T=s.Matrix(c['reconstruction_T']);Z=s.Matrix(c['cycle_constraints_Z']);B=s.Matrix(c['extra_constraints_B'])
need(E==s.Matrix(c['E']),'frozen incidence');need((E.rank(),A.rank(),Z.rank(),B.rank(),Z.col_join(B).rank())==(7,4,2,4,6),'frame ranks')
need(A==B*E and B==A*T,'decoration pullback')
g=s.Matrix([1]*4+[-1]*4);r0=s.zeros(1,8);r0[0,4]=1
need(T*E==s.eye(8)+g*r0 and E*g==s.zeros(9,1),'integral gauge inverse')
need(Z*E==s.zeros(2,8),'cycles')
# Build independent magic line constraints in the theorem's edge order.
grid=(6,7,3,2,1,5,4,0,8);lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
H=[]
for line in lines[1:]:
 row=[0]*9
 for k in line:row[grid[k]]+=1
 for k in lines[0]:row[grid[k]]-=1
 H.append(row)
H=s.Matrix(H);need(H.rank()==6 and H.col_join(Z).col_join(B).rank()==6,'exact magic rowspace')
need(len([p for p in it.permutations(range(8)) if p<p[::-1] and edges(p)<=set(U)])==4,'four frame paths')
routes=((1,5,2),(1,7,2),(1,6,3,4,2));need(set().union(*(edges(p) for p in routes))==set(U)-{(0,7)},'theta routes')
# Independent K33 enumeration and cover-partner components.
V=tuple(range(6));K=frozenset((i,j) for i in range(3) for j in range(3,6));ps=[p for p in it.permutations(V) if p<p[::-1] and edges(p)<=K]
pairs=[(i,j) for i,j in it.combinations(range(len(ps)),2) if edges(ps[i])|edges(ps[j])==K]
cycs={edges(p+(p[0],)) for p in it.permutations(V) if edges(p+(p[0],))<=K}
need((len(ps),len(cycs),len(pairs))==(36,6,args.expected_covers),'K33 enumeration')
need(all(sum(i in e for e in pairs)==2 for i in range(36)),'two partners')
need(sorted(map(len,components(range(36),pairs)))==[6]*6,'six partner hexagons')
for p in ps:
 need(sorted(map(len,components(V,K-edges(p))))==[2,4],'complement components')
P=(0,3,1,4,2,5);Q=(4,0,5,1,3,2)
need(edges(P)|edges(Q)==K and edges(P)&edges(Q)=={(1,3)},'explicit cover')
a=s.symbols('a1:4');b=s.symbols('b1:4');G=s.Matrix([[a[2]+b[2],a[1]+b[1],a[0]+b[0]],[a[0]+b[1],a[2]+b[0],a[1]+b[2]],[a[1]+b[0],a[0]+b[2],a[2]+b[1]]]);total=sum(a+b)
need(all(s.expand(sum(G.row(i))-total)==s.expand(sum(G.col(i))-total)==0 for i in range(3)),'semimagic')
need(s.expand(sum(G[i,i] for i in range(3))-total)==2*a[2]-a[0]-a[1],'A guard');need(s.expand(sum(G[i,2-i] for i in range(3))-total)==2*b[0]-b[1]-b[2],'B guard')
M=s.Matrix([[int(v in (i,3+j)) for v in V] for i in range(3) for j in range(3)]);R=[]
for i,j in it.product((1,2),repeat=2):
 row=[0]*9
 for k,v in ((3*i+j,1),(3*i,-1),(j,-1),(0,1)):row[k]+=v
 R.append(row)
R+= [[1,0,0,-2,0,0,1,0,0],[1,-2,1,0,0,0,0,0,0]]
need(M.rank()==5 and s.Matrix(R[:4]).rank()==4 and s.Matrix(R).rank()==6,'AP ranks')
# Symbolic positivity/edge-value identity; positivity uses C>x+y>0 and integer inputs.
C,x,y=s.symbols('C x y');AA=[1+4*i*x for i in range(3)];BB=[4*(C-x-y)-1+4*j*y for j in range(3)]
need(all(s.expand(AA[i]+BB[j]-4*(C+(i-1)*x+(j-1)*y))==0 for i,j in it.product(range(3),repeat=2)),'scaled AP labels')
# Exhaust ALL arrangements for a real example and the affine F3 plane, rather than only generated normal forms.
real=(-4,-3,-2,-1,0,1,2,3,4);nreal=0;n3=0
for p in it.permutations(range(9)):
 if all(sum(real[p[k]] for k in line)==sum(real[p[k]] for k in lines[0]) for line in lines[1:]):nreal+=1
 # F3^2 labels represented as base-3 digits; coordinatewise line sums.
 sums=[(sum(p[k]%3 for k in line)%3,sum(p[k]//3 for k in line)%3) for line in lines]
 if len(set(sums))==1:n3+=1
need(nreal==8 and n3==432,'exhaustive real and characteristic-three layouts')
# Direct residue and all-gauge test, not a search producer replay.
q=59;AA=[48,0,11];BB=[46,1,15];vals={(u+v)%q for u,v in it.product(AA,BB)};squares={u*u%q for u in range(1,q)}
good=[t for t in range(q) if len({(a+t)%q for a in AA}|{(b-t)%q for b in BB})==6]
need(len(vals)==9 and vals<=squares and len(good)==50,'F59 witness and gauges')
print(json.dumps({'status':'PASS','two_path_ranks':[7,4,2,4,6],'AP_ranks':[5,4,6],'K33_paths':len(ps),'K33_cycles':len(cycs),'covers':len(pairs),'partner_components':[6]*6,'real_layouts':nreal,'char3_layouts':n3,'F59_valid_gauges':len(good)},sort_keys=True,indent=2))
