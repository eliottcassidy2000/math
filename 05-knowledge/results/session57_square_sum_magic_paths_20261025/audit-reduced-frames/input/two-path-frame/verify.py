#!/usr/bin/env python3
"""Exact rational two-path-frame certificate; standard library, safe under -O."""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json,sys,hashlib
ROOT=Path(__file__).resolve().parent

def require(ok,message):
 if not ok: raise RuntimeError(message)
def rank(a):
 a=[[F(v) for v in row] for row in a];r=0
 for c in range(len(a[0])):
  k=next((k for k in range(r,len(a)) if a[k][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
  for k in range(len(a)):
   if k!=r:
    v=a[k][c];a[k]=[x-v*y for x,y in zip(a[k],a[r])]
  r+=1
  if r==len(a):break
 return r

def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def edges(p):return {tuple(sorted(e)) for e in zip(p,p[1:])}
def components(V,es):
 left=set(V);out=[]
 while left:
  seen={min(left)};front=list(seen)
  while front:
   u=front.pop()
   for a,b in es:
    if u not in (a,b):continue
    v=b if u==a else a
    if v not in seen:seen.add(v);front.append(v)
  out.append(sorted(seen));left-=seen
 return out

def main():
 P=[0,7,1,5,2,4,3,6];Q=[0,7,2,5,1,6,3,4]
 pe,qe=edges(P),edges(Q);U=sorted(pe|qe);shared=pe&qe
 require(U==[(0,7),(1,5),(1,6),(1,7),(2,4),(2,5),(2,7),(3,4),(3,6)],'union')
 require(len(shared)==5 and sorted(map(len,components(range(8),shared)))==[2,3,3],'shared forest')
 E=[[int(k in e) for k in range(8)] for e in U]
 A=[[-1,1,0,0,-1,1,0,0],[-1,0,1,0,-1,0,1,0],[-1,0,0,1,-1,0,0,1],[0,0,0,0,1,-1,-1,1]]
 # R0=0; reconstruct using the seven edges of P.
 T=[[1,1,0,-1,1,-1,0,0,0],[0,1,0,0,1,-1,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,0,1,0],[0]*9,[0,0,0,0,-1,1,0,0,0],[0,0,0,0,0,0,0,-1,1],[0,-1,0,1,-1,1,0,0,0]]
 Z=[[0,1,0,-1,0,-1,1,0,0],[0,1,-1,0,1,-1,0,-1,1]]
 B=mul(A,T);require(mul(B,E)==A,'A=BE')
 require(mul(Z,E)==[[0]*8]*2,'cycle identities')
 require([rank(a) for a in [E,A,Z,B,Z+B]]==[7,4,2,4,6],'ranks')
 gauge=[[1]]*4+[[-1]]*4
 require(mul(E,gauge)==[[0]]*9 and mul(A,gauge)==[[0]]*4,'gauge')
 TE=mul(T,E)
 # T E v = v + R0*(+1 on L,-1 on R).
 require(TE==[[int(i==j)+(1 if i<4 else -1)*int(j==4) for j in range(8)] for i in range(8)],'integral reconstruction')
 grid=[6,7,3,2,1,5,4,0,8]
 lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
 H=[]
 for line in lines[1:]:
  h=[0]*9
  for k in line:h[grid[k]]+=1
  for k in lines[0]:h[grid[k]]-=1
  H.append(h)
 require(rank(H)==rank(H+Z+B)==6,'magic row-space identity')
 # Independent six actual equal-line equations, recorded reproducibly.
 H6=[]
 for h in H:
  if rank(H6+[h])>len(H6):H6.append(h)
 require(len(H6)==6,'six magic equations')
 coeff=[[1,1,1],[1,0,0],[1,-1,1],[1,0,1],[1,0,-1],[1,1,-1],[1,1,0],[1,-1,-1],[1,-1,0]]
 require(mul(Z+B,coeff)==[[0]*3]*6 and rank(coeff)==3,'magic parameter plane')
 paths=[]
 for p in permutations(range(8)):
  if p[0]<p[-1] and edges(p)<=set(U):paths.append(list(p))
 require(len(paths)==4,'Hamiltonian count')
 # Theta branch vertices L1=1,R3=7? Compute and verify explicit routes below.
 deg=[sum(v in e for e in U) for v in range(8)]
 # Removing leaf L0 leaves theta with branch vertices L1,L2.
 routes=[[1,5,2],[1,7,2],[1,6,3,4,2]]
 core=set().union(*(edges(r) for r in routes))
 require(core==set(U)-{(0,7)} and [len(r)-1 for r in routes]==[2,2,4],'theta plus leaf')
 # A connected edge-potential counterexample, not magic.
 v=[[i] for i in range(8)];z=mul(E,v)
 require(mul(Z,z)==[[0],[0]] and any(t[0] for t in mul(B,z)),'undecorated negative control')
 source=ROOT/'frozen_minimal_union.json'
 frozen=json.loads(source.read_text())
 require([w['vertices'] for w in frozen['minimal_union_witness']]==[['L0','R3','L1','R1','L2','R0','L3','R2'],['L0','R3','L2','R1','L1','R2','L3','R0']],'frozen witness')
 return dict(vertex_order=['L0','L1','L2','L3','R0','R1','R2','R3'],paths=[P,Q],ordered_edges=U,shared_edges=sorted(shared),shared_components=components(range(8),shared),symmetric_difference_edges=sorted(pe^qe),E=E,A=A,reconstruction_T=T,extra_constraints_B=B,cycle_constraints_Z=Z,magic_grid_edge_indices=grid,magic_line_constraints_7=H,independent_magic_constraints_6=H6,edge_parameter_matrix=coeff,ranks=dict(E=7,A=4,Z=2,B=4,Z_and_B=6,magic=6),theta_routes=routes,leaf=[0,7],Hamiltonian_paths_up_to_reversal=paths,negative_control=dict(potentials=list(range(8)),edge_values=[x[0] for x in z],extra_residuals=[x[0] for x in mul(B,z)]),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
if __name__=='__main__':
 out=json.dumps(main(),indent=2,sort_keys=True)+'\n';p=ROOT/'certificate.json'
 if '--check' in sys.argv:require(p.read_text()==out,'stored certificate mismatch')
 else:p.write_text(out)
 print('PASS: exact ranks, integral reconstruction, row-space equality, shared forest, theta core, and four Hamiltonian paths')
