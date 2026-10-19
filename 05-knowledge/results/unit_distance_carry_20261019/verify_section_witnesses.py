#!/usr/bin/env python3
"""Independent exact verification of sharp section-defect witnesses.
No producer imports, no numerical libraries, no SAT/MIP claim.
"""
import json
from pathlib import Path
from itertools import product


def rank_mod(rows,p):
 a=[[x%p for x in row] for row in rows];rank=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(rank,len(a)) if a[i][j]),None)
  if pivot is None:continue
  a[rank],a[pivot]=a[pivot],a[rank]
  inv=pow(a[rank][j],-1,p);a[rank]=[x*inv%p for x in a[rank]]
  for i in range(rank+1,len(a)):
   c=a[i][j];a[i]=[(x-c*y)%p for x,y in zip(a[i],a[rank])]
  rank+=1
  if rank==len(a):break
 return rank


def verify(p):
 D=next(d for d in range(1,p) if pow(d,(p-1)//2,p)==p-1)
 V=list(product(range(p),repeat=2));m=p*p
 if p==3:
  # sqrt(2)'s residue is lifted as 7*sqrt(2), whose norm is -98=1 mod9.
  X=[(a,7*b%m) for a,b in V]
 else:
  z=json.loads(Path(__file__).with_name('defect-candidate-5.json').read_text())
  assert z['p']==p and z['vertices']==[list(v) for v in V]
  X=[((a+p*h[0])%m,(b+p*h[1])%m) for (a,b),h in zip(V,z['h'])]
 assert all((a%p,b%p)==v for v,(a,b) in zip(V,X))
 ix={v:i for i,v in enumerate(V)};E=[];bad=[];defect={}
 for i in range(len(V)):
  for j in range(i+1,len(V)):
   dx=V[j][0]-V[i][0];dy=V[j][1]-V[i][1]
   if (dx*dx-D*dy*dy)%p!=1:continue
   edge=frozenset((i,j));E.append(edge)
   dx=X[j][0]-X[i][0];dy=X[j][1]-X[i][1]
   err=(dx*dx-D*dy*dy-1)%m
   assert err%p==0
   defect[edge]=(err//p)%p
   if err:bad.append(edge)
 U=[v for v in V if (v[0]**2-D*v[1]**2)%p==1]
 dirs=[v for v in U if v<((-v[0])%p,(-v[1])%p)]
 cycles=[]
 for u in dirs:
  remain=set(V)
  while remain:
   a=min(remain)
   line=[((a[0]+t*u[0])%p,(a[1]+t*u[1])%p) for t in range(p)]
   remain.difference_update(line)
   C=[frozenset((ix[line[t]],ix[line[(t+1)%p]])) for t in range(p)]
   assert len(set(C))==p
   assert sum(bool(defect[e]) for e in C)==1
   assert sum(defect[e] for e in C)%p==(-2)%p
   cycles.append(C)
 flat=[e for C in cycles for e in C]
 assert len(flat)==len(set(flat))==len(E) and set(flat)==set(E)
 assert len(bad)==len(cycles)==len(E)//p
 rows=[]
 for edge in E:
  i,j=sorted(edge);dx=V[j][0]-V[i][0];dy=V[j][1]-V[i][1]
  row=[0]*(2*len(V))
  row[2*i]=-2*dx;row[2*j]=2*dx;row[2*i+1]=2*D*dy;row[2*j+1]=-2*D*dy
  rows.append(row)
 rank=rank_mod(rows,p);stress=len(E)-rank;extra=stress-len(cycles)
 assert rank=={3:12,5:47}[p] and extra=={3:0,5:13}[p]
 print(f'q={p}: {len(V)} vertices, {len(E)} edges; exactly {len(bad)} fail; every straight p-cycle has exactly one failure. OPTIMAL by cycle bound.')
 print(f'Rigidity rank {rank}; stress dimension {stress}, with {extra} beyond the {len(cycles)} disjoint-cycle stresses.')
 return {'q':p,'D':D,'vertices':V,'lifted_vertices':X,'edges':len(E),'failures':len(bad),'minimum_defect':len(bad),'cycles':len(cycles),'rigidity_rank':rank,'stress_dimension':stress,'cycle_stress_dimension':len(cycles),'additional_stress_dimension':extra}

if __name__=='__main__':
 result=[verify(3),verify(5)]
 Path(__file__).with_name('section-witnesses-verified.json').write_text(json.dumps(result,indent=2)+'\n')
 print('INDEPENDENT DIRECT NORM CHECKS PASSED')
