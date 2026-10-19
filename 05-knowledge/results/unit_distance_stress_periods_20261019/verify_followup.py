#!/usr/bin/env python3
"""Exact seam/cycle/wheel/global-period certificates. Standard library only."""
from itertools import product
from pathlib import Path
import json


def probe(p):
 V=list(product(range(p),repeat=2));dirs=[(1,0),(0,1),(-1,1)]
 add=lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p)
 edge=lambda a,b:tuple(sorted((a,b)))
 Q=lambda a:a[0]**2+a[0]*a[1]+a[1]**2
 selected=set()
 for d in dirs:
  for z in V:
   if (z[1]==0 if d==(0,1) else z[0]==0):selected.add(edge(z,add(z,d)))
 assert len(selected)==3*p
 cycles=[]
 for d in dirs:
  left=set(V)
  while left:
   v=min(left);line=[add(v,(t*d[0],t*d[1])) for t in range(p)];left.difference_update(line)
   C={edge(line[i],line[(i+1)%p]) for i in range(p)}
   assert len(C&selected)==1;cycles.append(C)
 hexdirs=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
 for v in V:
  ring=[add(v,d) for d in hexdirs]
  rim=sum(edge(ring[i],ring[(i+1)%6]) in selected for i in range(6))
  spokes=sum(edge(v,z) in selected for z in ring)
  assert rim==spokes
 T={}
 for x,y in V:
  v=(x,y);w=add(v,(-1,1));delta=(w[0]-x,w[1]-y)
  carry=((Q(delta)-1)//p)%p
  assert carry==(-int(x==0)-int(y==p-1))%p
  A=lambda z:-2*int(z!=0)
  known=(A(w[1])-A(y)-A(w[0])+A(x))%p
  residual=(-2*int(x==0)-carry-known)%p
  formula=(-3*int(x==0)+2*int(x==1)+2*int(y==0)-int(y==p-1))%p
  assert residual==formula
  T[v]=residual
 assert {sum(T[x,y] for x in range(p))%p for y in range(p)}=={p-1}
 assert {sum(T[x,y] for y in range(p))%p for x in range(p)}=={1}
 return {'p':p,'selected_edges':len(selected),'cycles':len(cycles),'integer_wheels_checked':len(V),'x_period_defect':p-1,'y_period_defect':1},selected


def full29(base):
 p=29;V=list(product(range(p),repeat=2));Q=lambda z:(z[0]**2+z[0]*z[1]+z[1]**2)%p
 mul=lambda a,b:((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0]+a[1]*b[1])%p)
 add=lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p)
 edge=lambda a,b:tuple(sorted((a,b)))
 U={v for v in V if Q(v)==1};assert len(U)==30
 omega=(0,1);hexunits=[];z=(1,0)
 for _ in range(6):hexunits.append(z);z=mul(z,omega)
 assert z==(1,0) and len(set(hexunits))==6
 remaining=set(U);reps=[];selected=set();allcycles=[]
 while remaining:
  u=min(remaining);orbit={mul(u,z) for z in hexunits};remaining.difference_update(orbit);reps.append(u)
  S={edge(mul(u,a),mul(u,b)) for a,b in base};assert len(S)==87
  assert not selected&S;selected.update(S)
  hd=[mul(u,z) for z in hexunits]
  for v in V:
   ring=[add(v,d) for d in hd]
   assert sum(edge(ring[i],ring[(i+1)%6]) in S for i in range(6))==sum(edge(v,z) in S for z in ring)
  for d in hd[:3]:
   left=set(V)
   while left:
    v=min(left);line=[add(v,(t*d[0],t*d[1])) for t in range(p)];left.difference_update(line)
    C={edge(line[i],line[(i+1)%p]) for i in range(p)}
    assert len(C&S)==1;allcycles.append(C)
 assert len(reps)==5 and len(selected)==435 and len(allcycles)==435
 assert len(set().union(*allcycles))==sum(map(len,allcycles))==12615
 return {'p':p,'triangular_family_representatives':reps,'selected_edges':435,'cycle_checks':435,'integer_wheel_checks':4205,'scope':'Explicit support passes these local checks but fails the proved global period certificate; it is not a section.'}


def rank_mod(rows,p):
 a=[[x%p for x in row] for row in rows];r=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];inv=pow(a[r][j],-1,p)
  a[r]=[x*inv%p for x in a[r]]
  for i in range(r+1,len(a)):
   if a[i][j]:
    k=a[i][j];a[i]=[(x-k*y)%p for x,y in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r


def stress_dimensions(p):
 V=list(product(range(p),repeat=2));vi={v:i for i,v in enumerate(V)}
 dirs=[(1,0),(0,1),(-1,1)]
 add=lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p)
 edge=lambda a,b:tuple(sorted((a,b)))
 E=[edge(v,add(v,d)) for d in dirs for v in V];ei={e:i for i,e in enumerate(E)}
 assert len(ei)==3*p*p
 R=[]
 for a,b in E:
  dx,dy=b[0]-a[0],b[1]-a[1];gx,gy=2*dx+dy,dx+2*dy
  row=[0]*(2*p*p)
  row[2*vi[a]]=-gx;row[2*vi[a]+1]=-gy;row[2*vi[b]]=gx;row[2*vi[b]+1]=gy
  R.append(row)
 wheels=[];hexdirs=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
 for v in V:
  ring=[add(v,d) for d in hexdirs];row=[0]*len(E)
  for i in range(6):row[ei[edge(ring[i],ring[(i+1)%6])]]+=1
  for z in ring:row[ei[edge(v,z)]]-=1
  wheels.append(row)
 cycles=[]
 for d in dirs:
  left=set(V)
  while left:
   v=min(left);line=[add(v,(t*d[0],t*d[1])) for t in range(p)];left.difference_update(line)
   row=[0]*len(E)
   for i in range(p):row[ei[edge(line[i],line[(i+1)%p])]]=1
   cycles.append(row)
 # Verify stresses against the rigidity matrix without assuming the proofs.
 for stress in wheels+cycles:
  total=[0]*(2*p*p)
  for a,row in zip(stress,R):
   if a:total=[x+a*y for x,y in zip(total,row)]
  assert all(x%p==0 for x in total)
 rankR=rank_mod(R,p);rankW=rank_mod(wheels,p);rankWC=rank_mod(wheels+cycles,p)
 assert rankR==2*p*p-3 and rankW==p*p-3 and rankWC==p*p+2
 return {'p':p,'rigidity_rank':rankR,'full_stress_dimension':len(E)-rankR,'wheel_span_dimension':rankW,'cycles_and_wheels_span_dimension':rankWC,'missing_period_dimensions':len(E)-rankR-rankWC}


if __name__=='__main__':
 out={'prime_controls':[]}
 for p in [5,11,17,23,29]:
  r,S=probe(p);out['prime_controls'].append(r)
  print(f'p={p}: {r["cycles"]} one-error cycles and {r["integer_wheels_checked"]} integer wheel checks pass; global period defects are -1,+1.')
 out['stress_dimension_controls']=[stress_dimensions(p) for p in [5,11]]
 for r in out['stress_dimension_controls']:print(f'p={r["p"]}: wheel rank {r["wheel_span_dimension"]}, cycles-plus-wheels rank {r["cycles_and_wheels_span_dimension"]}, full stress dimension {r["full_stress_dimension"]}; exactly one missing period.')
 out['full_q29_hostile']=full29(S)
 print('Full q29 hostile: 435 selected edges, all 435 cycle checks and 4205 integer wheel checks pass; no section realizes it.')
 Path(__file__).with_name('followup-exact-results.json').write_text(json.dumps(out,indent=2)+'\n')
