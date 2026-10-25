#!/usr/bin/env python3
"""Independent seven/nine root surface arrangement and invariants audit."""
from itertools import combinations,product
from functools import reduce
from math import gcd,comb
from collections import Counter
from pathlib import Path
import sympy as S
import json
checks=0
def check(v,m):
 global checks
 checks+=1
 if not v:raise RuntimeError(m)
def norm(v):
 v=tuple(map(int,v));g=reduce(gcd,(abs(x) for x in v));v=tuple(x//g for x in v)
 return tuple(-x for x in v) if next(x for x in v if x)!=abs(next(x for x in v if x)) else v
def modrank(mat,p):
 a=[[int(v)%p for v in row]for row in mat];r=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];inv=pow(a[r][j],-1,p);a[r]=[v*inv%p for v in a[r]]
  for i in range(len(a)):
   if i!=r:
    f=a[i][j];a[i]=[(v-f*w)%p for v,w in zip(a[i],a[r])]
  r+=1
 return r
rows=[]
for name,forms,expected_triples,expected_nodes,kcanonical,expected_pg in [
 ('H',[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)],3,24,1,7),
 ('M',[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1)],8,256,3,111)]:
 A=S.Matrix(forms);m=len(forms);points={}
 for i,j in combinations(range(m),2):
  p=norm(S.Matrix(forms[i]).cross(S.Matrix(forms[j])))
  van=tuple(k for k,f in enumerate(forms)if sum(x*y for x,y in zip(f,p))==0)
  points[p]=van
 triples={p:v for p,v in points.items()if len(v)==3}
 doubles={p:v for p,v in points.items()if len(v)==2}
 check(max(map(len,points.values()))==3,'maximum multiplicity')
 check(len(triples)==expected_triples,'triple count')
 check(len(doubles)==12,'double count')
 check(sum(comb(len(v),2)for v in points.values())==comb(m,2),'all pair incidences')
 relations=[]
 for p,v in triples.items():
  kernel=A.extract(v,range(3)).T.nullspace()
  check(len(kernel)==1,'triple relation unique')
  relation=kernel[0]
  check(all(c!=0 for c in relation),'nondegenerate local quadratic')
  relations.append([int(x) if x.q==1 else str(x) for x in relation])
 check(len(triples)*2**(m-4)==expected_nodes,'projective root sheet count')
 # Independent finite check over F17, where +/-1,+/-2 all square.
 pfield=17;base=[0,1,3];B=A.extract(base,range(3));C=[]
 for j in range(m):
  if j in base:continue
  coeff=A[j,:]*B.inv();row=[S.Rational(0)]*m;row[j]=1
  for i,c in zip(base,coeff):row[i]-=c
  check(all(S.Rational(x).q==1 for x in row),'integral CI elimination')
  C.append(list(map(int,row)))
 # Jacobian rank depends only on which roots vanish; test every arrangement point.
 for basepoint,zero_set in points.items():
  kept=[j for j in range(m)if j not in zero_set]
  restricted=S.Matrix(C).extract(range(m-3),kept)
  check(restricted.rank()==m-3-(len(zero_set)-2),'all crossing Jacobian ranks')
 for zero_set in [(),*((j,)for j in range(m))]:
  kept=[j for j in range(m)if j not in zero_set]
  check(S.Matrix(C).extract(range(m-3),kept).rank()==m-3,'smooth generic/simple branch')
 found=set()
 for p,v in triples.items():
  values=[sum(x*y for x,y in zip(f,p))%pfield for f in forms]
  choices=[[r for r in range(pfield)if r*r%pfield==z]for z in values]
  check(all(choices),'all branch root choices over F17')
  first=next(i for i,z in enumerate(values)if z)
  choices[first]=choices[first][:1]
  for roots in product(*choices):
   check(all(sum(c*r*r for c,r in zip(row,roots))%pfield==0 for row in C),'root equations')
   J=[[2*c*r%pfield for c,r in zip(row,roots)]for row in C]
   check(modrank(J,pfield)==m-4,'node Jacobian corank one')
   # Normalize geometric projective root vector for distinct-point test.
   inv=pow(roots[first],-1,pfield)
   found.add(tuple(r*inv%pfield for r in roots))
 check(len(found)==expected_nodes,'distinct rational nodes F17')
 # Complete intersection Hilbert coefficients via (1-T²)^(m-3)/(1-T)^m.
 degree=2**(m-3)
 pg=sum((-1)**j*comb(m-3,j)*comb(kcanonical-2*j+m-1,m-1) for j in range(kcanonical//2+1))
 check(2*(m-3)-m==kcanonical,'adjunction index')
 check(pg==expected_pg,'geometric genus Hilbert coefficient')
 rows.append(dict(surface=name,forms=forms,triples=[{'point':p,'vanishing_forms':v}for p,v in triples.items()],double_points=len(doubles),geometric_nodes=expected_nodes,verified_distinct_nodes_F17=len(found),canonical_index=kcanonical,degree=degree,K_squared=kcanonical*kcanonical*degree,p_g=pg,local_quadric_relations=relations))
# Exhaust every labelled subset of the nine forms independently.
allforms=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1),(1,1,-1),(1,-1,1)]
subset_counts={}
for n in range(4,10):
 counts=Counter()
 for indices in combinations(range(9),n):
  selected=[allforms[j]for j in indices]
  check(S.Matrix(selected).rank()==3,'all coordinate subsets span')
  crossings={}
  for i,j in combinations(range(n),2):
   point=norm(S.Matrix(selected[i]).cross(S.Matrix(selected[j])))
   multiplicity=sum(sum(x*y for x,y in zip(line,point))==0 for line in selected)
   crossings[point]=multiplicity
  check(max(crossings.values())<=3,'no fourfold subset crossing')
  counts[sum(z==3 for z in crossings.values())]+=1
 subset_counts[n]=dict(sorted(counts.items()))
expected={5:{0:28,1:76,2:22},6:{0:2,1:20,2:46,3:16},7:{2:4,3:18,4:12,5:2},8:{4:1,5:4,6:4},9:{8:1}}
for n,census in expected.items():check(subset_counts[n]==census,'complete labelled subset census')
# Stratified complex Euler characteristic, with reduced covers on branch strata.
euler_rows=[]
for n in range(4,10):
 k=n-6
 pg=0 if k<0 else sum((-1)**j*comb(n-3,j)*comb(k-2*j+n-1,n-1) for j in range(k//2+1))
 K2=k*k*2**(n-3)
 expected_euler=2**(n-4)*(n*n-9*n+24)
 check(expected_euler==12*(1+pg)-K2,'Noether consistency')
 for triples in subset_counts[n]:
  doubles=comb(n,2)-3*triples
  U=3-2*n+doubles+2*triples
  L=2*n-2*doubles-3*triples
  singular=2**(n-1)*U+2**(n-2)*L+2**(n-3)*doubles+2**(n-4)*triples
  resolved=singular+2**(n-4)*triples
  check(resolved==expected_euler,'Euler strata and node resolution')
 euler_rows.append({'n':n,'p_g':pg,'K_squared':K2,'resolved_euler':expected_euler})
out={'status':'PASS','always_active_checks':checks,'surfaces':rows,'all_coordinate_subset_triple_counts':subset_counts,'euler_crosschecks':euler_rows}
Path(__file__).with_name('torification_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
