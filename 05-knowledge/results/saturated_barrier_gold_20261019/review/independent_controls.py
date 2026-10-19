#!/usr/bin/env python3
"""Independent review: direct four-cell DP, set-relation construction, exterior controls.
Does not import packet code. Standard library only.
"""
import hashlib,json,random
from fractions import Fraction
from functools import lru_cache
from itertools import product,combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parent
Q=(0,0,0,0,7,13,14,12,206,5,685,535)
S=frozenset((1,4,5,6,7,9))
R=frozenset((i,j) for j,p in enumerate(Q) for i in range(12) if p&(1<<i))
checks=0

def check(test,label):
 global checks
 checks+=1
 if not test: raise RuntimeError(label)

def closure(vertices,edges):
 rr=set(edges)
 for k in vertices:
  rr.update((i,j) for i in vertices for j in vertices if (i,k) in rr and (k,j) in rr)
 if any((i,i) in rr for i in vertices): return None
 return frozenset(rr)

def incomparable(rr,a,b): return a!=b and (a,b) not in rr and (b,a) not in rr

def counts(vertices,rr):
 """One vector-valued DP, not four separate constrained-extension counts."""
 ix={q:i for i,q in enumerate(vertices)}
 predecessors=[sum(1<<ix[a] for a in vertices if (a,b) in rr) for b in vertices]
 final=(1<<len(vertices))-1
 @lru_cache(None)
 def rec(mask):
  if mask==final:return (1,0,0,0)
  ans=[0,0,0,0]
  for i,q in enumerate(vertices):
   if mask&(1<<i) or predecessors[i]&~mask:continue
   flag=0
   if q==9 and not mask&(1<<ix[7]):flag|=2
   if q==5 and not mask&(1<<ix[4]):flag|=1
   for c,v in enumerate(rec(mask|(1<<i))):ans[c|flag]+=v
  return tuple(ans)
 return rec(0)

def model(a,b,deleted=()):
 core=tuple(q for q in range(12) if q not in deleted)
 h=max((a[q] for q in core),default=0);j=max((b[q] for q in core),default=0)
 lower=tuple(range(12,12+h));upper=tuple(range(12+h,12+h+j));vertices=core+lower+upper
 edges={(x,y) for x,y in R if x in core and y in core}
 edges.update(zip(lower,lower[1:]));edges.update(zip(upper,upper[1:]))
 edges.update((x,y) for x in lower for y in upper)
 edges.update((x,q) for i,x in enumerate(lower) for q in core if i<h-a[q])
 edges.update((q,y) for i,y in enumerate(upper) for q in core if i>=b[q])
 rr=closure(vertices,edges)
 check(rr is not None,'chain construction has no cycle')
 check(frozenset((x,y) for x,y in rr if x in core and y in core)==frozenset((x,y) for x,y in R if x in core and y in core),'induced retained core')
 return vertices,rr

def profiles(a,b,c,f,g,h):return ((a,0,b,c,0,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,f,0,g,h))
def mm(cc):
 total=sum(cc)
 return tuple(total-cc[(r//2)*2]-cc[(r//2)*2+1]-cc[r] for r in range(4))

def modelcounts(a,b,deleted=()):return counts(*model(a,b,deleted))

A=sorted(profiles(a,b,c,0,0,0)[0] for a,b,c in product(range(2),range(5),range(2)))
B=sorted(profiles(0,0,0,f,g,h)[1] for f,g,h in product(range(2),range(3),range(2)))
minimum=[Fraction(1)]*4;first=[];bounded=[];atom_digest=hashlib.sha256()
for a,b in product(A,B):
 vertices,rr=model(a,b);cc=counts(vertices,rr);E=sum(cc)
 check(max(sum(incomparable(rr,x,y) for y in vertices) for x in vertices)<=7,'bounded degree')
 minimum=[min(x,Fraction(m,E)) for x,m in zip(minimum,mm(cc))]
 first.append(Fraction(cc[0]+cc[1],E));bounded.append([list(a),list(b),E,list(cc)])
 atom_digest.update((json.dumps(bounded[-1],separators=(',',':'))+'\n').encode())
check(list(map(str,minimum))==['78113/453252','25363/180473','50695/453252','4768/24065'],'four bounded minima')
check((str(min(first)),str(max(first)))==('66465/151084','84619/151084'),'first comparison endpoints')

free_min=[Fraction(1)]*4;free=[];free_digest=hashlib.sha256();generators={}
for a,c,f,h in product(range(2),repeat=4):
 for lower,upper in product(range(3),repeat=2):
  aa,bb=profiles(a,lower%2,c,f,upper%2,h)
  deleted=tuple(q for q,state in [(2,lower),(10,upper)] if state==2)
  cc=modelcounts(aa,bb,deleted);E=sum(cc)
  generators[a,c,f,h,lower,upper]=cc
  free_min=[min(x,Fraction(m,E)) for x,m in zip(free_min,mm(cc))]
  free.append([a,c,f,h,lower,upper,E,list(cc)])
  free_digest.update((json.dumps(free[-1],separators=(',',':'))+'\n').encode())
check(list(map(str,free_min))==['771/6245','3231/24464','386/6245','3116/16323'],'four free-port minima')

# All sixteen bounded coordinate choices, independently count larger free-cap models.
large_cap_tests=0
for a,c,f,h in product(range(2),repeat=4):
 for b,g in [(0,5),(5,0),(1,6),(6,1),(2,4),(4,3),(5,5)]:
  aa,bb=profiles(a,b,c,f,g,h);actual=modelcounts(aa,bb)
  termsL=[(b,1)] if b<=1 else [(1,1),(2,b-1)]
  termsU=[(g,1)] if g<=1 else [(1,1),(2,g-1)]
  expected=tuple(sum(x*y*generators[a,c,f,h,i,j][r] for i,x in termsL for j,y in termsU) for r in range(4))
  check(actual==expected,'all-cell free-cap cone identity beyond tabulated breakpoint');large_cap_tests+=1

# Derive all four polynomial coefficients from base and actual deletion models.
poly_cc=[]
c00=generators[1,0,1,0,1,1];c10=generators[1,0,1,0,2,1]
c01=generators[1,0,1,0,1,2];c11=generators[1,0,1,0,2,2]
for r in range(4):poly_cc.append([c11[r],c10[r]+c01[r]-2*c11[r],c00[r]-c10[r]-c01[r]+c11[r]])
Epoly=[sum(p[k] for p in poly_cc) for k in range(3)]
Mpoly=[Epoly[k]-2*poly_cc[2][k]-poly_cc[3][k] for k in range(3)]
check(Epoly==[12490,78784,122888],'derived E polynomial')
check(Mpoly==[772,8213,18243],'derived FT margin polynomial')
strict_diff=[6245*m-386*e for m,e in zip(Mpoly,Epoly)]
check(strict_diff[0]==0 and min(strict_diff[1:])>0,'symbolic strictness above infimum')

# Exhaust every induced one-point extension type of Q comparable with all S.
# No structural orientation is assumed before constructing the candidate relation.
one_point_types=0;thin_types=0
outside=tuple(q for q in range(12) if q not in S)
for directions in product((-1,1),repeat=len(S)):
 srel=dict(zip(sorted(S),directions))
 for rest in product((-1,0,1),repeat=len(outside)):
  d=srel|dict(zip(outside,rest));rr=set(R)
  rr.update((12,q) for q,v in d.items() if v==-1)
  rr.update((q,12) for q,v in d.items() if v==1)
  if closure(tuple(range(13)),rr)!=frozenset(rr):continue
  one_point_types+=1
  orient={d[q] for q in S};check(len(orient)==1,'connected barrier uniform orientation, exhaustive extensions')
  direction=next(iter(orient));check(all(v in (0,direction) for v in d.values()),'domination forces one-sided type, exhaustive extensions')
  if max(sum(incomparable(rr,x,y) for y in range(13)) for x in range(13))<=7:thin_types+=1

# Actual non-chain exteriors: generate varied A/B posets, close relations,
# condition on each exterior order, and test exact simultaneous mixture identity.
def extensions(vertices,rr):
 def go(word,left):
  if not left:
   yield tuple(word);return
  for q in sorted(left):
   if not any((p,q) in rr for p in left):yield from go(word+[q],left-{q})
 return list(go([],set(vertices)))

rng=random.Random(20261019);exterior_cases=[];nonchain=0
for trial in range(48):
 na=2+trial%3;nb=2+(trial//3)%3
 va=tuple(range(12,12+na));vb=tuple(range(12+na,12+na+nb));vs=tuple(range(12))+va+vb
 edges=set(R);edges.update((x,y) for x in va for y in vb)
 for part in (va,vb):
  edges.update((x,y) for i,x in enumerate(part) for y in part[i+1:] if rng.randrange(3)==0)
 for x in va:
  edges.update((x,q) for q in S|{8,10,11})
  edges.update((x,q) for q in (0,2,3) if rng.randrange(2)==0)
 for y in vb:
  edges.update((q,y) for q in S|{0,2,3})
  edges.update((q,y) for q in (8,10,11) if rng.randrange(2)==0)
 rr=closure(vs,edges);check(rr is not None and frozenset((x,y) for x,y in rr if x<12 and y<12)==R,'arbitrary exterior preserves induced Q')
 actual=counts(vs,rr);expected=[0]*4;orders=0
 for sigma,tau in product(extensions(va,rr),extensions(vb,rr)):
  aa=tuple(len(sigma)-1-max((i for i,x in enumerate(sigma) if (x,q) in rr),default=-1) for q in range(12))
  bb=tuple(min((i for i,y in enumerate(tau) if (q,y) in rr),default=len(tau)) for q in range(12))
  cc=modelcounts(aa,bb);orders+=1
  expected=[x+y for x,y in zip(expected,cc)]
 check(actual==tuple(expected),'conditioned mixture identity on actual non-chain exteriors')
 if any(incomparable(rr,x,y) for part in (va,vb) for x,y in combinations(part,2)):nonchain+=1
 exterior_cases.append({'vertices':len(vs),'exterior_order_pairs':orders,'cells':actual})

# End-to-end examples meeting the actual theorem hypotheses, with non-chain exteriors.
end_to_end=[]
for unbounded in (False,True):
 for trial in range(16):
  na=(5+trial%2) if unbounded else (2+trial%3)
  nb=(5+(trial//2)%2) if unbounded else (2+(trial//3)%3)
  va=tuple(range(12,12+na));vb=tuple(range(12+na,12+na+nb));vs=tuple(range(12))+va+vb
  edges=set(R);edges.update((x,y) for x in va for y in vb)
  for part in (va,vb):
   edges.update((x,y) for i,x in enumerate(part) for y in part[i+1:] if rng.randrange(4)==0)
  edges.update((x,q) for x in va for q in S|{8,10,11})
  edges.update((q,y) for y in vb for q in S|{0,2,3})
  for q in (0,3):
   missed=rng.choice(va+(-1,));edges.update((x,q) for x in va if x!=missed)
  for q in (8,11):
   missed=rng.choice(vb+(-1,));edges.update((q,y) for y in vb if y!=missed)
  if not unbounded:
   edges.update((x,2) for x in va if rng.randrange(3)==0)
   missed=set(rng.sample(vb,min(2,len(vb))))
   edges.update((10,y) for y in vb if y not in missed)
  rr=closure(vs,edges)
  check(rr is not None and frozenset((x,y) for x,y in rr if x<12 and y<12)==R,'end-to-end induced copy')
  check(all(not incomparable(rr,s,x) for s in S for x in va+vb),'end-to-end saturated barrier comparability')
  check(all(sum(incomparable(rr,q,x) for x in va+vb)<=1 for q in (0,3,8,11)),'bounded ports meet Theorem B')
  degree=max(sum(incomparable(rr,x,y) for y in vs) for x in vs)
  cc=counts(vs,rr);E=sum(cc);min_margin=min(mm(cc))
  if unbounded:
   check(degree>7,'free-port examples are outside seven-thin class')
   check(6245*min_margin>386*E,'end-to-end Theorem B sharp bound')
  else:
   check(degree<=7,'end-to-end Theorem A seven-thin hypothesis')
   check(453252*min_margin>=50695*E,'end-to-end Theorem A sharp bound')
  end_to_end.append({'family':'B beyond seven-thin' if unbounded else 'A seven-thin','order':len(vs),'maximum_incomparability_degree':degree,'cells':cc})

original=json.loads((ROOT/'original-normal.json').read_text())
check(atom_digest.hexdigest()==original['bounded_certificate']['count_table_sha256'],'entire independently rebuilt bounded count table matches packet hash')
check(free_digest.hexdigest()==original['unbounded_certificate']['count_table_sha256'],'entire independently rebuilt deletion count table matches packet hash')
result={'status':'PASS','checks':checks,'independent_method':'set-relation construction and one vector-valued ideal DP, no packet imports','bounded_atoms':len(bounded),'bounded_minima':list(map(str,minimum)),'first_probability_interval':[str(min(first)),str(max(first))],'bounded_table_sha256':atom_digest.hexdigest(),'deletion_generators':len(free),'free_minima':list(map(str,free_min)),'free_table_sha256':free_digest.hexdigest(),'large_free_cap_models':large_cap_tests,'sharpness_cell_polynomials_descending':poly_cc,'extension_polynomial_descending':Epoly,'FT_margin_polynomial_descending':Mpoly,'6245_times_margin_minus_386_times_E_descending':strict_diff,'exhaustive_valid_one_point_types_comparable_to_S':one_point_types,'seven_thin_one_point_types':thin_types,'actual_exterior_cases':len(exterior_cases),'actual_non_chain_cases':nonchain,'exterior_cases':exterior_cases,'end_to_end_theorem_cases':end_to_end}
(ROOT/'independent-results.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'independent-count-tables.json').write_text(json.dumps({'bounded':bounded,'free':free},indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('exterior_cases','end_to_end_theorem_cases')},indent=2))
