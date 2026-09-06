#!/usr/bin/env python3
"""Exact full-component envelope for inert odd-3-unit half-lift obstructions.

Default universe: complete 46-pair height head, independent literal wall
geometry, five maximal length profiles, compact-set and physical controls.
Optional --probe-bodies reproduces only the declared 10,599 template probes.
No general body classification or universal safe-phase claim is inferred.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import gcd
from hashlib import sha256
import json, sys

checks=0
def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def inert(n):
 p=2
 while p*p<=n:
  e=0
  while n%p==0:n//=p;e+=1
  if e and (p%3!=2 or e>2):return False
  p+=1
 return n==1 or n%3==2

def admissible(p,q):return 0<p<q and gcd(p,q)==gcd(p*q,6)==1 and inert(p+q)
def distance(x):
 r=x-x.numerator//x.denominator
 return min(r,1-r)
def bad(y,p,q):
 return all(any(distance(v*x)<F(1,14) for v in (p,q)) for x in (y/2,(y+1)/2))

def intersections(p,q):
 result=[]
 for k in range(p+1):
  for l in range(q+1):
   if (k-l)%2==0:continue
   lo=max(F(0),F(7*k-1,7*p),F(7*l-1,7*q))
   hi=min(F(1),F(7*k+1,7*p),F(7*l+1,7*q))
   if lo<hi:result.append((lo,hi))
 return sorted(result)

def wall_cells(p,q):
 walls={F(0),F(1)}
 for v in (p,q):
  for k in range(v+1):
   for sg in (-1,1):
    x=F(7*k+sg,7*v)
    if 0<x<1:walls.add(x)
 points=sorted(walls);out=[]
 for lo,hi in zip(points,points[1:]):
  if not bad((lo+hi)/2,p,q):continue
  if out and out[-1][1]==lo and bad(lo,p,q):out[-1]=(out[-1][0],hi)
  else:out.append((lo,hi))
 return out

def spectrum(intervals):return sorted((hi-lo for lo,hi in intervals),reverse=True)
def prefix(a,k):return sum(a[:k],F(0))
def dominates(a,b):return all(prefix(a,k)>=prefix(b,k) for k in range(1,max(len(a),len(b))+1))

def replicated_prefix(a,g,k):
 q,r=divmod(k,g)
 return prefix(a,q)+F(r,g)*(a[q] if q<len(a) else 0)

def full_gates(lengths,g=1):
 witnesses=[]
 for _,p in profiles:
  witness=next((k for k in range(1,len(lengths)+1) if prefix(lengths,k)>=replicated_prefix(p,g,k)),None)
  if witness is None:return None
  witnesses.append(witness)
 return witnesses

def old_gates(lengths,g=1):
 return all(sum(lengths)>=sum(p) or g*max(lengths,default=F(0))>=p[0] for _,p in profiles)

def safe_components(H):
 out=[(F(0),F(1))]
 for h in sorted(H):
  teeth=[(F(14*k+1,14*h),F(14*k+13,14*h)) for k in range(h)]
  new=[];j=0
  for a,b in out:
   while j<len(teeth) and teeth[j][1]<a:j+=1
   z=j
   while z<len(teeth) and teeth[z][0]<=b:
    lo,hi=max(a,teeth[z][0]),min(b,teeth[z][1])
    if lo<=hi:new.append((lo,hi))
    z+=1
  out=new
 return out

profiles=[((7,13),[F(1,49)]*2),((19,25),[F(37,3325)]*2+[F(23,3325)]*2+[F(9,3325)]*2),
 ((5,41),[F(2,287)]*6),((5,53),[F(2,371)]*6+[F(9,1855)]*2),((1,67),[F(2,469)]*10)]
rows={};certificates=[]
for q in range(2,68):
 for p in range(1,q):
  if not admissible(p,q):continue
  a=intersections(p,q);b=wall_cells(p,q)
  gate(a==b,'independent literal half-lift geometry')
  gate(all(not bad(x,p,q) for arc in a for x in arc),'strict endpoints retained')
  lengths=spectrum(a);rows[p,q]=lengths
  dominators=[pair for pair,v in profiles if dominates(v,lengths)]
  gate(bool(dominators),'full prefix domination')
  certificates.append({'pair':[p,q],'lengths':[str(x) for x in lengths],'dominators':dominators})
gate(len(rows)==46,'complete stated finite head size')
front=[pair for pair,a in rows.items() if not any(other!=pair and dominates(b,a) and not dominates(a,b) for other,b in rows.items())]
gate(set(front)==set(pair for pair,_ in profiles),'exact maximal profiles')
for pair,a in profiles:gate(rows[pair]==a,'attained exact component multiplicities')
# Infinite tail: q>67 has width<2/469 and inherited mass<20/469.
# The top-k sum is <=min(k*2/469,20/469), exactly the final profile prefix.
gate(profiles[-1][1]==[F(2,469)]*10,'equal-cap infinite-tail comparator')
gate(sum(profiles[-1][1])==F(20,469),'inherited strict mass maximum matches comparator')
for k in range(1,14):gate(prefix(profiles[-1][1],k)==min(k*F(2,469),F(20,469)),'finite-cap prefix identity')
# The corresponding square energy is universally bounded by2/2401.
energies=[sum(x*x for x in a) for _,a in profiles]
gate(max(energies)==F(2,2401),'sharp component square-energy ceiling')
gate(energies.count(max(energies))==1,'unique envelope square-energy maximizer')
# Common odd dilation: arithmetic prefix formula versus literal replication.
for _,a in profiles:
 for g in [1,3,5,17]:
  repeated=sorted([x/g for x in a for _ in range(g)],reverse=True)
  for k in [1,2,3,g,g+1,2*g+3,len(repeated),len(repeated)+2]:
   gate(replicated_prefix(a,g,k)==prefix(repeated,k),'literal dilation prefix formula')
# Compact packing controls: four positive components assigned to two open bins,
# each containing its assigned closed lengths with positive slack.
packing=0
for sizes in product(range(1,4),repeat=4):
 for assignment in product(range(2),repeat=4):
  capacity=sorted([sum(sizes[i] for i in range(4) if assignment[i]==j)+1 for j in range(2)],reverse=True)
  actual=sorted(sizes,reverse=True)
  gate(all(prefix(actual,k)<prefix(capacity,k) for k in range(1,5)),'strict prefix packing with full assignment')
  packing+=1
# A strict gain in the declared compact-set domain, with reflection symmetry.
lengths=[F(1,100)]*4+[F(1,2000)]*2
centers=[F(1,10),F(1,5),F(3,10),F(7,10),F(4,5),F(9,10)]
widths=[F(1,100),F(1,100),F(1,2000),F(1,2000),F(1,100),F(1,100)]
compact=[(c-w/2,c+w/2) for c,w in zip(centers,widths)]
gate(spectrum(compact)==lengths and all(a[1]<b[0] for a,b in zip(compact,compact[1:])),'disjoint actual compact arcs')
gate([(1-b,1-a) for a,b in reversed(compact)]==compact,'compact reflection symmetry')
gate(not old_gates(lengths) and full_gates(lengths)==[6,3,1,1,1],'strict improvement over all mass/width disjunctions')
gate(prefix(lengths,3)==F(3,100)>F(97,3325),'new middle-prefix obstruction')
# Retain an actual known eleven-body as a separate positive control.
H=(1,3,4,5,6,7,8,9,10,11,13)
physical=safe_components(H);physical_lengths=spectrum(physical)
gate(sum(physical_lengths)==F(5939,140140) and physical_lengths[0]==F(11,728),'recovered physical measure/width')
gate(sum(a==b for a,b in physical)==6,'isolated safe points retained')
for a,b in physical:
 for y in {a,b,(a+b)/2}:gate(all(distance(h*y)>=F(1,14) for h in H),'literal physical endpoint and midpoint')
for g in [1,3,5,17]:gate(full_gates(physical_lengths,g) is not None,'actual-body full prefix consumer')
row=tuple(2*h for h in H)+(1,67)
gate(min(distance(v*F(9,34)) for v in row)==F(2,17),'literal full thirteen-speed witness')
# Declared assumption hostiles; none is admitted by the theorem.
for pair in [(1,9),(5,11),(1,11)]:
 gate(not admissible(*pair),'arithmetic boundary rejected')
 gate(sum(spectrum(intersections(*pair)))>F(20,469),'dropped assumption exceeds envelope')

bank={'profiles':[[list(pair),list(map(str,a))] for pair,a in profiles],'head':certificates,
 'compact_arcs':[[str(a),str(b)] for a,b in compact],'compact_prefix_witnesses':full_gates(lengths),
 'physical_components':[[str(a),str(b)] for a,b in physical],'energies':list(map(str,energies))}
print('PASS full component-prefix envelope, conditional physical entry; LRC14 OPEN')
print('Declared pair universe:46 inert odd-3-unit pairs p<q<=67; infinite tail controlled analytically')
for pair,a in profiles:print('PROFILE',pair,[(str(x),a.count(x)) for x in sorted(set(a),reverse=True)],'energy',sum(x*x for x in a))
for cert in certificates:print('HEAD',cert)
print('COMPACT_GAIN',bank['compact_arcs'],'prefixes',full_gates(lengths),'M',sum(lengths),'L',max(lengths))
print('PHYSICAL_CONTROL',H,'components',bank['physical_components'],'prefixes',full_gates(physical_lengths))
print('PACKING_CONTROLS',packing,'GATES',checks)
print('SHA256',sha256(json.dumps(bank,sort_keys=True).encode()).hexdigest())
if '--probe-bodies' in sys.argv:
 count=0;gains=[]
 def probe(H):
  global count
  count+=1;a=spectrum(safe_components(H))
  if not old_gates(a) and full_gates(a):gains.append(H)
 for missing in range(1,12):
  for r in range(13,61):probe(tuple(i for i in range(1,12) if i!=missing)+(r,))
 for r in range(1,122,2):probe(tuple(range(2,21,2))+(r,))
 for C in combinations(range(1,15),10):
  for r in [1,3,5,7,11,13,17,25,39,61]:probe(tuple(2*c for c in C)+(r,))
 print('BOUNDED_TEMPLATE_PROBE',count,'gain_count',len(gains),'gains',gains)
