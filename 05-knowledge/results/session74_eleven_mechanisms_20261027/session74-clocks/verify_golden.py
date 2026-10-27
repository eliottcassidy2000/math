#!/usr/bin/env python3
"""Exact all-period carry certificate + finite positive and hostile controls."""
from collections import defaultdict
from fractions import Fraction
import hashlib,json
CHECKS=0
def check(v,msg):
 global CHECKS
 CHECKS+=1
 if not v: raise RuntimeError(msg)
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):
 a,b=x;c,d=y
 return (a*c+b*d,a*d+b*c+b*d)
def power(n):
 z=(1,0)
 for _ in range(n):z=mul(z,(0,1))
 return z
def sign_phi(x):
 a,b=x;q=2*a+b;r=b
 if not r:return (q>0)-(q<0)
 if not q:return (r>0)-(r<0)
 if q*r>0:return (q>0)-(q<0)
 s=q*q-5*r*r
 return ((q>0)-(q<0))*((s>0)-(s<0))
def le(x,y):return sign_phi((x[0]-y[0],x[1]-y[1]))<=0
def absle(x,y):return le(x,y) and le((-x[0],-x[1]),y)
def conjugate(x):return (x[0]+x[1],-x[1])
def ring(n):
 a,b=power(n);a-=1
 det=a*(a+b)-b*b
 return a,b,abs(det)
def key(x,n):
 a,b,d=ring(n);u,v=x
 return (((a+b)*u-b*v)%d,(-b*u+a*v)%d)
def reading(w):
 z=(0,0)
 for bit in w:z=add(mul(z,(0,1)),(bit,0))
 return z
def legal(w):return all(not(w[i] and w[(i+1)%len(w)]) for i in range(len(w)))
def words(n):
 def rec(w):
  if len(w)==n:
   if legal(w):yield tuple(w)
   return
  yield from rec(w+[0])
  if not w or not w[-1]:yield from rec(w+[1])
 yield from rec([])
def fixed(w):
 a,b,h=1,0,0
 for t in w:
  if t:a,b=3*a,3*b+2**h
  else:h+=1
 return Fraction(b,2**h-a)
def standard(x):return 3*x+1 if x.numerator%2 else x/2
# Coarse complete carry region; mathematical bound proves |b| <= 1 and |a| <= 3.
coarse=[(a,b) for a in range(-3,4) for b in range(-1,2) if absle((a,b),(0,1)) and absle(conjugate((a,b)),(1,1))]
expected=[(-2,1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(2,-1)]
check(coarse==expected,'nine coarse carries')
fine=[x for x in coarse if absle(x,(1,0))]
check(len(fine)==7,'seven sharp carries')
vertices=[(*z,p,q) for z in coarse for p in (0,1) for q in (0,1)]
vs=set(vertices);edges={x:[] for x in vertices}
for a,b,p,q in vertices:
 for u in (0,1):
  for v in (0,1):
   y=(b-u+v,a+b,u,v)
   if not p*u and not q*v and y in vs:edges[a,b,p,q].append(y)
# Reachability SCCs, rather than numerical eigenvalues or sampled periods.
reach={x:set(edges[x]) for x in vertices}
for _ in range(len(vertices)):
 changed=False
 for x in vertices:
  new=set().union(*(reach[y] for y in tuple(reach[x]))) if reach[x] else set()
  if not new <= reach[x]:reach[x]|=new;changed=True
 if not changed:break
cyclic=[];seen=set()
for x in vertices:
 if x in reach[x] and x not in seen:
  s=sorted(y for y in vertices if y in reach[x] and x in reach[y]);cyclic.append(s);seen.update(s)
cyclic=sorted(cyclic)
expected_scc=sorted([
 sorted([(0,0,0,0),(0,0,1,1)]),
 sorted([(1,0,0,0),(-1,1,1,0)]),
 sorted([(-1,0,0,0),(1,-1,0,1)]),
 sorted([(-2,1,1,0),(2,-1,0,1)])])
check(cyclic==expected_scc,'four exact recurrent SCCs')
# Independent finite word/ring census (not needed for the all-n proof).
census=[];rational_checks=0
for n in range(1,23):
 groups=defaultdict(list);total=0
 for w in words(n):
  total+=1;groups[key(reading(w),n)].append(w)
  # Direct source legality and return, exhaustively through length 16.
  if n<=16:
   x=fixed(w);y=x
   for bit in w:
    check(y.denominator%2==1 and y.numerator%2==bit,'rational parity')
    y=standard(y)
   check(y==x,'rational cycle');rational_checks+=1
   shifted=w[1:]+w[:1]
   check(key(reading(shifted),n)==key(mul(reading(w),(0,1)),n),'equivariance')
 expected_bad=[] if n%2 else sorted([tuple([0]*n),tuple([0,1]*(n//2)),tuple([1,0]*(n//2))])
 bad=[sorted(v) for v in groups.values() if len(v)>1]
 check(bad==([] if n%2 else [expected_bad]),'only even zero collision')
 check(len(groups)==ring(n)[2],'surjective exact quotient census')
 census.append({'n':n,'words':total,'ring_order':ring(n)[2],'largest_fibre':max(map(len,groups.values()))})
check(reading((1,0,0))==reading((0,1,1)),'internal-11 hostile')
check(not legal((0,1,1)),'internal hostile rejected')
check(key(reading((1,0,1)),3)==key(reading((0,1,0)),3),'cyclic seam hostile')
check(not legal((1,0,1)),'cyclic seam hostile rejected')
# All period-five sources and their source-preserving lifts.
base=list(words(5));base_labels=sorted((str(fixed(w)),(reading(w)[0]+4*reading(w)[1])%11) for w in base)
lifts=[]
for r in range(1,23):
 q=(0,0)
 for j in range(r):q=add(q,power(5*j))
 lifted=[key(reading(w*r),5*r) for w in base]
 projected=[(reading(w*r)[0]+4*reading(w*r)[1])%11 for w in base]
 check(len(set(lifted))==11,'lift injection of all eleven sources')
 for w in base:
  check(reading(w*r)==mul(reading(w),q),'exact repeated-word lift')
  check(key(reading(w*r),5)==key((r*reading(w)[0],r*reading(w)[1]),5),'return multiplies by r')
 masses=defaultdict(int)
 for i,label in enumerate(projected):masses[label]+=i+1
 check(sum(masses.values())==66,'source mass preserved')
 if r%11==0:check(masses=={0:66},'elevenfold lift has blind projection')
 lifts.append({'r':r,'period':5*r,'lifted_distinct':len(set(lifted)),'projected_distinct':len(set(projected)),'ambient_order':ring(5*r)[2],'projected_mass':dict(sorted(masses.items()))})
# Exact Terras clock and integer shadows, not standard-clock period 11.
def terras(x):return (3*x+1)//2 if x%2 else x//2
cycle=[];x=-17
while x not in cycle:cycle.append(x);x=terras(x)
check(x==-17 and len(cycle)==11 and sum(x%2 for x in cycle)==7,'Terras eleven cycle')
std=[];x=Fraction(-17)
while x not in std:std.append(x);x=standard(x)
check(x==-17 and len(std)==18,'standard eighteen cycle')
M=2**5*3**3*11**2
polys=[(1,0),(1,1),(1,-1),(3,5),(2,1)]
def val(n,p):
 n=abs(n)
 if not n:raise ValueError('valuation at zero')
 v=0
 while n%p==0:n//=p;v+=1
 return v
shadows=[]
for k in (1,2,3,11,22):
 n=M*2**(11*k)-17;y=M*3**(7*k)-17;x=n
 for j in range(11*k):
  check(x%2==cycle[j%11]%2,'integer shadow parity');x=terras(x)
 check(x==y,'integer shadow exact endpoint')
 for p in (2,3,11):
  for a,b in polys:check(val(a*n+b,p)==val(a*y+b,p),'endpoint polynomial valuations')
 check(val(n+17,2)-val(y+17,2)==11*k,'root-at-shadow hostile 2-adic')
 check(val(y+17,3)-val(n+17,3)==7*k,'root-at-shadow hostile 3-adic')
 shadows.append({'k':k,'steps':11*k,'n':str(n),'y':str(y)})
receipt={'status':'PASS','claim':'all-period golden fibre proof certificate; finite controls separately bounded','checks':CHECKS,'coarse_carries':coarse,'sharp_carries':fine,'graph_vertices':len(vertices),'graph_edges':sum(map(len,edges.values())),'cyclic_sccs':cyclic,'all_edges':[[list(x),[list(y) for y in edges[x]]] for x in vertices],'word_census':census,'rational_sources_checked':rational_checks,'period5_labels':base_labels,'lifts':lifts,'terras_minus17_cycle':cycle,'integer_shadows':shadows}
print(json.dumps(receipt,indent=2,sort_keys=True))
