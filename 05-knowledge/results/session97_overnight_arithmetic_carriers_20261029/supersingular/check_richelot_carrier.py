#!/usr/bin/env python3
"""Exact, stdlib-only certificates. Checks survive python -O."""
from itertools import combinations, permutations
from collections import Counter
from math import factorial
CHECKS=0
def check(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(msg)
def matchings(items):
 if not items: yield (); return
 a,*rest=items
 for b in rest:
  for tail in matchings([x for x in rest if x!=b]):
   yield tuple(sorted(((a,b),)+tail))
def permute_mask(x,p): return sum(1<<p[i] for i in range(len(p)) if x>>i&1)
def permute_matching(m,p): return tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in m))
def cycle_type(p):
 seen=set(); out=[]
 for i in range(len(p)):
  if i not in seen:
   j=i; n=0
   while j not in seen: seen.add(j); n+=1; j=p[j]
   out.append(n)
 return tuple(sorted(out))
def fixed_formula(ct):
 ans=1
 for ell,m in Counter(ct).items():
  if ell%2:
   if m%2:return 0
   j=m//2;ans*=factorial(m)//(2**j*factorial(j))*ell**j
  else:
   ans*=sum(factorial(m)//(factorial(m-2*j)*2**j*factorial(j))*ell**j for j in range(m//2+1))
 return ans
def canon(x,n):return min(x,x^((1<<n)-1))
def pairing(a,b):return (a&b).bit_count()%2
def lagrangians(n):
 # Enumerate isotropic subspaces without assumptions about matching origin.
 # EVEN subsets, with representatives whose last bit is zero.
 vs=[x for x in range(1<<(n-1)) if x.bit_count()%2==0]
 spaces={frozenset({0})}
 for d in range((n-2)//2):
  nxt=set()
  for u in spaces:
   for v in vs:
    if v not in u and all(pairing(v,w)==0 for w in u):
     nxt.add(u|frozenset(canon(v^w,n) for w in u))
  spaces=nxt
 return spaces
def matching_plane(m,n):
 u={0}
 for a,b in m:
  v=(1<<a)|(1<<b);u|={canon(w^v,n) for w in list(u)}
 return frozenset(u)
def q2_lines(u,n):
 return sum(min(x.bit_count(),n-x.bit_count())==2 for x in u if x)
def induced(u,p,n):return frozenset(canon(permute_mask(x,p),n) for x in u)
def adjacency(ms):return [{j for j,b in enumerate(ms) if i!=j and set(a)&set(b)} for i,a in enumerate(ms)]
M6=list(matchings(list(range(6))));L6=lagrangians(6)
check(len(M6)==len(L6)==15,'six count')
check({matching_plane(m,6) for m in M6}==L6,'six bijection')
rows={};group_images=set()
for p in permutations(range(6)):
 ct=cycle_type(p);fixed=sum(permute_matching(m,p)==m for m in M6)
 check(fixed==fixed_formula(ct),('fixed formula',p))
 rows[ct]=fixed
 group_images.add(tuple(M6.index(permute_matching(m,p)) for m in M6))
check(len(group_images)==720,'faithful Sp4(2)=S6 action')
G=adjacency(M6)
check(all(len(x)==6 for x in G),'valency')
for i,j in combinations(range(15),2):check(len(G[i]&G[j])==(1 if j in G[i] else 3),'SRG 15,6,1,3')
p=(1,0,3,2,5,4);F=[m for m in M6 if permute_matching(m,p)==m];A=adjacency(F)
check(len(F)==7 and sorted(map(len,A))==[2]*6+[6],'friendship degrees')
check(sum(map(len,A))//2==9,'friendship edges')
check(sum(all(j in A[i] for i,j in combinations(c,2)) for c in combinations(range(7),3))==3,'three triangles')
# Independent eight-point isotropic-subspace census.
M8=list(matchings(list(range(8))));L8=lagrangians(8)
check(len(M8)==105 and len(L8)==135,'eight counts')
matched={matching_plane(m,8) for m in M8};exotic=L8-matched
check(len(matched)==105 and len(exotic)==30,'eight split')
check(Counter(q2_lines(u,8) for u in L8)=={4:105,0:30},'weight-two orbit split')
# Every exotic gives extended Hamming [8,4,4], and a Fano plane at a named puncture.
fanos=set()
for u in exotic:
 code=u|frozenset(x^255 for x in u)
 check(Counter(x.bit_count() for x in code)=={0:1,4:14,8:1},'Hamming enumerator')
 triples=frozenset(x&127 for x in code if x.bit_count()==4 and x>>7&1)
 check(len(triples)==7,'Fano lines')
 for a,b in combinations(range(7),2):check(sum((x>>a&1) and (x>>b&1) for x in triples)==1,'Fano pair unique')
 fanos.add(triples)
check(len(fanos)==30,'puncture-bijection to 30 labelled Fano planes')
# Check general formula on every S8 class representative and all involutions via types.
rep8={}
for p in permutations(range(8)):rep8.setdefault(cycle_type(p),p)
rows8={}
for ct,p in sorted(rep8.items()):
 f=sum(permute_matching(m,p)==m for m in M8)
 check(f==fixed_formula(ct),('eight matching fixed',ct))
 e=sum(induced(u,p,8)==u for u in exotic)
 whole=sum(induced(u,p,8)==u for u in L8)
 check(whole==f+e,('eight sectors',ct));rows8[ct]=(f,e,whole)
# Exact product decomposition: 9 line-products + 6 symplectic graphs.
S3=list(permutations(range(3)))
def compose(p,q):return tuple(p[q[i]] for i in range(3))
reps=[(0,1,2),(1,0,2),(1,2,0)];table=[]
for u in reps:
 row=[]
 for v in reps:
  split=sum(u[i]==i for i in range(3))*sum(v[i]==i for i in range(3))
  mixed=sum(compose(v,a)==compose(a,u) for a in S3)
  row.append((split,mixed,split+mixed))
 table.append(row)
check(table==[[(9,6,15),(3,0,3),(0,0,0)],[(3,0,3),(1,2,3),(0,0,0)],[(0,0,0),(0,0,0),(0,3,3)]],'product Frob table')
def count_curve(p,a,b):return 1+sum(sum((y*y-x*x*x-a*x-b)%p==0 for y in range(p)) for x in range(p))
check(count_curve(11,-1,0)==count_curve(11,0,1)==12,'supersingular p11 traces')
check(sum((x*x*x-x)%11==0 for x in range(11))==3,'split E1 2torsion')
check(sum((x*x*x+1)%11==0 for x in range(11))==1,'nonsplit E2 2torsion')
# Pentads: outer six-set built from the same incidence, not arbitrary labels.
pentads=[]
for ids in combinations(range(15),5):
 edges=[e for i in ids for e in M6[i]]
 if len(set(edges))==15:pentads.append(frozenset(ids))
check(len(pentads)==6,'six one-factorizations')
new_duads=[frozenset(j for j,pd in enumerate(pentads) if i in pd) for i in range(15)]
check(len(set(new_duads))==15 and all(len(x)==2 for x in new_duads),'syntheme to new-duad bijection')
for i,j in combinations(range(15),2):check((j in G[i])==(not new_duads[i]&new_duads[j]),'outer graph duality')
tau=(1,0,2,3,4,5)
tau_m=[M6.index(permute_matching(m,tau)) for m in M6]
tau_pd=[pentads.index(frozenset(tau_m[i] for i in pd)) for pd in pentads]
check(cycle_type(tau_pd)==(2,2,2),'outer transposition becomes triple transposition')
# Frobenius (1,7) yields the two cyclic Fano planes and 21-element normalizers.
p7=(1,2,3,4,5,6,0,7)
fixed_exotic=[u for u in exotic if induced(u,p7,8)==u]
def puncture_fano(u):
 code=u|frozenset(x^255 for x in u)
 return frozenset(x&127 for x in code if x.bit_count()==4 and x>>7&1)
def develop(base):return frozenset(sum(1<<((a+i)%7) for a in base) for i in range(7))
check({puncture_fano(u) for u in fixed_exotic}=={develop((0,1,3)),develop((0,1,5))},'two cyclic Fano planes')
for u in fixed_exotic:
 f=puncture_fano(u)
 keep=[]
 for a in range(1,7):
  for b in range(7):
   perm=tuple((a*i+b)%7 for i in range(7))
   if frozenset(permute_mask(x,perm) for x in f)==f:keep.append((a,b))
 check(len(keep)==21 and {a for a,b in keep}=={1,2,4},'21 normalizer and 1-2-4 multipliers')
# Rabin irreducibility test for x^7+x+1 in F5[x]. Degree 7 is prime.
def rem(a,b,p):
 a=list(a)
 while len(a)>=len(b):
  t=a[-1]*pow(b[-1],-1,p)%p;d=len(a)-len(b)
  for j,x in enumerate(b):a[j+d]=(a[j+d]-t*x)%p
  while a and not a[-1]:a.pop()
 return a
def mul(a,b,f,p):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
 return rem(c,f,p)
def power(a,n,f,p):
 z=[1]
 while n:
  if n&1:z=mul(z,a,f,p)
  a=mul(a,a,f,p);n//=2
 return z
def gcd_poly(a,b,p):
 while b:a,b=b,rem(a,b,p)
 return a
f=[1,1,0,0,0,0,0,1]
check(power([0,1],5**7,f,5)==[0,1],'x^(5^7)=x')
check(len(gcd_poly(f,[0,4,0,0,0,1],5))==1,'gcd(f,x^5-x)=1')
check(fixed_formula((1,7))==0 and len(fixed_exotic)==2,'arithmetic no-matching/two-kernel hostile')
print('CHECKS',CHECKS)
print('S6 cycle type : rational matching kernels')
for ct,x in sorted(rows.items()):print(ct,x)
print('S8 cycle type : matching, exotic, total invariant Lagrangians')
for ct,x in sorted(rows8.items()):print(ct,x)
print('Product Frobenius table split,mixed,total:',table)
print('Fano puncture bijection: 30 exotic Lagrangians -> 30 labelled STS(7)')
