#!/usr/bin/env python3
"""Exact, bounded additive and quadratic-algebra controls; no assertions."""
from math import isqrt,gcd
from itertools import product
import json
from pathlib import Path

def ck(x,msg):
 if not x: raise RuntimeError(msg)
def b(n):
 m=n+1;return (3*m-isqrt(5*m*m)-1)//2
def Q(n):return (n-2*b(n),b(n))
def mul(x,y,s=1,u=1):
 a,b=x;c,d=y;return (a*c+u*b*d,a*d+b*c+s*b*d)
def norm(x,s=1,u=1):
 a,b=x;return a*a+s*a*b-u*b*b
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def support(n):
 fs=[1,2]
 while fs[-1]<=n:fs.append(sum(fs[-2:]))
 out=[]
 for i in range(len(fs)-1,-1,-1):
  if fs[i]<=n:out.append(i);n-=fs[i]
 return frozenset(out)
def color(n):
 z=min(support(n));return 'C' if z==0 else 'A' if z%2 else 'B'
def square(n):
 return any(isqrt(n-a*a)**2==n-a*a for a in range(isqrt(n)+1))
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
N=300
carry={}; overlap_zero=None; product_counter=None
for a in range(N+1):
 for c in range(N+1):
  delta=b(a)+b(c)-b(a+c)
  ck(delta in (-1,0,1),'carry bound')
  ck(sub(add(Q(a),Q(c)),Q(a+c))==(-2*delta,delta),'carry identity')
  overlap=bool(support(a)&support(c))
  ck(overlap or delta==0,'disjoint implication')
  if overlap and delta==0 and overlap_zero is None:overlap_zero=[a,c]
  carry.setdefault(str(delta),[a,c])
  x=mul(Q(a),Q(c));v=x[0]+2*x[1]
  ck(v==a*c-b(a)*b(c),'quadratic evaluation defect')
  ck(norm(x)==norm(Q(a))*norm(Q(c)),'norm multiplicativity')
  if x!=Q(a*c) and a>0 and c>0 and product_counter is None:product_counter=[a,c,list(x),list(Q(a*c))]
# Multiple summands and association independent carries, all triples <=30.
for a,c,d in product(range(31),repeat=3):
 D=b(a)+b(c)+b(d)-b(a+c+d)
 # ceil(2alpha)=1, hence -2 <=D <=1.
 ck(-2<=D<=1,'ternary carry bound')
 ck(D==(b(a)+b(c)-b(a+c))+(b(a+c)+b(d)-b(a+c+d)),'cocycle')
# Universal quadratic norms including Gaussian, golden, split, dual numbers.
for s,u in product(range(-2,3),repeat=2):
 for x,y in product(product(range(-2,3),repeat=2),repeat=2):
  ck(norm(mul(x,y,s,u),s,u)==norm(x,s,u)*norm(y,s,u),'universal norm')
# Mod2 nilpotents and units.
F=list(product(range(2),repeat=2))
def modmul(x,y,s,u):return tuple(z%2 for z in mul(x,y,s,u))
golden_nil=[x for x in F if x!=(0,0) and modmul(x,x,1,1)==(0,0)]
gaussian_nil=[x for x in F if x!=(0,0) and modmul(x,x,0,-1)==(0,0)]
ck(golden_nil==[] and gaussian_nil==[(1,1)],'mod2 obstruction')
# Positive other end: explicit golden-to-Gaussian algebra/norm map modulo3.
def map3(x):return ((x[0]+2*x[1])%3,x[1]%3)
for x,y in product(product(range(3),repeat=2),repeat=2):
 ck(map3(mul(x,y))==tuple(v%3 for v in mul(map3(x),map3(y),0,-1)), 'mod3 multiplication isomorphism')
 ck(norm(x)%3==norm(map3(x),0,-1)%3,'mod3 norm isometry')
ck(len({map3(x) for x in product(range(3),repeat=2)})==9,'mod3 bijection')
# Fibonacci/golden unit/Gaussian norm bridge.
f=[0,1]
for k in range(1,121):f.append(f[-1]+f[-2])
for k in range(1,60):
 ck(f[2*k+1]==f[k]**2+f[k+1]**2,'odd Fibonacci Gaussian norm')
 ck(gcd(f[k],f[k+1])==1,'primitive consecutive Fibonacci')
 raw=(f[k+1]**2-f[k]**2,2*f[k]*f[k+1],f[2*k+1])
 content=gcd(gcd(raw[0],raw[1]),raw[2])
 ck(content==(2 if k%3==1 else 1),'Euclid parity content')
 v=(1,0)
 for j in range(k):v=mul(v,(0,1))
 ck(norm(v)==(-1)**k,'golden norm unit')
# Same owner color and n mod4, opposite arithmetic membership.
hostiles={}
for predicate,name in [(prime,'prime'),(square,'two_squares')]:
 for c in 'CAB':
  found=None
  for a in range(2,501):
   if color(a)!=c:continue
   for d in range(a+1,501):
    if color(d)==c and a%4==d%4 and predicate(a)!=predicate(d):found=[a,d,predicate(a),predicate(d)];break
   if found:break
  ck(found is not None,'missing arithmetic hostile');hostiles[name+'_'+c]=found
# Finite, complete ordered Goldbach counts, partitioned by owner and carry.
limit=2000
ps=[p for p in range(2,limit+1) if prime(p)];P=set(ps)
goldbach={};all_nonempty=True
for n in range(4,limit+1,2):
 h={};total=0
 for p in ps:
  q=n-p
  if q in P:
   key=color(p)+color(q)+':'+str(b(p)+b(q)-b(n));h[key]=h.get(key,0)+1;total+=1
 ck(sum(h.values())==total,'partition count')
 all_nonempty &= total>0
 if n in (10,34,36,70,100,1000,2000):goldbach[str(n)]={'total':total,'bins':h}
# Exact multivariate transfer with a shared support state, nine digit positions.
def transfer(L,disjoint):
 states=[(0,0),(1,0),(0,1)] if disjoint else list(product(range(2),repeat=2))
 fs=[1,2]
 for i in range(2,L):fs.append(sum(fs[-2:]))
 D={((0,0),0,0):1}
 for weight in fs[:L]:
  E={}
  for (previous,a,c),v in D.items():
   for current in states:
    if any(x*y for x,y in zip(previous,current)):continue
    key=(current,a+weight*current[0],c+weight*current[1])
    E[key]=E.get(key,0)+v
  D=E
 out={}
 for (state,a,c),v in D.items():out[(a,c)]=out.get((a,c),0)+v
 return out
L=9;full=transfer(L,False);disjoint=transfer(L,True)
for (a,c),coefficient in full.items():
 ck(coefficient==1,'canonical pair uniqueness')
 ck(disjoint.get((a,c),0)==int(not (support(a)&support(c))),'shared-state support matrix')
ck(all(pair in full for pair in disjoint),'restricted support')
# Existing stripe branch is data, not extended using canonical candidate.
old='RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB'.translate(str.maketrans({'R':'C','K':'A'}))
new=''.join(color(i) for i in range(1,36))
ck([i+1 for i,(x,y) in enumerate(zip(old,new)) if x!=y]==[35],'branch discrepancy')
# On additive diagonal36 the change is genuine: unordered pair1+35 shifts CC to CB.
branch_pair={'sum':36,'pair':[1,35],'canonical':color(1)+color(35),'historical':old[0]+old[34]}
receipt={'status':'PASS','mod3_algebra_pairs':81,'transfer_positions':L,'full_transfer_monomials':len(full),'disjoint_transfer_monomials':len(disjoint),'pair_bound':N,'pair_tests':(N+1)**2,'triple_bound':30,'universal_quadratic_norm_tests':25*625,'carry_examples':carry,'overlap_zero_carry':overlap_zero,'product_hostile':product_counter,'golden_nonzero_square_zero':golden_nil,'gaussian_nonzero_square_zero':gaussian_nil,'arithmetic_hostiles_same_color_and_mod4':hostiles,'goldbach_even_bound':limit,'goldbach_no_missing_finite_only':all_nonempty,'selected_goldbach_partitions':goldbach,'historical_branch_pair':branch_pair}
print(json.dumps(receipt,indent=2,sort_keys=True))
