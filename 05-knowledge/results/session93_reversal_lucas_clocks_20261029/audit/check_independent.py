#!/usr/bin/env python3
"""Independent Session93 audit. No imports from producer code; active -O checks."""
from itertools import product
from math import gcd
import json
checks=0
def check(b):
 global checks
 checks+=1
 if not b: raise RuntimeError('audit check '+str(checks))
def mul(a,b):return tuple(a[i] for i in b)
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def power(a,n):
 r=tuple(range(len(a)))
 for _ in range(n):r=mul(a,r)
 return r
def conj(a,b):return mul(mul(a,b),inv(a))
def root(a):
 # Independent repeated powering, not cyclewise producer implementation.
 r=a
 while mul(r,r)!=a:r=mul(r,a)
 return r
def neg(k):
 return tuple(sum((-(x//3**i)%3)*3**i for i in range(k)) for x in range(3**k))
def rot(k,depth):
 # Rotate digit at depth from left only when prefix is all zero.
 place=3**(k-depth-1)
 return tuple(x if x//(3*place) else x+(1 if (x//place)%3<2 else -2)*place for x in range(3**k))
def embed(a,block):
 n=len(a);r=list(range(3*n));r[block*n:(block+1)*n]=[block*n+x for x in a];return tuple(r)
def triangle(s,t,u):return mul(inv(root(mul(u,s))),mul(root(mul(u,t)),root(mul(t,s))))
def group(gens,n):
 seen={tuple(range(n))};todo=list(seen)
 for a in todo:
  for g in gens:
   b=mul(g,a)
   if b not in seen:seen.add(b);todo.append(b)
 return seen
# Independent direct-label reconstruction and isolated-cycle certificates.
def fastpower(a,e):
 r=tuple(range(len(a)))
 while e:
  if e&1:r=mul(r,a)
  a=mul(a,a);e//=2
 return r
for n,p,alpha,exponents,exponent,length in [
 (7,29,24,[7,1,9,26,25,16,27,7,18,5],14,13),
 (11,199,62,[161],644028,23),
 (13,521,422,[168,442,105],1774410,7)]:
 labels={}
 for mask in range(1<<n):
  if any((mask>>i&1) and (mask>>((i+1)%n)&1) for i in range(n)):continue
  label=sum(pow(alpha,i,p) for i in range(n) if mask>>i&1)%p
  check(label not in labels);labels[label]=mask
 check(len(labels)==p)
 reflection=[]
 for x in range(p):
  mask=labels[x];reflection.append(sum(pow(alpha,(-i)%n,p) for i in range(n) if mask>>i&1)%p)
 reflection=tuple(reflection);g=tuple(range(p))
 for a in exponents:g=mul(reflection,mul(tuple((i+a)%p for i in range(p)),g))
 h=fastpower(g,exponent);moved={i for i in range(p) if h[i]!=i}
 check(len(moved)==length)
 x=next(iter(moved));orbit={x};y=h[x]
 while y!=x:orbit.add(y);y=h[y]
 check(orbit==moved);check(p-length>=3)

reports=[]
for k in range(2,5):
 n=3**k; alpha=neg(k);x=rot(k,0)
 for j in range(k-1):
  v=rot(k-1,j);y=embed(v,0)
  check(conj(alpha,y)==inv(y));check(power(x,3)==tuple(range(n)))
  s=alpha;t=mul(power(x,2),alpha);u=mul(power(y,2),alpha)
  for a in [s,t,u]:check(power(a,2)==tuple(range(n)))
  h=triangle(s,t,u)
  expected=mul(embed(v,1),embed(inv(v),2))
  check(h==expected);check(mul(h,alpha)==mul(alpha,h))
 # Exact group generation at k2/3 only; k4 enormous.
 if k<=3:
  gens=[]
  for d in range(1,k):
   for j in range(d):
    v=rot(d,j);h=mul(embed(v,1),embed(inv(v),2))
    while len(h)<n:h=embed(h,0)
    gens.append(h)
  C=group(gens,n);check(len(C)==3**((n-2*k-1)//4))
  reports.append({'k':k,'triangle_generators':len(gens),'generated_order':len(C)})
# Arithmetic observer, complete modulus range, including non-squarefree and ramified.
for m in range(2,41):
 fibres={}
 for u,v in product(range(m),repeat=2):
  key=((2*u+v)%m,(3*u+4*v)%m)
  fibres.setdefault(key,[]).append((u,v))
  x,y=u,v;hist=[]
  for k in range(6):
   hist.append((2*x+y)%m);x,y=(x+y)%m,(x+2*y)%m
  check(all((hist[i+2]-3*hist[i+1]+hist[i])%m==0 for i in range(4)))
 check(set(map(len,fibres.values()))=={gcd(m,5)})
 check(len(fibres)*gcd(m,5)==m*m)
# Sparse-clock joint information: independent direct iteration, all states.
for m in range(2,21):
 for strides in [(2,),(3,),(2,3),(2,4),(3,5)]:
  fibres={};fs=[0,1]
  for _ in range(2*max(strides)):fs.append(fs[-1]+fs[-2])
  expected=gcd(m,*[5*fs[2*d] for d in strides])
  for u,v in product(range(m),repeat=2):
   x,y=u,v;observed=[(2*u+v)%m]
   for d in range(1,max(strides)+1):
    x,y=(x+y)%m,(x+2*y)%m
    if d in strides:observed.append((2*x+y)%m)
   key=tuple(observed);fibres[key]=fibres.get(key,0)+1
  check(set(fibres.values())=={expected})
  if m==6 and strides==(2,3):check(len(fibres)==36)

# Named paths / blue quotient and Paley hostile.
def paths(rows):
 n=len(rows)
 def walk(p,left):
  if not left:yield tuple(p);return
  for v in left:
   if not p or v in rows[p[-1]]:yield from walk(p+[v],left-{v})
 return list(walk([],set(range(n))))
for name,n in [('Paley7',7),('W2',9)]:
 if n==7:
  rows=[{(x+d)%7 for d in (1,2,4)} for x in range(7)]
  A={tuple((a*x+b)%7 for x in range(7)) for a in (1,2,4) for b in range(7)};tau=tuple(-x%7 for x in range(7))
 else:
  rows=[{y for y in range(9) if y!=x and (((y//3-x//3)%3==1) if y//3!=x//3 else (y-x)%3==1)} for x in range(9)]
  A={tuple(3*((x//3+r)%3)+(x%3+b[x//3])%3 for x in range(9)) for r in range(3) for b in product(range(3),repeat=3)};tau=neg(2)
 C={a for a in A if mul(a,tau)==mul(tau,a)};I={conj(a,tau) for a in A}
 H=group({triangle(tau,t,u) for t,u in product(I,repeat=2)},n)
 X=paths(rows);F=[p for p in X if tuple(tau[v] for v in p[::-1])==p]
 check(all(tuple(h[v] for v in p)!=p for h in C if h!=tuple(range(n)) for p in F))
 check((len(X),len(F),len(C),len(H))==((189,9,3,1) if n==7 else (3159,33,3,3)))
 reports.append({'object':name,'H':len(X),'fixed_paths':len(F),'centralizer':len(C),'triangle_subgroup':len(H),'blue_presentations':len(F)//len(C)})
# Final independent port-orbit / base audit, literal C3 action.
from itertools import combinations
def digits(x,k):return tuple(x//3**j%3 for j in range(k-1,-1,-1))
def port_signature(xs,k):
 tags=[];folds=[]
 for x in xs:
  ds=digits(x,k);first=next((i for i,v in enumerate(ds) if v),None)
  if first is None:tags.append((-1,0));folds.append(())
  else:
   side=ds[first];tags.append((k-first-1,side))
   folds.append(tuple(v if side==1 else -v%3 for v in ds[first+1:]))
 pair=[]
 for a,b in combinations(range(len(xs)),2):
  if tags[a][0]!=tags[b][0] or tags[a][0]<0:continue
  u,v=folds[a],folds[b];j=next((i for i in range(len(u)) if u[i]!=v[i]),len(u))
  pair.append((a,b,j,0 if j==len(u) else (v[j]-u[j])%3))
 return tuple(tags),tuple(pair)
k=3;n=27;gens=[]
for d in range(1,k):
 for j in range(d):
  v=rot(d,j);h=mul(embed(v,1),embed(inv(v),2))
  while len(h)<n:h=embed(h,0)
  gens.append(h)
C=group(gens,n)
for arity in (1,2,3):
 remaining=set(product(range(n),repeat=arity));sig_to_orbit={};orbits=0
 while remaining:
  x=min(remaining);orbit={tuple(g[v] for v in x) for g in C};remaining-=orbit
  sig=port_signature(x,k);check(sig not in sig_to_orbit)
  sig_to_orbit[sig]=orbits;orbits+=1
  check(all(port_signature(y,k)==sig for y in orbit))
 if arity==1:check(orbits==2*k+1)
 if arity==2:check(orbits==8*k*k+1)
check(port_signature((3,6),2)!=port_signature((3,8),2))
# Exhaust all smaller candidate bases at k3, and verify a size4 base.
nonid=list(C-{tuple(range(n))});fix=[sum(1<<j for j,g in enumerate(nonid) if g[v]==v) for v in range(n)]
for size in range(4):
 for B in combinations(range(n),size):
  stabilizer=(1<<len(nonid))-1
  for v in B:stabilizer &=fix[v]
  check(stabilizer!=0)
stabilizer=(1<<len(nonid))-1
for v in (3,9,12,15):stabilizer &=fix[v]
check(stabilizer==0)
# Trace-conic substitution for all states, composite and ramified included.
for m in range(1,25):
 for d in range(1,8):
  r,s=1,0
  for _ in range(d):r,s=r+s,r+2*s
  L=2*r+s
  for u,v in product(range(m),repeat=2):
   t=(2*u+v)%m;z=((2*r+s)*u+(r+3*s)*v)%m
   check((t*t+z*z-L*t*z+5*s*s*(u*u+u*v-v*v))%m==0)

print(json.dumps({'status':'PASS','checks':checks,'reports':reports},sort_keys=True))
