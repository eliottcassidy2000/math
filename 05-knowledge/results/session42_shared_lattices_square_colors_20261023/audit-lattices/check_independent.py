#!/usr/bin/env python3
import itertools as it,json,math
checks=0
def need(v,why):
 global checks
 checks+=1
 if not v:raise RuntimeError(why)
def subgroups(elements,identity,mul):
 def close(H,g):
  gens=list(H)+[g];seen={identity};todo=[identity]
  for a in todo:
   for b in gens:
    c=mul(a,b)
    if c not in seen:seen.add(c);todo.append(c)
  return frozenset(seen)
 first=frozenset([identity]);seen={first};todo=[first]
 for H in todo:
  for g in elements:
   if g not in H:
    K=close(H,g)
    if K not in seen:seen.add(K);todo.append(K)
 return todo
cyclic=0
for m in range(1,7):
 for n in range(1,7):
  el=list(it.product(range(m),range(n)))
  hs=subgroups(el,(0,0),lambda a,b:((a[0]+b[0])%m,(a[1]+b[1])%n))
  for H in hs:
   h=next(t for t in range(1,m+1) if (t%m,0) in H)
   v=next(t for t in range(1,n+1) if any((a,t%n) in H for a in range(m)))
   s=next(t for t in range(h) if (t,v%n) in H)
   need(m%h==n%v==0 and (s*(n//v))%h==0,'periods')
   predicted={(a,b) for a,b in el if b%v==0 and (v*a-s*b)%(h*v)==0}
   need(predicted==set(H),'subgroup hnf')
   need(len(H)*h*v==m*n,'covolume')
   need(sum(b==0 for a,b in H)==m//h,'row kernel')
   need(sum(a==0 for a,b in H)==n*math.gcd(h,s)//(h*v),'column kernel')
   for a,b in el:
    actual={((a+x)%m,(b+y)%n) for x,y in H}
    want={(x,y) for x,y in el if (y-b)%v==0 and (v*(x-a)-s*(y-b))%(h*v)==0}
    need(actual==want,'coset hnf')
   cyclic+=1
# All subgroups of S3 x S3, with full nonabelian kernel quotient checks.
S=list(it.permutations(range(3)));e=tuple(range(3))
def mul(p,q):return tuple(p[q[i]] for i in range(3))
def inv(p):return tuple(p.index(i) for i in range(3))
el=list(it.product(S,S));hs=subgroups(el,(e,e),lambda a,b:(mul(a[0],b[0]),mul(a[1],b[1])))
for H in hs:
 P={a for a,b in H};Q={b for a,b in H};K={a for a,b in H if b==e};L={b for a,b in H if a==e}
 need(len(H)==len(P)*len(L)==len(Q)*len(K),'projection cardinalities')
 need(all(mul(mul(p,k),inv(p)) in K for p in P for k in K),'normal row kernel')
 need(all(mul(mul(q,l),inv(q)) in L for q in Q for l in L),'normal column kernel')
 need(all((mul(inv(a),c) in K)==(mul(inv(b),d) in L) for a,b in H for c,d in H),'quotient graph')
# Independent R9 module detection, exact full group and isolated base size.
def arc(x,y):
 a,u=x%3,x//3;b,v=y%3,y//3
 return x!=y and ((b-a)%3==1 if a!=b else (v-u)%3==1)
for x,y in it.combinations(range(9),2):
 common=sum(arc(x,z) and arc(y,z) for z in range(9))
 need(common==(3 if x%3==y%3 else 1),'intrinsic R9 modules')
autos=[]
for outer in range(3):
 for rotations in it.product(range(3),repeat=3):
  p=tuple((x%3+outer)%3+3*((x//3+rotations[x%3])%3) for x in range(9))
  need(all(arc(x,y)==arc(p[x],p[y]) for x in range(9) for y in range(9)),'full R9 chart')
  autos.append(p)
def base(group):
 for k in range(10):
  for marks in it.combinations(range(9),k):
   if sum(all(p[x]==x for x in marks) for p in group)==1:return k
need(base(autos)==3,'isolated full R9 base')
need(base([tuple((x+a)%9 for x in range(9)) for a in range(9)])==1,'isolated cyclic base')
print(json.dumps({'status':'PASS','checks':checks,'cyclic_subgroups_m_n_through6':cyclic,'S3xS3_subgroups':len(hs),'R9_module_and_isolated_base':'PASS'},indent=2))
