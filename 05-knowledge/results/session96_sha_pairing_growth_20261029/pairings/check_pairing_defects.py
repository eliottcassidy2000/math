#!/usr/bin/env python3
"""Finite pairing controls, not arithmetic Sha computations. Standard library only."""
from itertools import product
from math import gcd,lcm,isqrt
from collections import Counter
import json
checks=0
forms=0
perfect=0

def ck(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)

def partitions(n,lo=1):
 if not n:yield ()
 for i in range(lo,n+1):
  for tail in partitions(n-i,i):yield (i,)+tail

def group(ns):return list(product(*(range(n) for n in ns)))
def add(x,y,ns):return tuple((a+b)%n for a,b,n in zip(x,y,ns))
def scale(k,x,ns):return tuple(k*a%n for a,n in zip(x,ns))
def order(x,ns):return lcm(*(n//gcd(n,a) for a,n in zip(x,ns))) if ns else 1

def matrices(ns):
 L=lcm(*ns) if ns else 1;s=len(ns)
 coords=[(i,j) for i in range(s) for j in range(i,s)]
 vals=[]
 for i,j in coords:
  if i==j:vals.append([0,L//2] if ns[i]%2==0 else [0])
  else:vals.append([k*(L//gcd(ns[i],ns[j])) for k in range(gcd(ns[i],ns[j]))])
 for choices in product(*vals):
  B=[[0]*s for _ in range(s)]
  for (i,j),a in zip(coords,choices):B[i][j]=a;B[j][i]=(-a)%L if j!=i else a
  yield L,B

def inspect(ns,L,B):
 global forms,perfect
 forms+=1;G=group(ns);s=len(ns);zero=(0,)*s
 def b(x,y):return sum(x[i]*B[i][j]*y[j] for i in range(s) for j in range(s))%L
 signatures={tuple(sum(x[i]*B[i][j] for i in range(s))%L for j in range(s)) for x in G}
 if len(signatures)!=len(G):return False
 perfect+=1
 candidates=[y for y in G if all(sum(B[i][j]*y[j] for j in range(s))%L==B[i][i] for i in range(s))]
 ck(len(candidates)==1,('unique c',ns,B));c=candidates[0]
 ck(scale(2,c,ns)==zero,('two-torsion c',ns,B))
 delta={x:0 if b(x,x)==0 else 1 for x in G}
 ck(all(b(x,x) in {0,L//2} and b(x,c)==b(x,x) for x in G),('diagonal char',ns,B))
 eps=delta[c];f={x:add(x,scale(delta[x],c,ns),ns) for x in G}
 ck(all(f[f[x]]==(x if eps==0 else f[x]) for x in G),('repair square',ns,B))
 ck(all(f[add(x,y,ns)]==add(f[x],f[y],ns) for x in G for y in G),('repair homomorphism',ns,B))
 def beta(x,y):return b(x,f[y])
 ck(all(beta(x,x)==0 for x in G),('alternating repair',ns,B))
 radical=[x for x in G if all(beta(x,y)==0 for y in G)]
 ck(set(radical)==({zero} if eps==0 else {zero,c}),('repair radical',ns,B))
 domain=G if eps==0 else [x for x in G if b(x,c)==0]
 if eps:
  ck(len(domain)*2==len(G),('odd split',ns,B))
  ck(len({add(x,y,ns) for x in domain for y in [zero,c]})==len(G),('split bijection',ns,B))
 # Peeling uses beta when eps=0, original b on c-perp when eps=1.
 pairing=beta if eps==0 else b
 residual=domain;ms=[]
 while len(residual)>1:
  x=max(residual,key=lambda x:order(x,ns));m=order(x,ns)
  ys=[y for y in residual if pairing(x,y)==L//m]
  ck(bool(ys),('hyperbolic partner',ns,B,x,m));y=ys[0]
  H={add(scale(i,x,ns),scale(j,y,ns),ns) for i in range(m) for j in range(m)}
  perp=[z for z in residual if pairing(z,x)==pairing(z,y)==0]
  ck(len(H)==m*m and len(H)*len(perp)==len(residual),('peel size',ns,B))
  ck({add(v,w,ns) for v in H for w in perp}==set(residual),('peel reconstruct',ns,B))
  ms.append(m);residual=perp
 R=1
 for m in ms:R*=m
 ck(len(G)==2**eps*R*R,('square law',ns,B))
 ck(isqrt(len(G)//2**eps)**2==len(G)//2**eps,('integer square',ns,B))
 return True

families=[()]
for e in range(1,5):families.extend(tuple(2**j for j in part) for part in partitions(e))
families +=[(3,),(3,3),(5,5),(2,3),(2,2,3,3),(4,4,2),(8,8)]
summary=[]
for ns in families:
 before=(forms,perfect)
 for L,B in matrices(ns):inspect(ns,L,B)
 summary.append({'cyclic_orders':ns,'forms':forms-before[0],'perfect':perfect-before[1]})
# A square cardinality is not enough, and even defect does not imply alternating.
ck(not any(inspect((4,),L,B) for L,B in matrices((4,))),'C4 hostile')
ck(inspect((2,2),2,[[1,0],[0,1]]),'isotropic nonzero characteristic hostile')

def carry(a,b):
 e,R=a;f,S=b
 return e^f,2**(e*f)*R*S
for e,f,g in product(range(2),repeat=3):
 for R,S,T in product(range(1,5),repeat=3):
  a=(e,R);b=(f,S);c=(g,T)
  ck(carry(carry(a,b),c)==carry(a,carry(b,c)),('carry associative',a,b,c))
  out=carry(a,b);ck(2**out[0]*out[1]**2==(2**e*R*R)*(2**f*S*S),('carry order',a,b))
# Group-isotype inventory, not isometry counts.
D=20;A=[1]+[0]*D
for j in range(2,D+1,2):
 for n in range(j,D+1):A[n]+=A[n-j]
for n in range(D+1):
 literal=sum(all(v%2==0 for v in Counter(p).values()) for p in partitions(n))
 ck(A[n]==literal,('even multiplicity inventory',n))
 B=A[n]+(A[n-1] if n else 0)
 literal2=sum(all(v%2==0 for k,v in Counter(p).items() if k>1) for p in partitions(n))
 ck(B==literal2,('defect inventory',n))
print(json.dumps({'status':'PASS','checks':checks,'pairing_forms_examined':forms,'perfect_forms':perfect,'families':summary,'scope':'finite abelian pairings; no arithmetic Sha computation'},sort_keys=True))
