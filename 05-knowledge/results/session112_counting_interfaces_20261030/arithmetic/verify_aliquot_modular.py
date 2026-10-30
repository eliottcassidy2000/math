#!/usr/bin/env python3
import math,json

def need(b,msg):
 if not b:raise RuntimeError(msg)
def S(p,j):return (p**j-1)//(p-1)
def val(x,p,k):
 x%=p**k
 if x==0:return k
 z=0
 while x%p==0:z+=1;x//=p
 return z
def decode(p,t,k,A,B):
 m=p**k; h=min(val(B-A,p,k),k-1)
 return (S(p,h+1)*pow(S(p,t),-1,m)*(B-p**t*A)-A)%m
def encode(p,t,k,A,r):
 m=p**k;h=min(val(r,p,k),k-1)
 return (A+S(p,t)*pow(S(p,h+1),-1,m)*(r+p**(h+1)*A))%m
def sigma(n):
 z=0
 for d in range(1,math.isqrt(n)+1):
  if n%d==0:z+=d+(n//d if n//d!=d else 0)
 return z
arith=0;perms=0;metric=0;states=0;compat=0
for p in [2,3,5,7]:
 for t in [1,2,3]:
  for n in range(1,501):
   A=sigma(n)-n;B=sigma(p**t*n)-p**t*n
   for k in [1,2,3,4]:
    need(decode(p,t,k,A,B)==n%p**k,'actual divisor inverse');arith+=1
for p,k in [(2,1),(2,2),(2,3),(2,5),(2,7),(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2)]:
 m=p**k
 for t in [1,2,3]:
  for A in range(m):
   ds=[decode(p,t,k,A,B) for B in range(m)]
   need(len(set(ds))==m,'permutation');perms+=1
   for B,r in enumerate(ds):
    need(encode(p,t,k,A,r)==B,'inverse')
    if k>1:need(r%p**(k-1)==decode(p,t,k-1,A,B),'compatibility');compat+=1
    rr,aa=r,A
    for _ in range((k+t-1)//t):rr,aa=(p**t*rr)%m,encode(p,t,k,aa,rr)
    need(rr==0 and encode(p,t,k,aa,rr)==aa,'finite stabilization');states+=1
   if m<=32:
    for i in range(m):
     for j in range(m):need(val(ds[i]-ds[j],p,k)==val(i-j,p,k),'vertical isometry');metric+=1
latin=0;joint=0;conjugacy=0
for p,k in [(2,1),(2,2),(2,3),(2,5),(3,1),(3,2),(3,3),(5,1),(5,2)]:
 m=p**k
 for t in [1,2,3]:
  basins={c:0 for c in range(m)}
  depths={c:[0]*((k+t-1)//t+1) for c in range(m)}
  for B in range(m):
   column=[decode(p,t,k,A,B) for A in range(m)]
   need(len(set(column))==m,'Latin columns');latin+=1
   for A,r in enumerate(column):
    h=min(val(r,p,k),k-1)
    need((S(p,h+1)*B-S(p,t)*r)*pow(S(p,h+t+1),-1,m)%m==A,'first inverse')
  for r in range(m):
   h=min(val(r,p,k),k-1)
   for A in range(m):
    C=(A+r)*pow(1-p**(h+1),-1,m)%m
    rp=p**t*r%m;Bp=encode(p,t,k,A,r);hp=min(val(rp,p,k),k-1)
    Cp=(Bp+rp)*pow(1-p**(hp+1),-1,m)%m
    need(Cp==C,'normal form conjugacy');conjugacy+=1
    basins[C]+=1
    depth=0 if r==0 else (k-val(r,p,k)+t-1)//t
    depths[C][depth]+=1
  need(set(basins.values())=={m},'uniform basins')
  expect=[1]+[p**min(j*t,k)-p**min((j-1)*t,k) for j in range(1,(k+t-1)//t+1)]
  need(all(x==expect for x in depths.values()),'exact transient layers')
  if m<=9:
   for A in range(m):
    for B in range(m):
     for AA in range(m):
      for BB in range(m):
       old=min(val(A-AA,p,k),val(B-BB,p,k))
       new=min(val(A-AA,p,k),val(decode(p,t,k,A,B)-decode(p,t,k,AA,BB),p,k))
       need(old==new,'joint sup isometry');joint+=1
need(decode(2,1,2,0,1)==decode(2,1,2,1,0)==1,'nonperiodic observer hostile')
shifted=0;semigroup=0
for p in [2,3,5,7]:
 for k in range(2,7):
  m=p**k
  for j in range(1,k):
   n1=p**(k-j);n2=p**(k-j+1)
   need(val(n2-n1,p,k)==k-j,'sharp source difference')
   for ell in range(j,j+5):
    # Exact prime-power divisor sum, independent of decoder.
    a1=S(p,k-j+ell)-0
    a2=S(p,k-j+1+ell)-0
    need(a1%m==a2%m,'all later shifted collision');shifted+=1
for p,k in [(2,4),(3,3),(5,2)]:
 m=p**k
 for r in range(m):
  for A in range(m):
   for t,u in [(1,1),(1,2),(2,3)]:
    r1=p**t*r%m;A1=encode(p,t,k,A,r)
    r2=p**u*r1%m;A2=encode(p,u,k,A1,r1)
    need((r2,A2)==(p**(t+u)*r%m,encode(p,t+u,k,A,r)),'semigroup');semigroup+=1
print(json.dumps({'status':'PASS','direct_divisor_checks':arith,'vertical_permutations':perms,'state_inverse_and_stabilization_checks':states,'precision_compatibility_checks':compat,'metric_pair_checks':metric,'Latin_columns':latin,'joint_metric_checks':joint,'conjugacy_and_basin_checks':conjugacy,'shifted_precision_hostiles':shifted,'semigroup_checks':semigroup,'scope':'Universal claims proved in ALIQUOT_MODULAR.md; these finite checks use direct divisor enumeration and exhaustive residue planes.'},indent=2))
