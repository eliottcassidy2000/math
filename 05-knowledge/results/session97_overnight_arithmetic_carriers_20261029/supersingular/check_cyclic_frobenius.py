#!/usr/bin/env python3
"""Independent polynomial-factor and code certificates for (1,n) Frobenius."""
from itertools import product
from collections import Counter
from math import comb
from math import gcd
CHECKS=0
def check(x,label):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(label)
def deg(f):return f.bit_length()-1
def rem(f,g):
 while f and deg(f)>=deg(g):f^=g<<(deg(f)-deg(g))
 return f
def mul(f,g):
 z=0
 while g:
  if g&1:z^=f
  f<<=1;g>>=1
 return z
def divexact(f,g):
 q=0
 while f and deg(f)>=deg(g):
  i=deg(f)-deg(g);q^=1<<i;f^=g<<i
 check(f==0,'exact polynomial division')
 return q
def gcdpoly(f,g):
 while g:f,g=g,rem(f,g)
 return f
def nullspace(rows,n):
 rows=list(rows);piv=[];r=0
 for j in range(n):
  k=next((k for k in range(r,len(rows)) if rows[k]>>j&1),None)
  if k is None:continue
  rows[r],rows[k]=rows[k],rows[r]
  for k in range(len(rows)):
   if k!=r and rows[k]>>j&1:rows[k]^=rows[r]
  piv.append(j);r+=1
 out=[]
 for j in range(n):
  if j not in piv:
   v=1<<j
   for row,p in zip(rows,piv):
    if row>>j&1:v|=1<<p
   out.append(v)
 return out
def factor(f):
 n=deg(f)
 columns=[rem(1<<(2*j),f)^(1<<j) for j in range(n)]
 rows=[sum(((col>>i)&1)<<j for j,col in enumerate(columns)) for i in range(n)]
 basis=nullspace(rows,n);factors=[f]
 for v in basis:
  nxt=[]
  for h in factors:
   a=gcdpoly(h,v)
   if a!=1 and a!=h:nxt.extend((a,divexact(h,a)))
   else:nxt.append(h)
  factors=nxt
 check(len(factors)==len(basis),'Berlekamp complete factorization')
 z=1
 for h in factors:z=mul(z,h)
 check(z==f,'factor product')
 return sorted(factors)
def reciprocal(f):return int(bin(f)[2:][::-1],2)
def order2(n):
 a=2%n;k=1
 while a!=1:a=2*a%n;k+=1
 return k
def phi(n):return sum(gcd(i,n)==1 for i in range(1,n+1))
def primes(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   out.append(p)
   while n%p==0:n//=p
  p+=1
 if n>1:out.append(n)
 return out
def formula(n):
 if any(order2(p)%2==0 for p in primes(n)):return 0
 k=sum(phi(d)//order2(d) for d in range(3,n+1,2) if n%d==0)
 check(k%2==0,'reciprocal pairs integral')
 return 1<<(k//2)
def cyclotomic(n):
 unseen=set(range(1,n));cs=[]
 while unseen:
  a=min(unseen);c=set();b=a
  while b not in c:c.add(b);b=2*b%n
  unseen-=c;cs.append(frozenset(c))
 return cs
def code_basis(g,n):return [g<<i for i in range(n-deg(g))]
def rank(vs):
 piv={}
 for v in vs:
  while v:
   k=deg(v)
   if k in piv:v^=piv[k]
   else:piv[k]=v;break
 return len(piv)
rows=[]
weight_rows=[]
for n in range(3,128,2):
 fs=[f for f in factor((1<<n)|1) if f!=3]
 fixed=[f for f in fs if reciprocal(f)==f]
 check(all(reciprocal(f) in fs for f in fs),'reciprocal closure')
 a=0 if fixed else 1<<(len(fs)//2)
 b=formula(n);check(a==b,('support theorem',n))
 cs=cyclotomic(n)
 check(sorted(map(len,cs))==sorted(map(deg,fs)),('factor degrees',n))
 selfnegative=[c for c in cs if frozenset((-x)%n for x in c)==c]
 check(len(selfnegative)==len(fixed),('selfreciprocal cosets',n))
 if n in [3,5,7,11,15,21,23,31,49,71,73,89,127]:rows.append((n,tuple(map(deg,fs)),len(fixed),a))
 # Construct every positive kernel for selected n; no matching assumption.
 if n in [7,23,31,49,71,73,89,127]:
  pairs=[(f,reciprocal(f)) for f in fs if f<reciprocal(f)]
  distinct=set()
  for choice in product((0,1),repeat=len(pairs)):
   g=3 # x+1 enforces augmentation/even weight
   for pair,j in zip(pairs,choice):g=mul(g,pair[j])
   basis=code_basis(g,n)
   check(len(basis)==(n-1)//2,'Lagrangian dimension')
   check(rank(basis)==len(basis),'independence')
   for i,u in enumerate(basis):
    check(u.bit_count()%2==0,'augmentation')
    for v in basis[i:]:check((u&v).bit_count()%2==0,'isotropic')
   rotated=[((v<<1)&((1<<n)-1))|(v>>(n-1)) for v in basis]
   check(rank(basis+rotated)==len(basis),'Frobenius invariance')
   reflected=[sum(((v>>i)&1)<<((-i)%n) for i in range(n)) for v in basis]
   other=code_basis(reciprocal(g),n)
   check(rank(reflected+other)==len(basis),'reflection switches reciprocal choices')
   check(rank(reflected+basis)>len(basis),'no dihedral invariant kernel')
   distinct.add(g)
   if n in [7,23,31]:
    words=[0]
    for v in basis:words += [w^v for w in words]
    extended=words+[w^((1<<(n+1))-1) for w in words]
    weights=Counter(w.bit_count() for w in extended)
    weight_rows.append((n,g,dict(sorted(weights.items()))))
    if n==7:
     check(weights=={0:1,4:14,8:1},'extended Hamming enumerator')
     check(2**4*sum(comb(7,i) for i in range(2))==2**7,'perfect single error balls')
    if n==23:
     check(weights=={0:1,8:759,12:2576,16:759,24:1},'extended Golay enumerator')
     check(2**12*sum(comb(23,i) for i in range(4))==2**23,'perfect triple error balls')
  check(len(distinct)==a,('all constructed kernels',n))
print('CHECKS',CHECKS)
print('Complete odd n census: 3 <= n <= 127; independent Berlekamp and cyclotomic tests.')
for row in rows:print(row)
print('All positive kernels explicitly constructed for n=7,23,31,49,71,73,89,127.')

print('Extended-code weight enumerators (n,generator bitmask,weights):')
for row in weight_rows:print(row)
