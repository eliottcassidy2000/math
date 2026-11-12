#!/usr/bin/env python3
"""Independent exact controls. Uses no lane implementation or cached receipts."""
from collections import Counter
from fractions import Fraction
from itertools import product, permutations
from math import factorial
import json

def require(c,m):
 if not c: raise RuntimeError(m)

def canonical(xs):
 if len(xs)==1:return xs[0],0
 d=len(xs)//3
 p=[canonical(xs[i*d:(i+1)*d]) for i in range(3)]
 w=tuple(x[0] for x in p)
 return min(w,w[1:]+w[:1],w[2:]+w[:2]),sum(x[1] for x in p)+int(w[0]==w[1]==w[2])

def conv(a,b):
 c=Counter()
 for i,x in a.items():
  for j,y in b.items():c[i+j]+=x*y
 return c

def labeled_step(a,I):
 c=conv(conv(a,a),a)
 for s,x in a.items():
  v=x*3**(2*I-2*s)
  c[3*s]-=v;c[3*s+1]+=v
 return +c

def H(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   z=dp[mask][v]
   if not z:continue
   avail=T[v]&~mask
   while avail:
    b=avail&-avail;avail-=b
    dp[mask|b][b.bit_length()-1]+=z
 return sum(dp[-1])

def inflate(T,v,m):
 old=[x for x in range(len(T)) if x!=v]
 labels=old+[v]*m;n=len(labels);A=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   forward=(i<j) if labels[i]==labels[j] else (T[labels[i]]>>labels[j])&1
   if forward:A[i]|=1<<j
   else:A[j]|=1<<i
 return A

def skeleton(T,v):
 old=[x for x in range(len(T)) if x!=v];out=Counter()
 for r in range(1,len(T)+1):
  for w in set(permutations(old+[v]*r)):
   if any(a==b for a,b in zip(w,w[1:])):continue
   if all((T[a]>>b)&1 for a,b in zip(w,w[1:])):out[r]+=1
 return out

def stirling(m,r):
 return sum((-1)**(r-j)*__import__('math').comb(r,j)*j**m for j in range(r+1))

def main():
 checks={}
 for q in (2,3):
  a=Counter({0:q});I=0
  for e in (1,2):
   a=labeled_step(a,I);I=3*I+1
   brute=Counter(canonical(xs)[1] for xs in product(range(q),repeat=3**e))
   require(a==brute,f'color recurrence q{q} e{e}')
   checks[f'color_q{q}_e{e}']=dict(sorted(a.items()))
 # Every labelled core of order three, and several order four cores.
 cores=[]
 for n in (3,4):
  edges=[(i,j) for i in range(n) for j in range(i+1,n)]
  for bits in range(1<<len(edges)):
   if n==4 and bits not in (0,5,13,21,29,42,63):continue
   T=[0]*n
   for k,(i,j) in enumerate(edges):
    if bits>>k&1:T[i]|=1<<j
    else:T[j]|=1<<i
   cores.append(T)
 count=0
 for T in cores:
  for v in range(len(T)):
   a=skeleton(T,v)
   require(max(a)<=((len(T)-1)//2+1),'run rank bound')
   for m in range(1,5):
    pred=sum(c*stirling(m,r) for r,c in a.items())
    require(H(inflate(T,v,m))==pred,'physical inflation identity')
    count+=1
 checks['inflation_literal_cases']=count
 # C3 inflation: exact H=2^m+1, Hensel cells for p=3,5,7,11.
 for p in (3,5,7,11):
  q=((2**(p-1)-1)//p)%p
  require(q!=0,'ordinary Fermat quotient control')
  for a in (1,2,3):
   modulus=p**a
   values={pow(2,1+(p-1)*t,modulus)+1 for t in range(p**(a-1))}
   values={x%modulus for x in values}
   wanted={x for x in range(modulus) if x%p==3%p}
   require(values==wanted,'exact phase Hensel image')
 checks['hensel_control_primes']=[3,5,7,11]
 # A vanished first tangent does not eliminate local interior.
 p=1093;z=2**(p-1)-1;v=0
 while z%p==0:z//=p;v+=1
 require(v==2,'Wieferich higher tangent control')
 values={ (pow(2,1+(p-1)*t,p**3)+1)%(p**3) for t in range(p) }
 require(values=={3+p*p*t for t in range(p)},'higher tangent fills next digit')
 checks['Wieferich_1093_exact_log_valuation']=v

 # Independent rational boundary check, not importing the response verifier.
 q=Fraction(104,471);u=Fraction(13,471)
 A=Fraction(-34104960,367);B=Fraction(3945600,367);C=Fraction(-140160,367)
 D=Fraction(-27648000,471**2);E=Fraction(-829440,471)
 def res(x,m):return x.numerator*pow(x.denominator,-1,m)%m
 def phi(x,y):
  X=sum((q**i*z for i,z in enumerate(x)),Fraction(0));Y=sum((q**i*z for i,z in enumerate(y)),Fraction(0))
  return A*X+B*Y+D*X*Y
 for r in range(1,5):
  seen13=set();seen8=set();paired=set()
  for s in product((0,1),repeat=2*r):
   n=len(s)
   g=sum((A*q**i+B*q**(n-1-i)+C*q**(n-1))*s[i] for i in range(n))
   g+=sum((E*q**(n-2) if j==i+1 else D*q**(n-1-j+i))*s[i]*s[j] for i in range(n) for j in range(i+1,n))
   f=g+13*u**(n-1)*(sum(s)-9*sum(a*b for a,b in zip(s,s[1:])))
   z=phi(s[:r],s[r:][::-1])
   require(res(f-z,13**r)==0,'finite 13-adic precision')
   require(res((g-z)/128,8**r)==0,'finite 2-adic precision')
   a=res(z,13**r);b=res(z/128,8**r)
   seen13.add(a);seen8.add(b);paired.add((a,b))
  require(len(seen13)==len(seen8)==len(paired)==4**r,'shared address cardinality')
 checks['independent_reset_depths']=[1,2,3,4]
 # All-order additive proof reduces to these fixed branch tables.
 roots=[0,3,11,5]
 derivative=[{0,11,3,1},{0,2,3,5},{0,11,7,5},{0,2,7,9}]
 for target in range(13):
  good=[idx for idx in product(range(4),repeat=3) if sum(roots[i] for i in idx)%13==target and {sum(ds)%13 for ds in product(*(derivative[i] for i in idx))}==set(range(13))]
  require(bool(good),'persistent triple derivative coverage')
 require({(x+y)%13 for x in roots for y in roots}!=set(range(13)),'two13 summands insufficient')
 dyadic={0,5,7,4}
 require({(x+y)%8 for x in dyadic for y in dyadic}==set(range(8)),'two dyadic digits cover')
 checks['additive_basis_orders']={'2':2,'13':3}

 print(json.dumps(checks,sort_keys=True,indent=2))

if __name__=='__main__':main()
