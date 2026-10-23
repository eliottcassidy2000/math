"""Independent exact signed-Collatz controls; no repository imports."""
from itertools import product
from collections import Counter
from fractions import Fraction
import json

def need(c, label):
    if not c: raise ArithmeticError(label)
def v2(n):
    need(n != 0,'nonzero')
    return (abs(n)&-abs(n)).bit_length()-1
def stats(w):
    A=B=0
    for a in w:
        B=3*B+2**A
        A+=a
    return A,B,2**A-3**len(w)
def cylinder(w,s):
    A,B,D=stats(w); M=2**(A+1)
    return ((2**A-s*B)*pow(3**len(w),-1,M))%M,M
def actual(n,s,k):
    w=[]
    for _ in range(k):
        q=3*n+s;a=v2(q);w.append(a);n=q//2**a
    return tuple(w),n
def prime(n):
    return n>=2 and all(n%d for d in range(2,__import__('math').isqrt(n)+1))
c=Counter(); cycles=[]; examples=[]
for k in range(1,6):
 for w in product(range(1,5),repeat=k):
    A,B,D=stats(w)
    rp,M=cylinder(w,1);rm,_=cylinder(w,-1)
    need((rp+rm)%M==0,'signed reflection')
    for s,r in ((1,rp),(-1,rm)):
      for t in (0,1,2,11):
        n=r+M*t;aw,z=actual(n,s,k)
        need(aw==w,'exact cylinder')
        need(2**A*z==3**k*n+s*B,'affine endpoint')
        c['cylinder_endpoint_checks']+=1
      # Neighbor half-cylinder has incorrect final oddness.
      other=(r+M//2)%M
      need(actual(other,s,k)[0]!=w,'last oddness bit is necessary')
      c['last_bit_controls']+=1
      if B%D==0 and s*D>0:
        n=s*B//D;aw,z=actual(n,s,k)
        need(aw==w and z==n,'integral fixed point actual cycle')
        c['integral_cycle_checks']+=1
for w in ((1,), (1,2), (1,1,1,2,1,1,4)):
 A,B,D=stats(w);n=-B//D
 need(B%D==0 and n>0,'known minus cycle')
 aw,z=actual(n,-1,len(w));need(aw==w and z==n,'known cycle')
 cycles.append({'word':w,'A':A,'B':B,'D':D,'start':n})
 for s in (-1,1):
    r,M=cylinder(w,s)
    p=next(n for n in range(r,r+M*10000,M) if prime(n))
    need(actual(p,s,len(w))[0]==w,'prime witness')
    # Guaranteed composite of the same residue, including r=1.
    composite=r*(1+M) if r>1 else (1+M)**2
    need(composite%M==r and not prime(composite),'composite witness')
    examples.append({'word':w,'sign':s,'residue':r,'modulus':M,'prime':p,'composite':composite})
# Entire odd residue range independently enumerated for small words.
for k in range(1,4):
 for w in product(range(1,4),repeat=k):
  for s in (-1,1):
    r,M=cylinder(w,s)
    hits=[n for n in range(1,M,2) if actual(n,s,k)[0]==w]
    need(hits==[r],'unique exhaustive residue')
    c['exhaustive_residue_words']+=1
print(json.dumps({'counts':dict(c),'cycles':cycles,'prime_composite_examples':examples},indent=2))
