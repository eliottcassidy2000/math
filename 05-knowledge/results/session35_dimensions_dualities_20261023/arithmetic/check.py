#!/usr/bin/env python3
"""Exact bounded tests, standard library only; no probabilistic primality."""
import json, math
from collections import Counter
from pathlib import Path

def primes_upto(n):
 s=bytearray(b'\1')*(n+1); s[:2]=b'\0\0'
 for p in range(2,math.isqrt(n)+1):
  if s[p]:s[p*p::p]=b'\0'*(((n-p*p)//p)+1)
 return [p for p in range(2,n+1) if s[p]]
primes=primes_upto(10**6)
def factor(n):
 out={}
 for p in primes:
  if p*p>n:break
  while n%p==0:out[p]=out.get(p,0)+1;n//=p
 if n>1:out[n]=1
 return out

def roots(a,n):return [r for r in range(n) if r*r%n==a]
def leg(a,p):return 0 if a%p==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def crt(u,v,p,q):return (u+p*((v-u)*pow(p,-1,q)%q))%(p*q)
rows=[]; blums=[]
# Complete universe 2<=x<y and x^y+y^x <=10^12: monotonicity
# gives x<=10, and y<=39 since 2^40>10^12.
for x in range(2,40):
 for y in range(x+1,40):
  n=x**y+y**x
  if n>10**12:break
  f=factor(n); assert math.prod(p**e for p,e in f.items())==n
  blum=len(f)==2 and all(e==1 and p%4==3 for p,e in f.items())
  if blum:
   assert math.gcd(x,y)==1 and (x+y)%2==1
   even=x if x%2==0 else y
   assert math.isqrt(even)**2!=even
   assert all(leg(-even,p)==1 for p in f)
   blums.append([x,y,n,f])
  rows.append([x,y,n,blum])
# The square base exclusion and splitting property apply to every primitive odd row.
for x,y,n,_ in rows:
 if (x+y)%2 and math.gcd(x,y)==1:
  even=x if x%2==0 else y
  assert all(leg(-even,p)==1 for p in factor(n))
for y in range(3,2000,2):
 residue=(pow(2,y,9)+y*y)%9
 assert (residue==0)==(y%18 in (11,13))
 if y%3:assert residue%3==0
# Every residue of every p*q for distinct Blum primes p<q<50.
checked=0
ps=[p for p in primes if p<50 and p%4==3]
for i,p in enumerate(ps):
 for q in ps[i+1:]:
  n=p*q; buckets=[[] for _ in range(n)]
  for r in range(n):buckets[r*r%n].append(r)
  for a,rs in enumerate(buckets):
   if not rs:continue
   ap,aq=a%p,a%q
   bp,bq=pow(ap,(p+1)//4,p),pow(aq,(q+1)//4,q)
   assert leg(bp,p)>=0 and leg(bq,q)>=0
   generated={crt(sp*bp%p,sq*bq%q,p,q) for sp in (1,-1) for sq in (1,-1)}
   assert generated==set(rs)
   assert len(rs)==2**((ap!=0)+(aq!=0))
   # Local sign-flip operations commute and preserve target, not just count.
   for r in rs:
    fp=crt(-r%p,r%q,p,q);fq=crt(r%p,-r%q,p,q)
    assert fp in rs and fq in rs
    assert crt(-fq%p,fq%q,p,q)==crt(fp%p,-fp%q,p,q)==(-r)%n
   checked+=1
# Explicit hostiles: mod9 repeated prime breaks 4/2/1 law; mod65 sign cannot select QR root.
assert roots(0,9)==[0,3,6]
assert leg(2,17)==leg(-2,17)==1
assert roots(3,57)==[]
assert (13*13+2)%57==0
out={'status':'FINITE-EXACT','leyland_bound':10**12,'leyland_rows':len(rows),'blum_rows':blums,'root_fibres_checked':checked,'roots57':{str(a):roots(a,57) for a in (0,1,9)},'square_image57_histogram':dict(Counter(len(roots(a,57)) for a in range(57))), 'base2_mod9':{str(y):(pow(2,y,9)+y*y)%9 for y in range(1,19,2)}}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
