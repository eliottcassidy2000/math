#!/usr/bin/env python3
"""Independent finite controls of S192's exact digit/product identities.
All gates remain active under python -O. No third-party packages.
"""
import hashlib,itertools,json,math
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS={}
def require(ok,msg,key):
 if not ok:raise RuntimeError(msg)
 CHECKS[key]=CHECKS.get(key,0)+1

def count(w):
 h,A,p=1,0,0
 for b in w:
  if not b:h+=A+p
  A=2*A+b;p=b
 return h

def slot(w):
 if not w:return 1
 return w[0]+1-w[-1]+sum(a==0 and b==1 for a,b in zip(w,w[1:]))+sum(1<<(i-j-1) for j in range(len(w)) for i in range(j+1,len(w)) if w[j] and not w[i])

def graph(w):
 rows=[0]*(len(w)+1)
 for i,b in enumerate(w,1):
  if b:rows[0]|=1<<i
  else:rows[i]|=1
  for j in range(i+1,len(rows)):rows[i]|=1<<j
 return rows

def dp(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   ways=D[mask][v]
   if not ways:continue
   avail=rows[v]&~mask
   while avail:
    bit=avail&-avail;avail-=bit
    D[mask|bit][bit.bit_length()-1]+=ways
 return sum(D[-1])

def literal(rows):
 return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:])) for p in itertools.permutations(range(len(rows))))

def bits(A,n):return tuple((A>>j)&1 for j in range(n-1,-1,-1))
def integer(w):
 x=0
 for b in w:x=2*x+b
 return x

def gadget(k,j,e):return bits((1<<j)-2,k)+((1,0) if e else (0,1))+(1,)*k

def bank(k,e):return sum((gadget(k,j,e>>(j-1)&1) for j in range(1,k)),())
def floor_root_scaled(h,power):
 # floor(sqrt(h*2**power)), exact even if power is negative.
 return math.isqrt(h<<power) if power>=0 else math.isqrt(h>>(-power))
def ends(n,a,b,h):
 A=min((1<<a)-1,floor_root_scaled(h,2*a-n))
 B=min((1<<b)-1,floor_root_scaled(h,2*b-n))
 return bits(A,a),tuple(1-x for x in bits(B,b)[::-1])
def selector(n,k,r,h):
 L=2*k*k-2;a=(n-L)//2;b=n-L-a
 u,v=ends(n,a,b,h);base=u+bank(k,0)+v
 e=((r-count(base))% (1<<k))//2
 w=u+bank(k,e)+v
 E=(1<<(n-a))+(1<<(n-b))+n//2+1
 return w,E

def product_data(w):
 n=len(w);A=integer(w);B=integer(tuple(1-b for b in w[::-1]))
 # H*2**n = A*B + boundary*2**n - R*2**n.
 Rscaled=sum(1<<(n+i-j-1) for i in range(n) for j in range(i+1,n) if not w[i] and w[j])
 N=sum(a==0 and b==1 for a,b in zip(w,w[1:]))
 return A,B,Rscaled,N

for n in range(1,13):
 for w in itertools.product((0,1),repeat=n):
  h=count(w);A,B,R,N=product_data(w)
  require(h==slot(w),'slot mismatch','slot')
  require(h*(1<<n)==A*B+(w[0]+1-w[-1]+N)*(1<<n)-R,'product identity','product')
  require((R==0 if N==0 else 0<=R<N*(1<<n)),'reverse correction bound','reverse_bound')
  if w[0]==1 and w[-1]==0:
   D=h*(1<<n)-A*B
   require(2*(1<<n)<=D<=(n//2+1)*(1<<n),'strong product bound','strong_error')
  if n<=7:require(dp(graph(w))==h,'vertex DP','dp')
  if n<=5:require(literal(graph(w))==h,'literal permutations','permutations')

# Complete low-order deterministic selector tests; arbitrary middles separately.
for n in range(2,13):
 for h in range(1<<(n-2),(1<<n)+1):
  w,E=selector(n,1,1,h)
  require(len(w)==n and w[0]==1 and w[-1]==0,'selector type','selector_type')
  require(abs(count(w)-h)<=E,'selector magnitude','selector_magnitude')
for a,b,L in itertools.product(range(1,5),range(1,5),range(5)):
 n=a+b+L
 for h in (1<<(n-2),3*(1<<(n-2)),1<<n):
  u,v=ends(n,a,b,h)
  E=(1<<(n-a))+(1<<(n-b))+n//2+1
  for mid in itertools.product((0,1),repeat=L):
   require(abs(count(u+mid+v)-h)<=E,'arbitrary middle','middle')

# Each gadget in every short exterior; independent choices in full banks.
for k in range(2,9):
 for j in range(1,k):
  for ulen,vlen in itertools.product(range(4),repeat=2):
   for u in itertools.product((0,1),repeat=ulen):
    for v in itertools.product((0,1),repeat=vlen):
     delta=count(u+gadget(k,j,1)+v)-count(u+gadget(k,j,0)+v)
     require(delta%(1<<k)==1<<j,'insulated switch','switch')
 for u,v in [((1,),(0,)),((1,0,1),(0,1,0)),((0,1,1),(1,0,1))]:
  base=count(u+bank(k,0)+v);values=[]
  for e in range(1<<(k-1)):
   val=count(u+bank(k,e)+v)
   require((val-base)%(1<<k)==2*e,'whole bank','bank')
   values.append(val%(1<<k))
  require(sorted(values)==list(range(1,1<<k,2)),'bank onto','bank_onto')

simultaneous=[]
for k in range(1,8):
 L=2*k*k-2
 for n in (L+2,L+3,L+10,L+21):
  targets=[1<<(n-2),3*(1<<(n-2)),(1<<n)-1,1<<n]
  for h in targets:
   for r in range(1,1<<k,2):
    w,E=selector(n,k,r,h);value=count(w)
    require(len(w)==n and w[0]==1 and w[-1]==0,'joint strong type','joint_type')
    require(value%(1<<k)==r and abs(value-h)<=E,'joint residue magnitude','joint')
  simultaneous.append({'k':k,'chain_length':n,'middle_length':L,'error_bound':E})

moments=[]
for n in range(2,19):
 vals=[count((1,)+mid+(0,)) for mid in itertools.product((0,1),repeat=n-2)]
 mean=Fraction(sum(vals),len(vals));second=Fraction(sum(x*x for x in vals),len(vals))
 if n>=3:require(mean==9*Fraction(2)**(n-4)+Fraction(1,2),'exact finite mean','mean')
 unique=sorted(set(vals));gaps=[(y-x,x,y) for x,y in zip(unique,unique[1:])]
 gap=max(gaps,default=(0,None,None));E=(1<<(n//2))+(1<<((n+1)//2))+n//2+1
 require(gap[0]<=2*E,'maximum gap bound','gap')
 moments.append({'chain_length':n,'presentations':len(vals),'distinct_counts':len(unique),'mean':str(mean),'normalized_variance':str((second-mean*mean)/ (1<<(2*n))),'maximum_gap':gap})

# Exact independent hostile: Christoffel carry coprimality does not transfer to H.
w=(1,0,0);clock=2**3-3**1
require(count(w)==5 and math.gcd(count(w),clock)==5,'Christoffel hostile','hostile')
receipt={'status':'finite exact controls of elementary all-order proofs; no global spectrum claim','checks':CHECKS,'moments_and_gaps':moments,'joint_selector_controls':simultaneous,'limit_mean':'9/16','limit_variance':'55/2304','sources':{'S25':'../session25-spectrum/PROOF.md','S191':'../session191-arithmetic/PROOF.md'},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':CHECKS,'receipt_sha256':hashlib.sha256((HERE/'receipt.json').read_bytes()).hexdigest()},sort_keys=True))
