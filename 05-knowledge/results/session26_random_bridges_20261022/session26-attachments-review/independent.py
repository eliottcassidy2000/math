from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial,comb
from hashlib import sha256
from pathlib import Path
src=Path(__file__).resolve().parent.parent/'session26-attachments'
for name in ('proof.md','audit.py','audit.out'):
 print('SHA256',name,sha256((src/name).read_bytes()).hexdigest())
def require(c,msg):
 if not c: raise RuntimeError(msg)
def matrix(k,mask):
 a=[[False]*k for _ in range(k)]
 for i,(u,v) in enumerate(combinations(range(k),2)):
  a[u][v]=bool(mask>>i&1);a[v][u]=not a[u][v]
 return a
# Different profile mechanism: count paths on each subset, then ordered subset partitions.
def covers(a):
 k=len(a);full=(1<<k)-1;paths=[0]*(1<<k)
 for s in range(1,1<<k):
  vs=[v for v in range(k) if s>>v&1]
  paths[s]=sum(all(a[u][v] for u,v in zip(w,w[1:])) for w in permutations(vs))
 dp=[[0]*(1<<k) for _ in range(k+1)];dp[0][0]=1
 for r in range(1,k+1):
  for s in range(1,1<<k):
   t=s
   while t:
    dp[r][s]+=paths[t]*dp[r-1][s^t];t=(t-1)&s
 return [dp[r][full] for r in range(k+1)]
def spectrum(p):
 k=len(p)-1;out=[F(0)]*(k+2)
 # Expand equation (3), independently from the producer's surjection expansion.
 for r in range(1,k+1):
  for b in range(r):
   v=F(p[r]*comb(r-1,b)*(-1)**(r-1-b),4**r)
   for d,c in enumerate((1,2,1)):out[b+d]+=v*c
 return out
def value(c,m): return sum(v*j**m for j,v in enumerate(c))
# Enumerate actual vertex words, not masks, slots, runs, or Stirling numbers.
# Each core-core edge/base-base edge has deterministic legality; each cross edge contributes 1/2.
def wordmean(a,m):
 k=len(a);s=F(0)
 for w in permutations(range(k+m)):
  cross=0;good=True
  for u,v in zip(w,w[1:]):
   if u<k and v<k:
    if not a[u][v]:good=False;break
   elif u>=k and v>=k:
    if u>v:good=False;break
   else:cross+=1
  if good:s+=F(1,2**cross)
 return s
n=0;words=0;tomography=0
for k in range(1,6):
 for mask in range(1<<(k*(k-1)//2)):
  a=matrix(k,mask);p=covers(a);c=spectrum(p);n+=1
  require(c[-1]==F(factorial(k),4**k),'leading')
  require(c[-2]==F((k+1)*factorial(k),4**k),'second')
  cyc=sum(sum(a[u][v] for u in t for v in t if u!=v and u==t[0])==1 and all(sum(a[u][v] for v in t if u!=v)==1 for u in t) for t in combinations(range(k),3))
  if k>=3:require(p[k-2]==factorial(k-2)*(comb(k,3)+2*cyc+3*comb(k,4)),'triangle')
  for z in range(-k-4,k+4):
   poly=lambda x:sum(v*x**j for j,v in enumerate(c))
   require(poly(-2-z)==(-1)**(k+1)*poly(z),'reflection')
  # Independent inversion: solve low coefficients using the polynomial expansion of each basis profile.
  h=(k+1)//2;low=[0]*(k+1)
  for m in range(h):
   unit=[0]*(k+1);unit[m+1]=1
   coeff=value(spectrum(unit),m)
   require(coeff==F(factorial(m),4**m),'triangular diagonal')
   low[m+1]=(value(c,m)-value(spectrum(low),m))/coeff
  B=[F(0)]*k
  for b in range(h):
   B[b]=low[b+1]-sum(comb(k-1-j,b-j)*B[j] for j in range(b))
   B[k-1-b]=B[b]
  require([0]+[sum(comb(k-1-b,r-1-b)*B[b] for b in range(r)) for r in range(1,k+1)]==p,'tomography');tomography+=1
  ms=range(4) if k<=4 else range(2)
  for m in ms:
   require(wordmean(a,m)==value(c,m),f'word mean {k} {mask} {m}');words+=1
print('PASS subset-partition profiles, triangle and top channels, polynomial reflection:',n)
print('PASS independent weighted full-word expectation pairs:',words)
print('PASS half-prefix tomography:',tomography)
# Hankel determinant by permutation expansion (not Gaussian elimination or Vandermonde).
p=covers(matrix(4,63));c=spectrum(p)
require(c==[F(x,32) for x in (-5,-9,6,22,15,3)],'T4')
H=[[value(c,i+j+1) for j in range(5)] for i in range(5)]
d=F(0)
for w in permutations(range(5)):
 sign=(-1)**sum(w[i]>w[j] for i in range(5) for j in range(i+1,5));term=F(sign)
 for i,j in enumerate(w):term*=H[i][j]
 d+=term
require(d==F(-16238475,1024),'Hankel')
print('PASS T4 Hankel determinant by Leibniz expansion:',d)
# Fair mean versus real integer observations and conditioning hostile.
a=matrix(3,5);hist={};strongtotal=strongn=0
for mask in range(8):
 b=[row+[bool(mask>>i&1)] for i,row in enumerate(a)]+[[not bool(mask>>i&1) for i in range(3)]+[False]]
 H=sum(all(b[u][v] for u,v in zip(w,w[1:])) for w in permutations(range(4)))
 hist[H]=hist.get(H,0)+1
 reach=[{i}|{j for j in range(4) if b[i][j]} for i in range(4)]
 for t in range(4):
  for i in range(4):
   if t in reach[i]:reach[i]|=reach[t]
 if all(len(x)==4 for x in reach):strongtotal+=H;strongn+=1
require(hist=={3:2,5:6},'integer outcomes')
print('PASS C3+m1 distribution:',hist,'unconditional mean',F(sum(h*n for h,n in hist.items()),8),'strong-conditioned mean',F(strongtotal,strongn))
print('PASS all checks explicit exceptions, safe under python -O')
