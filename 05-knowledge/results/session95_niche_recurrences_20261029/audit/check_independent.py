"""Independent Session95 controls, no producer imports. Python standard library."""
from fractions import Fraction as Q
from itertools import product
from collections import defaultdict
from math import comb,gcd,lcm
import json
checks=0

def check(ok, msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)

def det(a):
 a=[[Q(v) for v in row] for row in a]; n=len(a); d=Q(1)
 for j in range(n):
  pivot=next((i for i in range(j,n) if a[i][j]),None)
  if pivot is None:return Q(0)
  if pivot!=j:a[j],a[pivot]=a[pivot],a[j];d=-d
  v=a[j][j];d*=v
  for i in range(j+1,n):
   t=a[i][j]/v
   for k in range(j+1,n):a[i][k]-=t*a[j][k]
 return d

def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

# Signed source windows: compare literal moving sums with boundary-prefix formula.
cases=0
for span in range(5):
 for cs in product((-1,0,1),repeat=span+1):
  if not cs[0] or not cs[-1]:continue
  S=sum(cs); pref=[sum(cs[:j+1]) for j in range(span)]
  for n in range(1,9):
   sites=defaultdict(int)
   for i,c in enumerate(cs):
    for j in range(n):sites[i+j]+=c
   sites={i:c for i,c in sites.items() if c};p=len(sites)
   if n>span:
    pred=(n-span)*(S!=0)+sum(x!=0 for x in pref)+sum(S-x!=0 for x in pref)
    check(p==pred,('pole-count',cs,n,p,pred))
   if p<=9:
    moments=[sum(c*i**r for i,c in sites.items()) for r in range(2*p+3)]
    annih=[1]
    for i in sites:annih=mul(annih,[-i,1])
    for r in range(len(moments)-p):check(sum(a*moments[r+j] for j,a in enumerate(annih))==0,('annihilator',cs,n,r))
    if p:
     target=1
     for c in sites.values():target*=c
     ss=sorted(sites)
     for i,x in enumerate(ss):
      for y in ss[i+1:]:target*=(y-x)**2
     check(det([[moments[i+j] for j in range(p)] for i in range(p)])==target,('Hankel',cs,n))
   # Independent time difference in n from direct moments.
  for r in range(7):
   vals=[sum(c*(i+j)**r for i,c in enumerate(cs) for j in range(n)) for n in range(r+4)]
   for _ in range(r+2):vals=[b-a for a,b in zip(vals,vals[1:])]
   check(not any(vals),('time-polynomial',cs,r))
  cases+=1
# The zero-site convention is essential: exact-prefix order1, eventual order0.
check([sum(c*u**r for u,c in {0:7}.items()) for r in range(4)]==[7,0,0,0],'zero-pole prefix')

# Generalized valuation blocks: literal compositions by division clock.
def compositions(total, length):
 if length==0:
  if total==0:yield ()
 elif length==1:
  if total>=1:yield (total,)
 else:
  for a in range(1,total-length+2):
   for rest in compositions(total-a,length-1):yield (a,)+rest

def stopped(word,h):
 d=2**h-1;live=1
 for j,a in enumerate(word):
  live+=d-1 if a%h==0 else -1
  if live==0:return j==len(word)-1
 return False

def coeff_formula(h,m,t):
 d=2**h-1;n=d*m+1;L=(d-1)*m+1
 poly=[1]
 for _ in range(L):poly=mul(poly,[0]+[1]*(h-1))
 s=0
 for r,c in enumerate(poly):
  rem=t-h*m-r
  if rem>=0 and rem%h==0:s+=c*comb(n+rem//h-1,rem//h)
 return comb(n,m)*s//n
for h,m,maximum in [(2,0,14),(2,1,13),(2,2,12),(3,0,14),(3,1,12),(4,0,14)]:
 n=(2**h-1)*m+1
 for t in range(n,maximum+1):
  raw=sum(stopped(w,h) for w in compositions(t,n))
  check(raw==coeff_formula(h,m,t),('tree-clock',h,m,t,raw))
for h in range(2,9):
 d=2**h-1
 probs=[Q(2**(h-r),d) for r in range(1,h+1)]
 check(sum(probs)==1 and probs[-1]*d==1,('critical',h))
 leafmean=sum(Q(r,2**r) for r in range(1,h))/sum(Q(1,2**r) for r in range(1,h))
 for m in range(10):
  n=d*m+1;L=(d-1)*m+1
  explicit=h*m+L*leafmean+Q(n*h,d)
  simple=2*n+Q(d*(2-h)-h,d*(d-1))
  check(explicit==simple,('clock-mean',h,m))
# Every periodic valuation event: divisibility is necessary and sufficient.
for h in range(1,11):
 M=2**h-1
 seen={sum(2**(h-r) for r in range(1,h+1) if mask>>(r-1)&1) for mask in range(2**h)}
 check(seen==set(range(M+1)),('binary-event bijection',h))
 for d in range(3,32):
  valid=M%d==0
  if valid:
   K=M//d;S={r for r in range(1,h+1) if K>>(h-r)&1}
   check(sum(Q(2**(h-r),M) for r in S)==Q(1,d),('event probability',h,d))
   check(d%2==1 and S and len(S)<h,('nonconstant odd event',h,d))
for d,h in [(3,2),(5,4),(7,3),(9,6),(11,10),(15,4)]:
 check(pow(2,h,d)==1 and all(pow(2,j,d)!=1 for j in range(1,h)),('minimal clock',d,h))
for d,h in [(3,2),(5,4),(7,3),(9,6),(11,10),(15,4)]:
 M=2**h-1;K=M//d;S={r for r in range(1,h+1) if K>>(h-r)&1}
 law={r:Q(2**(h-r),M) for r in range(1,h+1)}
 pi=Q(1,d);pl=1-pi
 mi=sum(r*p for r,p in law.items() if r in S)/pi
 ml=sum(r*p for r,p in law.items() if r not in S)/pl
 vi=sum((r-mi)**2*p for r,p in law.items() if r in S)/pi
 vl=sum((r-ml)**2*p for r,p in law.items() if r not in S)/pl
 delta=mi-ml
 for m in range(8):
  n=d*m+1;L=(d-1)*m+1
  check(m*mi+L*ml+Q(n*h,M)==2*n-delta/d,('general mean',d,m))
  direct=m*vi+L*vl+Q(n*h*h*2**h,M*M)
  reduced=n*(2-Q(d-1,d*d)*delta**2)+(vl-vi)/d
  check(direct==reduced,('general variance',d,m))
# Direct ordered valuation compositions, not a producer generating function.
for d,h,cap in [(5,4,12),(7,3,12),(9,6,12)]:
 K=(2**h-1)//d;S={r for r in range(1,h+1) if K>>(h-r)&1}
 n=d+1
 B=[int(r in S) for r in range(h+1)];A=[int(r not in S and r>0) for r in range(h+1)]
 poly=B
 for _ in range(d):poly=mul(poly,A)
 for t in range(n,cap+1):
  raw=0
  for word in compositions(t,n):
   live=1
   for j,a in enumerate(word):
    residue=(a-1)%h+1
    live+=d-1 if residue in S else -1
    if live==0:break
   raw+=live==0 and j==n-1
  predicted=0
  for power,c in enumerate(poly):
   rem=t-power
   if rem>=0 and rem%h==0:predicted+=c*comb(n+rem//h-1,rem//h)
  check(raw==predicted,('general event clock',d,h,t,raw,predicted))
# Cyclic translation sources: a clock needs one fundamental spatial period.
cyclic_cases=0
for modulus in range(2,6):
 for source in product(range(modulus),repeat=modulus):
  spatial=next(d for d in range(1,modulus+1) if modulus%d==0 and all(source[i]==source[(i+d)%modulus] for i in range(modulus)))
  total=sum(source[:spatial])%modulus
  expected=spatial*modulus//gcd(modulus,total)
  accumulated=[0]*modulus;shift=list(source);period=None
  for t in range(1,modulus*modulus+1):
   accumulated=[(x+y)%modulus for x,y in zip(accumulated,shift)]
   shift=shift[-1:]+shift[:-1]
   if not any(accumulated) and tuple(shift)==source:period=t;break
  check(period==expected,('cyclic source clock',modulus,source,period,expected));cyclic_cases+=1
# One geometric valuation drives both events; coprime clocks do not imply independence.
for h in range(2,7):
 for k in range(2,7):
  H=lcm(h,k);M=2**H-1
  joint=sum((Q(2**(H-r),M) for r in range(1,H+1) if r%h==0 and r%k==0),Q(0))
  check(joint==Q(1,2**H-1),('shared clock',h,k))
  covariance=(2**h-1)*(2**k-1)*joint-1
  if (h,k)==(2,3):check(covariance==Q(-2,3),'coprime covariance')
  if (h,k)==(2,4):check(covariance==2,'nested covariance')
# Actual path-native complement vs converse: distinct parent cycle types.
S=(2,12,9,1);P=(3,0,1,2)
path=set(zip(P,P[1:]));C=[0]*4
for i in range(4):
 for j in range(i+1,4):
  a,b=(i,j) if S[i]>>j&1 else (j,i)
  if (a,b) not in path:a,b=b,a
  C[a]|=1<<b
check(all((C[a]>>b)&1 for a,b in [(0,1),(1,2),(0,2)]),'native parent is transitive')
check(all((S[a]>>b)&1 for a,b in [(0,1),(1,2),(2,0)]),'original cyclic parent')
print(json.dumps({'status':'PASS','checks':checks,'signed_source_cases':cases,'cyclic_source_cases':cyclic_cases,'producer_imports':False},sort_keys=True))
