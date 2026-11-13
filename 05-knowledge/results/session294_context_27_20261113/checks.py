from itertools import combinations,product
from math import gcd,comb,isqrt,factorial
from collections import Counter,defaultdict
import json
from pathlib import Path

def require(x):
 if not x: raise RuntimeError('failed check')

def D(n,f):
 a=[f(i) for i in range(n+1)]
 b=[]
 for _ in range(n):
  a=[y-x for x,y in zip(a,a[1:])]
  b.append(a[0])
 return b

def v3(a):
 if not a: return float('inf')
 s=0
 while a%3==0:
  a//=3;s+=1
 return s

def G(F):
 return gcd(*F),gcd(*(f*factorial(k) for k,f in enumerate(F)))

def conv(p,q):
 r=defaultdict(int)
 for k,a in p.items():
  for j,b in q.items():r[k+j]+=a*b
 return dict(r)

def path(l):
 ps=[{0:1}]
 for i in range(1,l+1):
  d=Counter()
  for j in range(i):
   for k,v in ps[j].items():d[k+1]+=v*(1 if j==i-1 else 2)
  ps.append(d)
 return ps[l]

def forest(parts):
 p={0:1}
 for l in parts:p=conv(p,path(l))
 return G([p.get(i,0)*factorial(i) for i in range(1,sum(parts)+1)])

def edges(V):
 return {(a,b) for a,b in combinations(V,2) if isqrt(a+b)**2==a+b}

scan=json.loads(Path(__file__).with_name('scan.json').read_text())
coeff=Counter()
for r in scan['order8']['found']:
 H,F=r['H'],r['F']
 t=(F[1]//2-H-126)//6
 A=((F[2]-5796)//36-t)//10
 B=t-5*A; C=H-1-A-B
 coeff[A,B,C]+=1
 f=lambda x:x**8+A*x**6+B*x**4+C*x**2
 require(D(8,f)==F)
 require(all(f(x)%27==(18 if x%9 in (4,5) else 0) for x in range(27)))
 require(all(f(pow(8,j,27))%27==0 for j in range(6)))
 for r,s in product(range(1,7),range(7)):
  w=lambda x:f(x)**r*(x**3+2*x)**s
  got=tuple(map(v3,G(D(8*r+3*s,w))))
  want=(min(2*r+2*s,3*r+s),min(2*r+2*s+1,3*r+s))
  require(got==want)
 for m in range(12):
  g=lambda x:f(x)*x**m
  require(G(D(8+m,g))==(9,27))

S=[]
for n in range(1,8):
 hit=0
 for a in product(range(27),repeat=(n-1)//2 if n%2 else n//2-1):
  ev=lambda x:x**n+sum(v*x**(n-2*j) for j,v in enumerate(a,1))
  if all(ev(x)%27==0 for x in range(1,4)):
   if n==7:require(a[0]%9==a[1]%9==4 and a[2]%9==0)
   hit+=1
 S.append(hit)
require(S==[0]*6+[9])

require(edges([1,2,4,6,8,9,10,11])=={(1,8),(6,10)})
require(forest([4,2,2])==forest([4,3,2])==(48,96))
require(forest([2,2,1,1,1,1])==(2880,345600))
levels=sorted(a*9**j for j in range(6) for a in (2,14))
require(edges(levels)=={(2*9**j,14*9**j) for j in range(6)})

from scan import profile
F=profile([56,117,129,22,100,76,141,59])[1]
# An ordered path on the quotient cyclically visits its maximal runs.
words=Counter()
for start in range(3):
 a=[0]*3
 for size in range(1,25):
  a[(start+size-1)%3]+=1
  if max(a)<=8:words[tuple(a)]+=1
state={(0,0,0):1}; cyclic=[]
for k in range(1,25):
 nxt=defaultdict(int)
 for a,h in state.items():
  for b,j in words.items():
   c=tuple(x+y for x,y in zip(a,b))
   if max(c)<=8:nxt[c]+=h*j
 state=nxt
 cyclic.append(sum(h*F[a[0]-1]*F[a[1]-1]*F[a[2]-1] for a,h in state.items() if min(a)>0))
require(cyclic[0]==163567386770762025 and G(cyclic)==(6561,6561))
print(json.dumps(dict(abstract_residue_counts=S,order8_polynomial_multiplicities={str(k):v for k,v in sorted(coeff.items())},
                      cyclic_three_copies=dict(H=cyclic[0],gcds=G(cyclic)),checks='pass'),indent=2))
