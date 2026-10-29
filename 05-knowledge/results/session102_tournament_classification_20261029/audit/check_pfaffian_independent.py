"""Elimination-Pfaffian oracle independent of producer's recursive matching sum."""
from fractions import Fraction as F
from itertools import combinations,product
import random,json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def pf(a):
 a=[[F(x) for x in r] for r in a];ans=F(1)
 while a:
  n=len(a);need(n%2==0,'even pf domain')
  candidates=[j for j in range(1,n) if a[0][j]]
  if not candidates:return F(0)
  j=candidates[0]
  if j!=1:
   a[1],a[j]=a[j],a[1]
   for row in a:row[1],row[j]=row[j],row[1]
   ans=-ans
  c=a[0][1];ans*=c
  a=[[a[i][j]+(a[1][i]*a[0][j]-a[0][i]*a[1][j])/c for j in range(2,n)] for i in range(2,n)]
 return ans
def inverse(a):
 n=len(a);b=[[F(x) for x in r]+[F(i==j) for j in range(n)] for i,r in enumerate(a)]
 for i in range(n):
  j=next(j for j in range(i,n) if b[j][i]);b[i],b[j]=b[j],b[i];c=b[i][i];b[i]=[x/c for x in b[i]]
  for j in range(n):
   if j!=i:
    c=b[j][i];b[j]=[x-c*y for x,y in zip(b[j],b[i])]
 return [r[n:] for r in b]
def bil(x,I,y):return sum(x[i]*I[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))
def tournament(n,bits):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):a[i][j]=1 if bits>>k&1 else -1;a[j][i]=-a[i][j]
 return a
def combine(a,U,e):return [list(a[i])+list(U[i]) for i in range(len(a))]+[[-U[i][j] for i in range(len(a))]+list(e[j]) for j in range(len(e))]
def cofactor(a):return [(-1)**i*pf([[a[j][k] for k in range(len(a)) if k!=i] for j in range(len(a)) if j!=i]) for i in range(len(a))]
rng=random.Random(20261029);counts={'mixed_root':0,'even_transfer':0,'odd_pivots':0}
for m,n in [(2,1),(2,3),(4,1)]:
 for repeat in range(100):
  a=tournament(m,rng.randrange(1<<(m*(m-1)//2)));b=tournament(n,rng.randrange(1<<(n*(n-1)//2)))
  u=[rng.choice([-1,1]) for _ in range(m)];v=[rng.choice([-1,1]) for _ in range(n)];x=[rng.choice([-1,1]) for _ in range(m)];y=[rng.choice([-1,1]) for _ in range(n)]
  whole=combine(combine(a,[[s*t for t in v] for s in u],b),[[t] for t in x+y],[[0]])
  p=cofactor(b);expected=pf(a)*(sum(s*t for s,t in zip(p,y))+sum(s*t for s,t in zip(p,v))*bil(u,inverse(a),x))
  need(pf(whole)==expected,'mixed root sign');counts['mixed_root']+=1
for n in range(1,8):
 for k in range(1,5):
  if (n+k)%2:continue
  for repeat in range(20):
   a=tournament(n,rng.randrange(1<<(n*(n-1)//2)));e=tournament(k,rng.randrange(1<<(k*(k-1)//2)));U=[[rng.choice([-1,1]) for _ in range(k)] for _ in range(n)];actual=pf(combine(a,U,e))
   if n%2==0:
    I=inverse(a);cols=list(zip(*U));res=[[e[i][j]+bil(cols[i],I,cols[j]) for j in range(k)] for i in range(k)]
    need(actual==pf(a)*pf(res),'even multiport');counts['even_transfer']+=1
   else:
    for pivot in range(n):
     keep=[i for i in range(n) if i!=pivot];c=[[a[i][j] for j in keep] for i in keep];I=inverse(c);b=[a[i][pivot] for i in keep];cols=[[U[i][j] for i in keep] for j in range(k)]
     gamma=[U[pivot][j]+bil(b,I,cols[j]) for j in range(k)]
     res=[[0]+gamma]+[[-gamma[i]]+[e[i][j]+bil(cols[i],I,cols[j]) for j in range(k)] for i in range(k)]
     need(actual==(-1)**pivot*pf(c)*pf(res),'odd singular mode/pivot sign');counts['odd_pivots']+=1
print(json.dumps({'status':'FINITE-EXACT','seed':20261029,'oracle':'fractional two-by-two skew elimination, not recursive matching expansion','checks':counts},indent=2,sort_keys=True))
