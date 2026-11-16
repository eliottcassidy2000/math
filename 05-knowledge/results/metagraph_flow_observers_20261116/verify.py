from itertools import product, permutations, combinations
from functools import lru_cache
from math import gcd, factorial, isqrt

def dp(r,root=None):
 n=len(r);N=1<<n
 d=[[0]*n for _ in range(N)]
 for j in range(n) if root is None else (root,):d[1<<j][j]=1
 for i,row in enumerate(d):
  for j,v in enumerate(row):
   if not v:continue
   z=r[j]&~i
   while z:
    b=z&-z;z^=b
    d[i|b][b.bit_length()-1]+=v
 return [1]+[sum(t) for t in d[1:]]
def prof(t,x):
 d=dp(t);z=dp(t,x);N=len(d)-1
 ev,od=(sum(z[s]*d[N^s] for s in range(1,N+1) if s.bit_count()%2==k) for k in (0,1))
 return d[-1],d[N^(1<<x)],ev,od
def instr(t,x,w):
 for b in map(int,w):
  n=len(t);m=0 if b else 1<<x
  t=tuple(a|((not(m>>j&1))<<n) for j,a in enumerate(t))+(m,)
 return t
def subset(t,order):
 return tuple(sum(((t[u]>>v)&1)<<j for j,v in enumerate(order)) for u in order)
def op(t):return tuple(((1<<len(t))-1)^v^(1<<j) for j,v in enumerate(t))
def tt(n):return tuple((1<<n)-(1<<(j+1)) for j in range(n))
def tau(t):
 n=len(t)
 assert all(t[i]>>(i+1)&1 for i in range(n-1))
 return tuple(v^sum(1<<j for j in range(n) if abs(j-i)>1) for i,v in enumerate(t))
def join(a,b):
 n,m=len(a),len(b)
 return tuple(r|(((1<<m)-1)<<n) for r in a)+tuple(r<<n for r in b)
def perm_by_degree(t):
 groups=[tuple(i for i,r in enumerate(t) if r.bit_count()==d) for d in range(len(t))]
 for p in product(*(tuple(permutations(c)) for c in groups if c)):yield sum(p,())
def canon(t):return min(subset(t,p) for p in perm_by_degree(t))
def rigid(t):
 ref=subset(t,tuple(sorted(range(len(t)),key=lambda i:t[i].bit_count())))
 return sum(subset(t,p)==ref for p in perm_by_degree(t))==1

def metagraph():
 def c(t,m):return dp(tau(join(t,tt(m))))[-1]
 a,b=(10,4,9,2),(6,12,8,1)
 assert tau(a)==b and canon(a)==canon(b)
 for t in (a,b):assert subset(t,(3,2,1,0))==op(t)
 assert (c(a,1),c(b,1))==(11,13)
 for x in ((26,4,25,18,2),(22,4,24,19,2)):
  u=op(tau(x));k=len(x)-1
  h,g,r,s=prof(u,k);e=s-r;R=s+r
  assert dp(x)[-1]==9 and (h,e,R)==(13,2,18)
  a,b,d=(c(x,i) for i in (1,2,3))
  assert (a,b,d)==(29,49,91)
  assert e==(d-3*a)//2 and R==b-a-e and h==2*a+2*e-b
  assert dp(tau(subset(x,range(k))))[-1]==g
  assert prof(instr(u,k,'00'),k)[0]==(55 if g==5 else 59)
  F=[1,1,1]
  for i in range(3,8):F.append(sum(F[-3:]))
  Q=lambda m:sum(F[:m])
  for m in range(2,8):
   assert c(x,m)==h*F[m]-e*Q(m-2)+R*Q(m)
 S=(12,21,24,18,1);P=(0,2,3,1,4,5,6,7)
 A,B=[instr(S,0,w) for w in ('100','000')]
 assert prof(A,0)==(97,3,54,57) and prof(B,0)==(97,3,48,51)
 assert sorted(r.bit_count() for r in A)==sorted(r.bit_count() for r in B)
 assert canon(A) not in (canon(B),canon(op(B)))
 assert [dp(tau(subset(t,P)))[-1] for t in (A,B)]==[199,197]
 print('Blue loop lift 11/13; identical padding observables split at 55/59; H97 collision lifts to partner H199/197.')

def cover(t):
 n=len(t);H=dp(t);F=[[0]*(1<<n) for _ in range(n+1)];F[0][0]=1
 for k in range(1,n+1):
  for v in range(1,1<<n):
   s=v
   while s:
    F[k][v]+=H[s]*F[k-1][v^s]
    s=v&(s-1)
 return [F[k][-1] for k in range(1,n+1)]

def probe(t,m):
 n=len(t)
 a=[u|sum(1<<(n+i) for i in range(m)) for u in t]
 a += [(1<<(n+m+i))+sum(1<<(n+j) for j in range(i+1,m)) for i in range(m)]
 a += [((1<<n)-1)+sum(1<<(n+j) for j in range(m) if j!=i)+sum(1<<(n+m+j) for j in range(i+1,m)) for i in range(m)]
 return tuple(a)

def skeletons(host,k):
 # Replace the module by k ordered anchor tokens, and count words
 # in which all internal gaps are nonempty, with anchors in order.
 n=len(host);m=n-1
 def edge(u,v):
  if u>=m and v>=m:return u<v
  return (host[0 if u>=m else u+1]>>(0 if v>=m else v+1))&1
 count=0
 for w in permutations(range(m+k)):
  if tuple(i for i in w if i>=m)!=tuple(range(m,m+k)):continue
  if any(i>=m and j>=m for i,j in zip(w,w[1:])):continue
  count+=all(edge(i,j) for i,j in zip(w,w[1:]))
 return count

def modules():
 trap=(2,4,73,19,39,15,59); rel=(2,12,17,37,11,87,31)
 D=lambda t:[factorial(i)*v for i,v in enumerate(cover(t))]
 assert rigid(trap) and rigid(rel)
 assert gcd(*D(trap))==7 and gcd(*D(rel))==1
 assert [dp(probe(trap,m))[-1] for m in range(3)]==[49,861,15925]
 assert [dp(probe(rel,m))[-1] for m in range(3)]==[49,873,16069]
 assert all(16069%d for d in range(2,isqrt(16069)+1))
 for t in (rel,trap):
  q=cover(t)
  for m in (0,1,2):
   host=probe((0,),m)
   coef=[skeletons(host,k) for k in range(1,m+2)]
   assert coef[-1]==factorial(m)
   assert all(a%factorial(i)==0 for i,a in enumerate(coef))
   assert sum(a*b for a,b in zip(coef,q))==dp(probe(t,m))[-1]
 H63=(102,252,248,241,225,192,128,1)
 assert rigid(H63) and gcd(*D(H63))==9 and dp(probe(H63,1))[-1]==1935
 delay=(24,25,67,84,100,15,163,63)
 assert [dp(probe(delay,i))[-1] for i in range(3)]==[343,5999,112983]
 assert [j%7 for j in D(delay)][:3]==[0,0,3]
 # Ordered-cover coefficients also recover polynomial values at positive integers.
 from math import comb
 A=(2,20,9,3,45,79,31);B=(10,12,17,68,11,31,183,63)
 def P(t,k):return sum(comb(k,i+1)*a for i,a in enumerate(cover(t)) if i<k)
 assert rigid(A) and rigid(B)
 assert [P(A,k) for k in (1,2,3)]==[35,676,5649]
 assert [P(B,k) for k in (1,2,3)]==[119,2576,25191]
 assert [P(A,k)%7 for k in (1,2,3)]==[0,4,0]
 assert [P(B,k)%7 for k in (1,2,3)]==[0,0,5]
 print('Rigid H49 modules have contents 7 and 1; latter gives prime 16069; H63 content=9; 35/119 supports disjoint mod 7.')
if __name__=='__main__':
 metagraph();modules()
