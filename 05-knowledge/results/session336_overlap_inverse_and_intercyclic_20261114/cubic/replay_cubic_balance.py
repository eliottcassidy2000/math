import collections,itertools
from collections import Counter
from itertools import product, combinations
from random import Random


def arcs(word):
 n=len(word); rows=[0]*(n+2)
 for i,d in enumerate(word):
  rows[i]=sum(1<<j for j in range(i+1,n))
  for b in range(2):
   if d>>(1-b)&1: rows[n+b]|=1<<i
   else: rows[i]|=1<<(n+b)
 rows[n]=1<<n+1|rows[n]
 return rows

def frame_hist(rows):
 n=len(rows);s=n-2; hist=collections.Counter();ex={}
 for x,y in [(s,s+1),(s+1,s)]:
  for a in range(-1,s):
   if a>=0 and not(rows[a]>>x&1):continue
   for b,c in [(-1,-1)]+[(b,c) for b in range(s) for c in range(b,s)]:
    if b<0:
     if not(rows[x]>>y&1):continue
    elif not(rows[x]>>b&1 and rows[c]>>y&1):continue
    for d in range(-1,s):
     if d>=0 and not(rows[y]>>d&1):continue
     sets=[{t for t in g if t>=0} for g in [(a,), (b,c), (d,)]]
     if any(p&q for p,q in itertools.combinations(sets,2)):continue
     fre=[i for i in range(s) if all(i not in g for g in sets)];h=collections.Counter(); e={}
     for i in fre:
      elig=([0] if i<a else [])+([1] if b<i<c else [])+([2] if d>=0 and d<i else [])
      if not elig:break
      h[len(elig)]+=1;e[i]=elig
     else:
      hist[h[2],h[3]]+=1
      for i,L in e.items():
       if len(L)==3:ex[x,y,a,b,c,d,i]=e
 return hist,ex

def count_dp(rows):
 n=len(rows);W=1<<n;f=[[0]*n for j in range(W)]
 for x in range(n): f[1<<x][x]=1
 for M in range(1,W):
  for x,nx in enumerate(f[M]):
   if nx:
    m=rows[x]&~M
    while m:
     y=m&-m;m^=y
     f[M|y][y.bit_length()-1]+=nx
 return sum(f[-1])

def predicted(word):
 x=[(b,a) for b,a in combinations(range(len(word)),2) if word[b]&2 and not(word[a]&2)]
 y=[(d,c) for d,c in combinations(range(len(word)),2) if word[d]&1 and not(word[c]&1)]
 cnt=Counter()
 for b,a in x:
  for d,c in y:
   beta=min(a,c)-max(b,d)-1
   if beta>0 and len({a,b,c,d})==4:
    cnt[abs(b-d)+abs(a-c)-2,beta]+=2
 return cnt

def formula(a,b):
 return 1+(1<<b+3)+(1<<a)*((1<<b+3)+18*3**b-1)
def decode(N):
 if N<3 or not N%2:return []
 v=((N-1)&-(N-1)).bit_length()-1
 ans=[]
 lo=max(1,v-2);hi=N.bit_length()
 while lo<hi:
  mid=(lo+hi)//2
  if formula(v,mid)<N:lo=mid+1
  else:hi=mid
 if formula(v,lo)==N:ans.append((v,lo))
 if v>=4:
  b=v-3
  q,r=divmod(N-1-(1<<v),(1<<b+3)+18*3**b-1)
  a=q.bit_length()-1
  if q>0 and r==0 and q==1<<a and a>v:ans.append((a,b))
 if v>=5 and formula(v-1,v-4)==N:ans.append((v-1,v-4))
 return ans

def cubic(r):
 S=range(3*r);L=range(3*r,6*r);G=[set() for i in range(6*r)]
 for i,j,d in product(range(3),range(r),range(3)):
  a=i*r+j;b=3*r+d*r+(j-d*i)%r
  G[a].add(b);G[b].add(a)
 assert all(len(g)==3 for g in G)
 assert all(len(G[a]&G[b])<=1 for a,b in combinations(S,2))
 D=set(range(r))|set(range(3*r,4*r))
 assert len(D)*3==len(G) and all(i in D or G[i]&D for i in range(len(G)))
 for t in range(r+1):
  role=[0]*r+[1]*r+[2]*t+[0]*(r-t)
  c=Counter(len({role[j-3*r] for j in G[i]}) for i in S)
  a,b=3*r-3*t,3*t
  assert c==Counter({2:a,3:b})
  w=[2]+[0]*a+[1]+[0]*b+[0,0]
  rows=arcs(w)
  x,y=len(w),len(w)+1
  B,D,A,C=0,a+1,a+b+2,a+b+3
  assert rows[x]>>B&1 and rows[C]>>y&1 and rows[y]>>D&1 and rows[A]>>x&1
  eligible=lambda i: int(i<A)+int(B<i<C)+int(i>D)
  assert Counter(eligible(i) for i in range(len(w)) if i not in (B,D,A,C))==Counter({2:a,3:b})
 return [2**(3*r-3*t)*3**(3*t) for t in range(r+1)]

rng=Random(20261114)
examples=[list(w) for w in product(range(4),repeat=5)]
examples+=[[rng.randrange(4) for _ in range(rng.randint(4,10))] for i in range(128)]
for w in examples:
 h,_=frame_hist(arcs(w))
 assert Counter({q:m for q,m in h.items() if q[1]})==predicted(w)
for word,N in [('12111111',385),('11111120',573)]:
 w=list(map(int,word));h,_=frame_hist(arcs(w))
 assert count_dp(arcs(w))==N and not any(q[1] for q in h)
print(len(examples), 'exact triangle-interval/frame histogram comparisons; 385 and 573 remain beta-zero.')
for a,b in product((1,2,3),repeat=2):
 w=[2]+[0]*a+[1]+[0]*b+[0,0]
 assert count_dp(arcs(w))==formula(a,b)
for a,b in product(range(1,41),repeat=2):
 N=formula(a,b);sol=decode(N)
 assert (a,b) in sol and all(formula(*p)==N for p in sol)
print('Independent interpolation grid and 1600 three-branch inverse checks passed.')
print('36-vertex fixed cubic lift has seven frame multiplicities:',cubic(6))
