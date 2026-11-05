from pathlib import Path
from itertools import permutations,combinations,product
from functools import lru_cache
import json
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
@lru_cache(None)
def Hsub(a):
 n=len(a);need(all(not(a[i]>>i&1)for i in range(n))and all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n)for j in range(i+1,n)),'literal tournament');N=1<<n;dp=[[0]*n for _ in range(N)];H=[1]+[0]*(N-1)
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   R=S^(1<<v)
   if not R:dp[S][v]=1
   else:
    U=R&a[v]
    while U:
     b=U&-U;U-=b;dp[S][v]+=dp[R][b.bit_length()-1]
  H[S]=sum(dp[S])
 return tuple(H)
def count(a):return Hsub(tuple(a))[-1]
def b2(a):
 H=Hsub(tuple(a));F=len(H)-1;return sum(H[S]*H[F^S]for S in range(F+1))
def graphs(n):
 edges=list(combinations(range(n),2))
 for bits in range(1<<len(edges)):
  a=[0]*n
  for k,(u,v)in enumerate(edges):a[u if bits>>k&1 else v]|=1<<(v if bits>>k&1 else u)
  yield tuple(a)
def substitute(Q,r,T):
 blocks=[T if v==r else (0,)for v in range(len(Q))];starts=[];n=0
 for B in blocks:starts.append(n);n+=len(B)
 a=[0]*n
 for v,B in enumerate(blocks):
  for i,row in enumerate(B):a[starts[v]+i]|=row<<starts[v]
  for w,C in enumerate(blocks):
   if Q[v]>>w&1:
    for i in range(len(B)):a[starts[v]+i]|=((1<<len(C))-1)<<starts[w]
 return tuple(a)
def weight(Q,r):
 F=(1<<len(Q))-1;H=Hsub(Q);tri=[]
 for o in range(len(Q)):
  for i in range(len(Q)):
   if Q[r]>>o&1 and Q[o]>>i&1 and Q[i]>>r&1:tri.append((o,i,H[F^(1<<r)^(1<<o)^(1<<i)]))
 return H[F^(1<<r)],sum(t[2]for t in tri),tri
contexts=[(Q,r)for n in range(1,5)for Q in graphs(n)for r in range(n)]
inputs=[T for n in range(1,4)for T in graphs(n)]+[(6,12,24,17,3)]
weighted=0;strict=0
for Q,r in contexts:
 d,w,tri=weight(Q,r);need(w>=len(tri),'weighted root triangle count')
 for T in inputs:
  h=count(T);b=b2(T);out=count(substitute(Q,r,T));need(out>=d*h+w*b,'weighted physical inequality');strict+=out>d*h+w*b;weighted+=1
# A root triangle with an exterior C3 gives a strictly weighted coefficient3.
Q=(58,60,57,16,32,8);d,w,tri=weight(Q,0)
need(d==3 and w==3 and len(tri)==1,'exterior weighted triangle')
for T in inputs:need(count(substitute(Q,0,T))==3*(count(T)+b2(T)),'weighted equality family')
# Literal fixed-label cut injection for every labelled input through5.
cutgraphs=cuts=0
for n in range(1,6):
 for T in graphs(n):
  paths=[w for w in permutations(range(n))if all(T[u]>>v&1 for u,v in zip(w,w[1:]))];configs=set()
  for w in paths:
   for j in range(n+1):configs.add((w[:j],w[j:]))
  need(len(paths)==count(T) and len(configs)==(n+1)*len(paths),'injective labelled cuts')
  need(b2(T)>=len(configs) and b2(T)>=2**n,'both independent profile bounds');cutgraphs+=1;cuts+=len(configs)
# All choices in a small adaptive template bank, with actual internal padding.
bank=[((2,0),0),((2,4,1),0),((2,4,9,3),0),((14,4,8,2),0)]
adaptive=active=0
for T0 in [(0,),(2,0),(2,4,1)]:
 for indices in product(range(4),repeat=3):
  T=T0;h0=count(T);n0=len(T);lower=h0;j=0;hist=[]
  for k in indices:
   Q,r=bank[k];d,w,_=weight(Q,r);h=count(T);n=len(T);b=b2(T);U=substitute(Q,r,T);hu=count(U)
   need(hu>=h,'context monotonicity')
   if w:
    j+=1;lower*=n0+2*j;need(hu>=(n+2)*h and hu>=h+2**n,'local active bounds');hist.append((n,h,d,w,b));active+=1
   T=U
  M=count(T);need(M>=lower,'adaptive product')
  for n,h,d,w,b in hist:need(2**n<=M-h and w<=(M-d*h)//b,'internal and weighted pruning')
  adaptive+=1
# Integer target-depth boundaries, including no-active and M<h0 cases.
cutoffs=0
for n in range(1,7):
 for h in (1,3,5,9):
  for M in range(1001):
   if M<h:cutoffs+=1;continue
   z=h;s=0
   while z*(n+2*(s+1))<=M:s+=1;z*=n+2*s
   need(z<=M<z*(n+2*(s+1)),'exact largest depth');cutoffs+=1
# Transitive internal padding is retained, and its local budget is sharp.
transitive=[]
for n in range(1,10):
 T=tuple(sum(1<<j for j in range(i+1,n))for i in range(n));out=count(substitute((2,4,1),0,T));need(out==1+2**n,'padding obstruction');transitive.append(out)
print(json.dumps({'status':'PASS','rooted_contexts':len(contexts),'physical_weighted_comparisons':weighted,'strict_weighted_comparisons':strict,'weighted_exterior_triangle':[3,3],'labelled_cut_graphs':cutgraphs,'literal_labelled_cuts':cuts,'adaptive_sequences':adaptive,'active_steps':active,'integer_depth_boundaries':cutoffs,'transitive_wrapper_counts':transitive},sort_keys=True,indent=2))
