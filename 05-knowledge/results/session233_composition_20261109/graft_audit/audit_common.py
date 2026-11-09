import itertools


def need(x,m):
 if not x: raise RuntimeError(m)

def literal(adj,mask):
 n=len(adj); hs=[0]*(1<<n); ss=hs.copy(); es=hs.copy(); hs[0]=1
 # Independent ordered tuples of vertices, not endpoint recurrence.
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(adj[a]>>b&1 for a,b in zip(p,p[1:])):
    s=sum(1<<v for v in p);hs[s]+=1;ss[s]+=mask>>p[0]&1;es[s]+=mask>>p[-1]&1
 f=(1<<n)-1
 return (hs[f],sum(hs[s]*hs[f^s] for s in range(f+1)),sum(hs[s]*ss[f^s] for s in range(f+1)),sum(es[s]*hs[f^s] for s in range(f+1)),sum(es[s]*ss[f^s] for s in range(f+1)),ss[f],es[f])

def product(x,y):
 h,q,a,g,t,s,e=x;H,Q,A,G,T,S,E=y
 return (h*H,q*Q,a*Q+h*A,g*H+q*G,t*H+h*T+a*G+e*S,s*H,h*E)

def singleton(b):return (1,2,b,b,0,b,b)

def fold(xs):
 it=iter(xs);v=tuple(next(it))
 for x in it:v=product(v,x)
 return v

def word(w):return fold(singleton(int(b)) for b in w)

def ordinal(a,b):return [v|(((1<<len(b))-1)<<len(a)) for v in a]+[v<<len(a) for v in b]

def close(a,mask):
 n=len(a);return [v|(1<<n)|((not(mask>>i&1))<<(n+1)) for i,v in enumerate(a)]+[1<<(n+1),mask]

def count(adj):
 n=len(adj);N=1<<n;dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   z=dp[s][v]
   if z:
    nxt=adj[v]&(~s)&(N-1)
    while nxt:
     b=nxt&-nxt;nxt-=b;dp[s|b][b.bit_length()-1]+=z
 return sum(dp[-1])
