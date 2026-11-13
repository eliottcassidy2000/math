from itertools import permutations,product,combinations
from collections import Counter,defaultdict
from math import factorial

def make(m,eps,b):
 E=list(combinations(range(m),2));A=[[0]*m for _ in range(m)];G=[[0]*m for _ in range(m)]
 for k,(i,j) in enumerate(E):
  A[i][j]=b>>k&1;A[j][i]=1-A[i][j]
  G[i][j]=G[j][i]=b>>(len(E)+k)&1
 return A,G,[(b>>(2*len(E)+i))&1 for i in range(m)] if eps else []
def graph(D,eps):
 A,G,C=D;m=len(A);n=2*m+eps;T=[[0]*n for _ in range(n)]
 def add(x,y,p):T[x][y]=p;T[y][x]=1-p
 for i in range(m):
  add(2*i,2*i+1,1)
  if eps:add(2*i,2*m,C[i]);add(2*m,2*i+1,C[i])
 for i,j in combinations(range(m),2):
  add(2*i,2*j,A[i][j]);add(2*i+1,2*j+1,A[j][i])
  add(2*i,2*j+1,G[i][j]);add(2*j,2*i+1,G[i][j])
 return T
def flip(D,eps,L):
 A,G,C=([r[:] for r in D[0]],[r[:] for r in D[1]],D[2][:])
 for u,v in zip(L,L[1:]):
  i,j=u//2,v//2;K=G if (u^v)&1 else A
  K[i][j]^=1;K[j][i]^=1
 i=L[-1]//2
 if eps:C[i]^=1
 else:
  for j in range(len(A)):
   if j==i:continue
   x,y=A[i][j],G[i][j]
   A[i][j]=1-y;A[j][i]=y;G[i][j]=G[j][i]=1-x
 return A,G,C
def dp(T,m,eps):
 from functools import cache
 @cache
 def end(S,v):
  i=v//2
  if S==1<<i:return (int(v%2==0),int(v%2==1))
  C=[end(S^(1<<i),w) for w in range(2*m) if w//2!=i and S>>(w//2)&1 and T[w][v]]
  return (sum(z[0] for z in C),sum(z[1] for z in C))
 acc=[end((1<<m)-1,v) for v in range(2*m) if (T[v][2*m] if eps else v%2==0)]
 return (sum(z[0]+z[1] for z in acc),sum(z[1] for z in acc))
def counts(m,eps):
 E=list(combinations(range(m),2))
 def walk(p):
  seen=set();vs=[]
  for i in range(m):
   if i not in seen:
    j=i;k=0
    while j not in seen:seen.add(j);k+=1;j=p[j]
    vs.append(k)
  if any(k%2==0 for k in vs):return {}
  sizes=[];S=set()
  for i,j in E:
   if (i,j) in S:continue
   k=0
   while (i,j) not in S:S.add((i,j));k+=1;i,j=sorted((p[i],p[j]))
   sizes.append(k)
  d={(0,0):1<<len(sizes)}
  for x,y in [(k,0) for k in sizes]+([(0,k) for k in vs] if eps else []):
   new=Counter(d)
   for (a,b),v in d.items():new[a+x,b+y]+=v
   d=new
  return d
 agg=Counter()
 for p in permutations(range(m)):agg.update(walk(p))
 assert all(x%factorial(m)==0 for x in agg.values())
 return {k:v//factorial(m) for k,v in agg.items()}
def can(D):
 A,G,C=D;m=len(A);E=list(combinations(range(m),2))
 return min(tuple(A[p[i]][p[j]] for i,j in E)+tuple(G[p[i]][p[j]] for i,j in E)+tuple(C[i] for i in p if C) for p in permutations(range(m)))
for n in range(2,8):
 m,eps=divmod(n,2);E=list(combinations(range(m),2))
 classes=defaultdict(set);ports=0
 for b in range(1<<(2*len(E)+eps*m)):
  D=make(m,eps,b);T=graph(D,eps);A,G,C=D
  k=sum(G[i][j] for i,j in E),sum(C);classes[k].add(can(D))
  ts=closed=0
  for p in permutations(range(m)):
   for sig in product(range(2),repeat=m):
    L=[2*i+s for i,s in zip(p,sig)]
    P=L+([2*m] if eps else [])+[u^1 for u in L[::-1]]
    if not all(T[x][y] for x,y in zip(P,P[1:])):continue
    ts+=1;closed+=T[P[-1]][P[0]]
    U=[r[:] for r in T]
    for x,y in zip(P,P[1:]):U[x][y],U[y][x]=0,1
    if not eps:
     i=L[-1]//2
     perm=list(range(n));perm[2*i],perm[2*i+1]=perm[2*i+1],perm[2*i]
     U=[[U[x][y] for y in perm] for x in perm]
    assert U==graph(flip(D,eps,L),eps)
  assert (ts,closed)==dp(T,m,eps);ports+=ts
 assert {k:len(v) for k,v in classes.items()}==counts(m,eps)
 print(n,sum(map(len,classes.values())),ports)

def add(*P):
 z=Counter()
 for Q in P:z.update(Q)
 return z
def mul(P,Q):
 z=Counter()
 for a,s in P.items():
  for b,t in Q.items():z[tuple(x+y for x,y in zip(a,b))]+=s*t
 return z
def mon(*a):return {a:1}
def times(n,P):
 Q=mon(0,0,0,0,0)
 for _ in range(n):Q=mul(Q,P)
 return Q
one=mon(0,0,0,0,0);X=mon(1,0,0,0,0);Y=mon(0,1,0,0,0);V=mon(0,0,1,0,0);W=mon(0,0,0,1,0);U=mon(0,0,0,0,1)
twin=add(one,mul(X,Y))
M=[[twin,add(X,X)],[add(Y,Y),twin]]
for n in range(3,9):
 m,eps=divmod(n,2);E=list(combinations(range(m),2));expected=Counter()
 for b in range(1<<(2*len(E)+eps*m)):
  D=make(m,eps,b);T=graph(D,eps)
  for sig in product(range(2),repeat=m):
   L=[2*i+s for i,s in enumerate(sig)]
   if not all(T[x][y] for x,y in zip(L,L[1:])):continue
   if not (T[L[-1]][2*m] if eps else sig[-1]==0):continue
   U2=flip(D,eps,L)
   expected[sum(D[1][i][j] for i,j in E),sum(U2[1][i][j] for i,j in E),sum(D[2]),sum(U2[2]),sig[0]]+=1
 K=[one,U]
 for _ in range(m-2+eps):K=[add(*(mul(K[i],M[i][j]) for i in range(2))) for j in range(2)]
 if eps:
  q=add(mul(K[0],V),mul(K[1],W))
  q=mul(q,times(m-1,add(one,mul(V,W))))
  q=mul(q,times(len(E)-(m-1),add(twin,twin)))
 else:
  q=add(mul(K[0],add(one,X)),mul(K[1],add(one,Y)))
  q=mul(q,times(m-2,mul(add(one,X),add(one,Y))))
  q=mul(q,times((m-1)*(m-2)//2-(m-2),add(twin,twin)))
 assert q==expected
 print('joint blue',n,len(q),sum(q.values()))
def strong(T):
 n=len(T)
 for i in range(n):
  S={i}
  for _ in range(n):S|={j for j in range(n) if any(T[k][j] for k in S)}
  if len(S)!=n:return False
 return True
A=[1,1,1,2]  # A000568 at 0..3, immediate by enumeration
S=[]
for m in range(4):
 reps={}
 for b in range(1<<(m*m)):
  D=make(m,1,b)
  if strong(graph(D,1)):reps[can(D)]=sum(D[2])
 S.append(sum((-1)**v for v in reps.values()))
assert [sum(S[j]*(-1)**(i-j)*A[i-j] for j in range(i+1)) for i in range(4)]==[1,0,0,0]
print('strong signed reciprocal coefficients',S)
