"""Exact physical pair-entry and rooted one-seam constructions."""
from itertools import permutations

def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def tournament(r):
 n=len(r)
 return all(0<=x<1<<n and not x>>i&1 for i,x in enumerate(r)) and all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))

def inventories(r):
 require(tournament(r),'invalid tournament');n=len(r);N=1<<n
 end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
 for v in range(n):end[1<<v][v]=start[1<<v][v]=1
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    old=M^(1<<v)
    if old:
     end[M][v]=sum(end[old][u] for u in range(n) if old>>u&1 and r[u]>>v&1)
     start[M][v]=sum(start[old][u] for u in range(n) if old>>u&1 and r[v]>>u&1)
 H=[sum(x) for x in end];H[0]=1
 R=[sum(H[M]*start[(N-1)^M][a] for M in range(N)) for a in range(n)]
 L=[sum(end[M][a]*H[(N-1)^M] for M in range(N)) for a in range(n)]
 return dict(H=H[-1],e=end[-1],s=start[-1],R=R,L=L,allH=H,end=end,start=start)

def root_profile(r,a):
 p=inventories(r);return p['H'],p['e'][a],p['R'][a]

def cycles(r):
 if len(r)<3:return 0
 n=len(r);N=1<<n;D=[[0]*n for _ in range(N)];D[1][0]=1
 for M in range(1,N):
  for u in range(n):
   if D[M][u]:
    for v in range(n):
     if not M>>v&1 and r[u]>>v&1:D[M|1<<v][v]+=D[M][u]
 return sum(D[-1][u] for u in range(n) if r[u]&1)

def literal_h(r):
 return sum(all(r[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(r))))

def gadget(L):
 require(L>=1,'positive gadget index');bits=[] if L==1 else [1-int(x) for x in bin(L-1)[2:][::-1]]
 n=1+len(bits);r=[0]*n
 for i in range(1,n):
  r[0 if bits[i-1] else i]|=1<<(i if bits[i-1] else 0)
  for j in range(i+1,n):r[i]|=1<<j
 return r

def constants(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=constants(L//2);return g,Q,2*B
 ell=(L+1)//2;g,Q,B=constants(ell);return g+ell-ell%2,Q+B,2*B

def glue(T,entries,G):
 n=len(T);r=list(T)+[x<<n for x in G]
 for u in range(n):
  for v in range(len(G)):
   if v==0 and u in entries:r[n]|=1<<u
   else:r[u]|=1<<(n+v)
 return r,n

def pair_response(T,a,b):
 p=inventories(T);N=1<<len(T)
 K=sum(p['end'][M][a]*p['start'][(N-1)^M][b]+p['end'][M][b]*p['start'][(N-1)^M][a] for M in range(N))
 return p['H'],p['e'][a]+p['e'][b],p['R'][a]+p['R'][b],K

def pair_transfer(response,L):
 h,E,R,K=response;g,Q,B=constants(L)
 return g*h-(L%2)*E+L*R-K,h-int(L==1)*E,Q*h+B*R

def single_transfer(state,L):
 h,e,R=state;g,Q,B=constants(L)
 return g*h-(L%2)*e+L*R,h-int(L==1)*e,Q*h+B*R

def birth(T,ports,mask,drop):
 n=len(T);r=list(T)+[0]
 for u in range(n):
  if u in ports and mask>>ports.index(u)&1:r[n]|=1<<u
  else:r[u]|=1<<n
 p=list(ports)+[n];p.pop(drop)
 return r,p

def adjacency_map(A,B,f):
 return len(A)==len(B)==len(f) and sorted(f)==list(range(len(A))) and all(((A[u]>>v)&1)==((B[f[u]]>>f[v])&1) for u in range(len(A)) for v in range(len(A)))
