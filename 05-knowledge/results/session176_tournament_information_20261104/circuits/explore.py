from itertools import combinations

def tables(a):
 n=len(a); N=1<<n
 E=[[0]*n for _ in range(N)];B=[[0]*n for _ in range(N)]
 for u in range(n): E[1<<u][u]=B[1<<u][u]=1
 for S in range(1,N):
  for u in range(n):
   if not S>>u&1 or S==1<<u:continue
   R=S^(1<<u)
   E[S][u]=sum(E[R][v] for v in range(n) if a[v]>>u&1)
   B[S][u]=sum(B[R][v] for v in range(n) if a[u]>>v&1)
 return E,B

def kernel(a):
 n=len(a);E,B=tables(a);F=(1<<n)-1;K=[[0]*n for _ in a]
 for S in range(1,F):
  for u in range(n):
   for v in range(n): K[u][v]+=E[S][u]*B[F^S][v]
 return sum(E[F]),K

if __name__=='__main__':
 for n in range(2,6):
  pairs=list(combinations(range(n),2)); low=10**9; witness=None
  for bits in range(1<<len(pairs)):
   a=[0]*n
   for k,(u,v) in enumerate(pairs):
    if bits>>k&1:a[u]|=1<<v
    else:a[v]|=1<<u
   H,K=kernel(a)
   m=min(K[u][v]+K[v][u] for u,v in pairs)
   if m<low:low=m;witness=a
  print(n,low,witness)
