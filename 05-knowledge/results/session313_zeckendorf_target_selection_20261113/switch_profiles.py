from itertools import combinations, product
def mat(n,c):
 A=[0]*n
 for i,(u,v) in enumerate(combinations(range(n),2)):
  if c>>i&1:A[u]|=1<<v
  else:A[v]|=1<<u
 return A
def stats(A):
 n=len(A);N=(1<<n)-1
 s=[[0]*n for _ in range(N+1)];e=[[0]*n for _ in range(N+1)]
 for i in range(n):s[1<<i][i]=e[1<<i][i]=1
 for S in range(1,N+1):
  if S&(S-1)==0:continue
  for i in range(n):
   if S>>i&1:
    T=S^(1<<i)
    s[S][i]=sum(s[T][j] for j in range(n) if A[i]>>j&1)
    e[S][i]=sum(e[T][j] for j in range(n) if A[j]>>i&1)
 h=[1]+list(map(sum,s[1:]))
 ans={}
 for B in range(N+1):
  start=sum(s[N][i] for i in range(n) if B>>i&1)
  end=sum(e[N][i] for i in range(n) if not B>>i&1)
  sstar=sum(h[N^U]*sum(s[U][i] for i in range(n) if B>>i&1) for U in range(1,N+1))
  estar=sum(h[N^U]*sum(e[U][i] for i in range(n) if not B>>i&1) for U in range(1,N+1))
  split=sum(sum(s[U][i] for i in range(n) if B>>i&1)*sum(e[N^U][i] for i in range(n) if not B>>i&1) for U in range(1,N))
  ans[B]=((start,sstar),(end,estar,split))
 return h[N],sum(h[U]*h[N^U] for U in range(N+1)),ans
if __name__=='__main__':
 for n in range(1,7):
  for c in range(1<<(n*(n-1)//2)):
   A=mat(n,c)
   H,Q,ans=stats(A);seen={}
   for B,(x,y) in ans.items():
    if x in seen and y!=seen[x][0]:
     print(n,A,H,Q,B,seen[x],x,y);break
    seen[x]=(y,B)
   else:continue
   break
  else:continue
  break
