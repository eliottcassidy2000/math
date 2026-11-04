import itertools
def plus(a,b,scale=1):
 out=dict(a)
 for s,v in b.items():out[s]=out.get(s,0)+scale*v
 return {s:v for s,v in out.items()if v}
def mul(a,b):
 out={}
 for s,v in a.items():
  for t,w in b.items():
   if not(s&t):out[s|t]=out.get(s|t,0)+v*w
 return {s:v for s,v in out.items()if v}
def exterior(rows,into):
 n=len(rows);G=[[{}for _ in rows]for _ in rows]
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
    mask=sum(1<<i for i in p);g=G[p[0]][p[-1]];g[mask]=g.get(mask,0)+1
 S=[{}for _ in rows];E=[{}for _ in rows]
 for u in range(n):
  for v in range(n):S[u]=plus(S[u],G[u][v]);E[v]=plus(E[v],G[u][v])
 alpha={0:1};beta={0:1};tau={};ells=[];rhos=[]
 for v in range(n):
  if into>>v&1:alpha=plus(alpha,E[v])
  else:beta=plus(beta,S[v])
 for u in range(n):
  for v in range(n):
   if not(into>>u&1)and into>>v&1:tau=plus(tau,G[u][v])
 for i in range(n):
  ell={};rho={}
  for v in range(n):
   if into>>v&1:ell=plus(ell,G[i][v])
   else:rho=plus(rho,G[v][i])
  ells.append(ell);rhos.append(rho)
 return alpha,beta,tau,S,E,ells,rhos

def endpoints(rows,start=None):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)]
 if start is None:
  for i in range(n):dp[1<<i][i]=1
 else:dp[1<<start][start]=1
 for mask in range(1,N):
  if start is not None and not(mask>>start&1):continue
  for i,x in enumerate(dp[mask]):
   if not x:continue
   opts=rows[i]&~mask
   while opts:
    b=opts&-opts;opts-=b;dp[mask|b][b.bit_length()-1]+=x
 return dp

def capacities(rows):
 n=len(rows);full=(1<<n)-1;E=endpoints(rows);rev=[sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)];S=endpoints(rev)
 return [[sum(E[a][i]*S[full^a][j]for a in range(full+1))for j in range(n)]for i in range(n)]

