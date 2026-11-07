import itertools,json

def graphs(n):
 es=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(es)):
  a=[0]*n
  for k,(i,j) in enumerate(es):
   x,y=(i,j) if mask>>k&1 else (j,i)
   a[x]|=1<<y
  yield a

def profile(a,v):
 n=len(a); full=(1<<n)-1
 dp=[[0]*n for _ in range(1<<n)]
 H=[0]*(1<<n); H[0]=1
 for x in range(n): dp[1<<x][x]=1
 for S in range(1,1<<n):
  for x in range(n):
   if S>>x&1 and S!=(1<<x):
    dp[S][x]=sum(dp[S^(1<<x)][y] for y in range(n) if a[y]>>x&1)
  H[S]=sum(dp[S])
 # reverse endpoint DP supplies starts
 rev=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)]
 st=[[0]*n for _ in range(1<<n)]
 for x in range(n): st[1<<x][x]=1
 for S in range(1,1<<n):
  for x in range(n):
   if S>>x&1 and S!=(1<<x):
    st[S][x]=sum(st[S^(1<<x)][y] for y in range(n) if rev[y]>>x&1)
 L=sum(dp[S][v]*H[full^S] for S in range(1,1<<n))
 R=sum(st[S][v]*H[full^S] for S in range(1,1<<n))
 C=sum(H[S]*H[full^S] for S in range(1<<n))
 return dict(h=H[full],s=st[full][v],e=dp[full][v],L=L,R=R,C=C,d=H[full^(1<<v)])

def join(A,a,B,b):
 n=len(A);m=len(B)
 rows=[x|(((1<<m)-1)<<n) for x in A]+[x<<n for x in B]
 rows[a]^=1<<(n+b);rows[n+b]|=1<<a
 return rows

def main():
 inputs=[(a,v,profile(a,v)) for n in range(1,4) for a in graphs(n) for v in range(n)]
 count=0; zero=0
 for A,a,p in inputs:
  for B,b,q in inputs:
   x=profile(join(A,a,B,b),a)
   pred=dict(h=p['h']*q['h']+p['R']*q['L']-p['e']*q['s'],R=p['R']*q['C']-p['d']*q['R'],e=p['d']*q['e'],d=p['d']*q['h'])
   for k,value in pred.items():
    if x[k]!=value: raise RuntimeError((A,a,B,b,k,x[k],value))
   increment=x['h']-p['h']*q['h']
   sink=A[a]==0; source=B[b]==((1<<len(B))-1)^(1<<b)
   if (increment==0)!=(sink and source): raise RuntimeError('zero classification')
   if increment<0 or increment%2: raise RuntimeError('increment parity')
   count+=1;zero+=increment==0
 print(json.dumps(dict(status='FINITE-EXACT',universe='All labeled input tournaments of orders 1 through 3, every pair of roots',rooted_inputs=len(inputs),ordered_pairs=count,zero_increments=zero,coordinates=['h','R','e','d']),indent=2))
if __name__=='__main__':main()
