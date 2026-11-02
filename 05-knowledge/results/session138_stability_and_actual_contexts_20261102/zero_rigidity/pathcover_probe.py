import sys,itertools,json
from reference_paths import tournament,c3_blowup,extend,H,strong
from quotient_kernel import kernel

def induce(a,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def profile(a):
 n=len(a);N=1<<n;dp=[[None]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=[0,1]+[0]*(n-1)
 for S in range(1,N):
  for i,p in enumerate(dp[S]):
   if p is None:continue
   for j in range(n):
    if S>>j&1:continue
    D=dp[S|1<<j]
    if D[j] is None:D[j]=[0]*(n+1)
    q=D[j]
    for r in range(1,n):q[r+1]+=p[r]
    if a[i]>>j&1:
     for r in range(1,n+1):q[r]+=p[r]
 return [sum(dp[-1][i][r] for i in range(n)) for r in range(1,n+1)]

def zero_profiles(q,r):
 T=c3_blowup(q);V=[i for i in range(len(T)) if i!=3*r];B=induce(T,V);h,C=kernel(q,r);n=len(B);A=[sum(C[i][j] for j in range(n) if B[i]>>j&1) for i in range(n)];target=H(T);original=sum(1<<i for i,v in enumerate(V) if T[3*r]>>v&1);zeros=[]
 for t in itertools.product(range(4),repeat=len(q)):
  S=0
  for i,tt in enumerate(t):
   if i==r:
    for j in range(2):
     if tt>>j&1:S|=1<<V.index(3*i+j+1)
   else:
    for j in range(tt):S|=1<<V.index(3*i+j)
  val=h+sum(A[i] for i in range(n) if S>>i&1)-sum(C[i][j] for i in range(n) for j in range(i) if S>>i&1 and S>>j&1)
  if val==target and S!=original:zeros.append((S,t))
 return T,B,original,zeros

if __name__=='__main__':
 rows=[]
 for k,code in [(2,1),(3,0),(4,2)]:
  q=tournament(k,code)
  if strong(q):raise RuntimeError('nonstrong universe')
  T=c3_blowup(q);f=profile(T)
  for r in range(k):
   T,B,original,zeros=zero_profiles(q,r)
   for S,t in zeros:
    U=extend(B,S);ff=profile(U);d=[b-a for a,b in zip(f,ff)]
    row={'k':k,'code':code,'root':r,'cut':S,'profile_cut':t,'difference':d,'weakly_dominates':all(v>=0 for v in d),'strict_somewhere':any(v>0 for v in d)};rows.append(row);print(json.dumps(row,sort_keys=True),flush=True)
 print('PASS rows',len(rows),'all_weakly_dominate',all(row['weakly_dominates'] for row in rows))
