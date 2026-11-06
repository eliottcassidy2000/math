import sys,itertools,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import verify as v

def scalar(a,r):
 n=len(a);N=1<<n
 dp=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)]
 for j in range(n):dp[1<<j][j]=1
 sp[1<<r][r]=1
 for m in range(1,N):
  for j in range(n):
   if not dp[m][j]:continue
   q=a[j]&~m
   while q:
    b=q&-q;q-=b;k=b.bit_length()-1
    dp[m|b][k]+=dp[m][j];sp[m|b][k]+=sp[m][j]
 H=[sum(d) for d in dp];H[0]=1
 h=H[-1];s=sum(sp[-1]);L=sum(dp[m][r]*H[(N-1)^m] for m in range(1,N))
 return h+L-s,s
checks=0
for n in range(4):
 for w in itertools.product(range(2),repeat=n):
  rows=list(v.seed);ports=[0,1,2];x=v.x
  for bit in w:
   mask,drop,rem=v.ops[bit];k=len(rows)
   rows=v.append(rows,ports,mask);ports=[(ports+[k])[a] for a in rem];x=v.M[bit]*x
  root=ports[2]
  for _ in range(2):rows=v.append(rows,ports,3)
  for i in range(len(rows)):
   if rows[i]>>i&1:raise RuntimeError('loop')
   for j in range(i):
    if ((rows[i]>>j)&1)+((rows[j]>>i)&1)!=1:raise RuntimeError('orientation')
  if scalar(rows,root)!=tuple(v.C*x):raise RuntimeError(w)
  checks+=1
print(json.dumps({'status':'PASS','independent_subset_scalar_checks':checks,'max_word_length':3}))
