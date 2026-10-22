from pathlib import Path
from itertools import combinations
import numpy as np,json

def hp(n,flips):
 out=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 for i,j in flips:out[i]^=1<<j;out[j]^=1<<i
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   c=dp[mask][v]
   if c:
    rem=out[v]&~mask
    while rem:
     b=rem&-rem;rem-=b;dp[mask|b][b.bit_length()-1]+=c
 return sum(dp[-1])
rows=[]
for n in range(3,11):
 edges=[(i,j) for i in range(n) for j in range(i+2,n)];m=len(edges)
 one=[hp(n,[e]) for e in edges]; H=np.zeros((m,m),dtype='int64')
 for i in range(m):
  for j in range(i):H[i,j]=H[j,i]=hp(n,[edges[i],edges[j]])-one[i]-one[j]+1
 ev=np.linalg.eigvalsh(H.astype(float)); pos=sum(ev>1e-7);neg=sum(ev< -1e-7);zero=m-pos-neg
 rows.append({'n':n,'m':m,'signature_float':[int(pos),int(zero),int(neg)],'matrix':H.tolist()})
 print(n,pos,zero,neg,flush=True)
json.dump(rows,open(str(Path(__file__).with_name('transitive_hessian_probe.json')),'w'),indent=2)
