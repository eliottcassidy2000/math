"""Explore one-birth head lifts of the actual 35-grid s1 module family."""
from collections import defaultdict
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
B=(8,5,9,322,15,31,55,127,247)

def inflate(T,v,r):
 labels=list(range(v))+[v]*r+list(range(v+1,len(T)));R=[0]*len(labels)
 for i in range(len(labels)):
  for j in range(i+1,len(labels)):
   a,b=(i,j) if labels[i]==labels[j] or T[labels[i]]>>labels[j]&1 else(j,i)
   R[a]|=1<<b
 return R

def ends(T):
 n=len(T);M=(1<<n)-1;dp=[[0]*n for _ in range(M+1)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,M+1):
  for i,x in enumerate(dp[mask]):
   if not x:continue
   todo=T[i]&~mask&M
   while todo:
    b=todo&-todo;todo-=b;dp[mask|b][b.bit_length()-1]+=x
 return dp

def converse(T):
 return [sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T))]

def birth_counts(T):
 n=len(T);M=(1<<n)-1;e=ends(T);s=ends(converse(T));H=sum(e[M]);K=[[0]*n for _ in range(n)]
 for mask in range(1,M):
  left=[(i,x) for i,x in enumerate(e[mask]) if x];right=[(j,y) for j,y in enumerate(s[M^mask]) if y]
  for i,x in left:
   for j,y in right:K[i][j]+=x*y
 L=[s[M][i]-e[M][i]+sum(K[j][i] for j in range(n)) for i in range(n)]
 Q=[[K[i][j]+K[j][i] for j in range(n)] for i in range(n)]
 counts=[H]*(M+1)
 for mask in range(1,M+1):
  b=mask&-mask;i=b.bit_length()-1;rest=mask^b;cost=L[i];x=rest
  while x:
   c=x&-x;x-=c;cost-=Q[i][c.bit_length()-1]
  counts[mask]=counts[rest]+cost
 return H,counts

def unpack(mask,r):
 inds=list(range(4))+list(range(4+r,8+r))
 E=sum((mask>>j&1)<<i for i,j in enumerate(inds));W=(mask>>4)&((1<<r)-1)
 return E,W

def main():
 data=[]
 for r in range(1,8):
  T=inflate(B,4,r);h,counts=birth_counts(T);assert h==3+70*2**r
  hit=defaultdict(list)
  for mask,count in enumerate(counts):
   if count==2*h-1:
    E,W=unpack(mask,r);hit[E].append(W)
  row={'r':r,'source_H':h,'target':2*h-1,'hits':dict(hit)};data.append(row)
  print(json.dumps({'r':r,'target':2*h-1,'hits':sum(map(len,hit.values())),'externalcuts':len(hit),'data':dict(hit)}),flush=True)
 (HERE/'head_lift_receipt.json').write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':main()
