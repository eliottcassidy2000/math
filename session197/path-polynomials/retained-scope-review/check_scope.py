import json,hashlib
from pathlib import Path

def T(n,mask):
 rows=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   a,b=(i,j)if mask>>k&1 else(j,i);rows[a]|=1<<b;k+=1
 return rows

def direct(A,a,m,j):
 n=len(A);out=A+[0]*m
 for u in range(n):
  for v in range(n,n+m):
   if u==a and v==n+j:out[v]|=1<<u
   else:out[u]|=1<<v
 for u in range(n,n+m):
  for v in range(u+1,n+m):out[u]|=1<<v
 return out

def sequence(A,a,m,j):
 out=list(A)
 for i in range(m):
  n=len(out);out.append(0)
  for old in range(n):
   if i==j and old==a:out[n]|=1<<old
   else:out[old]|=1<<n
 return out
checks=0
for n in range(1,5):
 for mask in range(1<<(n*(n-1)//2)):
  A=T(n,mask)
  for a in range(n):
   for m in range(1,8):
    for j in range(m):
     if direct(A,a,m,j)!=sequence(A,a,m,j):raise RuntimeError('labelled stationary transitive decomposition')
     checks+=1
# Independent general coefficient vectors in original h,e,R,d coordinates:
# scalar(K)=scalar(A); scalar(KK)=scalar(KA); scalar(DK)=scalar(AA).
K=[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]]
D=[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]]
A=[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]
def mul(B,C):return [[sum(B[i][k]*C[k][j]for k in range(4))for j in range(4)]for i in range(4)]
if K[0]!=A[0]or mul(K,K)[0]!=mul(A,K)[0]or mul(K,D)[0]!=mul(A,A)[0]:raise RuntimeError('general continuation row identity')
P=Path(__file__).resolve().parent
print(json.dumps(dict(status='PASS',labelled_transitive_identifications=checks,general_row_identities=3,pinned_primary_sha256=hashlib.sha256((P/'S195_COMBINED_GROWTH.md').read_bytes()).hexdigest()),sort_keys=True,indent=2))
