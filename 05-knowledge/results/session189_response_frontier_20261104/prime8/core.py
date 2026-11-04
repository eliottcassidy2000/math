"""Standard-library exact tournament and repeated-run primitives."""
def converse(T):return tuple(((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T))
def endpoint(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,1<<n):
  for i in range(n):
   if m>>i&1:
    v=dp[m][i];todo=T[i]&~m
    while todo:
     b=todo&-todo;todo-=b;j=b.bit_length()-1;dp[m|b][j]+=v
 return dp[-1]
def repeated_endpoint(T,j):
 n=len(T);full=(1<<n)-1;dp={}
 for i in range(n):dp[(1<<i,0,i)]=1
 for length in range(1,n+1):
  nd={}
  for (m,b,i),v in dp.items():
   todo=T[i]&~m
   while todo:
    z=todo&-todo;todo-=z;k=z.bit_length()-1;key=(m|z,b,k);nd[key]=nd.get(key,0)+v
   if m>>j&1 and not b and T[i]>>j&1:
    key=(m,1,j);nd[key]=nd.get(key,0)+v
  dp=nd
 return [dp.get((full,1,i),0) for i in range(n)]
