from itertools import combinations
from collections import Counter,defaultdict

def adj_from(n,m):
 adj=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if m>>k&1:adj[j]|=1<<i
  else:adj[i]|=1<<j
 return adj

def dpstats(adj):
 n=len(adj);N=1<<n;dp=[[0]*n for i in range(N)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,N):
  for v,c in enumerate(dp[mask]):
   if c:
    todo=adj[v]&~mask
    while todo:
     bit=todo&-todo;todo-=bit;w=bit.bit_length()-1;dp[mask|bit][w]+=c
 H=sum(dp[-1]);dp=[[0]*n for i in range(N)];dp[1][0]=1
 for mask in range(1,N):
  for v,c in enumerate(dp[mask]):
   if c:
    todo=adj[v]&~mask
    while todo:
     bit=todo&-todo;todo-=bit;w=bit.bit_length()-1;dp[mask|bit][w]+=c
 C=sum(dp[-1][v] for v in range(n) if adj[v]&1)
 return H,C

if __name__=='__main__':
 for n in range(3,7):
  d=defaultdict(Counter);examples={}
  for m in range(1<<(n*(n-1)//2)):
   h,c=dpstats(adj_from(n,m));d[c][h]+=1;examples.setdefault((h,c),m)
  print(n,{c:sorted(d[c]) for c in [0,1,2]},flush=True)
  for h in sorted(d[1]):print('UNIQUE',n,h,examples[h,1],flush=True)
