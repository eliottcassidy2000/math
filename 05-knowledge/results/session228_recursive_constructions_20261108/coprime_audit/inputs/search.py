import itertools,json
from collections import Counter
from pathlib import Path

def stats(adj):
 n=len(adj); full=(1<<n)-1
 starts=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rem=mask^(1<<v)
    starts[mask][v]=sum(starts[rem][w] for w in range(n) if adj[v]>>w&1) if rem else 1
  h[mask]=sum(starts[mask])
 alpha=[sum(starts[mask][v]*h[full^mask] for mask in range(1,1<<n)) for v in range(n)]
 return h[full],sum(h[s]*h[full^s] for s in range(1<<n)),alpha

def run(n):
 bank={};counts=Counter()
 for code in range(1<<(n*(n-1)//2)):
  adj=[0]*n
  for i,(u,v) in enumerate(itertools.combinations(range(n),2)):
   if code>>i&1:adj[u]|=1<<v
   else:adj[v]|=1<<u
  H,Q,alpha=stats(adj)
  if Q%H or (Q//H)&(Q//H-1):continue
  q=Q//H;sol=[]
  for s in range(1,1<<n):
   a=sum(alpha[v] for v in range(n) if s>>v&1)
   if a%H==0 and 0<a//H<q-1:sol.append([s,a//H])
  key=(H,q,tuple(sorted(alpha)),tuple(sorted(set(m for s,m in sol))))
  counts[(H,q)]+=1
  if key not in bank:bank[key]=dict(n=n,code=code,adj=adj,H=H,Q=Q,alpha=alpha,sol=sol)
 out=list(bank.values());Path(__file__).with_name(f'bank{n}.json').write_text(json.dumps(out,indent=2));print(n,'counts',counts,'types',len(out),flush=True)
 for x in out:
  if x['H']%3:print(x,flush=True)
if __name__=='__main__':
 import sys
 for n in map(int,sys.argv[1:] or [5]):run(n)
