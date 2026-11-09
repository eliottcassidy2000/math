import random, math, json
from collections import defaultdict
from pathlib import Path

def data(a):
 n=len(a);N=1<<n;full=N-1;h=[1]+[0]*(N-1);en=[[0]*n for _ in range(N)];st=[[0]*n for _ in range(N)]
 for s in range(1,N):
  vs=[v for v in range(n) if s>>v&1]
  for v in vs:
   r=s^(1<<v)
   en[s][v]=sum(en[r][u] for u in vs if a[u]>>v&1) if r else 1
   st[s][v]=sum(st[r][u] for u in vs if a[v]>>u&1) if r else 1
  h[s]=sum(en[s])
 H=h[-1];Q=sum(h[s]*h[full^s] for s in range(N))
 aw=[sum(st[s][v]*h[full^s] for s in range(1,N)) for v in range(n)]
 return H,Q,aw,en,st,h

def inspect(a):
 H,Q,aw,en,st,h=data(a)
 if math.gcd(H,Q)!=1:return None
 n=len(a);N=1<<n;fib=defaultdict(list)
 for m in range(N):fib[sum(aw[i] for i in range(n) if m>>i&1)].append(m)
 mx=max(map(len,fib.values()))
 if mx<3:return None
 K=[[sum(en[s][i]*st[(N-1)^s][j] for s in range(1,N-1)) for j in range(n)] for i in range(n)]
 gw=[sum(en[s][v]*h[(N-1)^s] for s in range(1,N)) for v in range(n)]
 best=None
 for al,ms in fib.items():
  if len(ms)<3:continue
  states={}
  for m in ms:
   vs=[i for i in range(n) if m>>i&1]
   state=(sum(gw[i] for i in vs),sum(st[-1][i] for i in vs),sum(en[-1][i] for i in vs),sum(K[i][j] for i in vs for j in vs))
   states.setdefault(state,m)
  if best is None or len(states)>best['distinct']:
   best=dict(adj=a,H=H,Q=Q,alpha_weights=aw,alpha=al,masks=ms,distinct=len(states),states=[list(s)+[m] for s,m in states.items()],ratio=Q/H)
 return best

def universe():
 rng=random.Random(2291209)
 for n in range(7,13):
  for trial in range(180):
   if trial<120:
    nn=n+(n%2==0);a=[sum(1<<j for j in range(n) if 0<(j-i)%nn<=(nn-1)//2) for i in range(n)]
    for i,j in rng.sample([(i,j) for i in range(n) for j in range(i+1,n)],1+trial%4):a[i]^=1<<j;a[j]^=1<<i
   else:
    a=[0]*n
    for i in range(n):
     for j in range(i+1,n):
      if rng.randrange(2):a[i]|=1<<j
      else:a[j]|=1<<i
   yield n,trial,a
if __name__=='__main__':
 best=0;hits=[]
 for n,t,a in universe():
  b=inspect(a)
  if b:
   b.update(n=n,trial=t);hits.append(b)
   score=b['distinct']/b['ratio']
   if score>best:
    best=score;print(json.dumps(b),flush=True)
    Path(__file__).with_name('best.json').write_text(json.dumps(b,indent=2))
 Path(__file__).with_name('hits.json').write_text(json.dumps(hits,indent=2))
