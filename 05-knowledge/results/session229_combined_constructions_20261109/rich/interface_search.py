from search import *
rng=random.Random(2291211)
for t in range(4000):
 n=7+(t%3);a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if rng.randrange(2):a[i]|=1<<j
   else:a[j]|=1<<i
 H,Q,aw,en,st,h=data(a)
 if math.gcd(H,Q)!=1:continue
 N=1<<n;gw=[sum(en[s][v]*h[(N-1)^s] for s in range(1,N)) for v in range(n)];fs=defaultdict(list)
 for m in range(N):
  vs=[i for i in range(n) if m>>i&1];sig=tuple(sum(w[i] for i in vs) for w in [aw,gw,st[-1],en[-1]])
  fs[sig].append(m)
 dup=[(s,ms) for s,ms in fs.items() if len(ms)>1]
 if not dup:continue
 K=[[sum(en[s][i]*st[(N-1)^s][j] for s in range(1,N-1)) for j in range(n)] for i in range(n)]
 for sig,ms in dup:
  vals=[(sum(K[i][j] for i in range(n) for j in range(n) if m>>i&1 and m>>j&1),m) for m in ms]
  if len(set(v for v,m in vals))>1:print(json.dumps(dict(trial=t,adj=a,H=H,Q=Q,signature=sig,values=vals)),flush=True)
