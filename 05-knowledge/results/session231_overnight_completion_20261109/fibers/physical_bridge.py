from pathlib import Path
import json
root=Path(__file__).parent
scc={tuple(v[:2]) for c in json.loads((root/'NONZERO_SCC.json').read_text()) for v in c}
out=[]
for n in range(1,10):
 selected=[];h95=0;hits=[]
 for line in (root/f'classes9_{n}.txt').read_text().splitlines():
  if not line or line[0]=='#':continue
  words=line.split();code=int(words[0]);H=int(words[1])
  if H!=95:continue
  h95+=1;adj=[0]*n;bit=0
  for i in range(n):
   for j in range(i+1,n):
    if code>>bit&1:adj[i]|=1<<j
    else:adj[j]|=1<<i
    bit+=1
  N=1<<n;F=N-1;hp=[1]+[0]*F;en=[[0]*n for _ in range(N)];st=[[0]*n for _ in range(N)]
  for s in range(1,N):
   for v in range(n):
    if not(s>>v&1):continue
    r=s^(1<<v)
    en[s][v]=sum(en[r][u] for u in range(n) if adj[u]>>v&1) if r else 1
    st[s][v]=sum(st[r][u] for u in range(n) if adj[v]>>u&1) if r else 1
   hp[s]=sum(en[s])
  Q=sum(hp[s]*hp[F^s] for s in range(N))
  if Q!=1276:continue
  aw=[sum(st[s][v]*hp[F^s] for s in range(N)) for v in range(n)]
  gw=[sum(en[s][v]*hp[F^s] for s in range(N)) for v in range(n)]
  selected.append(code)
  for m in range(N):
   A=sum(aw[v] for v in range(n) if m>>v&1);G=sum(gw[v] for v in range(n) if m>>v&1)
   if (A,G) in scc:hits.append(dict(code=code,adj=adj,mask=m,A=A,G=G))
 out.append(dict(n=n,H95_classes=h95,H95_Q1276_classes=selected,hits=hits))
print(json.dumps(out,indent=2))
