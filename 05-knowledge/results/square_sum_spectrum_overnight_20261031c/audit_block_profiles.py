#!/usr/bin/env python3
"""All fixed-parent profiles from independent endpoint and cut DPs."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def need(x,m):
 if not x:raise RuntimeError(m)

def parent(k):
 if k==-1:return [8,1,3,6]
 n=k+8;out=[0]*n;D=[8,1,3,6]
 for u in range(n):
  for v in range(n):
   if u==v:continue
   if u<4 and v<4:b=D[u]>>v&1
   elif u>=n-4 and v>=n-4:b=D[u-n+4]>>(v-n+4)&1
   else:b=int(u<v)
   if b:out[u]|=1<<v
 return out

def direct(out):
 n=len(out);full=(1<<n)-1
 end=[[0]*n for _ in range(full+1)];start=[[0]*n for _ in range(full+1)]
 for v in range(n):end[1<<v][v]=start[1<<v][v]=1
 for mask in range(1,full+1):
  for v in range(n):
   if not(mask>>v&1):continue
   rest=mask^(1<<v)
   if rest:
    end[mask][v]=sum(end[rest][u] for u in range(n) if out[u]>>v&1)
    start[mask][v]=sum(start[rest][u] for u in range(n) if out[v]>>u&1)
 H=sum(end[full]);W=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  rest=full^mask
  for u in range(n):
   if end[mask][u]:
    for v in range(n):W[u][v]+=end[mask][u]*start[rest][v]
 E=2*H+sum(map(sum,W))
 b=[[(W[u][v]+W[v][u])//2 for v in range(n)] for u in range(n)]
 need(all(2*b[u][v]==W[u][v]+W[v][u] for u in range(n) for v in range(n)),'weight integrality')
 deg=[sum(b[u][v] for v in range(n) if out[u]>>v&1) for u in range(n)]
 pw=[end[full][u]+sum(W[u]) for u in range(n)]
 qw=[start[full][v]+sum(W[u][v] for u in range(n)) for v in range(n)]
 rows=[];F=[H]*(full+1)
 for mask in range(full+1):
  if mask:
   bit=mask&-mask;v=bit.bit_length()-1;rest=mask^bit
   F[mask]=F[rest]+2*deg[v]-2*sum(b[u][v] for u in range(n) if rest>>u&1)
  s=sum(start[full][v] for v in range(n) if mask>>v&1)
  t=sum(end[full][v] for v in range(n) if not(mask>>v&1))
  p=sum(pw[v] for v in range(n) if not(mask>>v&1));q=sum(qw[v] for v in range(n) if mask>>v&1)
  rows.append([H,E,p,q,s,t,F[mask],mask])
 return rows

receipts=[]
for k in range(-1,6):
 name='D' if k==-1 else 'DD' if k==0 else f'DT{k}D'
 actual=direct(parent(k));stored=[list(map(int,line.split())) for line in(HERE/f'{name}_profiles.tsv').read_text().splitlines()]
 need(actual==stored,'independent whole profile table '+name)
 receipts.append({'parent':name,'order':len(parent(k)),'cuts':len(actual),'H':actual[0][0],'E':actual[0][1]})
print(json.dumps({'status':'passed','profiles':sum(r['cuts'] for r in receipts),'parents':receipts,
                  'scope':'Every retained profile checked by independent full endpoint DP and weighted-cut identity; no producer or profile-composition import.'},indent=2,sort_keys=True))
