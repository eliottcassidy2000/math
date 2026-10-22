import json,hashlib,itertools
from collections import defaultdict,Counter
from pathlib import Path
root=Path(__file__).resolve().parent; source=root.parent/'session26-last-wiring'
data=json.loads((source/'catalogue.json').read_text()); masks=defaultdict(dict)
for line in (root/'selected.txt').read_text().splitlines():
 n,c,a=map(int,line.split());masks[n][c]=a

def matrix(n,c):
 a=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=1 if j==i+1 else (c>>k)&1
   if j!=i+1:k+=1
   a[j][i]=1-a[i][j]
 return a

def labels(a,aps):
 words=[]
 for p in range(len(a)):
  if not aps>>p&1:continue
  left=set(range(len(a)))-{p};order=[]
  while left:
   v=next((v for v in left if all(v==u or a[v][u] for u in left)),None)
   assert v is not None
   order.append(v);left.remove(v)
  words.append(''.join(str(a[p][v]) for v in order))
 assert words
 return min(words),min(words+[''.join('1' if x=='0' else '0' for x in w[::-1]) for w in words])

def paths(a):
 n=len(a);dp={(1<<v,v):1 for v in range(n)}
 for s in range(1,1<<n):
  for v in range(n):
   val=dp.get((s,v),0)
   if val:
    for u in range(n):
     if not s>>u&1 and a[v][u]:dp[s|1<<u,u]=dp.get((s|1<<u,u),0)+val
 return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

reports=[]
for d in data:
 n=d['n'];pairs=[(i,j) for i in range(n) for j in range(i+2,n)];idx={p:k for k,p in enumerate(pairs)};full=(1<<len(pairs))-1
 def r(c):return sum(((c>>idx[n-1-j,n-1-i])&1)<<k for k,(i,j) in enumerate(pairs))
 labs={c:labels(matrix(n,c),aps) for c,aps in masks[n].items()}
 ordinary=Counter(o for o,m in labs.values());merged=Counter(m for o,m in labs.values());internal={c for c in labs if c^full in labs}
 edges={}
 for c in internal:
  orb=sorted({c,c^full,r(c),r(c)^full});key=min(orb)
  edges[key]={'id':key,'endpoints':sorted([labs[c][1],labs[c^full][1]]),'color':'blue' if r(c)==c else 'black','tilings':orb}
 assert sorted(edges.values(),key=lambda x:(x['endpoints'],x['id']))==d['edges']
 touched={w for e in edges.values() for w in e['endpoints']}
 stats={'ordinary_classes':len(ordinary),'merged_nodes':len(merged),'selected_tilings':len(labs),'blue_incidences':sum(r(c)==c for c in labs),'internal_tilings':len(internal),'internal_edges':len(edges),'internal_blue_edges':sum(e['color']=='blue' for e in edges.values()),'internal_black_edges':sum(e['color']=='black' for e in edges.values()),'internal_loops':sum(e['endpoints'][0]==e['endpoints'][1] for e in edges.values()),'isolated_nodes':len(set(merged)-touched)}
 for k,v in stats.items():assert d[k]==v,(n,k,v,d[k])
 assert d['nodes']==[{'word':w,'ordinary_tiling_mass':ordinary[w],'merged_tiling_mass':merged[w]} for w in sorted(merged)]
 report={'n':n,**stats}
 if n==8:
  assert all(masks[n][c] in [1<<2,1<<5] and masks[n][c^full] != masks[n][c] for c in internal)
  assert Counter(masks[n][c] for c in internal)=={1<<2:12,1<<5:12}
  report['internal_masks']=[{'mask':c,'apex_position':masks[n][c].bit_length(),'ordinary_word':labs[c][0],'merged_word':labs[c][1]} for c in sorted(internal)]
  report['endpoint_H']={}
  for w in sorted(touched):
   c=next(c for c in internal if labs[c][0]==w);h=paths(matrix(n,c));assert h==ordinary[w];report['endpoint_H'][w]=h
  borbs={min(c,c^full,r(c),r(c)^full):('blue' if r(c)==c else 'black') for c in labs if c^full not in labs}
  report['boundary']=dict(Counter(borbs.values()));assert report['boundary']=={'blue':18,'black':1928}
 reports.append(report)
(root/'audit-results.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
