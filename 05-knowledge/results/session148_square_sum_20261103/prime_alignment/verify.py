import json, math, csv, functools, io
from pathlib import Path
D=Path(__file__).resolve().parent
paths=json.loads((D/'paths.json').read_text()); primes=[n for n in range(2,100) if all(n%d for d in range(2,math.isqrt(n)+1))]
orders=[15,16,17,23,*range(25,46)]
ogg={2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}
groups={'A':{2,3,5,7,11,13},'B':{17,19,23,31,47},'C':{29,41,59,71},'D':{37,43,53,61,71},'E':{73,79,83,89,97}}
def edge(x,y):return x!=y and math.isqrt(x+y)**2==x+y
def exists(n):
 a=[sum(1<<j for j in range(n) if edge(i+1,j+1)) for i in range(n)]
 @functools.lru_cache(None)
 def go(v,unseen):
  if not unseen:return True
  t=a[v]&unseen
  while t:
   b=t&-t;t-=b
   if go(b.bit_length()-1,unseen^b):return True
  return False
 full=(1<<n)-1
 return any(go(v,full^(1<<v)) for v in range(n))
def demand(x, message):
 if not x: raise ValueError(message)
def valid_path(p,n):
 return sorted(p)==list(range(1,n+1)) and all(edge(a,b) for a,b in zip(p,p[1:]))
for n in orders:
 p=paths[str(n)]
 demand(valid_path(p,n),n)
demand([n for n in range(2,25) if exists(n)]==[15,16,17,23], 'early exhaustive mismatch')
demand(len(primes)==len(orders)==25, 'cardinality mismatch')
def graph_stats(n):
 from collections import deque
 a=[[j for j in range(1,n+1) if edge(i,j)] for i in range(1,n+1)]
 seen=set(); components=0
 for start in range(1,n+1):
  if start in seen: continue
  components+=1; seen.add(start); q=[start]
  for v in q:
   for u in a[v-1]:
    if u not in seen: seen.add(u);q.append(u)
 best=None
 for start in range(1,n+1):
  dist={start:0}; par={start:0}; q=deque([start])
  while q:
   v=q.popleft()
   for u in a[v-1]:
    if u not in dist: dist[u]=dist[v]+1;par[u]=v;q.append(u)
    elif par[v]!=u:
     cycle=dist[v]+dist[u]+1
     best=cycle if best is None else min(best,cycle)
 return {'components':components,'edges':sum(map(len,a))//2,'girth':best}
G=[graph_stats(n) for n in range(1,98)]
with io.StringIO(newline='') as f:
 wr=csv.writer(f);wr.writerow(['i','prime','supersingular','group_as_supplied','rank_graph_order','prime_graph_order','admissible_graph_order','prime_graph_components','prime_graph_girth','admissible_graph_edges','admissible_graph_girth'])
 for i,(p,n) in enumerate(zip(primes,orders),1):
  gp=G[p-1]; ga=G[n-1];wr.writerow([i,p,int(p in ogg),'/'.join(k for k,v in groups.items() if p in v) or 'MISSING',i,p,n,gp['components'],gp['girth'] or '',ga['edges'],ga['girth'] or ''])
 demand(f.getvalue().encode('utf-8') == (D/'alignment.csv').read_bytes(), 'alignment.csv differs from regenerated table')
demand(sum(len(v) for v in groups.values())==25 and set.union(*groups.values())==set(primes)-{67}, 'group mismatch')
demand(not edge(1,2), 'nonedge negative control')
demand(not valid_path(paths['15'][:-1],15), 'missing-vertex control')
corrupt=paths['15'].copy(); corrupt[0],corrupt[1]=corrupt[1],corrupt[0]
demand(not valid_path(corrupt,15), 'bad-edge control')
print('PASS: 25 primes, 25 admissible-order witnesses, and exhaustive n=2..24; literal n=1 is separately admissible. stored CSV matches regeneration.')
