import json,itertools,math,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'collision'
receipt=json.loads((ROOT/'receipt.json').read_text())
def need(b,msg):
 if not b:raise RuntimeError(msg)
def literal(r):
 n=len(r);full=(1<<n)-1
 def extend(last,used):
  if used==full:return 1
  opts=r[last]&(full^used);total=0
  while opts:
   bit=opts&-opts;opts-=bit
   total+=extend(bit.bit_length()-1,used|bit)
  return total
 return sum(extend(u,1<<u) for u in range(n))
def extension(r):
 n=len(r);out=list(r)+[1<<(n+1),0]
 for v in range(n):
  if v==0:out[v]|=1<<n
  else:out[n]|=1<<v
  out[n+1]|=1<<v
 return out
def properties(r):
 n=len(r)
 for u in range(n):
  need(not(r[u]>>u&1),'loop')
  for v in range(u):need((r[u]>>v&1)+(r[v]>>u&1)==1,'tournament')
 reach=[[bool(r[u]>>v&1) or u==v for v in range(n)] for u in range(n)]
 for k in range(n):
  for u in range(n):
   for v in range(n):reach[u][v]|=reach[u][k] and reach[k][v]
 need(all(all(row) for row in reach),'strong')
 modules=[]
 for mask in range(1,(1<<n)-1):
  if mask.bit_count()<2:continue
  vs=[v for v in range(n) if mask>>v&1];out=((1<<n)-1)^mask
  if all((r[v]&out)==(r[vs[0]]&out) for v in vs):modules.append(mask)
 return modules
base=receipt['base_rows'];out=[]
need(literal(base)==3489,'base Hamilton')
need(literal(extension(base))==17055,'base query')
for h in receipt['scalar_hostiles']:
 r=h['rows'];need(properties(r)==[],'prime')
 need(literal(r)==3489,'hostile Hamilton');need(literal(extension(r))==15579,'hostile query')
 flips=[(i,j) for i in range(11) for j in range(i+1,11) if bool(base[i]>>j&1)!=bool(r[i]>>j&1)]
 need(flips==[tuple(x) for x in h['reversed_pairs']],'edges')
 out.append({'flips':flips,'H':3489,'response':15579,'prime':True,'strong':True})
lines=(ROOT/'perturb.out').read_text().splitlines();expected={}
for line in lines:
 if line.startswith('n='):
  d=dict(re.findall(r'(\w+)=(\d+)',line));key=tuple(int(d[k]) for k in ['n','seed','radius']);expected[key]=tuple(int(d[k]) for k in ['total','no_pair','same_H'])
actual={};matches=[]
for line in (HERE/'census.out').read_text().splitlines():
 fields=line.split()
 if fields[0]=='ROW':actual[tuple(map(int,fields[1:4]))]=tuple(map(int,fields[4:]))
 else:matches.append(tuple(map(int,fields[1:])))
need(len(actual)==27,'census finished');need(actual==expected,'independent complete census')
need(sum(v[0] for v in actual.values())==464709,'edit sets');need(sum(v[1] for v in actual.values())==28293,'pairfree');need(sum(v[2] for v in actual.values())==3,'scalar matches')
need(sorted(m[2:] for m in matches)==sorted(tuple(h['rows']) for h in receipt['scalar_hostiles']),'full match list')
for (n,j,k),v in actual.items():need(v[0]==math.comb(math.comb(n,2),k),'binomial universe')
print(json.dumps({'status':'ACCEPT','cases':464709,'pairfree':28293,'scalar_matches':3,'named_query_source':17055,'hostiles':out,'literal_count_method':'Uncached DFS over every valid directed Hamilton word','census_method':'Independently built inflations, incoming-based pair predicate, predecessor endpoint DP'},sort_keys=True,indent=2))
