import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import itertools,json,hashlib
from collections import defaultdict,deque
from core import states,compose
BANKS={'B29':[26,44,57,16,34,9],'B95':[64,121,35,37,13,81,28],'B389':[20,85,72,195,204,159,33,71]}
def need(x,msg):
 if not x:raise RuntimeError(msg)
def literal(adj):
 n=len(adj);h=0;q=0;aw=[0]*n;gw=[0]*n;sw=[0]*n;ew=[0]*n;K=[[0]*n for _ in adj]
 for p in itertools.permutations(range(n)):
  bad=[not(adj[p[i]]>>p[i+1]&1) for i in range(n-1)]
  if not any(bad):h+=1;sw[p[0]]+=1;ew[p[-1]]+=1
  for k in range(n+1):
   if any(bad[:max(k-1,0)]) or any(bad[k:]):continue
   q+=1
   if k<n:aw[p[k]]+=1
   if k:gw[p[k-1]]+=1
   if 0<k<n:K[p[k-1]][p[k]]+=1
 ans=[]
 for m in range(1<<n):
  vs=[v for v in range(n) if m>>v&1]
  ans.append((h,q,sum(aw[v] for v in vs),sum(gw[v] for v in vs),sum(K[u][v] for u in vs for v in vs),sum(sw[v] for v in vs),sum(ew[v] for v in vs)))
 return ans

def graph(ss):
 h,q=ss[0][:2];target=defaultdict(list)
 for m,s in enumerate(ss):target[s[2:4]].append(m)
 total=0;records=[];dig=hashlib.sha256()
 for src in sorted(target):
  seen={src:None};todo=deque([src]);hits=[]
  while todo:
   a,g=todo.popleft();x=g%2
   for y in (0,1):
    ap=2*a+h*y-q*x;gp=(g+q*y-h*x)//2;v=(ap,gp)
    if not(-h<ap<q and -h<gp<q):continue
    if v in target and not(src in [(0,0),(q-h,q-h)] and v==src):
     cur=(a,g);path=[(x,y)]
     while seen[cur] is not None:
      prev,px,py=seen[cur];path.append((px,py));cur=prev
     hits.append((target[v],path[::-1]))
    if v not in seen:seen[v]=((a,g),x,y);todo.append(v)
  dig.update(repr((src,sorted(seen))).encode())
  total+=len(seen)
  if hits:records.append(dict(source=target[src],states=len(seen),hits=hits))
 return total,dig.hexdigest(),records

def fold(parts):
 s=parts[0]
 for z in parts[1:]:s=compose(s,z)
 return s
one=[(1,2,b,b,0,b,b) for b in (0,1)]
if __name__=='__main__':
 out={}
 for name,adj in BANKS.items():
  ss=states(adj);need(ss==literal(adj),'literal mismatch '+name)
  total,digest,rs=graph(ss)
  fibers=defaultdict(set)
  for s in ss:fibers[s[2:4]].add(s[4:])
  out[name]=dict(h=ss[0][0],q=ss[0][1],AG_fibers=len(fibers),max_AG_fiber=max(map(len,fibers.values())),total_reachable=total,reachable_hash=digest,return_sources=len(rs),return_ingresses=sum(len(r['hits']) for r in rs),shortest_return=min((len(p) for r in rs for ms,p in r['hits']),default=None))
  if name in ('B29','B389'):need(not rs,'unexpected return')
  if name=='B95':
   rec=next(r for r in rs if r['source']==[14,22]);targets,path=next((ms,p) for ms,p in rec['hits'] if ms==[46,54])
   left=[fold([one[1],ss[m]]+[one[y] for x,y in path]+[one[1]]) for m in rec['source']]
   right=[fold([one[1]]+[one[x] for x,y in path]+[ss[m],one[1]]) for m in targets]
   zs=left+right;need(len({s[:4]+s[5:] for s in zs})==1,'interface mismatch')
   vals=sorted(s[4] for s in zs);need([v-vals[0] for v in vals]==[0,8,40,48],'wrong offsets')
   out['four_fiber']=dict(path=path,left_masks=rec['source'],right_masks=targets,states=zs,star_counts=sorted(2*s[2]+s[0]-s[4] for s in zs))
 # Empty contexts are omitted: the nonempty composition formula has no naive empty identity.
 need(states([0])==one,'singleton control')
 print(json.dumps(out,sort_keys=True,indent=2))
