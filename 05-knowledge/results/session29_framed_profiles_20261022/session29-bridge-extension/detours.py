from pathlib import Path
from itertools import combinations,permutations
from collections import deque
import json
import check as c
D=Path(__file__).parent
P={n:list(combinations(range(n),2)) for n in range(1,11)}
def adj(m,n):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(P[n]):a[i][j]=(m>>k)&1;a[j][i]=1-a[i][j]
 return a
def mask(a):return sum(a[i][j]<<k for k,(i,j) in enumerate(P[len(a)]))
def hp(a):
 n=len(a);out=[]
 def rec(p,left):
  if not left:out.append(p);return
  for x in left:
   if not p or a[p[-1]][x]:rec(p+[x],left-{x})
 rec([],set(range(n)));return out
def flip(a,p):
 n=len(a);b=[[int(i!=j and not a[i][j]) for j in range(n)] for i in range(n)]
 for u,v in zip(p,p[1:]):b[u][v]=1;b[v][u]=0
 return b
def scc(a):
 n=len(a);reach=[]
 for i in range(n):
  seen={i};todo=[i]
  while todo:
   u=todo.pop()
   for v in range(n):
    if a[u][v] and v not in seen:seen.add(v);todo.append(v)
  reach.append(seen)
 groups=[];left=set(range(n))
 while left:
  i=min(left);g={j for j in left if j in reach[i] and i in reach[j]};groups.append(g);left-=g
 return sorted(groups,key=lambda g:-len(reach[next(iter(g))]))
def cm(a):
 n=len(a)
 return min(mask([[a[i][j] for j in p] for i in p]) for p in permutations(range(n)))
def sourcekey(a,groups):return tuple((len(g),cm([[a[i][j] for j in sorted(g)] for i in sorted(g)])) for g in groups)
def rootcan(a,vertices,root,op=False):
 best=None;bw=None
 for p in permutations([v for v in vertices if v!=root]):
  p=p+(root,);m=mask([[(a[j][i] if op else a[i][j]) for j in p] for i in p])
  if best is None or m<best:best=m;bw=p
 return best,bw
def targetaddr(a):
 for aa in combinations(range(10),5):
  A=set(aa);B=set(range(10))-A;es=[(u,v) for u in A for v in B if a[u][v]]
  if len(es)==1:
   u,v=es[0];q,px=rootcan(a,A,u);r,py=rootcan(a,B,v,True)
   if q not in c.ci or r not in c.ci:return None
   conv=q>r;p=(py+px) if conv else (px+py);w=[p.index(i) for i in range(10)]
   key=tuple(sorted((q,r)))
   return key,w,conv
 return None
def bfs(start,target,ban):
 todo=deque([start]);prev={start:None}
 while todo and target not in prev:
  u=todo.popleft()
  for e,v in c.adj[u]:
   if e==ban or v in prev:continue
   prev[v]=(u,e);todo.append(v)
 if target not in prev:return None
 path=[];v=target
 while prev[v] is not None:v,e=prev[v];path.append(e)
 return list(reversed(path))
old=json.loads((D.parent/'session29-target-wiring/wiring-results.json').read_text());records=[]
for idx in (169,194):
 row=old['witness_rows'][idx];T=adj(row['target'],10);e=row['tiles'][0];seen=set();us=[]
 for p in hp(T):
  U=flip(T,p);g=scc(U);sizes=[len(x) for x in g];cum={sum(sizes[:i]) for i in range(1,len(sizes))}
  if 5 in cum or not ({4,6}&cum) or max(sizes)>6:continue
  key=sourcekey(U,g)
  if key not in seen:seen.add(key);us.append((U,p,key))
 print('remaining',e,'distinct alternative sources',len(us),flush=True)
 found=None;tested=0
 for U,p,key in us:
  for q in hp(U):
   tested+=1;V=flip(U,q);a=targetaddr(V)
   if a is None:continue
   addr,w,conv=a;dest=('Q',)+addr;path=bfs(('S',41,41),dest,(e,e))
   if path is None:continue
   found={'old_edge':idx,'tile':e,'target_mask':mask(T),'target_to_U_path':p,'U_mask':mask(U),'U_scc_key':key,'U_to_V_path':q,'V_mask':mask(V),'V_address':addr,'V_to_canonical_permutation':w,'V_converse':conv,'boundary_path_normalized':path,'boundary_path_actual':[[x,c.tau(y)] for x,y in path],'alternative_source_classes':len(us),'neighbor_paths_tested':tested};break
  if found:break
 records.append(found if found else {'old_edge':idx,'tile':e,'status':'unresolved','alternative_sources':len(us),'tested':tested})
 print(records[-1],flush=True)
(D/'detours-results.json').write_text(json.dumps(records,indent=2)+'\n')
