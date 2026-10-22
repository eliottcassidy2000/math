"""O(n^3) intrinsic bottleneck recognition with exact hostile comparisons."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
W=json.loads((HERE/'wiring-results.json').read_text())
def require(b,msg):
 if not b:raise RuntimeError(msg)
def adj(m,n):
 out=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   a,b=(i,j) if m>>k&1 else (j,i);out[a]|=1<<b;k+=1
 return out
def transpose(out):
 n=len(out);return [sum(((out[j]>>i)&1)<<j for j in range(n)) for i in range(n)]
def tree(out):
 seen={0};q=[0];es=[]
 for u in q:
  for v in range(len(out)):
   if out[u]>>v&1 and v not in seen:seen.add(v);q.append(v);es.append((u,v))
 require(len(seen)==len(out),'arborescence exists');return es
def scc(out):
 n=len(out);rev=transpose(out);seen=set();order=[]
 def dfs(u):
  seen.add(u)
  for v in range(n):
   if out[u]>>v&1 and v not in seen:dfs(v)
  order.append(u)
 for u in range(n):
  if u not in seen:dfs(u)
 seen=set();cs=[]
 def collect(u,c):
  seen.add(u);c.append(u)
  for v in range(n):
   if rev[u]>>v&1 and v not in seen:collect(v,c)
 for u in reversed(order):
  if u not in seen:c=[];collect(u,c);cs.append(c)
 return cs
def recognize(out):
 n=len(out);require(len(scc(out))==1,'input strong')
 candidates=set(tree(out))|{(v,u) for u,v in tree(transpose(out))};ans={};bridges=[]
 for u,v in sorted(candidates):
  a=out.copy();a[u]&=~(1<<v);cs=scc(a)
  if len(cs)==1:continue
  bridges.append((u,v));owner={x:i for i,c in enumerate(cs) for x in c};succ=[set() for c in cs];deg=[0]*len(cs)
  for x in range(n):
   for y in range(n):
    i,j=owner[x],owner[y]
    if i!=j and a[x]>>y&1 and j not in succ[i]:succ[i].add(j);deg[j]+=1
  seq=[];unused=set(range(len(cs)))
  while unused:
   zeros=[i for i in unused if deg[i]==0];require(len(zeros)==1,'total condensation order');i=zeros[0];seq.append(i);unused.remove(i)
   for j in succ[i]:deg[j]-=1
  require(owner[v]==seq[0] and owner[u]==seq[-1],'bridge returns sink to source')
  cut=0
  for i in reversed(seq[1:]):
   cut|=sum(1<<x for x in cs[i]);require(cut not in ans,'unique outgoing cut arc');ans[cut]=(u,v)
 require(len(candidates)<=2*n-2,'candidate budget');return ans,bridges,len(candidates)
def brute(out):
 n=len(out);ans={};full=(1<<n)-1
 for cut in range(1,full):
  cross=[(u,v) for u in range(n) if cut>>u&1 for v in range(n) if not(cut>>v&1) and out[u]>>v&1]
  if len(cross)==1:ans[cut]=cross[0]
 return ans
for row in W['witness_rows']:
 out=adj(row['target'],10);a,_,_=recognize(out);saved={sum(1<<v for v in c['vertices']):tuple(c['arc']) for c in row['one_arc_cuts']};require(a==saved,'all324 targets versus exhaustive cuts')
small=0
for n in (3,4,5):
 for m in range(1<<(n*(n-1)//2)):
  out=adj(m,n)
  if len(scc(out))>1:continue
  small+=1;require(recognize(out)[0]==brute(out),'exhaustive small recognition')
letters=sorted(int(t) for t in json.loads((HERE.parent/'session29-profiles/compressed-matrices.json').read_text()));flip=1023^sum(1<<k for k,(i,j) in enumerate(( (i,j) for i in range(5) for j in range(i+1,5))) if j==i+1)
large=[]
for r in (2,3,4,7,10,15,20):
 n=5*r;out=[0]*n
 for i in range(r):
  x=adj(letters[(7*i+r)%18]^flip,5)
  for u in range(5):out[5*i+u]=(x[u]<<(5*i))|((1<<(5*i))-1)
 for i in range(r-1):u=5*i+4;v=u+1;out[u]|=1<<v;out[v]&=~(1<<u)
 cuts,bs,c=recognize(out);require(all((1<<(5*j))-1 in cuts for j in range(1,r)),'all prescribed seams recovered')
 inner=sorted((x for x in cuts if 2<=x.bit_count()<=n-2),key=int.bit_count);require(all(x&y==x for x,y in zip(inner,inner[1:])),'interior cut chain')
 large.append({'r':r,'vertices':n,'candidates':c,'strong_bridges':len(bs),'one_arc_cuts':len(cuts)})
result={'status':'PASS','r2_targets_checked':324,'small_strong_tournaments_checked':small,'large_tests':large}
(HERE/'recognition-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
