"""Independent exact historical controls. No producer imports; explicit guards."""
from itertools import combinations
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def tournament(rows):
 n=len(rows)
 for i in range(n):
  need(not rows[i]>>i&1,'loop')
  for j in range(i):need(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,'not tournament')
 return rows
def substitute(children):
 offsets=[];n=0
 for c in children:offsets.append(n);n+=len(c)
 rows=[0]*n
 for i,c in enumerate(children):
  for u,r in enumerate(c):rows[offsets[i]+u]=r<<offsets[i]
  j=(i+1)%3
  for u in range(len(c)):
   for v in range(len(children[j])):rows[offsets[i]+u]|=1<<(offsets[j]+v)
 return tournament(rows)
def reverse(rows):
 n=len(rows);return [((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows)]
def strong(rows):
 n=len(rows)
 for root in range(n):
  reached=1<<root
  while True:
   nxt=reached
   for u in range(n):
    if reached>>u&1:nxt|=rows[u]
   if nxt==reached:break
   reached=nxt
  if reached!=(1<<n)-1:return False
 return True
def all_h(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)];out=[1]+[0]*((1<<n)-1)
 for s in range(1,1<<n):
  for v in range(n):
   if not s>>v&1:continue
   t=s^(1<<v)
   dp[s][v]=sum(dp[t][u] for u in range(n) if t>>u&1 and rows[u]>>v&1) if t else 1
  out[s]=sum(dp[s])
 return out
def autos(rows):
 n=len(rows);deg=[r.bit_count() for r in rows];p=[];used=set();answer=[]
 def rec():
  k=len(p)
  if k==n:answer.append(tuple(p));return
  for v in range(n):
   if v in used or deg[k]!=deg[v]:continue
   if any(((rows[k]>>i)&1)!=((rows[v]>>p[i])&1) for i in range(k)):continue
   p.append(v);used.add(v);rec();used.remove(v);p.pop()
 rec();return answer
def cyc(rows):return [s for s in combinations(range(len(rows)),3) if all(sum((rows[u]>>v)&1 for v in s)==1 for u in s)]
def orbit_count(G,n):return len({min(sum(1<<g[i] for i in range(n) if s>>i&1) for g in G) for s in range(1<<n)})
A0=tournament([14,4,8,2]);A1=reverse(A0)
T=substitute([A0]*3);U=substitute([A0,A0,A1]);GT=autos(T);GU=autos(U)
need((len(GT),len(GU))==(81,27),'81/27 full groups')
ht,hu=all_h(T),all_h(U);need(ht==hu,'all named subset H');need(ht[-1]==462915,'full H');need(cyc(T)==cyc(U) and len(cyc(T))==67,'67 cyclic triples')
for v in range(12):need({g for g in GT if g[v]==v}=={g for g in GU if g[v]==v},'literal point stabilizers')
need((orbit_count(GT,12),orbit_count(GU,12))==(176,512),'attachment workloads')
# All-strong run-order control: preserve exact singleton order, reverse only nonpath arcs.
A=tournament([10,12,9,16,7]);D=tournament([2,4,9,3]);W=substitute([A,D,D]);run=[]
need(strong(A) and strong(D) and strong(W),'all-strong source and children')
for middle,expected in [((0,1,2),2899675),((0,2,1),2906091)]:
 word=(3,)+middle+(4,);P=[]
 for i,v in enumerate(word):
  P.append(v)
  if i<4:P.extend((5+i,9+i))
 B=reverse(W)
 for u,v in zip(P,P[1:]):
  need(W[u]>>v&1,'Hamilton path condition');B[u]|=1<<v;B[v]&=~(1<<u)
 need(strong(B),'all-strong partner');h=all_h(B)[-1];need(h==expected,'all strong run-order H')
 run.append({'path':P,'H':h})
print(json.dumps({'status':'FINITE-EXACT','decorated_triangle':{'aut_orders':[len(GT),len(GU)],'cyclic_triples':67,'subsets':len(ht),'H':ht[-1],'attachment_orbits':[176,512],'point_stabilizers':'literally equal on all 12 labels'},'all_strong_run_order':run},sort_keys=True,indent=2))
