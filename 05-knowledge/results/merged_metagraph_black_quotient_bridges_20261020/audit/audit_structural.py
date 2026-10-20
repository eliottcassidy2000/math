#!/usr/bin/env python3
import itertools as it, json, math
from collections import Counter
from pathlib import Path
def adjacency(n,x,free):
    a=[[0]*n for _ in range(n)]
    for i in range(n-1):a[i][i+1]=1
    for k,(i,j) in enumerate(free):
        if x&(1<<k):a[i][j]=1
        else:a[j][i]=1
    return a

def classify(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    pairs=list(it.combinations(range(n),2))
    perms=list(it.permutations(range(n)))
    def code(a,p):
        return sum(a[p[i]][p[j]]*2**k for k,(i,j) in enumerate(pairs))
    mask=2**len(pairs)-1
    keys=[];selfconverse=[];As=[]
    for x in range(2**len(free)):
        a=adjacency(n,x,free); As.append(a)
        orbit=[code(a,p) for p in perms]
        ck=min(orbit);co=mask-max(orbit)
        keys.append(min(ck,co));selfconverse.append(ck==co)
    codes=sorted(set(keys));pi=[codes.index(k) for k in keys]
    return free,codes,pi,selfconverse,As
OUT=Path(__file__).resolve().parent

def fib(n):
 a=[1]+[0]*n
 for k in range(1,n+1):a[k]=a[k-1]+(a[k-3] if k>=3 else 0)
 p=a[n//2] if n%2==0 else a[n//2]+(a[n//2-1] if n>=3 else 0)
 return (a[n]-p)//2+(a[(n-4)//2] if n>=4 and n%2==0 else 0)

def run(n):
 free,codes,pi,sc,As=classify(n);N=len(pi);m=len(free);q=len(codes)
 reflperm=[free.index((n-1-j,n-1-i)) for i,j in free]
 rho=lambda x:sum(((x>>k)&1)<<p for k,p in enumerate(reflperm))
 masses=[pi.count(u) for u in range(q)]
 blue=[sum(pi[x]==u and rho(x)==x for x in range(N)) for u in range(q)]
 black=[[0]*q for _ in range(q)]
 for x in range(N):
  if rho(x)!=x:black[pi[x]][pi[x^(N-1)]]+=1
 assert all(black[u][v]%2==0 for u in range(q) for v in range(q))
 assert all(black[u][u]%4==0 for u in range(q))
 edges=[(u,v,black[u][v]//2) for u in range(q) for v in range(u+1,q) if black[u][v]]
 def comps(skip=None):
  adj=[set() for _ in range(q)]
  for u,v,k in edges:
   if (u,v)==skip:continue
   adj[u].add(v);adj[v].add(u)
  unseen=set(range(q));cc=0
  while unseen:
   todo=[unseen.pop()];cc+=1
   while todo:
    x=todo.pop();fresh=adj[x]&unseen;unseen-=fresh;todo.extend(fresh)
  return cc
 cc=comps();bridges=[(u,v) for u,v,k in edges if k==1 and comps((u,v))==cc+1]
 pure=[];signed=[]
 for u in range(q):
  a=As[pi.index(u)];score=list(map(sum,a));perms=list(it.permutations(range(n)))
  hp=[p for p in perms if all(a[i][j] for i,j in zip(p,p[1:]))]
  auto=sum(all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
  sign=lambda p:1-2*(sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2)
  R=sum(sign(p) for p in hp);g=math.prod(2*a[i][j]-1 for i in range(n) for j in range(i+1,n))
  z=sum(1-2*(x.bit_count()%2) for x in range(N) if pi[x]==u)
  epsilon=1 if sc[pi.index(u)] else 2
  assert z*auto==(-1)**m*epsilon*g*R
  signed.append({'code':codes[u],'z':z,'H':len(hp),'Aut':auto,'R':R,'g':g})
  if blue[u]!=masses[u]:continue
  # Strong iff all ordered pairs communicate.
  reach=[r[:] for r in a]
  for i in range(n):reach[i][i]=1
  for k in range(n):
   for i in range(n):
    for j in range(n):reach[i][j]|=reach[i][k]&reach[k][j]
  strong=all(all(r) for r in reach)
  if n%2==0:assert masses[u]==1 and not strong
  cycles=[p for p in hp if p[0]==0 and a[p[-1]][p[0]]]
  if strong:
   assert n%2==1 and len(set(score))==1
   for p in cycles:
    successor={p[i]:p[(i+1)%n] for i in range(n)}
    assert all(a[i][j]==a[successor[i]][successor[j]] for i in range(n) for j in range(n))
   assert len(cycles)<=auto-1
  pure.append({'code':codes[u],'mass':masses[u],'strong':strong,'HC':len(cycles),'Aut':auto})
 bound=fib(n)
 if n>=5:assert len(bridges)>=bound
 return {'n':n,'pure_blue':pure,'Q_bridges':bridges,'Q_bridge_lower_bound_formula':bound,'signed_bridge_nodes_verified':q,'signed_data':signed}

if __name__=='__main__':
 rows=[run(n) for n in range(3,7)]
 (OUT/'audit_structural.json').write_text(json.dumps(rows,indent=2)+'\n')
 for row in rows:print(json.dumps({k:v for k,v in row.items() if k!='signed_data'}))
