import itertools,json
from fractions import Fraction
from math import prod

def require(ok):
    if not ok: raise RuntimeError('check failed')
def parts(r,k):
    if r==0:return [()]
    return [a for a in itertools.product(range(k),repeat=r) if a[0]==0 and all(a[i]<=1+max(a[:i]) for i in range(1,r))]
def rank(a):
    a=[[Fraction(x) for x in row] for row in a]; t=0
    for c in range(len(a[0])):
        p=next((p for p in range(t,len(a)) if a[p][c]),None)
        if p is None:continue
        a[t],a[p]=a[p],a[t];v=a[t][c];a[t]=[x/v for x in a[t]]
        for p in range(t+1,len(a)):
            v=a[p][c]
            if v:a[p]=[x-v*y for x,y in zip(a[p],a[t])]
        t+=1
        if t==len(a):break
    return t
def count(n,edges,k,bound=()):
    return sum(all(c[u]!=c[v] for u,v in edges) for tail in itertools.product(range(k),repeat=n-len(bound)) for c in [bound+tail])
rows=[];checks=0
for r in range(5):
 for k in range(1,5):
    ps=parts(r,k);es=list(itertools.combinations(range(r),2));E=[tuple(e for i,e in enumerate(es) if m>>i&1) for m in range(1<<len(es))]
    a=[[int(all(p[u]!=p[v] for u,v in e)) for e in E] for p in ps]
    q=rank(a);require(q==len(ps));checks+=1;rows.append([r,k,len(ps),q])
# Two fragments on three common ports, each with one private vertex.
# Exhaust all 64 edge sets of each four-vertex fragment, palettes2..4.
es=list(itertools.combinations(range(4),2));families=[tuple(e for i,e in enumerate(es) if m>>i&1) for m in range(64)]
for k in range(2,5):
 ps=parts(3,k);fs=[[count(4,e,k,p) for p in ps] for e in families]
 for i,g in enumerate(families):
  for j,h in enumerate(families):
   edges=set(g)|{tuple(4 if x==3 else x for x in e) for e in h}
   expected=sum(prod(range(k-max(p),k+1))*fs[i][a]*fs[j][a] for a,p in enumerate(ps))
   require(expected==count(5,edges,k));checks+=1
for k in range(2,9):
 require(count(3,[(0,1)],k)==k*k*(k-1));require(count(3,[(0,1),(0,2)],k)==k*(k-1)**2);checks+=2
print(json.dumps({'status':'PASS','checks':checks,'partition_test_ranks':rows,'gluing_pairs':3*64*64,'hostile_palettes':list(range(2,9))},indent=2))
