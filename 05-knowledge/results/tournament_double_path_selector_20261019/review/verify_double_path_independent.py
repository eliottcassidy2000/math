#!/usr/bin/env python3
"""Independent literal-permutation and injective-search controls for D8."""
from itertools import permutations
from collections import Counter
import json
checks=0

def require(ok,context):
    global checks
    checks+=1
    if not ok:raise RuntimeError(context)

def transitive(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def tournament(out):
    n=len(out)
    require(all(not(out[i]>>i&1)for i in range(n)),('loops',out))
    require(all(((out[i]>>j)&1)+((out[j]>>i)&1)==1 for i in range(n)for j in range(i+1,n)),('tournament',out))
def histogram(edges,out):
    H=Counter()
    for f in permutations(range(len(out))):
        m=sum((not(out[f[a]]>>f[b]&1))<<j for j,(a,b)in enumerate(edges))
        H[m]+=1
    return H

def backtrack(edges,out,mask=0):
    n=len(out);directed=[(b,a)if(mask>>j&1)else(a,b)for j,(a,b)in enumerate(edges)]
    neighbors=[[]for _ in range(n)]
    for a,b in directed:neighbors[a].append((b,True));neighbors[b].append((a,False))
    order=sorted(range(n),key=lambda i:-len(neighbors[i]));f=[-1]*n
    def dfs(depth,used):
        if depth==n:return 1
        a=order[depth];total=0
        for x in range(n):
            if used>>x&1:continue
            if all(f[b]<0 or bool(out[x]>>f[b]&1)==forward for b,forward in neighbors[a]):
                f[a]=x;total+=dfs(depth+1,used|(1<<x));f[a]=-1
        return total
    return dfs(0,0)

D=[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7)]
hosts=[transitive(8),[210,188,169,81,164,201,22,72]]
records=[]
for idx,out in enumerate(hosts):
    tournament(out);H=histogram(D,out);F=H[0]+H[16]
    require(F==(480,617)[idx],('F count',idx,F))
    require(H[5]==(165,302)[idx],('D5 count',idx,H[5]))
    require(backtrack(D[:4]+D[5:],out)==F,('independent F',idx))
    require(backtrack(D,out,5)==H[5],('independent D5',idx))
    for m in range(128):
        epsilon=((m>>2)^(m>>4))&1
        require(H[m]%2==(1+epsilon*F)%2,('selector',idx,m))
        for e in range(7):
            require((H[m]+H[m^(1<<e)])%2==(F%2 if e in(2,4)else 0),('derivative',idx,m,e))
    records.append({'host':out,'F_count':F,'D_mask5_count':H[5],'orientation_counts':[H[m]for m in range(128)]})

undirected={frozenset(e)for e in D};auts=[]
for f in permutations(range(8)):
    if {frozenset((f[a],f[b]))for a,b in D}==undirected:auts.append(f)
require(len(auts)==2,('automorphism count',auts))
lookup={frozenset(e):j for j,e in enumerate(D)}
def image_mask(m,f):
    new=0
    for j,(a,b)in enumerate(D):
        if m>>j&1:a,b=b,a
        a,b=f[a],f[b];k=lookup[frozenset((a,b))]
        if(a,b)!=D[k]:new|=1<<k
    return new
orbits=set()
for m in range(128):
    orbit=tuple(sorted({image_mask(m,f)for f in auts}));orbits.add(orbit)
    require(len(orbit)==2,('free orbit',m,orbit))
    require(len({((v>>2)^(v>>4))&1 for v in orbit})==1,('selector invariance',orbit))
splits=Counter((((o[0]>>2)^(o[0]>>4))&1)for o in orbits)
require(splits=={0:32,1:32},('oriented classes',splits))

E7=[(0,1),(0,2),(0,3),(1,4),(2,5),(4,6)]
seven=[]
for idx,out in enumerate([transitive(7),[118,44,96,85,102,8,34]]):
    tournament(out);n=histogram(E7,out)[52]
    require(n==(66,81)[idx],('seven witness',idx,n))
    require(backtrack(E7,out,52)==n,('seven independent count',idx))
    if idx:require(all(x.bit_count()<6 for x in out),('source-free witness',out))
    seven.append({'host':out,'mask52_count':n})

result={'status':'PASS','checks':checks,'D8':records,'automorphisms':auts,'oriented_class_split':dict(splits),'sharp_seven_core':seven}
print(json.dumps(result,indent=2,sort_keys=True))
