"""Exact blue fold/observer-line census for C5 and one Petersen gauge."""
from itertools import permutations, combinations, product
from collections import Counter
from functools import cache


def observer(G, order):
    n=len(order)
    r=[0]*(2*n)
    for i,j in combinations(range(n),2):
        r[i]|=1<<j
        r[n+j]|=1<<(n+i)
    for i,j in product(range(n),repeat=2):
        if i==j or order[j] in G[order[i]]:r[i]|=1<<(n+j)
        else:r[n+j]|=1<<i
    return tuple(r)


def wholeheight(r):
    n=len(r)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s,row in enumerate(D):
        for i,v in enumerate(row):
            if not v:continue
            a=r[i]&~s
            while a:
                z=a&-a;a^=z;D[s|z][z.bit_length()-1]+=v
    return sum(D[-1])


def fold(r,closed):
    n=len(r)//2
    @cache
    def go(mask,last):
        if mask==(1<<n)-1:return r[last]>>((last+n)%(2*n))&1
        return sum(go(mask|1<<(j%n),j) for j in range(2*n) if not(mask>>(j%n)&1) and r[last]>>j&1)
    return sum(go(1<<(i%n),i) for i in range(2*n) if bool(r[(i+n)%(2*n)]>>i&1)==closed)


def bpaths(r):
    n=len(r)//2
    def rec(p,mask):
        if len(p)==n:
            if p[-1]<n:yield p+tuple((i+n)%(2*n) for i in reversed(p))
            return
        for i in range(2*n):
            if mask>>(i%n)&1:continue
            if not p or r[p[-1]]>>i&1:yield from rec(p+(i,),mask|1<<(i%n))
    return rec((),0)


def mate(r,p):
    s=list(r)  # use converse gauge: reverse only spine edges
    for i,j in zip(p,p[1:]):s[i]^=1<<j;s[j]^=1<<i
    return tuple(s)


def canon(r):
    groups=[]
    for d in sorted({x.bit_count() for x in r}):
        groups.append(tuple(i for i,x in enumerate(r) if x.bit_count()==d))
    return min(sum((r[i]>>j&1)<<k for k,(i,j) in enumerate(combinations(sum(pieces,()),2)))
               for pieces in product(*(tuple(permutations(g)) for g in groups)))


def cyc5():
    return {i:{(i+1)%5,(i-1)%5} for i in range(5)}

def petersen():
    V=list(combinations(range(5),2))
    return {v:{w for w in V if not(set(v)&set(w))} for v in V}


def petersen_chart(G):
    a=min(G);b=min(G[a]);F={a,b};D=F|G[a]|G[b]
    return sorted(set(G)-D)+sorted(D-F)+sorted(F)


G=cyc5();N=5
C={}
for p in permutations(G):
    if sum(p.index(x) for x in G[p[0]])==5:continue
    r=observer(G,p)
    k=canon(r)
    if k not in C:C[k]=(p,r)
assert len(C)==8
stats={}
nclasses=set(C)
global_neighbors=set()
inter=[]
for source,(order,r) in sorted(C.items()):
    h=wholeheight(r)
    O,Cl=fold(r,False),fold(r,True)
    stats[source]=(order,h,O,Cl)
    occur=Counter()
    closed=Counter()
    ds=Counter()
    for p in bpaths(r):
        mt=mate(r,p);t=canon(mt)
        occur[t]+=1
        closed[t]+=p[0]>=N
        if r[p[0]].bit_count()==N:ds[t]+=1
        if t in nclasses:inter.append((source,t,p))
    assert sum(occur.values())==O+Cl and sum(closed.values())==Cl
    global_neighbors.update(occur)
    print('C5', source,order,h,O,Cl,len(occur),dict(Counter(occur.values())),flush=True)
assert not inter and len(global_neighbors)==857
assert sum(O+Cl for order,h,O,Cl in stats.values())==980
def ordinary_paths(r):
    def rec(p,unused):
        if not unused:yield p;return
        options=unused & (r[p[-1]] if p else unused)
        while options:
            z=options&-options;options^=z
            yield from rec(p+(z.bit_length()-1,),unused^z)
    return rec((),(1<<len(r))-1)
nr=0
for _,r in C.values():
    d=[x.bit_count() for x in r]
    for P in ordinary_paths(r):
        if d[P[0]]!=d[P[-1]]+1:continue
        nr+=1
        assert canon(mate(r,P)) not in C
assert nr==10851
print('C5:',len(global_neighbors),'distinct external blue mates, all',nr,'score-compatible roots miss observer family')

G=petersen()
order=petersen_chart(G);r=observer(G,order)
print('Petersen order',order,'blue-open',fold(r,False),'blue-closed',fold(r,True),flush=True)
assert (fold(r,False),fold(r,True))==(960507,1278014)
P=tuple(range(10))+tuple(range(19,9,-1))
mt=mate(r,P)
print('Petersen canonical mate, for same anti:',fold(mt,False),fold(mt,True),flush=True)
def aut(r):
    groups=[tuple(i for i,x in enumerate(r) if x.bit_count()==v)
            for v in sorted({x.bit_count() for x in r})]
    c=0
    for target in product(*(tuple(permutations(g)) for g in groups)):
        P=[0]*len(r)
        for a,b in zip(sum(groups,()),sum(target,())):P[a]=b
        c+=all((r[i]>>j&1)==(r[P[i]]>>P[j]&1) for i,j in combinations(range(len(r)),2))
    return c
assert (aut(r),aut(mt))==(1,1)
assert (fold(mt,False),fold(mt,True))==(1195661,1270678)
