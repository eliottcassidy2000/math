from itertools import combinations, permutations
from collections import Counter

def order(a):
    return tuple(p for p in permutations(range(len(a)))
                 if all(a[u][v] for u,v in zip(p,p[1:])))

def scc(a):
    n=len(a)
    r=[[i==j or a[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n): r[i][j] |= r[i][k] and r[k][j]
    groups={frozenset(j for j in range(n) if r[i][j] and r[j][i])
            for i in range(n)}
    return sorted(groups,key=lambda g:-sum(r[next(iter(g))]))

def wrap(a):
    m=len(a)
    b,c,x=range(m,m+3)
    t=[[0]*(m+3) for _ in range(m+3)]
    for u in range(m):
        t[u][:m]=a[u]
        t[u][b]=t[u][x]=t[c][u]=1
    t[b][c]=t[c][x]=t[x][b]=1
    return t

def flip(t,u,v): t[u][v],t[v][u]=t[v][u],t[u][v]

def check():
    fans=0
    for m in range(1,5):
        es=list(combinations(range(m),2))
        for mask in range(1<<len(es)):
            a=[[0]*m for _ in range(m)]
            for j,(u,v) in enumerate(es):
                a[u][v]=mask>>j&1
                a[v][u]=1-a[u][v]
            paths=order(a)
            groups=scc(a)
            idx={v:j for j,g in enumerate(groups) for v in g}
            q=0
            for mask2 in range(1<<m):
                L=[u for u in range(m) if mask2>>u&1]
                R=[u for u in range(m) if not mask2>>u&1]
                q+=len(order([[a[u][v] for v in L] for u in L]))*len(order([[a[u][v] for v in R] for u in R]))
            t=wrap(a)
            assert len(order(t))==len(paths)+2*q
            b,c,x=range(m,m+3)
            for p in paths:
                fans+=1
                ring=(b,c,*p,x)
                edges=tuple(zip(ring,ring[1:]+ring[:1]))
                rev=[row[:] for row in t]
                for u,v in edges: flip(rev,u,v)
                boundary=Counter()
                for u,v in edges:
                    flip(rev,u,v)
                    sizes=tuple(map(len,scc(rev)))
                    flip(rev,u,v)
                    expect=(m+3,)
                    if (u,v)==(b,c): expect=(m+3,)
                    elif u==c: expect=(1,m+2)
                    elif v==b: expect=(m+2,1)
                    elif v==x: expect=(m+1,1,1) if m>1 else (1,1,1,1)
                    elif idx[u]!=idx[v]:
                        i=p.index(v)
                        expect=(i+1,m-i+2) if i>1 else (1,1,m+1)
                    assert sizes==expect, (m,p,u,v,sizes,expect)
                    if len(sizes)>1: boundary[sizes]+=1
                assert sum(boundary.values())==len(groups)+2
    print(f'{fans} fans: exact SCC compositions and count identity')

if __name__=='__main__': check()
