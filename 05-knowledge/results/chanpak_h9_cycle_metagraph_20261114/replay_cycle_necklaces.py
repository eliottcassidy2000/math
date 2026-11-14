from collections import Counter
from itertools import combinations, product
from math import comb

def require(b):
    if not b:
        raise AssertionError

def data(n):
    d=[x for x in range(1,n+1) if n%x==0]
    E={}; A={}
    for k in reversed(d):
        f=2**((comb(n,2)-n)//k) if k%2 else 0
        r=2**((comb(n,2)//k+n//2-(n//k)-(1 if n%2 else 2))//2) if k%2 else 0
        E[k]=f-sum(E[e] for e in E if e%k==0)
        A[k]=r-sum(A[e] for e in A if e%k==0)
    e=1 if n%2 else 2
    return E,A,sum(E[k]//(n//k) for k in d),sum(A[k]//e for k in d)

for n in range(3,8):
    edges=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(edges)}
    ring={(i,(i+1)%n) for i in range(n)}
    free=[i for i,e in enumerate(edges) if e not in ring and e[::-1] not in ring]
    base=sum(1<<i for i,e in enumerate(edges) if e in ring)
    def convert(r,mirror=False):
        return [(idx[tuple(sorted((a,b)))],int((a>b)!=mirror))
                for i,j in edges
                for a,b in [(r(i),r(j))]]
    maps=[convert(lambda i:(i+k)%n) for k in range(n)]
    s=convert(lambda i:n-1-i, True)
    def relabel(x,p):
        return sum((((x>>i)&1)^flip)<<j for i,(j,flip) in enumerate(p))
    def tri(x):
        score=[0]*n
        for i,(u,v) in enumerate(edges):
            score[u if (x>>i)&1 else v]+=1
        return comb(n,3)-sum(comb(a,2) for a in score)
    def strong(x,cut=-1):
        out=[0]*n; inc=[0]*n
        for i,(u,v) in enumerate(edges):
            if i==cut: continue
            if not (x>>i)&1: u,v=v,u
            out[u]|=1<<v; inc[v]|=1<<u
        for r in [out,inc]:
            reach=1; old=0
            while reach!=old:
                old=reach
                for i in range(n):
                    if (reach>>i)&1: reach|=r[i]
            if reach!=(1<<n)-1: return False
        return True
    er=Counter(); ar=Counter(); cyc=set(); sc=set(); fans=ns=0
    rmask=sum(1<<idx[tuple(sorted(e))] for e in ring)
    for mask in range(1<<len(free)):
        x=base|sum((1<<j) for i,j in enumerate(free) if (mask>>i)&1)
        rots={relabel(x,m) for m in maps}
        d=n//len(rots); er[d]+=1
        sx=relabel(x,s)
        if sx==x: ar[d]+=1
        if x==min(rots):
            cyc.add(x)
            if sx in rots: sc.add(x)
            # Lines indexed by a representative of each retained-edge orbit.
            yr=x^rmask; m=n//d
            changes=[tri(yr^(1<<idx[tuple(sorted((i,(i+1)%n)))]))-tri(x)
                     for i in range(m)]
            require(sum(changes)==-m)
            fans+=1
            if n<=6:
                bcount=sum(not strong(yr,j) for j in range(len(edges)))
                actual=sum(not strong(yr^(1<<idx[tuple(sorted((i,(i+1)%n)))]))
                           for i in range(m))
                require(actual*d==bcount)
                ns+=actual
    E,A,c,b=data(n)
    require(er==Counter({k:v for k,v in E.items() if v}))
    require(ar==Counter({k:v for k,v in A.items() if v}))
    require((len(cyc),len(sc))==(c,b))
    if n<=6:
        path={(i,i+1) for i in range(n-1)}
        pg=[i for i,e in enumerate(edges) if e not in path]
        tb=sum(1<<idx[e] for e in path)
        direct=sum(not strong(tb|sum((1<<j) for i,j in enumerate(pg) if (z>>i)&1))
                   for z in range(1<<len(pg)))
        require(ns==direct)
    print(n,len(cyc),len(sc),(len(cyc)+len(sc))//2,'fan checks',fans,
          'bridge boundary' if n<=6 else '',ns if n<=6 else '')
