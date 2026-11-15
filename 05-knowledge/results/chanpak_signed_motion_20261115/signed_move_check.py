#!/usr/bin/env python3
"""Independent subset path DP and explicit head tournament construction."""
from itertools import combinations


def tournament(n, code):
    g = [[0]*n for _ in range(n)]
    for k, (u, v) in enumerate(combinations(range(n), 2)):
        g[u][v] = (code >> k) & 1
        g[v][u] = 1-g[u][v]
    return g


def endpoint(g):
    n = len(g)
    d = [[0]*n for _ in range(1 << n)]
    for u in range(n):
        d[1 << u][u] = 1
    for s, row in enumerate(d):
        for u, h in enumerate(row):
            if h:
                for v in range(n):
                    if g[u][v] and not (s >> v) & 1:
                        d[s | (1 << v)][v] += h
    return d


def profile(g, m):
    d = endpoint(g)
    r = endpoint(tuple(zip(*g)))
    F = len(d)-1
    counts = [1]+list(map(sum, d[1:]))
    e = [sum(v for j,v in enumerate(R) if (m >> j) & 1) for R in d]
    l = [sum(v for j,v in enumerate(R) if (m >> j) & 1) for R in r]
    q = sum(counts[s]*counts[F^s] for s in range(F+1))
    E = sum(e[s]*counts[F^s] for s in range(F+1))
    L = sum(counts[s]*l[F^s] for s in range(F+1))
    J = sum(e[s]*l[F^s] for s in range(F+1))
    return counts[F],q,E,L,J


def build(g,m,a,b):
    n = len(g)
    U = a+n+b
    t,x = U,U+1
    P = [[0]*(U+2) for _ in range(U+2)]
    for u,v in combinations(range(U+2),2):
        if v == x:
            arc = (u == t) or (a <= u < a+n and (m >> (u-a)) & 1)
        elif a <= u < v < a+n:
            arc = g[u-a][v-a]
        else:
            arc = 1
        P[u][v]=int(arc)
        P[v][u]=1-int(arc)
    return P


def require(flag):
    if not flag:
        raise AssertionError("independent witness comparison failed")


def check(rows, expected, one_sided=False):
    for i,(n,code,m,u,v,p) in enumerate(rows):
        g=tournament(n,code)
        require(profile(g,m)==p)
        for s in (1,2,3):
            a,b=s+u,v if one_sided else s+v
            h,q,E,L,J=p
            lit=sum(endpoint(build(g,m,a,b))[-1])
            raw=(q << (a+b))+(E << a)-(L << b)+h-J
            require(lit==raw==expected(s,i))


if __name__=='__main__':
    check([(5,4,14,0,2,(5,64,40,22,12)),
           (5,4,13,0,2,(5,64,36,21,10)),
           (5,4,20,1,1,(5,64,17,41,8)),
           (5,4,24,1,1,(5,64,15,39,6))],
          lambda s,i:256*4**s-48*2**s-7+2*i)
    check([(4,4,9,0,1,(5,32,12,12,4)),
           (4,4,10,0,1,(5,32,14,13,2)),
           (4,4,4,0,1,(5,32,6,9,0))],
          lambda s,i:64*4**s-12*2**s+1+2*i)
    check([(6,32759,7,1,2,(9,168,35,112,20)),
           (6,32757,7,1,2,(9,168,31,110,18)),
           (6,32622,25,0,2,(33,336,150,132,40)),
           (6,32750,6,1,2,(9,168,29,109,14))],
          lambda s,i:1344*4**s-378*2**s-11+2*i)
    check([(6,32615,48,0,1,(37,376,138,96,20)),
           (6,32694,10,1,1,(13,208,29,88,10)),
           (6,32759,37,1,1,(9,168,109,71,38))],
          lambda s,i:890*2**s-175+2*i,True)
    print('14 independent profile checks and 42 explicit head tournament comparisons passed')


def full_masks(n,code):
    g=tournament(n,code)
    d=endpoint(g)
    r=endpoint(tuple(zip(*g)))
    F=(1<<n)-1
    counts=[1]+list(map(sum,d[1:]))
    q=sum(counts[z]*counts[F^z] for z in range(F+1))
    masks={}
    for m in range(1,F+1):
        ei=[sum(v for j,v in enumerate(row) if m>>j&1) for row in d]
        li=[sum(v for j,v in enumerate(row) if m>>j&1) for row in r]
        masks[m]=(sum(ei[z]*counts[F^z] for z in range(F+1)),
                  sum(counts[z]*li[F^z] for z in range(F+1)),
                  sum(ei[z]*li[F^z] for z in range(F+1)),ei[F],li[F])
    return counts[F],q,masks


def ordsum(g,r):
    n,m=len(g),len(r)
    return [list(row)+[1]*m for row in g]+[[0]*n+list(row) for row in r]


def halves(H,Q,ab,h,q,cd):
    E,L,J,T,S=ab
    e,l,j,t,s=cd
    return h*E+Q*e,q*L+H*l,H*h-h*J-H*j-T*s-L*e


def check_pair(n1,c1,n2,c2,ms,expected,shelves=None):
    h,q,A=full_masks(n1,c1)
    k,r,B=full_masks(n2,c2)
    require(q*r==expected[0])
    for ix,(m1,m2) in enumerate(ms):
        u,v,w=halves(h,q,A[m1],k,r,B[m2])
        require((u,v,w)==(expected[1],expected[2],expected[3]+2*ix))
        g=ordsum(tournament(n1,c1),tournament(n2,c2))
        # Compare the literal joined head at the least guard lengths.
        mask=m1|(m2<<n1)
        lit=sum(endpoint(build(g,mask,1,1))[-1])
        require(lit==4*q*r+2*(u-v)+w)
    shelves_all={}
    for a in A.values():
        for b in B.values():
            x,y,z=halves(h,q,a,k,r,b)
            shelves_all.setdefault((x,y),set()).add(z)
    require(all(not(d+2 in pool and d+4 in pool) for pool in shelves_all.values() for d in pool))
    if shelves is not None:
        for table,maskset in zip((A,B),shelves):
            seed=next(iter(maskset))
            a=table[seed]
            require({m for m,v in table.items() if a[:2]+a[3:]==v[:2]+v[3:]}==maskset)


if __name__=='__main__':
    check_pair(7,1464080,8,45442625,[(25,57),(35,101)],
               (5059200,2554038,1952430,-942835),({25,35},{57,101}))
    check_pair(5,955,5,955,[(21,21),(7,28)],(10816,6786,6786,-3935))
    G=ordsum(tournament(5,955),tournament(5,955))
    for mark,J in ((21|(21<<5),4104),(7|(28<<5),4102)):
        require(profile(G,mark)==(169,10816,6786,6786,J))
        require(sum(v for i,v in enumerate(endpoint(G)[-1]) if mark>>i&1)==104)
        require(sum(v for i,v in enumerate(endpoint(tuple(zip(*G)))[-1]) if mark>>i&1)==104)
    print('Both coupled repairs and their fixed-core mask ceilings verified.')
