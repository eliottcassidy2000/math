#!/usr/bin/env python3
"""Erasure controls for the edge-rooted Moore route verifier."""
from itertools import combinations, product


def instance(n):
    if n==2:
        V=list(range(10))
        edges=[(i,(i+1)%5) for i in range(5)]+[(i+5,(i+2)%5+5) for i in range(5)]+[(i,i+5) for i in range(5)]
        u,v=0,1
    else:
        V=[(s,i,j) for s in range(2) for i in range(5) for j in range(5)]
        edges=[((s,i,j),(s,i,(j+1+s)%5)) for s,i,j in V]+[
          ((0,i,x),(1,j,(x+i*j)%5)) for i,x,j in product(range(5),repeat=3)]
        u,v=(0,0,0),(0,0,1)
    E=set(map(frozenset,edges))
    edge=lambda x,y:frozenset((x,y)) in E
    heads=lambda u,v:[x for x in V if x!=v and edge(u,x)]
    grid=[[next(x for x in V if edge(x,r) and edge(x,c)) for c in heads(v,u)] for r in heads(u,v)]
    return {(i,j):tuple(next(y for y in range(n) if edge(grid[i][x],grid[j][y])) for x in range(n))
            for i in range(n) for j in range(n) if i!=j}


def inv(p):
    return tuple(p.index(x) for x in range(len(p)))


def residual(P,n,K,i,j):
    if i==j:
        return tuple(frozenset((x,)) for x in range(n))
    out=[]
    for x in range(n):
        seen=[x,P[i,j][x]]+[P[t,j][P[i,t][x]] for t in K if t not in (i,j)]
        assert len(seen)==len(set(seen))
        out.append(frozenset(range(n))-set(seen))
    return tuple(out)


def check(P,n):
    K=range(n)
    assert all(sorted(v)==list(K) and inv(v)==P[j,i] for (i,j),v in P.items())
    assert all(not any(residual(P,n,K,i,j)) for i,j in combinations(K,2))


def repair_one(P,n,r):
    K=[i for i in range(n) if i!=r]
    a=K[0]
    F={i:tuple(next(iter(x)) if len(x)==1 else (_ for _ in ()).throw(ValueError('not a singleton'))
          for x in residual(P,n,K,i,a)) for i in K}
    if any(sorted(p)!=list(range(n)) for p in F.values()):
        raise ValueError('not a permutation')
    t=[]
    for z in range(n):
        used={F[i].index(z) for i in K}
        if len(used)!=n-1:
            raise ValueError('repeated old neighbor column')
        t.append(next(x for x in range(n) if x not in used))
    if sorted(t)!=list(range(n)):
        raise ValueError('repeated new column')
    R=dict(P)
    for i in K:
        R[i,r]=tuple(t[z] for z in F[i])
        R[r,i]=inv(R[i,r])
    check(R,n)
    return R


def split(res,first):
    """Two matchings of the residual, oriented by one bit per cycle."""
    n=len(res)
    assert all(len(row)==2 for row in res)
    col=[{i for i in range(n) if y in res[i]} for y in range(n)]
    assert all(len(c)==2 for c in col)
    p=[None]*n
    rowbits=[None]*n
    colbits=[None]*n
    for x in range(n):
        if rowbits[x] is not None:
            continue
        i,y=x,min(res[x])
        while rowbits[i] is None:
            rowbits[i]=colbits[y]=first
            z=next(iter(res[i]-{y}))
            p[i]=(y,z)
            i=next(iter(col[z]-{i}))
            y=z
        assert i==x
        first+=1
    rev=tuple(inv(tuple(p[i][b] for i in range(n))) for b in (0,1))
    return p,rev,rowbits,colbits,first


def factor_rules(P,n,K):
    a=K[0]
    models={}
    N=0
    for i in K[1:]:
        data=split(residual(P,n,K,i,a),N)
        models[i]=data[:4]
        N=data[4]
    rules=set()
    for i,j in combinations(K[1:],2):
        I,_,ib,_=models[i]
        J,Jr,jb,jc=models[j]
        allowed=residual(P,n,K,i,j)
        for x in range(n):
            for t,phase,u in product(range(2),repeat=3):
                y=I[x][phase^t]
                dest=Jr[phase^u][y]
                if dest not in allowed[x]:
                    rules.add((ib[x],t,jc[y],u))
            for z in range(n):
                if set(I[x])==set(J[z]):
                    for t,u in product(range(2),repeat=2):
                        if I[x][t]==J[z][u]:
                            rules.add((ib[x],t,jb[z],u))
    return models,N,rules


def model_ok(rules,bits):
    return all(bits[i]!=a or bits[j]!=b for i,a,j,b in rules)


def direct_residual_check(P,n,K,models,bits):
    for i,j in combinations(K[1:],2):
        I,_,ib,_=models[i]
        _,Z,_,zb=models[j]
        target=residual(P,n,K,i,j)
        for x in range(n):
            seen=[]
            for b in range(2):
                y=I[x][b^bits[ib[x]]]
                seen.append(Z[b^bits[zb[y]]][y])
            if set(seen)!=target[x]:
                return False
    return True


if __name__=='__main__':
    for n in (2,6):
        P=instance(n);check(P,n)
        for r in range(n):
            visible={p:x for p,x in P.items() if r not in p}
            assert repair_one(visible,n,r)==P
        print(n,'all single-row erasures restored')
    n=6
    P=instance(n)
    certs=0
    for u,v in combinations(range(n),2):
        K=[i for i in range(n) if i not in (u,v)]
        models,N,rules=factor_rules(P,n,K)
        actual=[None]*N
        for i,(opts,_,bits,_) in models.items():
            for x in range(n):
                b=opts[x].index(P[u,K[0]][P[i,u][x]])
                assert actual[bits[x]] in (None,b)
                actual[bits[x]]=b
        assert model_ok(rules,actual)
        for bits in product(range(2),repeat=N):
            assert model_ok(rules,bits)==direct_residual_check(P,n,K,models,bits)
            certs+=1
    print('HS, all 15 double erasures; Boolean factorizer matches',certs,'assignments')
