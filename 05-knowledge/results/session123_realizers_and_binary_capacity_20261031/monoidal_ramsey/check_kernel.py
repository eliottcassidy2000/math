#!/usr/bin/env python3
"""Finite controls, not a substitute for the uniform proof. No assert guards."""
import json
from functools import lru_cache
from fractions import Fraction
from random import Random

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def kernel(a,b):
    if a<0 or b<0:
        raise ValueError('nonnegative budgets required')
    if not a or not b:
        return 0
    aa,bb=a-1,b-1
    states={(True,True):1}
    for bit in range(max(aa.bit_length(),bb.bit_length())-1,-1,-1):
        na,nb=(aa>>bit)&1,(bb>>bit)&1
        nxt={}
        for (ta,tb),weight in states.items():
            for u in range(2):
                for v in range(2):
                    if (ta and u>na) or (tb and v>nb):
                        continue
                    key=(ta and u==na,tb and v==nb)
                    nxt[key]=nxt.get(key,0)+weight*(2 if u and v else 1)
        states=nxt
    return sum(states.values())

@lru_cache(None)
def witness(a,b):
    """Memoized construction DAG, not a graph of falsely rounded order."""
    if a<1 or b<1:
        raise ValueError('positive graph budgets required')
    if a==1:
        return ('empty_edges',b)
    if b==1:
        return ('clique',a)
    return ('C5',tuple(witness(a//2+int(a%2 and i in (0,2)),
                                 b//2+int(b%2 and i in (0,1))) for i in range(5)))

@lru_cache(None)
def stats(w):
    if w[0]=='empty_edges':
        return w[1],1,w[1]
    if w[0]=='clique':
        return w[1],w[1],1
    s=[stats(c) for c in w[1]]
    return (sum(t[0] for t in s),
            max(s[i][1]+s[(i+1)%5][1] for i in range(5)),
            max(s[i][2]+s[(i+2)%5][2] for i in range(5)))

def expand(w):
    if w[0]=='empty_edges':
        return [0]*w[1]
    if w[0]=='clique':
        n=w[1]
        return [((1<<n)-1)^(1<<i) for i in range(n)]
    parts=[expand(c) for c in w[1]]
    offsets=[]
    n=0
    for p in parts:
        offsets.append(n)
        n+=len(p)
    masks=[((1<<len(p))-1)<<o for p,o in zip(parts,offsets)]
    result=[]
    for i,(part,o) in enumerate(zip(parts,offsets)):
        outside=masks[(i-1)%5]|masks[(i+1)%5]
        result.extend((row<<o)|outside for row in part)
    return result

def brute_stats(adj):
    n=len(adj)
    clique=indep=0
    for subset in range(1,1<<n):
        size=subset.bit_count()
        if size<=min(clique,indep):
            continue
        vertices=[i for i in range(n) if subset>>i&1]
        if all((adj[i]&subset)==(subset^(1<<i)) for i in vertices):
            clique=max(clique,size)
        if all((adj[i]&subset)==0 for i in vertices):
            indep=max(indep,size)
    return n,clique,indep

def run():
    for a in range(65):
        for b in range(65):
            direct=sum(1<<((i&j).bit_count()) for i in range(a) for j in range(b))
            require(kernel(a,b)==direct,('digit DP',a,b))
    actual=0
    for a in range(1,65):
        for b in range(1,65):
            require(stats(witness(a,b))==(kernel(a,b),a,b),('constructor',a,b))
            if kernel(a,b)<=16:
                require(brute_stats(expand(witness(a,b)))==(kernel(a,b),a,b),('actual graph',a,b))
                actual+=1
    rng=Random(123)
    for _ in range(10000):
        a=[rng.randrange(65) for i in range(5)]
        b=[rng.randrange(65) for i in range(5)]
        A=max(a[i]+a[(i+1)%5] for i in range(5))
        B=max(b[i]+b[(i+2)%5] for i in range(5))
        require(sum(kernel(x,y) for x,y in zip(a,b))<=kernel(A,B),('weighted C5',a,b))
    for a in range(1,65):
        L=1<<(a-1).bit_length()
        slope=sum((Fraction(3,2)**i.bit_count() for i in range(a)),Fraction(0))
        for b in range(129):
            require(kernel(a,b+L)-kernel(a,b)==L*slope,('translation',a,b))
    huge=kernel(10**100+7,10**101+13)
    require(huge==kernel(10**101+13,10**100+7),'huge symmetry')
    return {'status':'PASS','direct_grid':[0,64],'symbolic_witness_grid':[1,64],
            'actual_graphs_enumerated':actual,'actual_graph_vertex_limit':16,
            'random_weighted_controls':10000,'translation_grid':[64,128],
            'huge_evaluation_decimal_digits':len(str(huge))}
if __name__=='__main__':
    print(json.dumps(run(),sort_keys=True))
