"""Literal merged fixed-Hamilton-path tiling line adjacency, then exact detours."""
from pathlib import Path
from itertools import permutations,product
from functools import lru_cache
from collections import deque
import sys
from atlas import reps,rep8,canon_lane,H,restrict

@lru_cache(None)
def canon(T):
    n=len(T);f=(1<<n)-1
    best=None
    for A in (T,tuple(f^r^(1<<i) for i,r in enumerate(T))):
        B={}
        for i,r in enumerate(A):B.setdefault(r.bit_count(),[]).append(i)
        for J in product(*(permutations(B[q]) for q in sorted(B))):
            C=restrict(A,tuple(i for b in J for i in b))
            if best is None or C<best:best=C
    return best

def frames(T):
    n=len(T);f=(1<<n)-1
    def rec(P,M):
        if not M:yield P;return
        S=T[P[-1]]&M if P else M
        while S:
            a=S&-S;S-=a
            yield from rec(P+(a.bit_length()-1,),M^a)
    return rec((),f)

def line(T,F):
    n=len(T);bits=(1<<n)-1
    U=[bits^r^(1<<i) for i,r in enumerate(T)]
    for i,j in zip(F,F[1:]):
        U[i]^=1<<j;U[j]^=1<<i
    return tuple(U)

def blue(T,F):
    R=dict(zip(F,reversed(F)))
    return all((T[i]>>j&1)==(T[R[j]]>>R[i]&1) for i in F for j in F)

def graphs(path,path8=None):
    T=reps(path)
    maxn,cap=7,43
    if path8 is not None:
        T.update(rep8(path8));maxn,cap=8,73
    E={canon(t) for t in T}
    for n in range(3,maxn+1):
        V={T:H(T) for T in E if len(T)==n}
        G={T:{} for T in V}
        for T in V:
            if V[T]>cap:continue
            for F in frames(T):
                D=line(T,F);U=canon(D)
                if V[U]>cap:continue
                G[T].setdefault(U,(T,F,D,blue(T,F)))
        yield n,V,G

def solve(p,q,GG,cap=None):
    if cap is None:cap=q
    ans=None
    for n,V,G in GG:
        Q=deque(T for T in V if V[T]==p)
        prev={T:None for T in Q}
        while Q:
            T=Q.popleft()
            if V[T]==q:
                L=[];R=T
                while prev[R] is not None:
                    S=prev[R];L.append(G[S][R]);R=S
                L=L[::-1]
                if ans is None or (len(L),n)<(len(ans[1]),ans[0]):ans=n,L
                break
            for U in G[T]:
                if V[U]<=cap and U not in prev:
                    prev[U]=T;Q.append(U)
    return ans

def run():
    GG=list(graphs(Path(sys.argv[1]),Path(sys.argv[2]) if len(sys.argv)>2 else None))
    T=range(3,75 if len(sys.argv)>2 else 45,2)
    distances={}
    for q in T:
        if q in (7,21):continue
        p=canon_lane(q)
        if p in (7,21):continue
        r=solve(p,q,GG)
        print(q,p, None if r is None else (r[0],len(r[1]),[H(z[2]) for z in r[1]]))
        distances[q]=None if r is None else len(r[1])
        if r and len(r[1])>1:
            print(' detour masks',r)
    assert distances[17] is None and solve(13,17,GG,19)[0]==6
    assert solve(1,17,GG)[1][0][-1]
    if len(sys.argv)>2:
        assert len(distances)==31 and tuple(sorted(q for q,l in distances.items() if l is None))==(17,63)
        assert sum(x==1 for x in distances.values())==25
        assert solve(47,63,GG,73) is None
        assert len(solve(29,63,GG)[1])==1
        for q in range(3,74,2):
            if q in (7,21):continue
            assert any(V[T]<q and any(V[u]==q for u in edges)
                       for n,V,G in GG for T,edges in G.items())
if __name__=='__main__':run()
