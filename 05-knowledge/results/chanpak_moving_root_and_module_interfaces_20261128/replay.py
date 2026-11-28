"""Verify rooted extension maps and module probe examples by independent DP."""
from itertools import product, combinations
from pathlib import Path
from collections import Counter
from math import factorial, comb
import sys

def arcs(n,mask):
    A=[0]*n
    for k,(x,y) in enumerate(combinations(range(n),2)):
        if not mask>>k&1:x,y=y,x
        A[x]|=1<<y
    return A

def dp(A):
    n=len(A); limit=1<<n
    incoming=[sum(1<<i for i in range(n) if A[i]>>j&1) for j in range(n)]
    H=[1]+[0]*(limit-1)
    D=[[0]*n for _ in range(limit)]
    for S in range(1,limit):
        if S.bit_count()==1:D[S][S.bit_length()-1]=1
        else:
            W=S
            while W:
                b=W&-W;j=b.bit_length()-1;W-=b
                U=(S^b)&incoming[j]
                while U:
                    x=U&-U;U-=x
                    D[S][j]+=D[S^b][x.bit_length()-1]
        H[S]=sum(D[S])
    return H,D

def state(A,a):
    H,D=dp(A)
    rev=[sum(1<<i for i in range(len(A)) if A[i]>>j&1) for j in range(len(A))]
    _,S=dp(rev)
    all_=len(H)-1
    return H[all_],D[all_][a],sum(S[i][a]*H[all_^i] for i in range(len(H))),H[all_^(1<<a)]

def move(v,op):
    h,e,u,g=v
    if op=="P":return h,0,2*u,g
    if op=="K":return h-e+u,g,2*u-g,g
    if op=="M":return h-e+u,h-e,h+u,h
    raise ValueError(op)

def make(word):
    A=[0];a=0;v=(1,1,1,1)
    for op in word:
        b=len(A)
        A=[r|((1<<b) if i!=a or op=="P" else 0) for i,r in enumerate(A)]
        A.append((1<<a) if op!="P" else 0)
        if op=="M":a=b
        v=move(v,op)
    return A,a,v

def probes():
    for k in range(7):
        for ops in product("KMP",repeat=k):
            A,a,v=make(ops)
            assert state(A,a)==v,(ops,state(A,a),v)
    print("Independent full register checks",sum(3**k for k in range(7)))

def moments(A):
    H,_=dp(A); U=len(H)-1
    return H[U],sum(h*H[U^S] for S,h in enumerate(H))
def triangle(X, m=1):
    n=len(X); A=[0]*(n+2*m)
    for x in range(n):
        A[x]=X[x]|sum(1<<y for y in range(n,n+m))
    for y in range(n,n+m):
        A[y]=sum(1<<i for i in range(y+1,n+2*m))
    for z in range(n+m,n+2*m):
        A[z]=(1<<n)-1|sum(1<<i for i in range(z+1,n+2*m))
    return A

def same_count():
    for n in range(3,7):
        cache={}
        for mask in range(1<<(n*(n-1)//2)):
            X=arcs(n,mask)
            h,j=moments(X)
            if h in cache and cache[h][1]!=j:
                prev,J=cache[h]
                assert dp(triangle(X))[0][-1]==h+j
                assert dp(triangle(arcs(n,prev)))[0][-1]==h+J
                print("First equal count with unequal triangle substitution:",
                      n,h,"masks, colored moments, outputs",
                      (prev,J,h+J),(mask,j,h+j))
                return
            cache[h]=(mask,j)

def residues(p):
    def reduce(v):
        return tuple(x%p for x in v)
    init=reduce((1,1,1,1)); seen={init}; Q=[init]
    for v in Q:
        for op in "KMP":
            w=reduce(move(v,op))
            if w not in seen:seen.add(w);Q.append(w)
    pool=set()
    for h in range(p):
        for e in range(p):
            for g in range(p):
                pool.add((h,e,(h+g-e)%p,g))
            for u in range(p):
                pool.add((h,e,u,e));pool.add((h,0,u,e))
    assert seen<=pool
    return len(seen),len(pool),len({v[0] for v in seen}),len(pool-seen-{(0,0,0,0)})

if __name__=="__main__":
    probes()
    same_count()
    for p in (3,5,7,9,11,13,17,19,23,25):
        print("Mod",p,residues(p))
    seeds=Path(__file__).with_name('seeds.txt')
    for l in seeds.read_text().splitlines():
        h,n,m=map(int,l.split())
        assert dp(arcs(n,m))[0][-1]==h
    print("Independently checked fallback seed masks",len(seeds.read_text().splitlines()))
    if len(sys.argv)>1:
        path=Path(sys.argv[1]);seen=set()
        for line in path.read_text().splitlines():
            h,*p=line.split();h=int(h)
            A,a,v=make(p[0] if p else '')
            assert h==v[0] and h not in seen
            seen.add(h)
        missing=set(range(1,max(seen)+1,2))-seen
        assert missing=={7,21}|{int(l.split()[0]) for l in seeds.read_text().splitlines()}
        print("Bank max",max(seen),"distinct",len(seen),"missing",sorted(missing))
