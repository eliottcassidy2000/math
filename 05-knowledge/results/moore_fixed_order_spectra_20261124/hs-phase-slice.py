#!/usr/bin/env python3
from itertools import permutations
from collections import Counter

def make(u,v,guard,selector=None):
    # Each site is a tagged (side,cohort,coordinate).
    X=[(0,i,0) for i in range(1,4)]
    if guard: X.append((1,0,1))
    X.extend([(1,1,1),(1,2,1),(1,1,2),(1,3,1),(1,1,3)])
    if selector is not None: X=X[:3]+[X[i] for i in selector]
    def edge(i,j):
        side,a,x=X[i];bside,b,y=X[j]
        if side==bside:
            return a==b and ((y-x)%5 in ((1,4) if side==0 else (2,3)))
        if side:
            a,b,x,y=b,a,y,x
        return y==(x+u[a]*v[b])%5
    mask=0
    n=len(X)
    for i in range(n):
        for j in range(i+1,n):
            t=i*n-i*(i+1)//2+j-i-1
            # Lexicographic pair bit index.
            assert t==sum(n-1-z for z in range(i))+j-i-1
            if edge(i,j): mask|=1<<t
    return X,mask

def read(X,mask):
    n=len(X)
    def e(a,b):
        i,j=sorted([X.index(a),X.index(b)])
        t=i*n-i*(i+1)//2+j-i-1
        return (mask>>t)&1
    def head(j,y):
        hits=[i for i in (1,2,3) if e((0,i,0),(1,j,y))]
        assert len(hits)<=1
        return hits[0] if hits else 4
    u=[0]*5
    for y in (1,2,3):
        i=head(1,y)
        assert u[i]==0
        u[i]=y
    u[u.index(0,1)]=4
    v=[0,1,*[pow(u[head(j,1)],-1,5) for j in (2,3)],0]
    v[4]=next(x for x in (2,3,4) if x not in v)
    assert set(v)==set(u)==set(range(5))
    return tuple(u),tuple(v)

def count(n,mask):
    R=[0]*n
    k=0
    for i in range(n):
        for j in range(i+1,n):
            if mask>>k&1:R[j]|=1<<i
            else:R[i]|=1<<j
            k+=1
    dp=[[0]*n for i in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    allbits=(1<<n)-1
    for s in range(1,(1<<n)):
        for i,x in enumerate(dp[s]):
            if not x: continue
            r=R[i]&(allbits^s)
            while r:
                b=r&-r;r-=b
                dp[s|b][b.bit_length()-1]+=x
    return sum(dp[-1])

def full_graph(u,v):
    nodes=[(side,i,x) for side in (0,1) for i in range(5) for x in range(5)]
    def e(s,t):
        side,a,x=s;tside,b,y=t
        if side==tside:return a==b and (y-x)%5 in ((1,4) if side==0 else (2,3))
        if side:a,b,x,y=b,a,y,x
        return y==(x+u[a]*v[b])%5
    adj=[{j for j,t in enumerate(nodes) if e(s,t)} for s in nodes]
    assert all(len(S)==7 for S in adj)
    for i in range(50):
        for j in range(i+1,50):
            assert len(adj[i]&adj[j])+(j in adj[i])==1

for g in (False,True):
    source=Counter()
    plain=Counter()
    mate=Counter()
    pair=set()
    masks=set()
    for u0 in permutations(range(1,5)):
        for v0 in permutations((2,3,4)):
            u=(0,*u0);v=(0,1,*v0)
            X,mask=make(u,v,g);n=len(X)
            assert read(X,mask)==(u,v) and mask not in masks
            masks.add(mask)
            full_graph(u,v)
            k=0;spine=0
            for i in range(n):
                for j in range(i+1,n):
                    if j==i+1:spine|=1<<k
                    k+=1
            width=(1<<k)-1
            plain[count(n,mask)]+=1
            if mask&spine:continue
            H=count(n,mask);K=count(n,mask^(width^spine))
            source[H]+=1;mate[K]+=1;pair.add((H,K))
            if u==v==tuple(range(5)):print("standard",n,mask,mask^(width^spine),H,K)
    print(n,"all src",len(plain), min(plain),max(plain))
    print(n,"legal",sum(source.values()),"source", len(source),min(source),max(source),
         "mate",len(mate),min(mate),max(mate),"different height pairs",len(pair))

# Two legal local orders suffice for all 144 tables on only eight sites.
def spine_mask(n):
    return sum(1<<(i*n-i*(i+1)//2) for i in range(n-1))

n=8
S=spine_mask(n)
alts=(5,3,4,6,7)
pair_to_phase={}
sources=set();mates=set()
for u0 in permutations(range(1,5)):
    for v0 in permutations((2,3,4)):
        u=(0,*u0);v=(0,1,*v0)
        X,mask=make(u,v,False)
        if mask&S: X,mask=make(u,v,False,alts)
        assert not (mask&S) and read(X,mask)==(u,v)
        C=count(n,mask);D=count(n,mask^(((1<<28)-1)^S))
        assert (C,D) not in pair_to_phase
        pair_to_phase[(C,D)]=(u,v)
        sources.add(C);mates.add(D)
assert len(pair_to_phase)==144
print("adaptive eight:",len(sources),min(sources),max(sources),
      len(mates),min(mates),max(mates),len(pair_to_phase))


# A nonpayload W-W arc can be a literal Chan--Pak +2 move.
u=v=tuple(range(5))
X,C=make(u,v,False)
assert (X[5],X[7])==((1,1,2),(1,1,3))
D=C^(1<<26)
assert not ((C|D)&S) and read(X,C)==read(X,D)==(u,v)
assert [count(8,k) for k in (C,D)]==[317,319]
print("physical nonpayload arc flip:",C,D,317,319,
      "tile mate after flip",count(8,D^(((1<<28)-1)^S)))
