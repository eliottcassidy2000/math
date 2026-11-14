from itertools import permutations, product
from math import comb
from random import Random
A=(2,12,17,37,75,23,47)
B=(2,12,17,37,43,71,31)
def H(A):
    n=len(A);F=(1<<n)-1
    s=[[0]*n for _ in range(F+1)]
    for i in range(n):s[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(s[S]):
            if not x:continue
            r=A[i]&~S
            while r:
                b=r&-r;r^=b
                s[S|b][b.bit_length()-1]+=x
    return sum(s[F]) if n else 1
def inflate(A,K):
    mats=[K.get(i,(0,)) for i in range(len(A))]
    shift=0
    blocks=[]
    for M in mats:
        blocks.append((((1<<len(M))-1)<<shift,shift));shift+=len(M)
    return tuple((r<<blocks[i][1]) | sum(F for j,(F,s) in enumerate(blocks) if A[i]>>j&1)
                 for i,M in enumerate(mats) for r in M)
def paint(A,G,gm,W,wm):
    E=[(i,j) for i in range(2,7) for j in range(2,7) if A[i]>>j&1]
    m=len(G);n=len(W)
    E=[(i-2,j-2) for i,j in E]
    for M,p,offset,isW in [(G,gm,5,False),(W,wm,5+m,True)]:
        for i,r in enumerate(M):
            E.extend((offset+i,offset+j) for j in range(len(M)) if r>>j&1)
            for hub in (0,1):E.append((offset+i,hub) if p>>i&1 else (hub,offset+i))
            if isW:E.extend((offset+i,j) for j in range(2,5+m))
            else:E.extend((j,offset+i) for j in range(2,5))
    P=[0]*(5+m+n)
    for i,j in E:P[i]|=1<<j
    return tuple(P)
def C2(U):
    ans=0;n=len(U);F=(1<<n)-1
    def sub(S):
        ix=[j for j in range(n) if S>>j&1]
        return tuple(sum(((U[i]>>v)&1)<<j for j,v in enumerate(ix)) for i in ix)
    return sum(H(sub(S))*H(sub(F^S)) for S in range(1,F))
def tt(n):return tuple(((1<<n)-1)^((2<<i)-1) for i in range(n))
def visits(A,v):
    B=[i for i in range(len(A)) if i!=v];out=[0]*len(A)
    def b(i,j):return A[i]>>j&1
    for p in permutations(B):
        r=0;g=b(v,p[0])+b(p[-1],v);ok=True
        for i,j in zip(p,p[1:]):
            a=b(i,v)&b(v,j);bad=1-b(i,j)
            if bad and not a:ok=False;break
            g+=a;r+=bad
        if not ok:continue
        for k in range(max(1,r),g+1):out[k]+=comb(g-r,k-r)
    return out
rng=Random(6769)
M=[(0,),tt(2),tt(3),(2,4,1)]
assert (H(A),H(B))==(67,69)
for X,Y,Z in product(M,repeat=3):
    a,b=[H(inflate(Q,{0:X,1:Y,6:Z})) for Q in (A,B)]
    expr=2*H(X)*H(Y)*H(Z)-C2(Z)*(2*H(X)+C2(X))*(H(Y)+C2(Y))
    assert b-a==expr
for i in range(120):
    G,W=rng.choice(M),rng.choice(M)
    gm,wm=rng.randrange(1<<len(G)),rng.randrange(1<<len(W))
    a,b=[H(paint(Q,G,gm,W,wm)) for Q in (A,B)]
    assert b-a==2*H(G)*H(W)
    assert H(paint(B,G,gm,(),0))-H(paint(A,G,gm,(),0))==2*H(G)
for n in range(1,6):
    a,b=[H(paint(Q,tt(2),2,tt(n),0)) for Q in (A,B)]
    assert a==124*3**n-21*2**n-53 and b==a+2
assert visits(paint(A,tt(2),2,tt(1),0),7)[1:4]==[277,351,124]
assert visits(paint(B,tt(2),2,tt(1),0),7)[1:4]==[279,351,124]
print('Literal 67/69, 64 arbitrary three-module substitutions, 120 two-bank examples, spectral clock certified.')

def four_tapes(U,C,V,G,W,um,vm,gm,wm,reverse=False):
    blocks=[U,C,V,(0,),(0,),G,W]
    ofs=[0]
    for L in blocks:ofs.append(ofs[-1]+len(L))
    T=[0]*ofs[-1]
    def add(u,v):T[u]|=1<<v
    for i,R in enumerate(blocks):
        for u,r in enumerate(R):T[ofs[i]+u]|=r<<ofs[i]
    for i,j in [(0,1),(0,2),(1,2),(0,5),(1,5),(2,5),(6,0),(6,1),(6,2),(6,5)]:
        for u in range(ofs[i],ofs[i+1]):
            for v in range(ofs[j],ofs[j+1]):add(u,v)
    p,q=ofs[3],ofs[4];a,b=ofs[1],ofs[1]+1
    add(*((q,p) if reverse else (p,q)))
    add(a,p);add(q,a);add(p,b);add(b,q)
    for j,mask in [(0,um),(1,0),(2,vm),(5,gm),(6,wm)]:
        for v in range(ofs[j],ofs[j+1]):
            if v in (a,b):continue
            for u in (p,q):add(*((v,u) if mask>>(v-ofs[j])&1 else (u,v)))
    return tuple(T)
assert [H(four_tapes((),(6,12,24,17,3),(),(),(),0,0,0,0,j)) for j in [0,1]]==[89,91]
for _ in range(30):
    U,C,V,G,W=rng.choice(M[:2]),(6,12,24,17,3),rng.choice(M[:2]),rng.choice(M[:2]),rng.choice(M[:2])
    masks=[rng.randrange(1<<len(R)) for R in [U,V,G,W]]
    old,new=[H(four_tapes(U,C,V,G,W,*masks,k)) for k in (0,1)]
    assert new-old==2*H(U)*H(V)*H(G)*H(W)
print('Literal 89/91 and 30 independent four-tape extensions certified.')
