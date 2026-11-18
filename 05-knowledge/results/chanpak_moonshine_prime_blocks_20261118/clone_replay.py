"""Literal and compressed-slot replay, including exact shared-ray collisions."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from math import comb, prod
from random import Random
from frame_replay import WITNESSES, TRIPLE, add, cut, join
from gap_replay import rand_t
from padding_replay import slots
from replay import profile

def inflate(A, blocks):
    n = len(A); pieces = [blocks.get(i, (0,)) for i in range(n)]
    rep = [0]*n; off = 0
    for i,T in enumerate(pieces):
        rep[i] = ((1<<len(T))-1)<<off
        off += len(T)
    result=[];off=0
    for i,T in enumerate(pieces):
        exterior=sum(rep[j] for j in range(n) if A[i]>>j&1)
        result += [(a<<off)|exterior for a in T]
        off += len(T)
    return tuple(result)

def vectors(A, marks, q, cap):
    n=len(A);f=(1<<n)-1
    reuse=sum(1<<i for i in marks)
    @lru_cache(None)
    def go(todo, last, slots_left, counts):
        if not slots_left:
            return Counter({counts:1}) if not todo and all(counts) else Counter()
        result=go(todo,n,slots_left-1,counts).copy()
        available=f if last==n else A[last]
        for i in range(n):
            if not available>>i&1:continue
            if i in marks:
                j=marks.index(i)
                if counts[j]<cap:
                    nxt=counts[:j]+(counts[j]+1,)+counts[j+1:]
                    result.update(go(todo,i,slots_left,nxt))
            elif todo>>i&1:
                result.update(go(todo^(1<<i),i,slots_left,counts))
        return result
    return go(f^reuse,n,q,(0,)*len(marks))

def read(T,k):
    return slots(T,k)[k][-1]

def TT(d):return tuple((1<<d)-(1<<(i+1)) for i in range(d))

def P(d):
    A=list(TT(d));A[0]^=1<<(d-1);A[-1]=1
    return tuple(A)

def fast(d,j):
    return j**d+j*((j+1)**(d-2)-(j-1)**(d-2))

def nonempty(T,r):
    return sum((-1)**(r-j)*comb(r,j)*read(T,j) for j in range(1,r+1))

def run():
    rng=Random(908141)
    for _ in range(70):
        n=rng.randrange(2,6);A=rand_t(n,rng);s=rng.randrange(1,3)
        m=tuple(rng.sample(range(n),s));Q=[rand_t(rng.randrange(1,4),rng) for _ in m]
        for q in (1,2):
            C=vectors(A,m,q,q+n-s)
            T=inflate(A,dict(zip(m,Q)))
            p=0
            vals=[{r:nonempty(V,r) for r in range(1,q+n-s+1)} for V in Q]
            for t,c in C.items():p+=c*prod(vals[i][r] for i,r in enumerate(t))
            assert p==read(T,q)
    print('140 literal/compressed checks, including two marks')

    # Cyclic support is infinite formally, but finite physical blocks
    # truncate each color at its number of vertices.
    A=join((2,4,1),(0,))
    m=(0,1,2)
    for q in (1,2):
        C=vectors(A,m,q,4)
        assert C[(4,4,4)]>0
        for k in ((1,2,2),(2,2,3),(3,2,2)):
            Q=tuple(TT(i) for i in k)
            vals=[{r:nonempty(V,r) for r in range(1,5)} for V in Q]
            predicted=sum(c*prod(vals[i][r] for i,r in enumerate(a)) for a,c in C.items())
            assert predicted==read(inflate(A,dict(zip(m,Q))),q)
    print('six cyclic-mark literal checks and a repeated triangle run')

    expect=((-2808,1696,-3728,-25456,-8904,1096,5040,34144,40),
            (-128324,-16292,-155716,15828,160748,24092,115340,3412,-14220,156))
    for (U,I,_,_,_),e in zip(WITNESSES,expect):
        A=add(U,I);B=cut(A,1<<len(U))
        assert tuple(profile(inflate(B,{i:TT(2)}))[1]-profile(inflate(A,{i:TT(2)}))[1]
                     for i in range(len(A)))==e
    T=[cut(TRIPLE,b) for b in (185,251,164)]
    for i,j in ((0,1),(0,2),(1,2)):
        for x in range(len(TRIPLE)):
            R=min(len(TRIPLE)+1,2+max(min(A[x].bit_count(),len(A)-1-A[x].bit_count()) for A in (T[i],T[j])))
            assert any(profile(inflate(T[i],{x:TT(k)}))[1]!=profile(inflate(T[j],{x:TT(k)}))[1]
                       for k in range(2,R+1))
    print('19 width-two defects and all 27 triple/mark tests reproduce')

    for d,j in product(range(3,9),range(1,7)):assert read(P(d),j)==fast(d,j)
    assert 2*6==3*4
    assert fast(4,2)*fast(4,6)==46080 and fast(4,3)*fast(4,4)==37440
    signatures=[(i,j,k) for i in range(1,7) for j in range(i,7) for k in range(j,7)]
    d=4
    while len({prod(fast(d,z) for z in u) for u in signatures})<len(signatures):d+=2
    print('shared-block test separates all',len(signatures),'cubic multisets on slots 1..6, even block order',d)

if __name__=='__main__':run()
