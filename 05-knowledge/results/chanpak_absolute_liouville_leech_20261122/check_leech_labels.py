"""Exact profile and literal Hamilton DP for the two missing Leech labels."""
from functools import lru_cache
from itertools import chain


def require(ok, context):
    if not ok:
        raise RuntimeError(context)


def from_string(s):
    return tuple(map(int, s.split(",")))


rows = (
    (9, 1, "136,17,3,22,101,15,47,126", 81, (277,1152,276)),
    (9, 7, "4,1,10,67,15,31,55", 96, (29,63,22)),
    (25,1, "2,48,3,7,13,29", 12, (23,100,22)),
    (25,7, "32,1,3,7,15,158,63,95", 144, (51,175,44)),
)
starts = {
    19:"2,40,3,5,15,29",
    37:"32,17,3,7,13,30",
    51:"64,9,3,5,15,31,62",
    73:"4,65,34,7,15,27,61",
    101:"2,112,67,7,13,29,57",
    145:"66,48,67,7,13,29,58",
    289:"192,5,49,7,11,27,62,126",
    577:"130,196,97,135,15,27,57,116",
    1153:"144,321,131,7,14,31,61,122,253",
}


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length()-1
        mask ^= low


def paths(A):
    n = len(A)
    require(all((A[i]>>j&1)+(A[j]>>i&1)==(i!=j) for i in range(n) for j in range(n)), A)
    incoming = tuple(sum((1<<i) for i in range(n) if A[i]>>j&1) for j in range(n))

    @lru_cache(None)
    def end(s,j):
        t = s ^ (1<<j)
        return 1 if not t else sum(end(t,i) for i in bits(t&incoming[j]))

    @lru_cache(None)
    def start(s,j):
        t = s ^ (1<<j)
        return 1 if not t else sum(start(t,i) for i in bits(t&A[j]))

    @lru_cache(None)
    def count(s):
        return 1 if not s else sum(end(s,j) for j in bits(s))

    return (1<<n)-1, end,start,count


def reaches(A, target):
    seen = target
    for _ in A:
        for i,adj in enumerate(A):
            if adj&seen:
                seen |= 1<<i
    return seen == (1<<len(A))-1


def strong(A):
    return all(reaches(A,1<<i) for i in range(len(A)))


def assemble(D,I,X):
    n,w=len(D),len(X)
    b=n+w
    dmask=(1<<n)-1
    xmask=((1<<w)-1)<<n
    return tuple(chain(
        (D[i] | (((I>>i)&1)<<b) for i in range(n)),
        ((a<<n) | dmask for a in X),
        (xmask | (dmask^I),),
    ))


def tt(w):
    return tuple(((1<<w)-1)^((1<<(i+1))-1) for i in range(w))


for a,scale,masks,I,profile in rows:
    D = from_string(masks)
    whole, end,start,H=paths(D)
    G=DII=0
    for S in range(1,whole+1):
        T=whole^S
        left=sum(end(S,i) for i in bits(S&I))
        G += left * H(T)
        DII += left * sum(start(T,j) for j in bits(T&I))
    d=H(whole)
    require((d,G,DII)==profile and d-DII==scale and reaches(D,I), (a,scale,d,G,DII))
    for X in (tt(1),tt(2),(2,4,1)):
        U,_,_,hx=paths(X)
        cover=sum(hx(s)*hx(U^s) for s in range(U+1))
        Y=assemble(D,I,X)
        W,_,_,hy=paths(Y)
        require(strong(Y) and hy(W)==G*cover+scale*hx(U), (a,scale,X))

for h,masks in starts.items():
    Y=from_string(masks)
    W,_,_,H=paths(Y)
    require(strong(Y) and H(W)==h,(h,H(W)))

print("4 new profiles, 12 strong whole-host replays, 9 strong finite starts: OK")
