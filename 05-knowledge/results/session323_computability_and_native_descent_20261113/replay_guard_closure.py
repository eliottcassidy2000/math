from itertools import combinations, product
from math import lcm, prod

def ck(ok,msg):
    if not ok: raise AssertionError(msg)

def maskset(mask,n):
    return {i for i in range(n) if mask>>i&1}

def test_closures():
    n=3
    top=2**n-1
    primes=(2,3,5,7,11,13,17,19)
    count=0
    for bits in range(1<<top):
        C={top}|maskset(bits,top)
        if any((a&b) not in C for a in C for b in C):
            continue
        d=[prod(primes[c] for c in C if not c>>i&1) for i in range(n)]
        achieved=set()
        for S in range(top+1):
            N=lcm(*(d[i] for i in maskset(S,n)))
            act=sum(1<<i for i in range(n) if N%d[i]==0)
            meet=top
            for c in C:
                if S&c==S: meet&=c
            ck(act==meet,(C,S))
            achieved.add(act)
        ck(achieved==C,C)
        count+=1
    print("closure systems on three guards:",count)

def source(S,R):
    return next((i for i in S if all(i==j or (i,j) in R for j in S)),None)

def test_shields():
    total=0
    for n in range(1,6):
        pairs=list(combinations(range(n),2))
        for bits in product(range(2),repeat=len(pairs)):
            R={e if b else e[::-1] for e,b in zip(pairs,bits)}
            triangles=[frozenset(c) for c in combinations(range(n),3)
                       if source(c,R) is None]
            for m in range(1,1<<n):
                S=maskset(m,n)
                hits=[t for t in triangles if t<=S]
                ck(bool(hits) or source(S,R) is not None,(R,S))
                if len(S)==3 and source(S,R) is None:
                    ck(hits==[frozenset(S)],(R,S))
            total+=1
    print("labelled original priority tournaments through order 5:",total)

def test_cycle():
    R={(0,1),(1,2),(2,0),(3,0),(3,1),(3,2)}
    D=(2,3,5,30)
    winners=set()
    for n in range(1,31):
        S={i for i,d in enumerate(D) if n%d==0}
        w=source(S,R)
        ck(not S or w is not None,n)
        if w is not None:
            winners.add(w)
            ck(n//D[w]<n,n)
    ck(winners==set(range(4)),winners)
    exposed={(source({i,j},R),j if source({i,j},R)==i else i)
             for i,j in combinations(range(3),2)}
    ck(exposed=={(0,1),(1,2),(2,0)},exposed)
    print("cyclic guarded selector: four reachable winners, total decreasing updates")

if __name__=='__main__':
    test_closures();test_shields();test_cycle()
