"""Exact two-contact bridge, targeted near-boundary seeds, and inverse."""
from functools import cache
from itertools import combinations
from math import isqrt
from bridge import table,bridge,rooted_profile
from physical_prefix import row,greedy,prime
from physical_seed import decode
from construct_target import actual,strong,enclosed_prefix,target_line

def profile(T,U):
    U=frozenset(U)
    n=len(T); full=(1<<n)-1
    h=table(T)
    @cache
    def endpoints(S,a,z):
        if a==z:return int(S==(1<<a))
        X=S^(1<<z)
        total=0
        for j in range(n):
            if X>>j&1 and T[j]>>z&1:total+=endpoints(X,a,j)
        return total
    R=E=0
    for S in range(1,1<<n):
        for a in U:
            if not S>>a&1:continue
            for z in range(n):
                if not S>>z&1:continue
                w=endpoints(S,a,z)*h[full^S]
                R+=w
                if z in U and not S.bit_count()&1:E+=w
    return h[-1],R,E

def multi(T,U,r):
    assert r>=1
    n=len(T); m=n+r+1
    assert all(u in range(n) for u in U)
    later=((1<<m)-1)^((1<<n)-1)
    out=[a|later for a in T]
    for u in U:out[u]^=1<<(m-1)
    for j in range(n,m-1):out.append(((1<<m)-1)^((1<<(j+1))-1))
    out.append(sum(1<<u for u in U))
    return out

def power2(x):
    return x.bit_length()-1 if x>=2 and not x&(x-1) else None

def inverse_pair(a,b,record):
    h,R,E=record
    S=a+2*E
    if R<=0 or S<=h or (S-h)%R or b%S:return None
    r,s=power2((S-h)//R),power2(b//S)
    return None if r is None or s is None else (r,s)

def inverse_total(N,record):
    h,R,E=record
    S=N+2*E
    if S<=0 or R<=0:return ()
    solutions=set()
    for d in range(1,isqrt(S)+1):
        if S%d:continue
        for e in (d,S//d):
            y=e-1
            s=power2(y)
            x=S//e-h
            r=power2(x//R) if x>0 and not x%R else None
            if r is not None and s is not None:solutions.add((r,s))
    return tuple(sorted(solutions))

def ready(T,U,r,s):
    W=multi(T,U,r)
    return bridge(W,len(W)-1,s)

def near_source(k,record):
    h,R,E=record
    u,v=row(k);A=u+2*v
    r=1
    while h+(R<<r)-2*E<A:r+=1
    a=h+(R<<r)-2*E
    assert A<=a<2*A
    b=2*(a+2*E)
    return u*a+v*b,(a,b),r

def solve_line(p,bank):
    z=target_line(p)
    if z is None:return None
    k,points=z
    known=[(T,U,profile(T,U)) for T,U in bank]
    for a,b in points:
        for T,U,D in known:
            width=inverse_pair(a,b,D)
            if width is not None:
                P=ready(T,U,*width)
                return k,(a,b),tuple(P),tuple(enclosed_prefix(P,k))
    return None

if __name__=='__main__':
    T=(24,17,3,6,12)
    U,V=(0,2),(0,1)
    assert profile(T,U)==(15,42,4) and profile(T,V)==(15,42,5)
    family=[(T,U,4,1,2053),(T,V,4,1,2051),
            ((20,1,2,7,14),(0,4),2,2,679),
            ((20,1,2,7,14),(1,3),2,2,677)]
    for G,C,r,s,N in family:
        P=ready(G,C,r,s)
        assert actual(P) and strong(P) and table(P)[-1]==N
        h,R,E=profile(G,C); S=h+(R<<r)
        assert (S<<s)+S-2*E==N
        assert (r,s) in inverse_total(N,(h,R,E))
        assert inverse_pair(S-2*E,S<<s,(h,R,E))==(r,s)
    assert prime(2053) and prime(293) and prime(677) and prime(97)
    p,z,r=near_source(1,profile(T,V))
    assert (p,z,r)==(1229903,(1349,2718),5) and prime(p)
    assert greedy(p)==(1,77,3802) and decode(*z) is None
    P=ready(T,V,r,1)
    W=enclosed_prefix(P,1)
    assert rooted_profile(P,len(P)-1)==(4067,4067,-1349)
    assert solve_line(p,((T,U),(T,V)))==(1,z,tuple(P),tuple(W))
    assert actual(W) and strong(W) and p*7==8609321
    print('literal target source',p,z,tuple(P),len(W))
    # All possible old tournaments through four vertices and two-contact sets.
    for n in range(2,5):
        pairs=tuple(combinations(range(n),2))
        for bits in range(1<<len(pairs)):
            T=[0]*n
            for k,(i,j) in enumerate(pairs):
                if bits>>k&1:T[i]|=1<<j
                else:T[j]|=1<<i
            for U in pairs:
                h,R,E=profile(T,U)
                v=multi(T,U,1)
                assert rooted_profile(v,len(v)-1)[:2]==(h+2*R-2*E,h+2*R)
    print('two-contact transfer, bivariate target inverse, and both prime/7p adjacent pairs passed')
