"""Exact dihedral tournament heights and target decoding. Stdlib only."""
from math import comb,isqrt
from functools import lru_cache
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import TT,inflate,vectors,read
from replay import circ,dp,profile,reads,wrap,pad

C=circ(7,(1,2,3))
DIA=(20,13,8,17,6)
AFF=(16,1,3,7,14)
LOW=((190,124,120,240,96,64,129,54),19)
P31=((80,41,67,5,14,29,58),29)
C11=(16,17,3,7,12)

def count(T):
    return sum(dp(T)[-1])

def polynomial(Q,marks):
    C7=[0]+[read(C,k)//7 for k in range(1,len(Q)+1)]
    V=vectors(Q,marks,1,len(Q)-len(marks)+1)
    return {e:sum(c*sum((-1)**(r-j)*comb(r,j)*C7[j] for j in range(r+1))
                  for (r,*s),c in V.items() if tuple(s)==e)
            for e in {i[1:] for i in V}}

def decode_width(p):
    a,b,c=4573,7392,2919
    D=b*b+4*c*(p-a)
    if D<0:return None
    d=isqrt(D)
    if d*d!=D:return None
    v,r=divmod(d-b,2*c)
    if r or v<0:return None
    w=v+2
    return w.bit_length()-1 if w&(w-1)==0 else None

def make(k):
    U=TT(k); W=tuple(((1<<i)-1) for i in range(k)); n=7+2*k
    T=inflate(DIA,{0:C,1:U,2:W})
    sigma=tuple((-i)%7 for i in range(7))+tuple(range(7+k,7+2*k))+tuple(range(7,7+k))+(n+1,n)
    rot=tuple((i+1)%7 for i in range(7))+tuple(range(7,n+2))
    return T,sigma,rot

def fold(T,sigma):
    n=len(T); c,=(i for i in range(n) if sigma[i]==i)
    K=tuple((i,sigma[i]) for i in range(n) if i<sigma[i])
    @lru_cache(None)
    def cont(last,m):
        if not m:return int(bool(T[last]>>c&1))
        ways=0
        for i,Q in enumerate(K):
            if not m>>i&1:continue
            for y in Q:
                if last==n or T[last]>>y&1:ways+=cont(y,m^(1<<i))
        return ways
    return cont(n,(1<<len(K))-1)

def fixed_profile():
    # Each copy pair is represented by the marks 7 and 8 of make(1).
    A,s,_=make(1); c=0; marks=(7,8); cap=10
    other=tuple((i,s[i]) for i in range(len(A)) if i<s[i] and i!=7)
    from itertools import product
    from collections import Counter
    result=Counter()
    for fixed in product(*other):
        @lru_cache(None)
        def dp(m,prev,t):
            z=Counter({t:1}) if not m and prev is not None and A[prev]>>c&1 else Counter()
            for i in m:
                if prev is None or A[prev]>>i&1:
                    z.update(dp(tuple(j for j in m if j!=i),i,t))
            if t<cap:
                for i in marks:
                    if prev is None or A[prev]>>i&1:z.update(dp(m,i,t+1))
            return z
        result.update(dp(fixed,None,0))
    return result

def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def run():
    assert polynomial(DIA,(0,1,2))=={(1,1):4573,(1,2):3696,(2,1):3696,(2,2):2919}
    weights=fixed_profile(); print('fixed run coefficients',sorted(weights.items()))
    for p,k in ((31033,2),(47278009,7)):
        assert prime(p) and decode_width(p)==k
        T,s,r=make(k)
        for i in range(len(T)):
            for j in range(len(T)):
                assert T[i]>>j&1==T[r[i]]>>r[j]&1==T[s[j]]>>s[i]&1
            assert s[r[s[i]]]==(r.index(i))
        z=fold(T,s)
        assert z==sum(t*sum((-1)**(j-a)*comb(j,a)*a**k for a in range(j+1)) for j,t in weights.items())
        assert (p-z)%2==0
        if k==2:
            assert count(T)==7*p
            assert count(wrap(*P31,TT(5)))==p
        print('dihedral',p,k,7*p,'fixed normalized',z,'double orbits',(p-z)//2)
    assert polynomial(AFF,(0,1))=={(1,):1863,(2,):1036}
    d,g,e=reads(*LOW)
    assert (g,d+e-g)==(1036,-209)
    assert profile(C11)==(11,100)
    for U,k,p in ((TT(1),7,265007),(C11,7,13258501)):
        h,b=profile(U); assert p==1036*b*2**k-209*h and prime(p)
        v,r=divmod(p+209*h,1036*b)
        assert r==0 and v==1<<k
        S=pad(U,k)
        if p<300000:
            assert count(wrap(*LOW,S))==p
            assert count(inflate(AFF,{0:C,1:S}))==7*p
        print('paired affine',p,7*p,'width',k)
if __name__=='__main__':run()
