#!/usr/bin/env python3
"""Construct one strong one-apex tournament for any prescribed odd congruence."""
from itertools import combinations, product
from math import gcd, isqrt


def power_order(m):
    if m == 1:
        return 1
    v, L = 2 % m, 1
    while v != 1:
        v, L = 2*v % m, L+1
    return L


def fan_count(pos, mod=None):
    def two(k):
        return 1 << k if mod is None else pow(2, k, mod)
    h = 1 + sum(two(t-1) for t in pos)
    h -= sum(two(max(t-s-1, 1)) for s, t in combinations(pos, 2))
    return h if mod is None else h % mod


def linear_fan_count(pos, m):
    R, h, prev = 0, 1, 0
    for t in pos:
        if prev:
            d = t-prev
            R = (pow(2,d,m)*R+pow(2,d-1,m)) % m
        h = (h+pow(2,t-1,m)-R-(prev == t-1)) % m
        prev = t
    return h


def marked_source(M, r, end=None, *, start=0, suffix=None):
    """Return last spine position n and ordered positions beating apex 0."""
    if M < 1 or (M % 2 == 0 and r % 2 == 0):
        raise ValueError("target not an odd class or modulus not positive")
    s = M & -M
    a, m = s.bit_length()-1, M//s
    L = power_order(m)
    D = max(2, a+1)
    g = gcd(3, m)
    q = m//g
    k = (q-1).bit_length()
    phases = 2 if g == 3 else 0
    marks = 2*max(0,a-1)+3*phases+6*k
    # Fixed last losing mark, common to every target at this modulus.
    bound = start+(marks+1)*(D+L)
    if suffix is None:
        if end is None:
            end = bound
        suffix = [end]
    elif end is not None:
        raise ValueError("supply an end or a suffix, not both")
    if not suffix or suffix[0] < bound or list(suffix) != sorted(set(suffix)):
        raise ValueError("fixed suffix does not leave room for the controls")
    end = suffix[-1]
    delta = (linear_fan_count(suffix,s)-2*phases-r) % s

    def compile(which, bits):
        ps = []
        def long(j):
            lo = (ps[-1] if ps else start)+D
            ps.append(lo+(j-lo) % L)
        def short(d):
            ps.append(ps[-1]+d)
        for i in range(1, a):
            long(0)
            if delta >> i & 1:
                short(i+1)
            else:
                long(ps[-1]+i+1)
        for j in range(phases):
            long(0)
            if j < which:
                short(1)
                long(2)
            else:
                short(2)
                long(1)
        def pair(d, swap):
            long(d if swap else 0)
            long(0 if swap else d)
        for d in range(k):
            active = (bits>>(k-1-d))&1
            pair(d+1, active)
            pair(d, not active)
            pair(d, not active)
        assert suffix[0]-(ps[-1] if ps else start) >= D
        return ps+list(suffix)

    which = ((linear_fan_count(compile(0,0),m)-r) % 3) if phases else 0
    h0 = linear_fan_count(compile(which,0), m)
    step = 3*pow(2,-k-1,m)%m if k else 0
    assert (r-h0) % g == 0
    I = ((r-h0)//g * pow(step//g,-1,q)) % q if q > 1 else 0
    ps = compile(which, I)
    assert linear_fan_count(ps,M) == r % M
    return end, ps


def approximate_with_residue(N, M, r):
    n=N.bit_length()
    start=n//2
    C,_=marked_source(M,r)
    u=start+C
    if u>=n:
        raise ValueError("target too small for leading-bit reservation")
    suffix=[i+1 for i in range(u,n) if N>>i&1]
    return marked_source(M,r,start=start,suffix=suffix)


def rows(n, pos):
    """Fan on vertices 0,1,...,n, outgoing bit rows."""
    rt = [((1<<(n+1))-1)^((1<<(i+1))-1) for i in range(n+1)]
    for t in pos:
        rt[0] ^= 1<<t
        rt[t] ^= 1
    return rt


def count(rt):
    n = len(rt)
    dp = [[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i] = 1
    for mask, row in enumerate(dp):
        for u, ways in enumerate(row):
            if not ways:
                continue
            left = rt[u] & ~mask
            while left:
                b = left & -left
                left ^= b
                dp[mask|b][b.bit_length()-1] += ways
    return sum(dp[-1])


def subgraph(rt, vs):
    return [sum(((rt[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs]


def word_count(s):
    if not s:
        return 1
    h,x,prev = 1,0,0
    for b in s:
        if b=='0':
            h += x+prev
        prev=int(b)
        x=2*x+prev
    return h


def next_deviations(s):
    h=word_count(s)
    return {word_count(s[:j]+b+s[j:])-2*h for j in range(len(s)+1) for b in '01'} & {-1,1}


def fork_check():
    large=single=none=0
    for n in range(2,14):
        for w in product('01',repeat=n-2):
            s='1'+''.join(w)+'0'
            h=word_count(s)
            for j in range(1,n):
                U,V=s[:j],s[j:]
                X,Z=int(U,2),int(''.join('1' if b=='0' else '0' for b in V[::-1]),2)
                K=word_count(U)+word_count(V)-1+(U[-1]+V[0]=='10')
                e0,e1=U[-1]+V[0]=='11',U[-1]+V[0]=='00'
                w0,w1=U+'0'+V,U+'1'+V
                s0,s1=word_count(w0)-2*h,word_count(w1)-2*h
                assert h==K+X*Z and (s0,s1)==(X+e0-K,Z+e1-K)
                if {s0,s1}!={-1,1}:
                    continue
                assert h in (K*K,K*K-2,K*K+K-1)
                if h<63:
                    continue
                large+=1
                for child,dev in ((w0,s0-2*e0),(w1,s1-2*e1)):
                    expected={dev}&{-1,1}
                    assert next_deviations(child)==expected
                    single+=bool(expected); none+=not expected
    print('Large fork presentations, next-single, next-none:',large,single,none)
    assert next_deviations('101000010')=={-1} and next_deviations('101010010')=={-1}


def main():
    for n in range(2, 9):
        for bits in product((0,1), repeat=n-1):
            ps = [i for i,b in enumerate(bits,2) if b]
            assert fan_count(ps) == count(rows(n,ps))
    total=0
    for M in range(1,141):
        signature = None
        for r in range(M):
            if not M%2 and not r%2: continue
            n,ps=marked_source(M,r)
            total+=1
            assert ps[-1]==n and ps[0]>=2 and fan_count(ps,M)==r%M
            assert marked_source(M,r,end=n+3)[0]==n+3
            if signature is None: signature=(n,len(ps))
            assert signature==(n,len(ps))
    for M in (97,301,2048,720,1024*27,8191,3**6*5**3,7**4,81,2**6*9991):
        for r in (7,21,63,1234567):
            n,ps=marked_source(M,r)
            assert fan_count(ps,M)==r%M
    n,ps=marked_source(8191,7)
    assert (n,len(ps))==(1185,79)
    print('Full residue replays:',total,'; modulus 8191:',n,'spine positions with',len(ps),'marks.')
    fork_check()


if __name__ == '__main__':
    main()
