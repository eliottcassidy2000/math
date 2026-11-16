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


def inv_word(x, y, I):
    assert 0 <= I <= x*y
    word = []
    while x+y:
        if y and I >= x:
            word.append(1)
            I -= x
            y -= 1
        else:
            word.append(0)
            x -= 1
    assert I == 0
    return word


def marked_source(M, r, end=None, *, start=0, suffix=None):
    """Return last spine position n and ordered positions beating apex 0."""
    if M < 1 or (M % 2 == 0 and r % 2 == 0):
        raise ValueError("target not an odd class or modulus not positive")
    s = M & -M
    a, m = s.bit_length()-1, M//s
    L = power_order(m)
    D = max(2, a+1)
    g = gcd(3, m)
    x = isqrt((m//g)-1)+1 if m > 1 else 0
    phases = 2 if g == 3 else 0
    marks = 2*max(0,a-1)+3*phases+2*x
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

    def compile(which, word):
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
        for bit in word:
            long(bit)
        assert suffix[0]-(ps[-1] if ps else start) >= D
        return ps+list(suffix)

    word = [0]*x+[1]*x
    which = ((linear_fan_count(compile(0,word),m)-r) % 3) if phases else 0
    h0 = linear_fan_count(compile(which,word), m)
    step = 3*pow(4,-1,m) % m if m > 1 else 0
    assert (r-h0) % g == 0
    q = m//g
    I = ((r-h0)//g * pow(step//g,-1,q)) % q if q > 1 else 0
    ps = compile(which, inv_word(x,x,I))
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


def main():
    for n in range(2, 9):
        for bits in product((0,1), repeat=n-1):
            ps = [i for i,b in enumerate(bits,2) if b]
            assert fan_count(ps) == count(rows(n,ps))
    for M in range(1, 151):
        signature = None
        for r in range(M):
            if not M%2 and not r%2:
                continue
            n, ps = marked_source(M,r)
            assert ps[-1] == n and ps[0]>=2
            sig = (n,len(ps))
            assert fan_count(ps,M) == r%M
            assert marked_source(M,r,end=n+1)[0] == n+1
            if signature is None:
                signature = sig
            assert sig == signature
    assert marked_source(2**7*3**3*5**2*11**2*13, 7)
    assert marked_source(2**7*3**3*5**2*11**2*13, 21)
    for M in (1,2,7,16,33,72):
        r=7
        C,_=marked_source(M,r)
        for j in (3,9,17):
            N=(1<<(2*C+j))+314159265358979323846
            _,ps=approximate_with_residue(N,M,r)
            H=fan_count(ps)
            assert H%M==r%M
            assert abs(H-N)<(2**(C+5))*isqrt(N)
    # Literal child, and 63's induced spectra using the THM-344 two-class census.
    parent, child, fan63 = rows(6,[2,4,6]), rows(7,[2,5,7]), rows(7,[3,4,7])
    assert [count(r) for r in (parent,child,fan63)] == [31,61,63]
    assert subgraph(child,[0,1,2,3,5,6,7]) == parent
    assert sorted(map(int.bit_count,fan63)) == [1,1,2,4,4,5,5,6]
    for k in range(1,8):
        vals = sorted({count(subgraph(fan63,vs)) for vs in combinations(range(8),k)})
        assert 31 not in vals
        print("63-fan induced",k,vals)
    assert count(subgraph(fan63,[0,1,2,4,5,6,7])) == 33
    prior = subgraph(fan63,[0,1,3,2,4,6,7])
    prior[2] ^= 1<<3
    prior[3] ^= 1<<2
    assert prior == parent
    mixed = {
        497:[224,213,9,179,69,22,172,52],
        511:[166,244,88,35,169,68,25,108],
        553:[40,197,41,50,71,210,141,29],
        581:[204,49,138,178,101,5,46,114],
    }
    for h, r in mixed.items():
        assert all((r[i]>>j&1)+(r[j]>>i&1)==1 for i,j in combinations(range(8),2))
        assert all(not (v>>i&1) and v < 256 for i,v in enumerate(r))
        assert count(r)==h
    for k in range(4,100):
        a,b = k//2, k+1-k//2
        h = fan_count([a,b,a+b])
        assert h == (1<<k)-1
        assert fan_count([a,b+1,a+b+1]) == 2*h+(1 if k%2==0 else -1)
        assert fan_count([a+1,b+1,a+b+2]) == 4*h+3
    def win_word(s):
        assert s[0]=='1'
        return fan_count([i for i,b in enumerate(s,1) if b=='0'])
    for m in range(1,70):
        A,B='1'+'0'*m+'1','1'*m+'0'
        h=(2**(m+2)+3)**2
        assert [win_word(A+x+B) for x in ['000','0000','0100']]==[h,2*h-1,2*h+1]
    from itertools import permutations
    circ = []
    for t in product((1,-1),repeat=3):
        ds = {k*s%7 for k,s in enumerate(t,1)}
        circ.append(sum(all((v-u)%7 in ds for u,v in zip(p,p[1:]))
                        for p in permutations(range(7))))
    assert sorted(circ) == [175]*6+[189]*2
    print("All modular families, induced counts and sevenfold symmetry checks passed.")


if __name__ == '__main__':
    main()
