"""
opus-2026-07-19-S403 (H5 / HYP-7980): THE SINGER POLE -- M of the PG(2,13)
difference set as a speed family (boxeph-S110's unworked metagraph-transport
lead, quantitative first data).

183 = 14^2 - 14 + 1 = |PG(2,13)|.  The deep well {1..12,182} (M = 14/183,
the covering-min extremal) is the 'transitive pole' mod 183; the Singer
perfect difference set (14 residues, every nonzero difference exactly once)
is the 'regular pole'.  Hypothesis: Singer speed families are conspicuously
LOOSE (M large) -- the anti-extremal twin.

Construction: GF(13^3) = F13[x]/(x^3 - 2)  (2 is a non-cube mod 13).
Frobenius: x -> 3x (since x^13 = (x^3)^4 * x = 2^4 x = 16x = 3x), so
Tr(a + bx + cx^2) = 3a: zero-trace <=> constant coefficient 0.
Find primitive g (order 2196); D = {i mod 183 : g^i has zero constant coeff}
(well-defined: zero-const is F13*-invariant; |D| = 14).  Verify perfect
difference set.  Lift to integer speeds, compute exact M (q <= 400), compare
with deep well 14/183, AP 1/14, and random 13-subsets of [1,182].
"""
import random
from math import gcd
from fractions import Fraction

P = 13
def mul(u, v):
    # (a,b,c) * (d,e,f) mod x^3 = 2
    a,b,c = u; d,e,f = v
    # product coeffs before reduction: x^0..x^4
    c0 = a*d
    c1 = a*e + b*d
    c2 = a*f + b*e + c*d
    c3 = b*f + c*e
    c4 = c*f
    # x^3 = 2, x^4 = 2x
    return ((c0 + 2*c3) % P, (c1 + 2*c4) % P, c2 % P)

def powe(u, n):
    r = (1,0,0)
    while n:
        if n & 1: r = mul(r, u)
        u = mul(u, u); n >>= 1
    return r

def order_ok(u):
    # primitive iff u^(2196/p) != 1 for p in {2,3,61}
    for p in (2,3,61):
        if powe(u, 2196//p) == (1,0,0): return False
    return True

def exact_max(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q:
                bg, bq, wit = m, q, (a,q)
    return Fraction(bg,bq), wit

if __name__ == "__main__":
    random.seed(183)
    g = None
    for b in range(1, P):
        cand = (b,1,0)  # x + b
        if order_ok(cand): g = cand; break
    assert g, "no primitive found among x+b"
    print(f"primitive element: x + {g[0]}")
    D = set()
    u = (1,0,0)
    for i in range(2196):
        if u[0] == 0:
            D.add(i % 183)
        u = mul(u, g)
    D = sorted(D)
    print(f"|D mod 183| = {len(D)}; D = {D}")
    # verify perfect difference set (lambda = 1)
    from collections import Counter
    diffs = Counter((a-b) % 183 for a in D for b in D if a != b)
    ok = len(diffs) == 182 and all(v == 1 for v in diffs.values())
    print(f"perfect difference set (every nonzero diff exactly once): {ok}")
    # speeds: nonzero residues of D (drop 0 if present)
    speeds = [d for d in D if d != 0]
    print(f"speeds ({len(speeds)}): {speeds}")
    M, wit = exact_max(speeds, 400)
    print(f"M(Singer) = {M} = {float(M):.5f} at t* = {wit[0]}/{wit[1]}")
    print(f"compare: AP 1/14 = {1/14:.5f}; deep well 14/183 = {14/183:.5f}")
    # random controls: 13-subsets of [1,182]
    vals = []
    for _ in range(25):
        Vr = random.sample(range(1,183), 13)
        m,_ = exact_max(Vr, 200)
        vals.append(m)
    vals.sort()
    print(f"25 random 13-subsets of [1,182]: min {vals[0]} = {float(vals[0]):.4f}, "
          f"median {vals[12]} = {float(vals[12]):.4f}, max {vals[-1]} = {float(vals[-1]):.4f}")
    # extra: the deep well itself as sanity
    Mdw, witdw = exact_max(list(range(1,13))+[182], 400)
    print(f"sanity deep well {{1..12,182}}: M = {Mdw} at {witdw} (expect 14/183)")
