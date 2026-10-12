"""Audit B (2026-10-12) of THM-4623(2): count the maps h : anc_D(u) u {u} -> N with h(u) = v and
h(T n) = T(h n), on the depth-D truncation of the ancestor tree of u = 2m, v = (m-1)/3, m = 4 mod 6.
THM-4623 claims 2^{aleph_0} such h for EVERY m = 4 (mod 6), m != 4.  Counterexample family: when
v = 0 (mod 3) (i.e. m = 10 mod 18) every value of h is a multiple of 3 and T^{-1}(3j) = {6j}, so h is
unique.  Also verifies the explicit construction of THM-4623(2) (values along the doubling chain plus
h(n) = 2^{depth(n)-k} t_k off the chain) is an endomorphism on a finite window.
"""
import sys
from functools import lru_cache
sys.setrecursionlimit(10000)

def T(n): return n // 2 if n % 2 == 0 else 3 * n + 1
def preimages(n):
    out = [2 * n]
    if n % 6 == 4 and (n - 1) // 3 >= 1: out.append((n - 1) // 3)
    return out

@lru_cache(maxsize=None)
def count(n, t, d):
    """number of h on the depth-d truncated ancestor tree of n with h(n) = t."""
    if d == 0: return 1
    total = 1
    for c in preimages(n):
        s = sum(count(c, tp, d - 1) for tp in preimages(t))
        total *= s
    return total

D = 22
print("m     u     v    v mod 3   #h on depth-%d truncation" % D)
for m in [10, 16, 22, 28, 34, 40, 46, 52, 58, 64, 70]:
    assert m % 6 == 4 and m != 4
    u, v = 2 * m, (m - 1) // 3
    print("%-5d %-5d %-4d %-9d %d" % (m, u, v, v % 3, count(u, v, D)))

# the construction of THM-4623(2) for m = 16: t_0 = v, odd preimage chosen whenever it is != 0 mod 3
def build_h(m, K):
    u, v = 2 * m, (m - 1) // 3
    t = [v]
    for k in range(K):
        tk = t[-1]
        if tk % 6 == 4 and ((tk - 1) // 3) % 3 != 0 and (tk - 1) // 3 >= 1: t.append((tk - 1) // 3)
        else: t.append(2 * tk)
    return u, v, t

def check_endomorphism(m, K, W):
    u, v, t = build_h(m, K)
    chain = {u * 2 ** k: k for k in range(K + 1)}
    # phi on n <= W: follow the forward orbit until it hits the chain or leaves anc(u)
    def phi(n):
        path = 0; x = n
        while x not in chain:
            if x == 1 or x > 2 ** 60: return n     # not an ancestor of u (orbit reaches 1 first, or escapes window)
            if x < n and x not in chain and False: pass
            x = T(x); path += 1
            if path > 2000: return n
        k = chain[x]
        if k == K and path > 0: return None       # beyond the truncation: h undefined
        return 2 ** path * t[k]
    bad = 0; tested = 0
    for n in range(1, W + 1):
        a, b = phi(n), phi(T(n))
        if a is None or b is None: continue
        tested += 1
        if T(a) != b: bad += 1
    return tested, bad

for m in (16, 22):
    tested, bad = check_endomorphism(m, 40, 200000)
    print("m=%d: construction checked at %d vertices n <= 200000 (phi defined at n and T n); violations: %d" % (m, tested, bad))
