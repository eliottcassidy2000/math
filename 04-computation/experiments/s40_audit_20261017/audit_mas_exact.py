"""S40 audit: exact maximum acyclic subgraph MAS(P_p), independent of the session's code.

Reductions (different from the session's 'vertex 0 first'):
  * vertex 0 LAST (translation is an automorphism);
  * for p = 31 also: the vertex just before 0 is an in-neighbour of 0 in every optimal order (else swapping the two
    gains an arc), i.e. a non-residue; multiplication by squares fixes 0 and is transitive on non-residues, so WLOG it
    is p - 1 = -1.
Pull DP: dp[S] = max_(v in S) dp[S - v] + #(in-neighbours of v inside S - v). Brute force over all orders with 0 last
for p <= 11. Natural order forward count by direct counting.
Usage: python audit_mas_exact.py [PMAX_FULL_DP=23] [--p31]
"""
import sys, os, json, time
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_common import legendre


@njit(cache=True)
def popcount(x):
    x = x - ((x >> 1) & 0x5555555555555555)
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F
    return (x * 0x0101010101010101) >> 56 & 0xFF


def beats(p):
    """B[u][v] = 1 iff u -> v in P_p (v - u a nonzero square)."""
    B = np.zeros((p, p), dtype=np.int64)
    for u in range(p):
        for v in range(p):
            if u != v and legendre(v - u, p) == 1:
                B[u, v] = 1
    return B


@njit(cache=True)
def dp_pull(inmask, m):
    """inmask[i] = bitmask (over the m free vertices) of free vertices that beat free vertex i."""
    size = 1 << m
    dp = np.zeros(size, dtype=np.int16)
    for S in range(1, size):
        best = -1
        for i in range(m):
            low = 1 << i
            if S & low:
                R = S ^ low
                val = dp[R] + popcount(inmask[i] & R)
                if val > best:
                    best = val
        dp[S] = best
    return dp[size - 1]


@njit(cache=True)
def brute_last0(B, p):
    """All (p-1)! orders of 1..p-1 followed by 0; max forward arcs."""
    a = np.arange(1, p)
    m = p - 1
    best = -1
    while True:
        f = 0
        for i in range(m):
            for j in range(i + 1, m):
                f += B[a[i], a[j]]
        for i in range(m):
            f += B[a[i], 0]
        if f > best:
            best = f
        # next permutation
        k = m - 2
        while k >= 0 and a[k] >= a[k + 1]:
            k -= 1
        if k < 0:
            break
        l = m - 1
        while a[l] <= a[k]:
            l -= 1
        a[k], a[l] = a[l], a[k]
        i, j = k + 1, m - 1
        while i < j:
            a[i], a[j] = a[j], a[i]
            i += 1; j -= 1
    return best


def mas_last0(p, B):
    free = list(range(1, p))
    m = len(free)
    inm = np.zeros(m, dtype=np.int64)
    for i, v in enumerate(free):
        msk = 0
        for j, u in enumerate(free):
            if B[u, v]:
                msk |= 1 << j
        inm[i] = msk
    into0 = int(sum(B[u, 0] for u in free))
    return int(dp_pull(inm, m)) + into0


def mas_last0_pen(p, B):
    pen = p - 1
    free = list(range(1, p - 1))
    m = len(free)
    inm = np.zeros(m, dtype=np.int64)
    for i, v in enumerate(free):
        msk = 0
        for j, u in enumerate(free):
            if B[u, v]:
                msk |= 1 << j
        inm[i] = msk
    assert B[pen, 0] == 1
    tail = int(sum(B[u, pen] for u in free)) + int(sum(B[u, 0] for u in range(1, p)))
    return int(dp_pull(inm, m)) + tail


def main():
    pmax = int(sys.argv[1]) if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else 23
    res = {}
    for p in [3, 7, 11, 19, 23]:
        if p > pmax:
            break
        B = beats(p)
        nat = int(sum(B[i, j] for i in range(p) for j in range(i + 1, p)))
        t = time.time()
        m1 = mas_last0(p, B)
        r = dict(p=p, natural_forward=nat, MAS_dp_last0=m1, seconds=round(time.time() - t, 1))
        if p <= 11:
            r['MAS_brute_last0'] = int(brute_last0(B, p))
        if p >= 7:
            r['MAS_dp_last0_pen'] = mas_last0_pen(p, B)
        res[p] = r
        print(r, flush=True)
    if '--p31' in sys.argv:
        p = 31
        B = beats(p)
        nat = int(sum(B[i, j] for i in range(p) for j in range(i + 1, p)))
        t = time.time()
        m2 = mas_last0_pen(p, B)
        r = dict(p=p, natural_forward=nat, MAS_dp_last0_pen=m2, seconds=round(time.time() - t, 1))
        res[p] = r
        print(r, flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    tag = '_p31' if '--p31' in sys.argv else ''
    json.dump(res, open(os.path.join(here, 'audit_mas_exact%s.json' % tag), 'w'), indent=1)


if __name__ == '__main__':
    main()
