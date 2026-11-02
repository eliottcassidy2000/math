#!/usr/bin/env python3
"""Independent literal-word/cycle count and Bareiss determinant referee."""
from fractions import Fraction as F
from itertools import combinations, permutations
from random import Random
from hashlib import sha256
from pathlib import Path

def require(condition, detail='check failed'):
    if not condition:
        raise AssertionError(detail)

def det(a):
    a = [r[:] for r in a]
    n = len(a)
    if not n:
        return 1
    sign = 1
    old = 1
    for k in range(n - 1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[k], a[p] = (a[p], a[k])
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                z = a[i][j] * pivot - a[i][k] * a[k][j]
                require(z % old == 0)
                a[i][j] = z // old
        for i in range(k + 1, n):
            a[i][k] = 0
        old = pivot
    return sign * a[-1][-1]

def measure(b, q=1):
    n = len(b)
    full = (1 << n) - 1
    stats = {}
    cycles = []
    for mask in range(1 << n):
        vs = [i for i in range(n) if mask >> i & 1]
        k = len(vs)
        hnum = 0
        for p in permutations(vs):
            w = 1
            for i, j in zip(p, p[1:]):
                w *= b[i][j]
            hnum += w
        h = F(hnum, q ** max(0, k - 1))
        a = [[b[i][j] for j in vs] for i in vs]
        u = F(det([[a[i][j] + q * (i == j) for j in range(k)] for i in range(k)]), q ** k)
        v = F((-1) ** k * det(a), q ** k)
        stats[mask] = (h, u, v)
        if k >= 2:
            for tail in permutations(vs[1:]):
                p = (vs[0],) + tail
                w = 1
                for i, j in zip(p, p[1:] + p[:1]):
                    w *= b[i][j]
                if w:
                    cycles.append((mask, k, F(w, q ** k)))
    h, u, v = stats[full]
    if n:
        require(h >= u + abs(v), (b, h, u, v))
    for mask, (hh, uu, vv) in stats.items():
        require(hh >= uu)
        if mask:
            require(uu + vv > 0 and uu - vv > 0)
        else:
            require((hh, uu, vv) == (1, 1, 1))
    for x in range(n):
        old = full ^ 1 << x
        hh, uu, vv = stats[old]
        pu = uu
        pv = F(0)
        ph = hh
        for cm, k, w in cycles:
            if not cm >> x & 1:
                continue
            rh, ru, rv = stats[full ^ cm]
            pu += (-1) ** (k - 1) * w * ru
            pv -= w * rv
            if k % 2:
                ph += 2 * w * rh
        require((ph, pu, pv) == (h, u, v), (n, x, (ph, pu, pv), (h, u, v)))
    odd = [c for c in cycles if c[1] % 2]

    def ocf(start, used):
        out = F(1)
        for t in range(start, len(odd)):
            cm, k, w = odd[t]
            if not cm & used:
                out += 2 * w * ocf(t + 1, used | cm)
        return out
    require(ocf(0, 0) == h)
    if q == 1:
        call = sum((w for _, _, w in cycles))
        cnon = sum((w for cm, _, w in cycles if cm != full))
        tris = [cm for cm, k, w in cycles if k == 3]
        pair = any((not a & b for i, a in enumerate(tris) for b in tris[i + 1:]))
        require(h - u >= call)
        if n:
            gap = h - u - abs(v) - cnon
            require(gap >= 0)
            require((gap == 0) == (not pair))
            if pair:
                require(gap >= 2)
        nontrans = bool(tris)
        if n:
            require((h == u + abs(v)) == (not nontrans or n == 3))
            require((h == u + v) == (not nontrans))
    if q != 1 and n:
        fractional = any((b[i][j] not in (0, q) for i, j in combinations(range(n), 2)))
        if n <= 2:
            require(h == u + abs(v))
        elif fractional:
            require(h > u + abs(v))
        if fractional:
            require(h > u + v)
    return (h, u, v)

def tournament(n, mask):
    b = [[0] * n for _ in range(n)]
    for t, (i, j) in enumerate(combinations(range(n), 2)):
        b[i][j] = mask >> t & 1
        b[j][i] = 1 - b[i][j]
    return b

def main():
    count = 0
    for n in range(6):
        for m in range(1 << n * (n - 1) // 2):
            measure(tournament(n, m))
            count += 1
        print('literal exhaustive unweighted n=%d passed (%d cumulative)' % (n, count), flush=True)
    rng = Random(137902)
    for n, nt in [(6, 100), (7, 20), (8, 3)]:
        for _ in range(nt):
            measure(tournament(n, rng.randrange(1 << n * (n - 1) // 2)))
        print('literal random unweighted n=%d count=%d passed' % (n, nt), flush=True)
    for n in [6, 7, 8]:
        b = tournament(n, 0)
        for group in [(0, 1, 2), (3, 4, 5)]:
            for i, j in zip(group, group[1:] + group[:1]):
                b[i][j] = 1
                b[j][i] = 0
        print('two C3 hostile n=%d (H,u,v)=%s' % (n, measure(b)), flush=True)
    weighted = 0
    for n in range(1, 7):
        for _ in range(30):
            q = 4
            b = [[0] * n for _ in range(n)]
            for i, j in combinations(range(n), 2):
                b[i][j] = rng.randrange(q + 1)
                b[j][i] = q - b[i][j]
            measure(b, q)
            weighted += 1
        print('exact rational weighted n=%d passed (%d cumulative)' % (n, weighted), flush=True)
    print('PASS all identities, bounds, equality classifications, empty bases; no floating point', flush=True)
    print('script_sha256', sha256(Path(__file__).read_bytes()).hexdigest(), flush=True)
if __name__ == '__main__':
    main()
