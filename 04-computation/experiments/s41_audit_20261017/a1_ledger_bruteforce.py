#!/usr/bin/env python3
"""S41 audit A1 (independent referee code): brute-force check of THM-4667 statements 1, 2, 3, 4, 8 (INTRA/INTER split)
and 9 on small N, with an own implementation of the clock parity f, the weights w and chi, R(N), the threshold
coefficients K_N, the seed map, the odd transport P with sign (-1)^a, the partial currents and the canonical depth.

Nothing is imported from the session's scripts. Exact Python integers throughout.
Usage: python a1_ledger_bruteforce.py
"""
import sys, random, json, os, time
sys.dont_write_bytecode = True

MEMO = {1: 0}


def clock(n):
    """Number of shortcut steps T(n) = n/2 or (3n+1)/2 until the first hit of 1 (clock(1) = 0)."""
    path = []
    x = n
    while x not in MEMO:
        path.append(x)
        x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    c = MEMO[x]
    for y in reversed(path):
        c += 1
        MEMO[y] = c
    return MEMO[n]


def f(n):
    return 1 - 2 * (clock(n) & 1)


def w(n):
    return 2 if n % 3 == 2 else -1


def chi(n):
    return (0, 1, -1)[n % 3]


def seed(m):
    while m % 8 == 5:
        m = (m - 1) // 4
    return m


def P(s):
    y = 3 * s + 1
    a = 0
    while y % 2 == 0:
        y //= 2
        a += 1
    assert a in (1, 2), (s, a)
    return seed(y), a


def R_direct(N):
    return sum(w(n) * f(n) for n in range(N + 1, 4 * N + 1))


def initial(N):
    c = {}
    for m in range(1, 4 * N + 1, 2):
        k = -3 * chi(m) if m < 2 * N else w(m)
        if k:
            s = seed(m)
            c[s] = c.get(s, 0) + k
    return {s: v for s, v in c.items() if v}


def block_id(s):
    """Sibling block of a seed: heads s = 1 (8) and their tails 2s+1 = 3 (16) share the head; others are singletons."""
    if s % 8 == 1:
        return ('pair', s)
    if s % 16 == 3:
        return ('pair', (s - 1) // 2)
    return ('single', s)


DEPTH = {1: 0}


def depth(x):
    path = []
    y = x
    while y not in DEPTH:
        path.append(y)
        y = P(y)[0]
    dd = DEPTH[y]
    for z in reversed(path):
        dd += 1
        DEPTH[z] = dd
    return DEPTH[x]


def check_N(N, stats):
    R = R_direct(N)
    c0 = initial(N)
    # source tracking: each source x carries (current image, accumulated sign)
    src = {x: (x, 1) for x in c0}
    c = dict(c0)
    d = 0
    prevE = None
    supE = 0
    while True:
        E = sum(v * v for v in c.values())
        supE = max(supE, E)
        # statement 1: invariant
        inv = sum(v * f(u) for u, v in c.items())
        assert inv == R, ('invariant', N, d, inv, R)
        # statement 2: partial currents
        Ru = {}
        for x, (u, sg) in src.items():
            Ru[u] = Ru.get(u, 0) + c0[x] * f(x)
            assert sg == f(x) * f(u), ('sigma_d = f(x) f(u)', N, d, x)
        for u in set(Ru) | set(c):
            assert Ru.get(u, 0) == f(u) * c.get(u, 0), ('R_u = f(u) c_d(u)', N, d, u)
        assert E == sum(v * v for v in Ru.values()), ('E_d = sum R_u^2', N, d)
        # statement 9: ROOT = shell partial sum
        shell = sum(c0[x] * f(x) for x in c0 if depth(x) <= d)
        assert c.get(1, 0) == shell, ('ROOT shell identity', N, d)
        stats['depth_checks'] += 1
        if len(c) <= 1 and (not c or 1 in c):
            break
        # one transport step, with the merge ledger
        groups = {}
        for s, v in c.items():
            u, a = P(s)
            groups.setdefault(u, []).append((s, a, v))
        nxt = {}
        cross = 0
        intra = 0
        inter = 0
        for u, lst in groups.items():
            tot = 0
            for s, a, v in lst:
                tot += (1 if a == 2 else -1) * v
            if tot:
                nxt[u] = tot
            for i in range(len(lst)):
                for k in range(i + 1, len(lst)):
                    s1, a1, v1 = lst[i]
                    s2, a2, v2 = lst[k]
                    t = 2 * (1 if a1 == 2 else -1) * (1 if a2 == 2 else -1) * v1 * v2
                    cross += t
                    if block_id(s1) == block_id(s2):
                        intra += t
                        lo, hi = min(s1, s2), max(s1, s2)
                        assert hi == 2 * lo + 1 and lo % 8 == 1, ('pair shape', s1, s2)
                    else:
                        inter += t
                        if v1 * v2 != 0:
                            ratio = max(s1, s2) / min(s1, s2)
                            stats['min_inter_ratio'] = min(stats['min_inter_ratio'], ratio)
        # INTRA as stated: -2 sum_(s = 1 mod 8) c(s) c(2s+1)
        intra_formula = -2 * sum(v * c.get(2 * s + 1, 0) for s, v in c.items() if s % 8 == 1)
        assert intra_formula == intra, ('INTRA formula', N, d)
        newE = sum(v * v for v in nxt.values())
        assert newE - E == cross, ('statement 3 merge identity', N, d)
        assert newE - E == intra + inter, ('statement 8 split', N, d)
        # move the sources
        for x in src:
            u, sg = src[x]
            v, a = P(u)
            src[x] = (v, sg * (1 if a == 2 else -1))
        c = nxt
        d += 1
    # statement 4: end state
    assert (not c and R == 0) or (list(c) == [1] and c[1] == R), ('end state', N)
    stats['N'] += 1
    return supE, R


def main():
    t0 = time.time()
    random.seed(41)
    Ns = list(range(1, 401)) + sorted(random.sample(range(401, 4096), 60)) + [7910, 10870, 8191, 8192]
    stats = dict(N=0, depth_checks=0, min_inter_ratio=1e9)
    rows = {}
    for N in Ns:
        supE, R = check_N(N, stats)
        rows[N] = (supE / N, R)
    best = max(rows, key=lambda n: rows[n][0])
    print('A1: statements 1, 2, 3, 4, 8 (INTRA/INTER, pair shape), 9 verified exactly for %d values of N (%d depth checks)'
          % (stats['N'], stats['depth_checks']))
    print('    smallest scale ratio between different blocks of one fibre carrying a nonzero cross term: %.4f'
          % stats['min_inter_ratio'])
    print('    largest sup_d E_d(N)/N among the tested N: %.4f at N = %d (R = %d); N = 1: %.4f; N = 4: %.4f; N = 7910: %.4f'
          % (rows[best][0], best, rows[best][1], rows[1][0], rows[4][0], rows[7910][0]))
    print('    [%.0f s]' % (time.time() - t0))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a1_ledger_bruteforce.json')
    json.dump(dict(stats=stats, supE_over_N={str(k): v[0] for k, v in rows.items()}, R={str(k): v[1] for k, v in rows.items()}),
              open(out, 'w'), indent=0)


if __name__ == '__main__':
    main()
