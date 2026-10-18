#!/usr/bin/env python3
"""S43b referee: independent recomputation of the half-period twin census (written from scratch).

Library of 10 binomial factors in (n, k); all products of 3 or 4 factors (with repetition); u(n) = sum_{k=0}^{n}.
- identity of sequences decided on n = 0..NID (default 40), compared with the session's n = 0..12;
- twin relation u((p-1)/2) == eps v((p-1)/2) (mod p^2) for ALL primes 5 <= p <= PMAX (default 400),
  eps in {+1, (-1)^n}; classes = connected components of the relation (full pairwise, no single-base screening);
- the trivial class (u == 0 mod p^2 at every half-period) is reported separately;
- robustness: the same with eps in {+-1, +-(-1)^n, +-(2/p), +-(-2/p)}, and with scale factors 4^n, 16^n allowed.
Usage: python audit_census.py [PMAX] [NID]
"""
import sys
import itertools
from math import comb
import numpy as np


def C(a, b):
    if a < 0 or b < 0 or b > a:
        return 0
    return comb(a, b)


LIB = [
    ('C(n,k)', lambda n, k: C(n, k)),
    ('C(n+k,k)', lambda n, k: C(n + k, k)),
    ('C(2k,k)', lambda n, k: C(2 * k, k)),
    ('C(2k,n)', lambda n, k: C(2 * k, n)),
    ('C(2n-2k,n-k)', lambda n, k: C(2 * n - 2 * k, n - k)),
    ('C(k,n-k)', lambda n, k: C(k, n - k)),
    ('C(n,2k)', lambda n, k: C(n, 2 * k)),
    ('C(2k,n-k)', lambda n, k: C(2 * k, n - k)),
    ('C(n+k,2k)', lambda n, k: C(n + k, 2 * k)),
    ('C(2n,n+k)', lambda n, k: C(2 * n, n + k)),
]


def primes_in(a, b):
    return [p for p in range(max(a, 2), b + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.p[max(a, b)] = min(a, b)


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    NID = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    L = len(LIB)
    combos = [c for r in (3, 4) for c in itertools.combinations_with_replacement(range(L), r)]
    print('products generated: %d (3-factor %d, 4-factor %d)' % (len(combos), sum(len(c) == 3 for c in combos),
                                                                  sum(len(c) == 4 for c in combos)))
    # exact small values
    fac = [[[LIB[i][1](n, k) for k in range(n + 1)] for n in range(NID + 1)] for i in range(L)]

    def exact_vals(c, N):
        out = []
        for n in range(N + 1):
            tot = 0
            for k in range(n + 1):
                v = 1
                for i in c:
                    v *= fac[i][n][k]
                    if not v:
                        break
                tot += v
            out.append(tot)
        return tuple(out)

    vals = {c: exact_vals(c, NID) for c in combos}
    groups40 = {}
    for c, v in vals.items():
        groups40.setdefault(v, []).append(c)
    groups12 = {}
    for c, v in vals.items():
        groups12.setdefault(v[:13], []).append(c)
    nz40 = {v: cs for v, cs in groups40.items() if any(v)}
    nz12 = {v: cs for v, cs in groups12.items() if any(v)}
    zero40 = [cs for v, cs in groups40.items() if not any(v)]
    print('distinct nonzero sequences: %d by n <= %d, %d by n <= 12 (session); identically-zero products: %d'
          % (len(nz40), NID, len(nz12), sum(len(cs) for cs in zero40)))
    # representatives (one per distinct sequence, n <= NID)
    reps = [cs[0] for cs in nz40.values()]
    name = lambda c: ' '.join(LIB[i][0] for i in c)
    primes = primes_in(5, PMAX)
    R = np.zeros((len(reps), len(primes)), dtype=np.int64)
    for j, p in enumerate(primes):
        n = (p - 1) // 2
        pp = p * p
        tab = [np.array([LIB[i][1](n, k) % pp for k in range(n + 1)], dtype=np.int64) for i in range(L)]
        for r, c in enumerate(reps):
            acc = tab[c[0]].copy()
            for i in c[1:]:
                acc = (acc * tab[i]) % pp
            R[r, j] = int(acc.sum() % pp)
    pp_arr = np.array([p * p for p in primes], dtype=np.int64)
    sgn_n = np.array([(-1) ** ((p - 1) // 2) for p in primes], dtype=np.int64)
    chi2 = np.array([legendre(2, p) for p in primes], dtype=np.int64)
    chim2 = np.array([legendre(-2, p) for p in primes], dtype=np.int64)
    four_n = np.array([pow(4, (p - 1) // 2, p * p) for p in primes], dtype=np.int64)
    sixteen_n = np.array([pow(16, (p - 1) // 2, p * p) for p in primes], dtype=np.int64)
    zero_rows = [r for r in range(len(reps)) if not R[r].any()]

    def classes(eps_list):
        dsu = DSU(len(reps))
        key = {}
        for r in range(len(reps)):
            if r in zero_rows:
                continue
            key.setdefault(tuple(R[r]), []).append(r)
        for r in range(len(reps)):
            if r in zero_rows:
                continue
            for eps in eps_list:
                t = tuple(((eps * R[r]) % pp_arr).tolist())
                for s in key.get(t, []):
                    if s != r:
                        dsu.union(r, s)
        comp = {}
        for r in range(len(reps)):
            if r in zero_rows:
                continue
            comp.setdefault(dsu.find(r), []).append(r)
        return [sorted(v) for v in comp.values() if len(v) > 1]

    base = classes([np.ones_like(sgn_n), sgn_n])
    print('\n[B1] twin classes, eps in {1, (-1)^n}, all %d primes 5 <= p <= %d: %d classes' % (len(primes), PMAX, len(base)))
    for cl in sorted(base, key=lambda v: v[0]):
        print('   class:')
        for r in cl:
            al = nz40[vals[reps[r]]]
            print('      %s%s' % (name(reps[r]), '   (+%d identical rewritings)' % (len(al) - 1) if len(al) > 1 else ''))
    print('[B2] trivial class (u == 0 mod p^2 at every tested half-period): %d distinct sequences' % len(zero_rows))
    for r in zero_rows[:12]:
        print('      %s' % name(reps[r]))
    ext = classes([s * e for s in (np.ones_like(sgn_n), -np.ones_like(sgn_n)) for e in (np.ones_like(sgn_n), sgn_n, chi2, chim2)])
    print('[B3] with eps in {+-1, +-(-1)^n, +-(2/p), +-(-2/p)}: %d classes (sizes %s)' % (len(ext), sorted(len(c) for c in ext)))
    new = [c for c in ext if c not in base]
    for cl in new:
        print('      new/merged class: %s' % [name(reps[r]) for r in cl])
    four_mn = np.array([pow(4, -((p - 1) // 2), p * p) for p in primes], dtype=np.int64)
    sixteen_mn = np.array([pow(16, -((p - 1) // 2), p * p) for p in primes], dtype=np.int64)
    scales = [np.ones_like(sgn_n), four_n, sixteen_n, four_mn, sixteen_mn]
    sc = classes([(e * s) % pp_arr for s in (np.ones_like(sgn_n), sgn_n) for e in scales])
    print('[B4] allowing the scale factors 4^(+-n), 16^(+-n) (and (-1)^n): %d classes (sizes %s)'
          % (len(sc), sorted(len(c) for c in sc)))
    for cl in sc:
        if cl not in base:
            print('      merged class: %s' % [name(reps[r]) for r in cl])
    # false-merge check: groups identical at n <= 12 but different at n <= NID
    split = 0
    for v12, cs in groups12.items():
        if len(set(vals[c] for c in cs)) > 1:
            split += 1
    print('[B5] sequences identical for n <= 12 but different for n <= %d: %d groups' % (NID, split))


if __name__ == '__main__':
    main()
