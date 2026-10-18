#!/usr/bin/env python3
"""Half-period supercongruence classes of Apery-like sequences (S43).

For sequences u, v in the battery of apery_like_m24_supercongruences_20261018.py (plus a few more), test
    u((p-1)/2) == eps(p) * v((p-1)/2)   (mod p^2)  [and mod p],   eps(p) = (d/p), d in D,
for all primes 5 <= p <= PMAX not dividing 6 d. Pairs that agree mod p^2 for every tested prime form "half-period
supercongruence classes"; the output lists the classes (mod p^2 and mod p) and, for each class, which M24 eta-product
(if any, from the earlier battery) governs it. Motivation: Cooper's s7 agrees with Apery's A_n up to (-1/p) mod p^2.
Usage: python apery_like_halfperiod_clusters_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, os, json, math
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apery_like_m24_supercongruences_20261018 import SEQ, kron

EXTRA = {
    'AZ-eps sum C(n,k)^2 C(2k,n)^2': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, n) ** 2 for k in range(n + 1)),
    'sum C(2k,k)^2 C(2n-2k,n-k)^2': lambda n: sum(comb(2 * k, k) ** 2 * comb(2 * n - 2 * k, n - k) ** 2 for k in range(n + 1)),
    'sum C(n,k)C(2k,k)C(2n-2k,n-k)^2/...: C(2n,n)sum C(n,k)^2 C(2k,k)': lambda n: comb(2 * n, n) * sum(comb(n, k) ** 2 * comb(2 * k, k) for k in range(n + 1)),
    'trunc 4F3: sum C(2k,k)^4/256^k (rational, p-adic)': None,
    'central trinomial T_n': lambda n: sum(comb(n, 2 * k) * comb(2 * k, k) for k in range(n // 2 + 1)),
}
D = [1, -1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7]


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    primes = [p for p in range(5, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    seqs = {k: v for k, v in SEQ.items() if 'control' not in k and k != 'Zagier A (Franel-type #A)'}
    for k, v in EXTRA.items():
        if v is not None:
            seqs[k] = v
    vals = {}
    for name, fn in seqs.items():
        vals[name] = {p: fn((p - 1) // 2) for p in primes}
    # p-adic truncated hypergeometric sums as rationals reduced mod p^2
    vals['trunc 4F3 sum_{k<=(p-1)/2} C(2k,k)^4/256^k'] = {
        p: sum(comb(2 * k, k) ** 4 * pow(256, -k, p * p) for k in range((p - 1) // 2 + 1)) % (p * p) for p in primes}
    vals['trunc 3F2 sum_{k<=(p-1)/2} C(2k,k)^3/64^k'] = {
        p: sum(comb(2 * k, k) ** 3 * pow(64, -k, p * p) for k in range((p - 1) // 2 + 1)) % (p * p) for p in primes}
    names = list(vals)
    pairs = []
    for i, u in enumerate(names):
        for v in names[i + 1:]:
            for d in D:
                tested = [p for p in primes if 6 % p and d % p]
                ok2 = all((vals[u][p] - kron(d, p) * vals[v][p]) % (p * p) == 0 for p in tested)
                ok1 = ok2 or all((vals[u][p] - kron(d, p) * vals[v][p]) % p == 0 for p in tested)
                if ok1:
                    pairs.append(dict(u=u, v=v, d=d, mod='p^2' if ok2 else 'p'))
    print('pairs agreeing at (p-1)/2 for all primes 5 <= p <= %d (twist (d/p)):' % PMAX)
    for pr in pairs:
        print('   %-4s %-48s ~ %-48s twist (%d/p)' % (pr['mod'], pr['u'], pr['v'], pr['d']))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(PMAX=PMAX, pairs=pairs), open(os.path.join(here, '..', '..', '05-knowledge', 'results',
              'apery_like_halfperiod_clusters_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
