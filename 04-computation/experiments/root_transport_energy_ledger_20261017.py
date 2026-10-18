#!/usr/bin/env python3
"""Energy ledger of the exact ROOT transport (S41): is the signed odd transport (nearly) orthogonal?

Sources: primitive odd s (s != 5 mod 8) with coefficient H_N(s) = sum_(j: B^j(s) <= 4N) K_N(B^j(s)), where
K_N(m) = w(m)[m <= 4N] + (-3 chi(m) - w(m))[m < 2N] (clock_duality_oct15's two threshold events, R1) and B(x) = 4x + 1.
One odd transport: s -> seed((3s+1)/2^a), a = nu_2(3s+1) in {1, 2}, weight (-1)^a; seed strips B. The signed sum
sum_s c_d(s) f(s) = R(N) = H(4N) - H(N) is invariant (f the Collatz depth-parity sign), and at full depth every source
sits at ROOT (seed 1) with coefficient R(N).

For each depth d the script records the l1 mass, the energy E = sum c^2, the ROOT coefficient, the number of terminals,
and the energy change. If merges were orthogonal (random relative signs), the energy would be conserved and R(N)^2 would
equal the initial energy: the sqrt(N) law as energy conservation. The cross terms are the correlations of the current
between merging fibres, i.e. across scale factors 4^j.
Usage: python root_transport_energy_ledger_20261017.py k1 k2 ...   (N = 2^k)
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np


def w(m):
    return 2 if m % 3 == 2 else -1


def chi(m):
    r = m % 3
    return 0 if r == 0 else (1 if r == 1 else -1)


def seed(u):
    while u % 8 == 5:
        u = (u - 1) // 4
    return u


def parity_clock(n):
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        c += 1
    return c & 1


def ledger(N, dmax=400):
    # initial coefficients on primitive seeds
    c = {}
    for m in range(1, 4 * N + 1, 2):
        k = w(m) * (m <= 4 * N) + (-3 * chi(m) - w(m)) * (m < 2 * N)
        if k:
            s = seed(m)
            c[s] = c.get(s, 0) + k
    c = {s: v for s, v in c.items() if v}
    R = sum(v * (1 - 2 * parity_clock(s)) for s, v in c.items())
    rows = []
    for d in range(dmax):
        l1 = sum(abs(v) for v in c.values()); E = sum(v * v for v in c.values())
        root = c.get(1, 0)
        rows.append(dict(d=d, l1=l1, energy=E, root=root, terminals=len(c)))
        if len(c) == 1 and 1 in c:
            break
        nxt = {}
        for s, v in c.items():
            if s == 1:
                nxt[1] = nxt.get(1, 0) + v
                continue
            u = 3 * s + 1
            a = 0
            while u % 2 == 0:
                u //= 2; a += 1
            assert a in (1, 2), (s, a)
            t = seed(u)
            nxt[t] = nxt.get(t, 0) + (v if a == 2 else -v)
        c = {s: v for s, v in nxt.items() if v}
        if not c:
            rows.append(dict(d=d + 1, l1=0, energy=0, root=0, terminals=0)); break
    return R, rows


def main():
    out = []
    for k in [int(a) for a in sys.argv[1:]]:
        N = 2 ** k
        R, rows = ledger(N)
        E0 = rows[0]['energy']
        print('N = 2^%d: R(N) = %d, R^2/N = %.3f; initial l1/N = %.3f, initial energy/N = %.3f; depth to ROOT = %d'
              % (k, R, R * R / N, rows[0]['l1'] / N, E0 / N, rows[-1]['d']))
        for r in rows[:12] + rows[12::6]:
            print('   d = %3d  l1/N = %.4f  energy/N = %.4f  root = %8d  terminals = %8d' % (r['d'], r['l1'] / N, r['energy'] / N, r['root'], r['terminals']))
        print('   final: root = %d (R = %d), energy/N = %.4f' % (rows[-1]['root'], R, rows[-1]['energy'] / N), flush=True)
        out.append(dict(k=k, N=N, R=R, rows=rows))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_energy_ledger_20261017.json'), 'w'))


if __name__ == '__main__':
    main()
