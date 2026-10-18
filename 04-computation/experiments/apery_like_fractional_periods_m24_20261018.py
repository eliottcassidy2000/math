#!/usr/bin/env python3
"""Apery-like sequences at the fractional periods n == -s (mod p), s in {1/3, 2/3, 1/4, 3/4, 1/6, 5/6}, against the
26 M24 eta-products and the truncated hypergeometric sums (S43, exploring around THM-4669).

For each prime p >= 7 and each s, n_s(p) is the least n >= 0 with n == -s (mod p) (so n < p). We test, for the
sequences u of apery_like_m24_supercongruences_20261018.py:
  (a) u(n_s(p)) == (d/p) a_g(p) (mod p^2) and (mod p), for all 26 Mason eta-products g and twists d;
  (b) A(n) == sum_{k<=n} b_k^2 and s_7(n) == (-1)^n sum_k b_k b_(n-k) with b_k = (s)_k (1-s)_k/k!^2, mod p^2 and mod p
      (the THM-4669 mechanism; at s != 1/2 a first-order correction is expected to spoil mod p^2).
Usage: python apery_like_fractional_periods_m24_20261018.py [PMAX]
Session opus-2026-10-17-S43.
"""
import sys, os, json, math
from math import comb
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m24_paley_golay_eta_20261018 import SHAPES, eta_product
from apery_like_m24_supercongruences_20261018 import SEQ, kron, D


def red(x, m):
    return x.numerator * pow(x.denominator, -1, m) % m


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    primes = [p for p in range(7, PMAX + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    N = PMAX + 2
    forms = {}
    for name, (sh, cen) in SHAPES.items():
        b = eta_product(sh, N)
        a = [0] * (N + 1)
        for n in range(N):
            a[n + 1] = b[n]
        forms[name] = (a, math.lcm(*sh.keys()) * min(sh.keys()))
    S = [Fr(1, 3), Fr(2, 3), Fr(1, 4), Fr(3, 4), Fr(1, 6), Fr(5, 6)]
    seqs = {k: v for k, v in SEQ.items() if 'control' not in k and 'Zagier A' not in k}
    out = {}
    for s in S:
        print('s = %s' % s)
        nsp = {}
        for p in primes:
            # least n >= 0 with n == -s mod p
            n = (-s.numerator * pow(s.denominator, -1, p)) % p
            nsp[p] = n
        found = []
        vals = {name: {p: fn(nsp[p]) for p in primes} for name, fn in seqs.items()}
        for name in seqs:
            for g, (a, level) in forms.items():
                for d in D:
                    tested = [p for p in primes if level % p and d % p]
                    if len(tested) < 10:
                        continue
                    ok2 = all((vals[name][p] - kron(d, p) * a[p]) % (p * p) == 0 for p in tested)
                    ok1 = ok2 or all((vals[name][p] - kron(d, p) * a[p]) % p == 0 for p in tested)
                    if ok1:
                        found.append((name, g, d, 'p^2' if ok2 else 'p'))
        for f in found:
            print('   %-4s %-22s ~ %-5s * (%d/p)' % (f[3], f[0], f[1], f[2]))
        # mechanism checks
        mech = {}
        for test in ('A vs sum b^2', 's7 vs conv'):
            res = {}
            for mod in (2, 1):
                ok = True
                for p in primes:
                    n = nsp[p]
                    m = p ** mod
                    bk = []
                    acc = Fr(1)
                    for k in range(n + 1):
                        if k > 0:
                            acc *= Fr((s + k - 1) * (1 - s + k - 1), k * k)
                        bk.append(acc)
                    if test == 'A vs sum b^2':
                        lhs = vals['Apery A (zeta3)'][p]
                        rhs = red(sum(x * x for x in bk), m)
                    else:
                        lhs = vals['Cooper s7'][p]
                        rhs = (-1) ** n * red(sum(bk[k] * bk[n - k] for k in range(n + 1)), m)
                    if (lhs - rhs) % m:
                        ok = False
                        break
                res['p^%d' % mod] = ok
            mech[test] = res
        print('   mechanism:', mech)
        out[str(s)] = dict(found=found, mechanism=mech)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'apery_like_fractional_periods_m24_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
