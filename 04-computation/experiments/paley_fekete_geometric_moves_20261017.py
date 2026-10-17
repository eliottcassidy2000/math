#!/usr/bin/env python3
"""Fekete positivity = local optimality of the natural Paley order against geometrically randomised insertions (S40).

Move vertex 0 of the natural order of P_p forward by a random distance T with P(T = t) = (1 - x) x^(t-1),
t = 1, ..., p - 1 (with probability x^(p-1) it moves to the end, a change of 0 since S(p - 1) = 0). The expected change
of the forward-minus-backward excess is

    E[Delta_x] = sum_t (1 - x) x^(t-1) (-2 S(t)) = -(2/x) f_p(x),     f_p(x) = sum_(t=1)^(p-1) (t/p) x^t,

so f_p > 0 on (0, 1) (which excludes real zeros of L(s, chi_(-p)), Fekete) iff every such randomised move strictly
lowers the expected excess; L+ (all S(t) >= 0, THM-4666 (5)) is the deterministic version. This script (a) checks the
identity in exact rational arithmetic, (b) counts, for p = 3 (mod 4) prime up to PMAX, the primes in L+ and the primes
whose Fekete polynomial has no zero in (0, 1) (NUMERICAL: sign of g(x) = f_p(x)/(1 - x) = sum_t S(t) x^t on a grid of
4000 points dense near 1, float64 Horner; a double root between grid points could be missed).
Usage: python paley_fekete_geometric_moves_20261017.py PMAX
Session opus-2026-10-17-S40.
"""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import is_prime


def legendre_table(p):
    chi = [-1] * p
    chi[0] = 0
    for x in range(1, p):
        chi[(x * x) % p] = 1
    return chi


def walk(p):
    chi = legendre_table(p)
    S = [0] * p
    s = 0
    for t in range(1, p):
        s += chi[t]
        S[t] = s
    return chi, S


def exact_check(p, x):
    chi, S = walk(p)
    f = sum(Fr(chi[t]) * x ** t for t in range(1, p))
    E = sum((1 - x) * x ** (t - 1) * (-2 * S[t]) for t in range(1, p))
    return E == -2 * f / x


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    ok = all(exact_check(p, x) for p in (7, 11, 19, 23, 43, 163) for x in (Fr(1, 3), Fr(7, 10), Fr(99, 100)))
    print('exact identity E[Delta_x] = -(2/x) f_p(x) at p in {7, 11, 19, 23, 43, 163}, x in {1/3, 7/10, 99/100}:', ok)
    u = np.concatenate([np.linspace(0.001, 0.9, 1000), 1 - np.logspace(-1, -7, 3000)])
    rows = []
    for p in range(7, PMAX + 1, 4):
        if not is_prime(p):
            continue
        chi, S = walk(p)
        Lplus = min(S[1:]) >= 0
        coeff = np.array(S[1:][::-1], dtype=np.float64)          # g(x) = sum_t S(t) x^t, Horner from the top
        g = np.zeros_like(u)
        for c in coeff:
            g = g * u + c
        g = g * u
        fek = bool(np.all(g > 0))
        rows.append((p, Lplus, fek))
    n = len(rows); nl = sum(r[1] for r in rows); nf = sum(r[2] for r in rows)
    print('p = 3 (mod 4) primes 7..%d: %d; in L+: %d (%.3f); Fekete-positive on (0,1) (numerical): %d (%.3f); L+ but not Fekete-positive: %d'
          % (PMAX, n, nl, nl / n, nf, nf / n, sum(r[1] and not r[2] for r in rows)))
    print('Fekete-positive but not in L+ (first 20):', [r[0] for r in rows if r[2] and not r[1]][:20])
    for lo, hi in [(7, 500), (500, 1000), (1000, 2000), (2000, 3000)]:
        sel = [r for r in rows if lo <= r[0] < hi]
        if sel:
            print('  [%d, %d): %d primes, L+ %.3f, Fekete-positive %.3f' % (lo, hi, len(sel), sum(r[1] for r in sel) / len(sel),
                                                                        sum(r[2] for r in sel) / len(sel)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(PMAX=PMAX, rows=rows), open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_fekete_geometric_moves_20261017.json'), 'w'))


if __name__ == '__main__':
    main()
