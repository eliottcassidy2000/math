#!/usr/bin/env python3
"""Where does a genuine Riemann-hypothesis (Weil/Hasse) term enter the Paley Hamiltonian-path count?

R(q) = sum over configurations of path-step runs; a run of length L contributes A_L/(q)_(L+1) with
A_L = sum over L+1 distinct points v_0..v_L of prod chi(v_i - v_(i-1)) (chi the quadratic character, q = 3 mod 4).
A_2 = q(q-1), A_3 = 0 (odd runs vanish), A_4 = 2q(q-1)(q-3) are exact polynomials in q (elementary Jacobi-sum
identities; checked for 11 primes).  This script computes A_5 (must be 0) and A_6 exactly for several primes and tests
whether A_6 is a polynomial in q (interpolate through some primes, predict the others).  A failure of polynomiality
would locate the first place where a curve over F_q (and hence the Hasse-Weil bound, the Riemann hypothesis for
curves) enters the count.
Session opus-2026-10-14-S34.
"""
import sys
import numpy as np
import numba as nb
from fractions import Fraction

QS = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [11, 19, 23, 31, 43, 47]


@nb.njit(cache=True)
def run_sum(chi, q, L):
    """A_L for L <= 6 by nested loops, v_0 = 0 (times q)"""
    tot = 0
    for v1 in range(1, q):
        s1 = chi[v1]
        for v2 in range(1, q):
            if v2 == v1:
                continue
            s2 = s1 * chi[(v2 - v1) % q]
            if L == 2:
                tot += s2
                continue
            for v3 in range(1, q):
                if v3 == v1 or v3 == v2:
                    continue
                s3 = s2 * chi[(v3 - v2) % q]
                if L == 3:
                    tot += s3
                    continue
                for v4 in range(1, q):
                    if v4 == v1 or v4 == v2 or v4 == v3:
                        continue
                    s4 = s3 * chi[(v4 - v3) % q]
                    if L == 4:
                        tot += s4
                        continue
                    for v5 in range(1, q):
                        if v5 == v1 or v5 == v2 or v5 == v3 or v5 == v4:
                            continue
                        s5 = s4 * chi[(v5 - v4) % q]
                        if L == 5:
                            tot += s5
                            continue
                        for v6 in range(1, q):
                            if v6 == v1 or v6 == v2 or v6 == v3 or v6 == v4 or v6 == v5:
                                continue
                            tot += s5 * chi[(v6 - v5) % q]
    return tot * q


def interp_predict(xs, ys, x):
    """Lagrange interpolation through (xs, ys), evaluated at x (exact rationals)"""
    tot = Fraction(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        term = Fraction(yi)
        for j, xj in enumerate(xs):
            if j != i:
                term *= Fraction(x - xj, xi - xj)
        tot += term
    return tot


def main():
    vals = {}
    for q in QS:
        chi = np.zeros(q, dtype=np.int64)
        sq = set((x * x) % q for x in range(1, q))
        for x in range(1, q):
            chi[x] = 1 if x in sq else -1
        A5 = run_sum(chi, q, 5)
        A6 = run_sum(chi, q, 6)
        vals[q] = A6
        print(f'   q = {q}: A_5 = {A5}, A_6 = {A6}, A_6/q^4 = {A6 / q ** 4:.5f}, (A_6 - 5 q^4)/q^3 = {(A6 - 5 * q ** 4) / q ** 3:.5f}', flush=True)
    # A_6 is divisible by q(q-1) (translations and the multiplier group x -> a x, a a square, act; also checks)
    for q in QS:
        print(f'   q = {q}: A_6 / (q (q-1)) = {Fraction(vals[q], q * (q - 1))}')
    # polynomiality: A_6 has degree <= 5 in q if polynomial (leading 5 q^4 by THM-438); fit degree-4 polynomial in q
    # for A_6/q (which is then of degree <= 4) through 5 primes and predict the 6th
    qs = sorted(vals)
    if len(qs) >= 6:
        for leave in qs:
            rest = [q for q in qs if q != leave][:5]
            pred = interp_predict(rest, [Fraction(vals[q], q) for q in rest], leave)
            print(f'   degree-4 interpolation of A_6/q through {rest} predicts q = {leave}: {float(pred * leave):.1f} vs exact {vals[leave]};  '
                  f'difference {float(pred * leave - vals[leave]):+.1f}')


if __name__ == '__main__':
    main()
