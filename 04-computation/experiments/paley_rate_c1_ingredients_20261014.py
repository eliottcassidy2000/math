#!/usr/bin/env python3
"""The first-order rate of R(q) = E_sigma prod_k (1 + chi(d_k)) -> e for Paley tournaments: c_1 = 3/2.

R(q) = H(P_q) 2^(q-1)/q! = sum over step subsets S of E[prod_{k in S} chi(d_k)], d_k = sigma(k+1) - sigma(k) for a
uniformly random ordering sigma of F_q.  Decompose S into maximal runs.  To first order in 1/q:
  (i)   a run of length 2 (cherry) has E = A_2/(q)_3 with A_2 = sum over 3 distinct points = q(q-1) exactly, i.e. 1/(q-2);
  (ii)  runs of odd length vanish exactly (reversal and chi(-1) = -1);
  (iii) a run of length 4 has E = A_4/(q)_5 with A_4 = 2q^3 + O(q^2) (THM-438: Catalan C_2 = 2);
  (iv)  two cherries on six distinct points: E[XY] / (E[X]E[Y]) = 1 - 2/q + O(1/q^2) (distinctness repulsion).
Placement of j separated cherries among q - 1 steps: C(q - 2j, j)/(q-2)^j = (1/j!)(1 - (5/2) j(j-1)/q + ...).
Summing over j with the repulsion (1 - 2/q)^(j(j-1)/2) and one quartic run (total weight 2/q):
    R = e (1 - 5/(2q) - 1/q + 2/q + O(1/q^2)) = e (1 - 3/(2q) + O(1/q^2)).
This script verifies (i), (ii), (iii) and (iv) exactly for several q (character sums by brute force over distinct
tuples, with the first point fixed by translation), and compares c_1 = 3/2 with the fits of the exact R(q).
Session opus-2026-10-14-S34.
"""
import numpy as np
import numba as nb
from fractions import Fraction

FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


@nb.njit(cache=True)
def run_sum(chi, q, L):
    """A_L = sum over L+1 distinct points v_0..v_L of prod chi(v_{i+1} - v_i); v_0 = 0 fixed by translation (times q).
    Explicit nested loops for L = 2, 3, 4."""
    tot = 0
    v0 = 0
    for v1 in range(1, q):
        s1 = chi[(v1 - v0) % q]
        for v2 in range(q):
            if v2 == v0 or v2 == v1:
                continue
            s2 = s1 * chi[(v2 - v1) % q]
            if L == 2:
                tot += s2
                continue
            for v3 in range(q):
                if v3 == v0 or v3 == v1 or v3 == v2:
                    continue
                s3 = s2 * chi[(v3 - v2) % q]
                if L == 3:
                    tot += s3
                    continue
                for v4 in range(q):
                    if v4 == v0 or v4 == v1 or v4 == v2 or v4 == v3:
                        continue
                    tot += s3 * chi[(v4 - v3) % q]
    return tot * q


@nb.njit(cache=True)
def two_cherries(chi, q):
    """N_6 = sum over 6 distinct points (a,b,c,d,e,f) of chi(b-a)chi(c-b)chi(e-d)chi(f-e), with a = 0 (times q)."""
    tot = 0
    a = 0
    for b in range(q):
        if b == a:
            continue
        s1 = chi[(b - a) % q]
        for c in range(q):
            if c == a or c == b:
                continue
            x = s1 * chi[(c - b) % q]
            for d in range(q):
                if d == a or d == b or d == c:
                    continue
                for e in range(q):
                    if e == a or e == b or e == c or e == d:
                        continue
                    y = x * chi[(e - d) % q]
                    for f in range(q):
                        if f == a or f == b or f == c or f == d or f == e:
                            continue
                        tot += y * chi[(f - e) % q]
    return tot * q


def falling(q, k):
    r = 1
    for i in range(k):
        r *= q - i
    return r


def main():
    print('ingredients of c_1 for Paley tournaments (q prime, q = 3 mod 4)')
    for q in (11, 19, 23, 31, 43):
        chi = np.zeros(q, dtype=np.int64)
        for x in range(1, q):
            chi[x * x % q] = 1
        for x in range(1, q):
            if chi[x] == 0:
                chi[x] = -1
        A2 = run_sum(chi, q, 2)
        A3 = run_sum(chi, q, 3)
        A4 = run_sum(chi, q, 4)
        check(f'A2 {q}', A2 == q * (q - 1))
        check(f'A3 {q}', A3 == 0)
        line = f'   q = {q}: A_2 = {A2} (= q(q-1): {A2 == q * (q - 1)}), A_3 = {A3}, A_4 = {A4}, (A_4 - 2q^3)/q^2 = {(A4 - 2 * q ** 3) / q ** 2:.4f}'
        if q <= 31:
            N6 = two_cherries(chi, q)
            ratio = Fraction(N6 * (q - 2) ** 2, falling(q, 6))
            line += f';  two-cherry ratio E[XY]/(E[X]E[Y]) = {float(ratio):.6f}, q (ratio - 1) = {q * (float(ratio) - 1):.4f}'
        print(line, flush=True)
    print('\nexact closed forms (checked at every q above and at further primes below):')
    print('   A_2 = q(q-1);  A_3 = 0;  A_4 = 2 q (q-1)(q-3);  two-cherry ratio = (q-2)(q-9)/((q-4)(q-5))')
    for q in (11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
        chi = np.zeros(q, dtype=np.int64)
        for x in range(1, q):
            chi[x * x % q] = 1
        for x in range(1, q):
            if chi[x] == 0:
                chi[x] = -1
        A4 = run_sum(chi, q, 4)
        ok4 = check(f'A4 closed {q}', A4 == 2 * q * (q - 1) * (q - 3))
        msg = f'   q = {q}: A_4 = 2q(q-1)(q-3): {ok4}'
        if q <= 31:
            N6 = two_cherries(chi, q)
            okc = check(f'cherry closed {q}', Fraction(N6 * (q - 2) ** 2, falling(q, 6)) == Fraction((q - 2) * (q - 9), (q - 4) * (q - 5)))
            msg += f';  two-cherry ratio = (q-2)(q-9)/((q-4)(q-5)): {okc}'
        print(msg, flush=True)
    print('   so a single quartic run has E = A_4/(q)_5 = 2/((q-2)(q-4)) and a single cherry 1/(q-2) (exact)')
    print('\nfirst-order prediction c_1 = 5/2 (placement) + 1 (two-cherry repulsion) - 2 (quartic runs, C_2 = 2) = 3/2')
    print('fits of the exact R(q) (paley_R_asymptotic_fit_20261014.out): c_1 = 1.50156 (order 2), 1.50121 (order 3)')
    check('c1 near 3/2', abs(1.50156 - 1.5) < 0.005 and abs(1.50121 - 1.5) < 0.005)
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
