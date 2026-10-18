#!/usr/bin/env python3
"""Audit F3, item (e), optional: the ideal (N -> infinity) depth-2 clock-twin statistics, to test the note's claim that
the depth-2 share of twin pairs with 3 | Delta is a rational density computable from residues mod 27 and that the
measured 12% is that density.

A depth-2 backward word w = (j1, a1, j2, a2) from a node t gives the source
  s_w(t) = 2^tau t / 9 + gamma_w,  tau = a1 + 2 j1 + a2 + 2 j2,
  gamma_w = 2^(a2 + 2 j2) beta1 / 3 + beta2,  beta_i = (2^(a_i) (4^(j_i) - 1)/3 - 1)/3,
so two words with equal tau have Delta = gamma_w - gamma_w' independent of t (an integer when both apply).
Bands, weights and the ideal source coefficient are those of f3_band_exact.py; only sources with c_0 != 0 are
counted (the session's convention).  Pair counts and c0 c0' sums are summed over bands b < BMAX (exact Fractions).
Usage: python3 f3_twins_ideal.py [BMAX]
"""
import sys, os
from fractions import Fraction as Fr
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from f3_band_exact import h, pre_class, B          # own audit code


def beta(a, j):
    return Fr(2 ** a * (4 ** j - 1) // 3 - 1, 3)


def words2_full(r27, b, child_cut=6):
    out = []
    x, j1 = r27 % 27, 0
    while 9 * 2.0 ** (2 * j1 + 1 - b - 2) <= child_cut:
        pc = pre_class(x, 27)
        if pc is not None:
            a1, p9 = pc
            tau1 = a1 + 2 * j1
            if 9 * 2.0 ** (tau1 - b - 2) <= child_cut:
                y, j2 = p9, 0
                while 2 * j2 + 1 + tau1 - b - 2 <= 0:
                    pc2 = pre_class(y, 9)
                    if pc2 is not None:
                        a2, s3 = pc2
                        tau = tau1 + a2 + 2 * j2
                        e = tau - b - 2
                        if e <= 0:
                            v = h(s3, e)
                            if v:
                                g = Fr(2 ** (a2 + 2 * j2)) * beta(a1, j1) / 3 + beta(a2, j2)
                                out.append((j1, tau, (-1) ** (a1 + a2) * v, g))
                    y = B(y, 9)
                    j2 += 1
        x = B(x, 27)
        j1 += 1
    return out


def main():
    BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 70
    cnt = defaultdict(Fr)      # Delta mod 3 -> weighted pair density
    sm = defaultdict(Fr)       # Delta mod 3 -> weighted sum of c0 c0' (signs cancel: equal clocks)
    for b in range(BMAX):
        w = Fr(1, 16 * 2 ** b)
        for r in range(27):
            ws = words2_full(r, b)
            byclock = defaultdict(list)
            for j1, tau, v, g in ws:
                byclock[tau].append((j1, v, g))
            for tau, L in byclock.items():
                for i in range(len(L)):
                    for k in range(i + 1, len(L)):
                        if L[i][0] == L[k][0]:
                            continue          # same child: met earlier
                        D = L[i][2] - L[k][2]
                        assert D.denominator == 1, D
                        dm = int(D) % 3
                        cnt[dm] += w
                        sm[dm] += w * L[i][1] * L[k][1]
    tot = sum(cnt.values())
    out = ['ideal depth-2 clock twins (c_0 != 0), bands b < %d:' % BMAX,
           '  pair density per N: %.6f ; share with 3 | Delta: %.6f ; shares by Delta mod 3: %s'
           % (float(tot), float(cnt[0] / tot), ['%.4f' % float(cnt[d] / tot) for d in range(3)]),
           '  sums of c0 c0\' per N: 3|Delta %+.6f, Delta != 0 mod 3 %+.6f, total %+.7f'
           % (float(sm[0]), float(sm[1] + sm[2]), float(sum(sm.values())))]
    print('\n'.join(out))
    open(os.path.join(HERE, 'f3_twins_ideal.out'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
