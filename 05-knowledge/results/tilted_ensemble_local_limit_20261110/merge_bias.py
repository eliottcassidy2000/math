#!/usr/bin/env python3
"""Sum the one-line totals of biasweights.c (several residue-class parts per order are allowed) and print, per order m
(labelled weights m!/|Aut|, exact):
  E_R[b_R^2]            second moment of the mod-4 bias of a vertex cube = deletion-pair correlation at order m+2
  rho_m = 2^m E_R[b_R^2] ratio to the twin term (= value for a uniformly random quadratic form)
  E_R|b_R|               the bound of Theorem C (4) for the mod-4 discrepancy at order m+1
  D_(m+1)                the exact labelled discrepancy P(H = 1 mod 4) - P(H = 3 mod 4) at order m+1
  nonzero-bias fraction by corank (model: 2^-corank)
Usage: python3 merge_bias.py file [file ...]"""
import sys, re
from fractions import Fraction
from collections import defaultdict
tot = defaultdict(lambda: dict(classes=0, weight=0, fails=[0, 0, 0, 0], S=0, U=0, A=0, r=defaultdict(lambda: [0, 0, 0, 0])))
for path in sys.argv[1:]:
    for line in open(path):
        g = re.match(r'order (\d+) classes (\d+) weight (\d+) failures\(mu odd, a odd, C3, aut\) (\d+) (\d+) (\d+) (\d+) S (-?\d+) U (-?\d+) A (\d+)(.*)', line)
        if not g: continue
        m = int(g.group(1)); t = tot[m]
        t['classes'] += int(g.group(2)); t['weight'] += int(g.group(3))
        for k in range(4): t['fails'][k] += int(g.group(4 + k))
        t['S'] += int(g.group(8)); t['U'] += int(g.group(9)); t['A'] += int(g.group(10))
        for h in re.finditer(r' r(\d+) (\d+) (\d+) (\d+) (\d+)', g.group(11)):
            r = int(h.group(1))
            for k in range(4): t['r'][r][k] += int(h.group(2 + k))
for m in sorted(tot):
    t = tot[m]; W = t['weight']
    complete = (W == 2 ** (m * (m - 1) // 2))
    print("order m = %d: classes %d, labelled weight %d (%s 2^C(m,2)); failures (mu odd, a odd, Theorem C (3), aut) = %s"
          % (m, t['classes'], W, "=" if complete else "NOT", t['fails']))
    e2 = sum(Fraction(v[3], 2 ** r) for r, v in t['r'].items()) / W
    e1 = Fraction(t['A'], W * 2 ** m)
    eh = sum(Fraction(v[2] + v[3]) * Fraction(1, 2) ** (Fraction(r, 2)) if r % 2 == 0 else 0 for r, v in t['r'].items()) / W
    D = Fraction(t['S'], W * 2 ** m)
    print("  E_R[b_R^2] = %s = %.9f ; rho_m = 2^m E_R[b_R^2] = %.6f" % (e2, float(e2), float(e2 * 2 ** m)))
    print("  E_R|b_R| = %.6f ; E_R[2^(-rank/2)] = %.6f ; sqrt(E_R[b_R^2]) = %.6f" % (float(e1), float(eh), float(e2) ** 0.5))
    print("  D_%d = P(H = 1 mod 4) - P(H = 3 mod 4) at order %d = %s = %.9f" % (m + 1, m + 1, D, float(D)))
    for r in sorted(t['r']):
        v = t['r'][r]
        if v[0] + v[1] == 0: continue
        print("    rank %2d (corank %d): classes %d zero-bias + %d ; labelled fraction %.7f ; nonzero-bias fraction (labelled) %.5f (model %.5f)"
              % (r, m - r, v[0], v[1], (v[2] + v[3]) / W, v[3] / (v[2] + v[3]), 2.0 ** (r - m)))
