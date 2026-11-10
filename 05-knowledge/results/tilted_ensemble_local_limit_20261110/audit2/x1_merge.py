#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  Add the raw one-line totals printed by x1_cube.c (several parts per order allowed) and
print, per order m, the exact labelled statistics of the vertex cubes:
  E_R[b_R^2] = Q / (4^m W), rho_m = 2^m E_R[b_R^2], E_R|b_R| = A / (2^m W), E_R b_R = U / (2^m W),
  D_(m+1) = S1 / (2^m W), the labelled rank law, the non-zero-bias frequency per corank, the decomposition
  rho_m - 1 = sum_c P(corank c) (2^c f_c - 1), the antipodal term P(all beta-degrees even), all failure counters.
Usage: python3 x1_merge.py raw_file [raw_file ...]"""
import sys, re
from fractions import Fraction
from collections import defaultdict

tot = {}
for path in sys.argv[1:]:
    for line in open(path):
        g = re.match(r'order (\d+) mode (\S+) classes (\d+) W (-?\d+) S1 (-?\d+) U (-?\d+) A (-?\d+) Q (\d+) wEven0 (\d+) wEven1 (\d+) minmu (\d+) maxmu (\d+)(?: affdimB (\d+) affdimLB (\d+))?', line)
        if not g: continue
        key = (int(g.group(1)), g.group(2))
        t = tot.setdefault(key, dict(classes=0, W=0, S1=0, U=0, A=0, Q=0, wE0=0, wE1=0, minmu=10 ** 30, maxmu=0, parts=0,
                                     fails=defaultdict(int), counts=defaultdict(int), mod8=[0] * 8, rank=defaultdict(lambda: [0, 0, 0, 0]),
                                     edges=defaultdict(int), deg={1: defaultdict(int), 2: defaultdict(int), 3: defaultdict(int)}, aff=[]))
        v = [int(x) if x is not None else -1 for x in g.groups()[2:]]
        t['classes'] += v[0]; t['W'] += v[1]; t['S1'] += v[2]; t['U'] += v[3]; t['A'] += v[4]; t['Q'] += v[5]; t['wE0'] += v[6]; t['wE1'] += v[7]
        t['minmu'] = min(t['minmu'], v[8]); t['maxmu'] = max(t['maxmu'], v[9]); t['parts'] += 1; t['aff'].append((v[10], v[11]))
        for h in re.finditer(r' (fail_\w+) (\d+)', line): t['fails'][h.group(1)] += int(h.group(2))
        for h in re.finditer(r' (n_\w+) (\d+)', line): t['counts'][h.group(1)] += int(h.group(2))
        h = re.search(r' mod8((?: -?\d+){8})', line)
        for k, x in enumerate(h.group(1).split()): t['mod8'][k] += int(x)
        for h in re.finditer(r' rank (\d+) (\d+) (\d+) (\d+) (\d+)', line):
            for k in range(4): t['rank'][int(h.group(1))][k] += int(h.group(2 + k))
        h = re.search(r' edges((?: \d+:\d+)*)', line)
        for x in h.group(1).split(): e, c = x.split(':'); t['edges'][int(e)] += int(c)
        for kk in (1, 2, 3):
            h = re.search(r' deg%d((?: -?\d+:\d+)*)' % kk, line)
            if h:
                for x in h.group(1).split(): d, c = x.rsplit(':', 1); t['deg'][kk][int(d)] += int(c)

for (m, mode) in sorted(tot):
    t = tot[(m, mode)]; W = t['W']; full = 2 ** (m * (m - 1) // 2)
    print("order m = %d (%s; %d part%s): %d %s, total weight %d (%s 2^C(m,2) = %d)"
          % (m, mode, t['parts'], '' if t['parts'] == 1 else 's', t['classes'], 'labelled tournaments' if mode == 'labelled' else 'classes', W, '=' if W == full else 'NOT', full))
    bad = {k: v for k, v in t['fails'].items() if v}
    print("  failures: %s ; controls: %s ; mu in [%d, %d]" % (bad if bad else 'none (all %d counters zero)' % len(t['fails']), dict(t['counts']), t['minmu'], t['maxmu']))
    e2 = Fraction(t['Q'], 4 ** m * W); e1 = Fraction(t['A'], 2 ** m * W); eb = Fraction(t['U'], 2 ** m * W); D = Fraction(t['S1'], 2 ** m * W)
    print("  E_R[b_R^2] = %s = %.9f ; rho_m = %.6f ; E_R|b_R| = %s = %.6f ; E_R b_R = %.6f" % (e2, float(e2), float(e2 * 2 ** m), e1, float(e1), float(eb)))
    print("  D_%d = %s = %.9f ; D^2 = %.3e <= (E|b|)^2 = %.3e <= E[b^2] = %.3e : %s"
          % (m + 1, D, float(D), float(D * D), float(e1 * e1), float(e2), D * D <= e1 * e1 <= e2))
    r8 = t['mod8']; s8 = sum(r8)
    print("  labelled residues of H mod 8 at order %d (1,3,5,7): %s ; even residues: %d" % (m + 1, " ".join("%.5f" % (r8[k] / s8) for k in (1, 3, 5, 7)), r8[0] + r8[2] + r8[4] + r8[6]))
    chk = Fraction(0); excess = Fraction(0)
    for r in sorted(t['rank']):
        v = t['rank'][r]; lab = v[2] + v[3]; c = m - r
        f = Fraction(v[3], lab) if lab else Fraction(0)
        chk += Fraction(v[3], W) / 2 ** r; excess += Fraction(lab, W) * (2 ** c * f - 1)
        print("    rank %2d (corank %d): classes %d zero-bias + %d non-zero ; labelled fraction %.7f (%s) ; non-zero-bias frequency %.5f (uniform-form model %.5f) ; contribution to rho - 1: %+.6f"
              % (r, c, v[0], v[1], lab / W, Fraction(lab, W), float(f), 2.0 ** (-c), float(Fraction(lab, W) * (2 ** c * f - 1))))
    assert chk == e2, "E[b^2] from ranks disagrees with the directly accumulated sum of squares"
    print("  rho - 1 = %.6f = sum of the corank contributions (%.6f) ; P(all beta-degrees even) = %s = %.6f (the antipodal pair X, V - X contributes exactly this to rho)"
          % (float(e2 * 2 ** m - 1), float(excess), Fraction(t['wE0'] + t['wE1'], W), (t['wE0'] + t['wE1']) / W))
    print("  number of beta-edges (classes): %s" % " ".join("%d:%d" % kv for kv in sorted(t['edges'].items())))
    if any(t['deg'][k] for k in (1, 2, 3)):
        print("  F_2-degree of digit k of X -> H(R + v_X) (%s): %s" % ('labelled' if mode == 'labelled' else 'classes',
              " ; ".join("digit %d: %s" % (k, dict(sorted(t['deg'][k].items()))) for k in (1, 2, 3))))
    if mode == 'labelled' and t['parts'] == 1:
        print("  affine span over F_2: beta_R dimension %d (C(m,2) - 1 = %d) ; (ell_R, beta_R) dimension %d (m + C(m,2) - 2 = %d)"
              % (t['aff'][0][0], m * (m - 1) // 2 - 1, t['aff'][0][1], m + m * (m - 1) // 2 - 2))
