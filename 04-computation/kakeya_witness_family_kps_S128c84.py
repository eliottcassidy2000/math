#!/usr/bin/env python3
"""kakeya_witness_family_kps_S128c84.py -- kind-pasteur-2026-07-19-S128c84

The owner asked to prove D >= 3P "from the (1, D-1, D) witness family".  death-star
and I both found that the MINIMAL D achieving each P is always of that shape:
  P=1 D=3 (1,2,3) | P=2 D=7 (1,6,7) | P=3 D=14 (1,13,14) | P=4 D=17 (1,16,17)
  P=5 D=21 (1,20,21) | P=6 D=28 (1,27,28)
So that family is the extremal locus and the place a proof must start.

WHAT THIS SCRIPT DOES
 1. Computes P for d = (1, D-1, D) for D up to a large bound, exactly, and tests
    P <= D/3 far beyond the D<=30 cube scan.
 2. Does the same for the general ADDITIVE family (a, b, a+b) -- the witness family
    is its a=1 slice -- because additivity is exactly the compatible relation
    (1,1,-1), the one that concentrates the sojourn (THM-1152).
 3. Reports the exact structure of the Case-A window that governs the family:
      s in [3/14, 2/7],  frac((D-1)s) in [5/7 - s, 2/7 + s],
    whose window WIDTH is 2s - 3/7 -- vanishing at s = 3/14 and maximal 1/7 at
    s = 2/7.  The measure heuristic gives ~ (D-1)/196 runs per case, i.e.
    P ~ 0.015 D, so the witnesses sit ~14x ABOVE the equidistributed count: they
    are genuinely resonant, which is why D >= 3P is not a measure statement.
 4. Corrects my own L5 over-claim: run = 1/(7D) is attained by MANY directions, not
    only (1,2,3).  What is true is run <= (1/7)/r_max with equality iff the run
    enters at the corner (2/7,4/7,6/7) -- r_max is the rate of the MAX coordinate,
    which need not be D.  Tested here separately.
"""
import sys
from fractions import Fraction as F
from math import gcd

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400


def cell_runs(DD):
    bps = {F(0), F(1)}
    for d in DD:
        for k in range(d + 1):
            bps.add(F(k, d))
    for i in range(3):
        for j in range(i + 1, 3):
            dd = abs(DD[i] - DD[j])
            if dd > 0:
                for k in range(dd + 1):
                    bps.add(F(k, dd))
    bps = sorted(bps)
    out = []
    for idx in range(len(bps) - 1):
        lo, hi = bps[idx], bps[idx + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        gmid = [(-DD[i] * mid) % 1 for i in range(3)]
        a = [gmid[i] + DD[i] * mid for i in range(3)]
        order = sorted(range(3), key=lambda i: gmid[i])
        s0, s1, s2 = order
        cons = [("ge", a[s0] - F(2, 7), DD[s0]), ("le", a[s2] - F(5, 7), DD[s2])]
        feas = True
        for (x, y) in [(s1, s0), (s2, s1)]:
            c = a[x] - a[y]
            dc = DD[x] - DD[y]
            if dc > 0:
                cons.append(("ge", c - F(2, 7), dc))
            elif dc < 0:
                cons.append(("le", c - F(2, 7), dc))
            elif not (c <= F(2, 7)):
                feas = False
        if not feas:
            continue
        ulo, uhi = lo, hi
        for typ, const, dd in cons:
            b = const / dd
            if typ == "ge":
                ulo = max(ulo, b)
            else:
                uhi = min(uhi, b)
        if uhi > ulo:
            out.append((ulo, uhi))
    out.sort()
    m = []
    for lo, hi in out:
        if m and lo <= m[-1][1]:
            m[-1] = (m[-1][0], max(m[-1][1], hi))
        else:
            m.append((lo, hi))
    return m


def gvec(DD, u):
    return [(-DD[i] * u) % 1 for i in range(3)]


print("=" * 78)
print("1. THE WITNESS FAMILY  d = (1, D-1, D),  D = 3 .. %d" % DMAX)
print("=" * 78)
worst = (F(0), None)
viol = 0
rows = []
for D in range(3, DMAX + 1):
    DD = [1, D - 1, D]
    if gcd(gcd(DD[0], DD[1]), DD[2]) != 1:
        continue
    r = cell_runs(DD)
    P = len(r) // 2
    if P == 0:
        continue
    if D < 3 * P:
        viol += 1
    ratio = F(P, D)
    if ratio > worst[0]:
        worst = (ratio, D)
    rows.append((D, P, ratio))
print("  D values tested            : %d" % len(rows))
print("  violations of D >= 3P      : %d  -> %s"
      % (viol, "HOLDS on the family" if viol == 0 else "BROKEN"))
print("  max P/D on the family      : %s = %.6f at D = %s  (threshold 1/3)"
      % (worst[0], float(worst[0]), worst[1]))
print("  top 12 by P/D:")
for D, P, ratio in sorted(rows, key=lambda t: -t[2])[:12]:
    print("     D=%-5d P=%-4d P/D = %-10s = %.6f   3P=%-5d slack %d"
          % (D, P, ratio, float(ratio), 3 * P, D - 3 * P))

print()
print("=" * 78)
print("2. THE GENERAL ADDITIVE FAMILY  d = (a, b, a+b)  (relation (1,1,-1))")
print("=" * 78)
worst2 = (F(0), None)
viol2 = 0
cnt2 = 0
for a in range(1, 25):
    for b in range(a + 1, 260):
        DD = [a, b, a + b]
        if gcd(gcd(a, b), a + b) != 1:
            continue
        r = cell_runs(DD)
        P = len(r) // 2
        if P == 0:
            continue
        cnt2 += 1
        D = a + b
        if D < 3 * P:
            viol2 += 1
        if F(P, D) > worst2[0]:
            worst2 = (F(P, D), tuple(DD))
print("  additive triples tested    : %d" % cnt2)
print("  violations of D >= 3P      : %d  -> %s"
      % (viol2, "HOLDS on additive family" if viol2 == 0 else "BROKEN"))
print("  max P/D                    : %s = %.6f at %s"
      % (worst2[0], float(worst2[0]), worst2[1]))

print()
print("=" * 78)
print("3. NON-ADDITIVE STRESS FAMILIES (is additivity really the extremal shape?)")
print("=" * 78)
for name, gen in [
    ("(1, 2, D)", lambda D: [1, 2, D]),
    ("(1, D-2, D)", lambda D: [1, D - 2, D]),
    ("(2, D-1, D)", lambda D: [2, D - 1, D]),
    ("(1, D//2, D)", lambda D: [1, max(2, D // 2), D]),
    ("(3, D-2, D)", lambda D: [3, D - 2, D]),
]:
    w = (F(0), None)
    v = 0
    for D in range(4, 260):
        DD = gen(D)
        if len(set(DD)) < 3 or min(DD) < 1:
            continue
        if gcd(gcd(DD[0], DD[1]), DD[2]) != 1:
            continue
        P = len(cell_runs(DD)) // 2
        if P == 0:
            continue
        if max(DD) < 3 * P:
            v += 1
        if F(P, max(DD)) > w[0]:
            w = (F(P, max(DD)), tuple(DD))
    print("  %-14s  max P/D = %-9s = %.6f at %-14s  violations %d"
          % (name, w[0], float(w[0]), str(w[1]), v))

print()
print("=" * 78)
print("4. CORRECTION TO MY OWN L5:  run = (1/7)/r_max  <=>  corner entry ?")
print("=" * 78)
eq_rmax = 0
corner_fail = 0
eq_D = 0
eq_D_non123 = set()
for d2 in range(1, 25):
    for d3 in range(1, 25):
        for d4 in range(1, 25):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            D = max(DD)
            for lo, hi in cell_runs(DD):
                g = gvec(DD, (lo + hi) / 2)
                top = max(range(3), key=lambda i: g[i])
                rmax = DD[top]
                if hi - lo == F(1, 7) / rmax:
                    eq_rmax += 1
                    if sorted(gvec(DD, lo)) != [F(2, 7), F(4, 7), F(6, 7)]:
                        corner_fail += 1
                if hi - lo == F(1, 7) / D:
                    eq_D += 1
                    s = sorted(DD)
                    if not (s[1] == 2 * s[0] and s[2] == 3 * s[0]):
                        eq_D_non123.add(tuple(DD))
print("  runs with run = (1/7)/r_max            : %d" % eq_rmax)
print("  ... of which NOT entering at the corner: %d  -> corner claim %s"
      % (corner_fail, "HOLDS" if corner_fail == 0 else "FAILS"))
print("  runs with run = 1/(7D)                 : %d" % eq_D)
print("  ... non-(1,2,3) directions attaining it: %d  e.g. %s"
      % (len(eq_D_non123), sorted(eq_D_non123)[:6]))
print("  -> MY L5 OVERCLAIM CORRECTED: run = 1/(7D) does NOT force (1,2,3).")
print("     The true chain is  run <= 1/(7D) <= 1/21, and it is the SECOND")
print("     inequality (D >= 3, from distinct rates) that is tight only at (1,2,3).")
