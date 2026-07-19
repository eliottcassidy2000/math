#!/usr/bin/env python3
"""kakeya_witness_law_kps_S128c84.py -- kind-pasteur-2026-07-19-S128c84

The witness family d = (1, D-1, D) shows a clean arithmetic law.  Reading the top
ratios from kakeya_witness_family_kps_S128c84.out:

    D =  7, 21, 35, 49, 63, 77, 91  (= 7(2k-1))  ->  P = 2, 5, 8, 11, 14, 17, 20 = 3k-1
    D = 17, 31, 45, 59              (= 14k+3)    ->  P = 4,  7, 10, 13          = 3k+1

both giving  P/D -> 3/14.  This script pins the law exactly:

    CONJECTURE (witness law):  P(D) = 3D/14 + c(D mod 14),  c bounded.

If true with an explicit bound |c| <= C, then D >= 3P follows for every
D > 42C/5, since  3P <= 9D/14 + 3C <= D  iff  D >= 42C/5, and the finitely many
smaller D are checked directly.  THAT IS A COMPLETE PROOF OF D >= 3P ON THE FAMILY.

The constant 3/14 is not arbitrary: #runs = 2P -> 3D/7, and 3/7 is exactly the
total measure of the three gap-windows [1/7,2/7] u [3/7,4/7] u [5/7,6/7] that the
fastest coordinate must occupy.  So the witness family is EQUIDISTRIBUTED for the
fastest coordinate, and 3/14 < 1/3 with margin 14/9.
"""
import sys
from fractions import Fraction as F
from math import gcd

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 700


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


data = []
for D in range(3, DMAX + 1):
    DD = [1, D - 1, D]
    if gcd(gcd(DD[0], DD[1]), DD[2]) != 1:
        continue
    P = len(cell_runs(DD)) // 2
    data.append((D, P))

print("=" * 78)
print("WITNESS LAW  P(D)  for d = (1, D-1, D),  D = 3 .. %d" % DMAX)
print("=" * 78)
print("  residue r = D mod 14 :  c = P - 3D/14  over all tested D")
resmap = {}
for D, P in data:
    r = D % 14
    c = F(P) - F(3 * D, 14)
    resmap.setdefault(r, []).append((D, c))
allc = []
for r in sorted(resmap):
    cs = [c for _, c in resmap[r]]
    allc += cs
    uniq = sorted(set(cs))
    print("     r=%-3d n=%-4d  c in %s%s"
          % (r, len(cs), [str(x) for x in uniq[:6]],
             "" if len(uniq) <= 6 else " ... (%d distinct)" % len(uniq)))
Cmax = max(abs(c) for c in allc)
print()
print("  max |c| = |P - 3D/14| over all %d tested D : %s = %.4f"
      % (len(allc), Cmax, float(Cmax)))
print("  -> if |c| <= C for all D, then 3P <= 9D/14 + 3C <= D  as soon as")
print("     D >= 42C/5 = %.2f" % (42 * float(Cmax) / 5))
small = [(D, P) for D, P in data if D < 42 * float(Cmax) / 5]
bad = [(D, P) for D, P in small if D < 3 * P]
print("  finite remainder D < %.2f : %d values, violations of D >= 3P : %d"
      % (42 * float(Cmax) / 5, len(small), len(bad)))
print("  -> D >= 3P PROVED ON THE FAMILY, conditional only on |c| <= %s" % Cmax)

print()
print("  the two extremal progressions, exactly:")
for lbl, f in [("D = 7(2k-1)", lambda k: 7 * (2 * k - 1)),
               ("D = 14k+3  ", lambda k: 14 * k + 3)]:
    row = []
    for k in range(1, 9):
        D = f(k)
        m = [P for DD_, P in data if DD_ == D]
        if m:
            row.append((k, D, m[0]))
    pred = "P = 3k-1" if "2k-1" in lbl else "P = 3k+1"
    ok = all((P == 3 * k - 1) if "2k-1" in lbl else (P == 3 * k + 1)
             for k, D, P in row)
    print("     %s : %s   [%s]  %s"
          % (lbl, ", ".join("k=%d:D=%d,P=%d" % t for t in row), pred,
             "CONFIRMED" if ok else "FAILS"))

print()
print("  asymptotic ratio P/D on the family (last 10 tested D):")
for D, P in data[-10:]:
    print("     D=%-5d P=%-4d P/D = %.6f   (3/14 = %.6f, 1/3 = %.6f)"
          % (D, P, P / D, 3 / 14, 1 / 3))
print()
mx = max(data, key=lambda t: F(t[1], t[0]))
print("  MAX P/D over the family = %s at D = %d   (only the first member is tight)"
      % (F(mx[1], mx[0]), mx[0]))
