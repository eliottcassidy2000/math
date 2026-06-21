#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_L7_gap_vs_ratio_THREADA_opus.py   (THREAD A, opus)

THE STRUCTURAL READ of L7.  The tight-2-cluster scan showed: the only HIGH measS7
values come from shapes that are SECRETLY one consecutive run (gap=1, f2=f1+s ->
E = {f1,...,f1+2s-1}).  A *genuinely separated* 2-cluster (a real GAP remaining in
the offset multiset) has measS7 collapse to ~0.12-0.16, far below cap.

This script makes that DICHOTOMY exact:

 (A) "MERGED" = the offset set E is an arithmetic progression (after primitivization),
     i.e. it is L1/consec in disguise.  These are NOT new -- already bounded by THM-534.
 (B) "GENUINELY SEPARATED" = E has an internal gap (consecutive differences include a
     value >= 2 somewhere strictly inside).  We measure max measS7 over (B) only.

We sweep tight 2-clusters and CLASSIFY each by its internal gap g = f2-(f1+s-1) AND by
whether the union is an AP.  Then:
  - max measS7 over genuinely-separated (g>=2, or g==1 but not an AP) shapes;
  - the margin of THAT max to cap;
  - confirm: separated worst margin >> consec margin, i.e. the binding case is consec.

Also: vary cluster sizes (s1,s2) unequal a bit (relaxing "balanced") to see if imbalance
helps the adversary.

EXACT Fractions.  cap_k as before.
"""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(int(e) for e in E if e != 0))
    bps = set([F(0), F(1)])
    for e in E:
        ae = abs(e)
        for m in range(0, 7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = set(int(((e * xm) % 1) * 7) for e in E)
        if len(secs) == 7: total += x1 - x0
    return total

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91),
       11: F(66, 91), 12: F(6, 7), 13: F(1)}

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

def is_AP(E):
    E = sorted(E)
    if len(E) < 2: return True
    d = E[1] - E[0]
    return all(E[i+1]-E[i] == d for i in range(len(E)-1))

def internal_max_gap(E):
    E = sorted(E)
    return max(E[i+1]-E[i] for i in range(len(E)-1))

if __name__ == "__main__":
    F1MAX = 40
    for (s1, s2) in [(4, 4), (5, 5), (3, 5), (4, 6)]:
        k = s1 + s2
        if k not in CAP: continue
        c = CAP[k]
        cm = measS7(list(range(1, k + 1)))
        print("=" * 78)
        print(f"  TWO TIGHT CLUSTERS sizes ({s1},{s2}), k={k}; cap={float(c):.5f} ({c}); "
              f"consec measS7={float(cm):.5f}")
        print("=" * 78)
        all_rows = []; sep_rows = []   # sep = genuinely separated (not an AP)
        for f1 in range(1, F1MAX + 1):
            for f2 in range(f1 + s1, int((f1 + s1) * 3.3) + 2):
                E = tuple(sorted(set(range(f1, f1 + s1)) | set(range(f2, f2 + s2))))
                if len(E) != k: continue
                if not primitive(E): continue
                m = measS7(E)
                ap = is_AP(E)
                g = internal_max_gap(E)
                c1 = F(2*f1 + s1 - 1, 2); c2 = F(2*f2 + s2 - 1, 2)
                rho = c2 / c1
                row = (m, rho, g, ap, f1, f2, E)
                all_rows.append(row)
                if not ap:
                    sep_rows.append(row)
        # overall worst (incl merged AP)
        bm_all = max(all_rows, key=lambda r: r[0])
        # genuinely separated worst
        bm_sep = max(sep_rows, key=lambda r: r[0]) if sep_rows else None
        print(f"  ALL shapes ({len(all_rows)}): worst measS7={float(bm_all[0]):.6f} "
              f"rho={float(bm_all[1]):.3f} AP={bm_all[3]} E={list(bm_all[6])} "
              f"margin={float(c-bm_all[0]):+.5f}")
        if bm_sep:
            m,rho,g,ap,f1,f2,E = bm_sep
            print(f"  SEPARATED (non-AP) shapes ({len(sep_rows)}): worst measS7={float(m):.6f} ({m})")
            print(f"     rho={float(rho):.4f}  internal_max_gap={g}  f1={f1},f2={f2}  E={list(E)}")
            print(f"     cap={float(c):.6f}  MARGIN={float(c-m):+.6f} ({c-m})  "
                  f"{'BREACH' if m>c else 'OK'}")
            print(f"     ==> separated margin {float(c-m):.4f} vs consec margin {float(c-cm):.4f}: "
                  f"separated is {'SAFER' if (c-m)>(c-cm) else 'TIGHTER'}")
        # top separated
        sep_rows.sort(reverse=True)
        print("  TOP 6 genuinely-separated:")
        for (m, rho, g, ap, f1, f2, E) in sep_rows[:6]:
            print(f"     measS7={float(m):.5f} rho={float(rho):.3f} gap={g} E={list(E)}")
        print()
