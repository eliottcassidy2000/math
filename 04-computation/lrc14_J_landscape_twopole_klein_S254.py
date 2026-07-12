# -*- coding: utf-8 -*-
# klein-2026-07-11-S254: THE J-FUNCTIONAL LANDSCAPE + TWO-POLE STRUCTURE (exact).
#
# The moment-ladder base rows are:
#   k=9 (deg-2):  J := E[N(7-N)] = 6 m1 - m2  >=  432/91   for ALL 9-cores  (THM-711).
#   k=8 (deg-3):  the cubic requirement (THM-712/713), min margin at consec.
# where N(x) = # empty sectors among the 7 arcs [s/7,(s+1)/7), for the set of phases
# {frac(e x): e in E}.  mac-mini cont.37 found TWO extremal poles: consec (three-gap
# bunching) and mod-7-aligned {1,8,15,...} (7-sector resonance); consec wins the J-race.
#
# THIS SCRIPT (exact rationals, breakpoint sweep):
#  (A) J-LANDSCAPE: J for consec, mod-7-aligned, dilated-AP, dissociated, random, and
#      HYBRIDS (interpolations consec<->mod7) at k=8,9 -- is consec the UNIQUE global J-min?
#  (B) DECORRELATION: J-margin (J - 432/91) vs spread -- does it grow (cores far from the
#      poles are SAFE), giving a threshold beyond which no core is near-extremal?
#  (C) EXACT FORM: is J an exact affine function of a pair-correlation profile
#      P(d) = sum over pairs with (e_i + e_j, e_i - e_j) structure?  Regress with exact
#      rationals; the two poles should be the two extreme points of that profile.
#
# Conventions == mac-mini/kps: sectors [c/7,(c+1)/7); e=0 (if present) occupies sector 0;
# N = 7 - #distinct occupied sectors; m_r = E[N(N-1)...(N-r+1)]; J = 6 m1 - m2.

import sys
from math import gcd
from fractions import Fraction as F
from itertools import combinations

def lcm(a, b):
    return a // gcd(a, b) * b

def moments(E):
    """Exact (m1, m2, m3, J) for core E via all-integer breakpoint sweep."""
    nz = [abs(e) for e in E if e != 0]
    has0 = len(nz) < len(E)
    L = 1
    for e in nz:
        L = lcm(L, e)
    D = 7 * L
    pts = set([0, D])
    for e in nz:
        step = L // e
        pts.update(range(0, D + 1, step))
    pts = sorted(pts)
    pn = [0] * 8
    base_hit = 1 if has0 else 0
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2
        hit = base_hit
        for e in nz:
            c = (7 * e * s // (2 * D)) % 7
            hit |= 1 << c
        n = 7 - bin(hit).count("1")
        pn[n] += t2 - t1
    p = [F(x, D) for x in pn]
    m1 = sum(F(j) * p[j] for j in range(8))
    m2 = sum(F(j * (j - 1)) * p[j] for j in range(8))
    m3 = sum(F(j * (j - 1) * (j - 2)) * p[j] for j in range(8))
    J = 6 * m1 - m2
    return m1, m2, m3, J

def norm(E):
    nz = [e for e in E if e != 0]
    g = 0
    for e in nz:
        g = gcd(g, e)
    return tuple(e // g for e in E) if g > 1 else tuple(E)

THRESH9 = F(432, 91)   # k=9 J-floor

def families(k):
    """Structured + hybrid families of size k."""
    out = {}
    out['consec'] = tuple(range(1, k + 1))
    out['mod7-align'] = tuple(1 + 7 * i for i in range(k))       # {1,8,15,...}
    out['dilate2'] = tuple(2 * i for i in range(1, k + 1))
    out['AP-step3'] = tuple(1 + 3 * i for i in range(k))
    # hybrids: first j elements consec, rest mod-7-shifted
    for j in range(1, k):
        hyb = list(range(1, j + 1)) + [j + 7 * i for i in range(1, k - j + 1)]
        out[f'hyb{j}'] = tuple(hyb)
    return out

def landscape(k):
    print(f"=== (A) J-LANDSCAPE at k={k} (floor 432/91 = {float(THRESH9):.4f} is the k=9 target) ===")
    fams = families(k)
    rows = []
    for name, E in fams.items():
        m1, m2, m3, J = moments(E)
        rows.append((float(J), name, E, J))
    rows.sort()
    for Jf, name, E, J in rows:
        marg = J - THRESH9
        print(f"  {name:12s} J = {J}  ~ {Jf:.5f}   margin(vs432/91) = {float(marg):+.5f}   {E}")
    print(f"  MIN J family: {rows[0][1]}  (consec is min: {rows[0][1]=='consec'})")

def decorrelation(k, W):
    """(B) exhaustive: for each min-spread threshold, the worst (lowest) J over cores with
    that spread. Does the worst-J RISE as min-spread grows (=> far cores safe)?"""
    print(f"\n=== (B) DECORRELATION at k={k}: worst J by spread band, box [1..{W}] ===")
    seen = set()
    # bucket by spread = max - min
    buckets = {}
    n = 0
    for combo in combinations(range(1, W + 1), k):
        key = norm(combo)
        if key in seen:
            continue
        seen.add(key)
        n += 1
        spread = combo[-1] - combo[0]
        _, _, _, J = moments(combo)
        b = spread // 5 * 5
        if b not in buckets or J < buckets[b][0]:
            buckets[b] = (J, combo)
    print(f"  {n} normalized cores; worst (min) J per spread-band:")
    for b in sorted(buckets):
        J, combo = buckets[b]
        print(f"    spread [{b},{b+5}): min J = {float(J):.5f} (margin {float(J-THRESH9):+.5f}) at {combo}")

def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    landscape(k)
    decorrelation(k, W)

if __name__ == "__main__":
    main()
