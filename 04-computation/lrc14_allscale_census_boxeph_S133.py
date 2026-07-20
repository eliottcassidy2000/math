#!/usr/bin/env python3
"""
lrc14_allscale_census_boxeph_S133.py  (HYP-7930; owner-directed)

PART A — the S132b exhaustive census extended to ALL NINE scales of the q=38
in-band slice (one per straddle pair (v, 38-v)):

  scale m (odd unit, m <= 17, plus m=1): allowed speeds A_m = {c in [1,37] :
  c*m mod 38 in [3,35]} (|A_m| = 33); straddlers s1 = 3*m^{-1} mod 38 and
  s2 = 38 - s1 forced; enumerate ALL 12-sets {s1,s2} <= C <= A_m with
  covering{2..12}, passing all five spread rungs {13,17,19,23,25}; exact M
  (all q <= 74 complete by Pinch, speeds <= 37) for every survivor.

Statement target: the min-M over the ENTIRE q=38 in-band universe (all scales,
all straddle pairs) — S132b gave 3/29 at m=1; does the global min stay on a
rung, and does any scale dip into the gap (1/13, 2/25)?  (Killed families have
M >= 2/25 by the rungs, so gap-emptiness of the slice = min over survivors > gap.)

Usage: python3 ... A <m1> [m2] ...   (run one or more scales)
boxeph-2026-07-19-S133.  Pure Python, exact integers.
"""

import sys
from math import gcd

RUNGS = [13, 17, 19, 23, 25]
DUTIES = [7, 8, 9, 10, 11, 12]

def upairs(p):
    return [j for j in range(1, (p + 1) // 2) if gcd(j, p) == 1 and gcd(p - j, p) == 1]

def covering(sp):
    return all(any(v % q == 0 for v in sp) for q in range(2, 13))

def M_int(sp, cap2):
    bn, bd, qstar = 0, 1, None
    for q in range(2, cap2 + 1):
        for b in range(1, q // 2 + 1):
            md = q
            for c in sp:
                r = (c * b) % q
                d = r if r <= q - r else q - r
                if d < md:
                    md = d
                    if md * bd <= bn * q:
                        break
            if md * bd > bn * q:
                bn, bd, qstar = md, q, q
    g = gcd(bn, bd)
    return bn // g, bd // g, qstar

def census38(m):
    band = {r for r in range(38) if 3 <= min(r, 38 - r) <= 19 and r != 0}
    band = {r for r in range(38) if min(r, 38 - r) >= 3}
    Am = [c for c in range(1, 38) if (c * m) % 38 in band]
    minv = pow(m, -1, 38)
    s1 = (3 * minv) % 38
    s2 = 38 - s1
    assert s1 in Am and s2 in Am and len(Am) == 33
    VALS = [c for c in Am if c not in (s1, s2)]
    # duty availability
    for q in DUTIES:
        if not any(c % q == 0 for c in Am):
            print("m=%2d straddle (%d,%d): duty %d EMPTY in A_m -> 0 survivors" % (m, s1, s2, q))
            return []
    dbit = {q: 1 << k for k, q in enumerate(DUTIES)}
    DFULL = (1 << len(DUTIES)) - 1
    off, width = {}, 0
    for p in RUNGS:
        off[p] = width; width += len(upairs(p)) + 1
    SEG = {p: (1 << len(upairs(p))) - 1 for p in RUNGS}

    def contrib(v):
        d = 0
        for q in DUTIES:
            if v % q == 0: d |= dbit[q]
        s = 0
        for p in RUNGS:
            if v % p == 0:
                s |= 1 << (off[p] + len(upairs(p)))
            elif gcd(v % p, p) == 1:
                s |= 1 << (off[p] + upairs(p).index(min(v % p, p - v % p)))
        return d, s

    D0 = contrib(s1)[0] | contrib(s2)[0]
    S0 = contrib(s1)[1] | contrib(s2)[1]
    DV = [contrib(v) for v in VALS]
    n = len(VALS)
    sufD = [0] * (n + 1); sufS = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        sufD[i] = sufD[i + 1] | DV[i][0]
        sufS[i] = sufS[i + 1] | DV[i][1]

    def rungs_ok(s):
        for p in RUNGS:
            seg = s >> off[p]
            L = len(upairs(p))
            if not ((seg >> L) & 1 or (seg & SEG[p]) == SEG[p]):
                return False
        return True

    survivors = []
    stack = []

    def rec(start, got_d, got_s):
        slots = 10 - len(stack)
        if slots == 0:
            if got_d != DFULL: return
            if rungs_ok(got_s):
                survivors.append(tuple(sorted([s1, s2] + stack)))
            return
        if (DFULL & ~(got_d | sufD[start])): return
        if not rungs_ok(got_s | sufS[start]): return
        if n - start < slots: return
        for i in range(start, n):
            stack.append(VALS[i])
            rec(i + 1, got_d | DV[i][0], got_s | DV[i][1])
            stack.pop()

    rec(0, D0, S0)
    # exact M for all survivors
    blocker_hist = {}
    best = (1, 1, None, None)   # M as (num, den, q, C)
    for C in survivors:
        blk = tuple(p for p in RUNGS if any(v % p == 0 for v in C))
        blocker_hist[blk] = blocker_hist.get(blk, 0) + 1
        nu, de, qs = M_int(list(C), 74)
        if nu * best[1] < best[0] * de:
            best = (nu, de, qs, C)
    assert all(covering(list(C)) for C in survivors[:20])
    top = sorted(blocker_hist.items(), key=lambda x: -x[1])[:3]
    print("m=%2d straddle (%2d,%2d): survivors=%6d  minM=%d/%d=%.4f (q=%s)  top blockers %s"
          % (m, s1, s2, len(survivors), best[0], best[1], best[0] / best[1], best[2],
             ["%s:%d" % (b, c) for b, c in top]))
    print("      minM family: %s" % (best[3],))
    return [(best[0], best[1], best[3])]

if __name__ == "__main__":
    mode = sys.argv[1]
    assert mode == "A"
    for marg in sys.argv[2:]:
        census38(int(marg))
