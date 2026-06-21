#!/usr/bin/env python3
"""
ROUTE 5, part 6 -- WHICH inclusion-exclusion level drives consec's
measS7-maximality?  Systematic sweep over ALL full-residue primitive k=8
shapes with span<=14, EXACT.

For each shape compute:
  measS7, and the level masses MISS_k (k=1..6).
Then test:
  (Q1) Does consec minimize MISS_1 (best single-sector marginal coverage)?
  (Q2) Is measS7-rank predicted by MISS_1-rank alone? (low-order dominance)
  (Q3) Is there a high-order (k>=2) sign-flip mechanism like the tournament?
       i.e. does any shape BEAT consec at some MISS level it loses at another?
  (Q4) Decompose measS7 = MISS_0 - MISS_1 + MISS_2 - ... ; attribute the
       consec gap over the runner-up to each level.

This pins the STRUCTURAL DIFFERENCE between the two extremalities:
  tournament = HIGH-ORDER (packing j>=2) AP advantage with a level sign-flip;
  LRC        = LOW-ORDER  (single-sector MISS_1) consec advantage, no flip.

Author: mac-mini-2026-06-21-S20 (ROUTE 5)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.path.insert(0, '04-computation')
from lrc14_route1_conductance_minimax_opus_0621 import (
    sec, breakpoints, is_full_residue, primitive, consec, measS7,
)
sys.stdout.reconfigure(line_buffering=True)

HALF = F(1, 14)


def occ_set(E, a, y):
    return set(sec(e, a, y) for e in E)


def miss_levels(E):
    exact = defaultdict(F)
    for a in range(1, 7):
        bps = breakpoints(E, a, HALF)
        for lo, hi in zip(bps, bps[1:]):
            if hi <= lo:
                continue
            mid = (lo + hi) / 2
            occ = occ_set(E, a, mid)
            empty = frozenset(s for s in range(1, 7) if s not in occ)
            exact[empty] += hi - lo
    miss_level = defaultdict(F)
    allsec = list(range(1, 7))
    for k in range(0, 7):
        tot = F(0)
        for S in itertools.combinations(allsec, k):
            Sset = frozenset(S)
            for E_set, v in exact.items():
                if Sset <= E_set:
                    tot += v
        miss_level[k] = tot
    return miss_level


def gen_shapes(k=8, span_max=14):
    """All full-residue primitive shapes {0=...} of size k, max element <= span_max,
    containing 0 (anchor), distinct, full residue mod 7."""
    out = []
    # 0 fixed; choose k-1 from 1..span_max
    for combo in itertools.combinations(range(1, span_max + 1), k - 1):
        E = [0] + list(combo)
        if not is_full_residue(E):
            continue
        if not primitive(E):
            continue
        out.append(E)
    return out


def main():
    print("ROUTE 5 part 6: LRC inclusion-exclusion level dominance sweep")
    print("=" * 70)
    shapes = gen_shapes(8, 14)
    print(f"full-residue primitive k=8 span<=14 shapes: {len(shapes)}")

    data = []
    for E in shapes:
        ml = miss_levels(E)
        ms = measS7(E)
        data.append((tuple(E), ms, ml))

    # sort by measS7 desc
    data.sort(key=lambda r: -r[1])
    consec_E = tuple(consec(8))

    # rank of consec by measS7
    ms_rank = next(i for i, r in enumerate(data) if r[0] == consec_E) + 1
    print(f"consec measS7-rank: {ms_rank} / {len(data)} "
          f"(measS7={float(data[0][1] if data[0][0]==consec_E else next(r[1] for r in data if r[0]==consec_E)):.6f})")

    # Q1: does consec minimize MISS_1?
    by_miss1 = sorted(data, key=lambda r: r[2][1])
    m1_rank = next(i for i, r in enumerate(by_miss1) if r[0] == consec_E) + 1
    print(f"consec MISS_1-rank (lower=better): {m1_rank} / {len(data)}")
    print(f"  consec MISS_1 = {float(next(r[2][1] for r in data if r[0]==consec_E)):.6f}")
    print(f"  min MISS_1 over all = {float(by_miss1[0][2][1]):.6f} "
          f"(shape {by_miss1[0][0]})")

    # Q2: correlation between measS7 and -MISS_1
    import statistics
    ms_vals = [float(r[1]) for r in data]
    m1_vals = [float(r[2][1]) for r in data]
    # Spearman via ranks
    def ranks(xs):
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        rk = [0]*len(xs)
        for pos, i in enumerate(order):
            rk[i] = pos
        return rk
    rms = ranks(ms_vals); rm1 = ranks([-v for v in m1_vals])
    n = len(ms_vals)
    d2 = sum((rms[i]-rm1[i])**2 for i in range(n))
    spear = 1 - 6*d2/(n*(n*n-1))
    print(f"Spearman(measS7, -MISS_1) = {spear:.4f}")

    # also each higher level
    for k in range(2, 6):
        mk = [float(r[2][k]) for r in data]
        rmk = ranks([-v for v in mk])
        d2 = sum((rms[i]-rmk[i])**2 for i in range(n))
        sp = 1 - 6*d2/(n*(n*n-1))
        print(f"Spearman(measS7, -MISS_{k}) = {sp:.4f}")

    # Q3: high-order sign flip? For the top runner-up, show per-level diffs.
    print("\nTop 6 by measS7 with level masses:")
    print(f"  {'shape':>26} {'measS7':>10} {'M1':>8} {'M2':>8} {'M3':>8} {'M4':>8} {'M5':>8}")
    for E, ms, ml in data[:6]:
        tag = " <-consec" if E == consec_E else ""
        print(f"  {str(E):>26} {float(ms):>10.6f} "
              f"{float(ml[1]):>8.4f} {float(ml[2]):>8.4f} {float(ml[3]):>8.4f} "
              f"{float(ml[4]):>8.4f} {float(ml[5]):>8.4f}{tag}")

    # Q4: attribution -- consec vs runner-up, per (-1)^k level
    cons = next(r for r in data if r[0] == consec_E)
    runner = next(r for r in data if r[0] != consec_E)  # top non-consec
    print(f"\nAttribution consec vs runner-up {runner[0]} (measS7 gap "
          f"{float(cons[1]-runner[1]):+.6f}):")
    netchk = F(0)
    for k in range(1, 7):
        contrib = ((-1)**k) * (cons[2][k] - runner[2][k])  # contribution to measS7 diff
        # measS7 = sum (-1)^k MISS_k, so diff = sum (-1)^k (MISS_k(cons)-MISS_k(run))
        # consec advantage = -[that]; let's report signed contribution to (cons-run)
        signed = ((-1)**k) * (cons[2][k] - runner[2][k])
        netchk += signed
        print(f"  level {k}: (-1)^{k} (M_k(cons)-M_k(run)) = {float(signed):+.6f}")
    print(f"  net (should equal measS7 gap): {float(netchk):+.6f}")

    # count how many shapes beat consec at SOME single MISS level (k>=2)
    flips = 0
    for E, ms, ml in data:
        if E == consec_E:
            continue
        # does this shape have LOWER (better) MISS_k than consec at some k>=2
        # while consec has lower at k=1? (the tournament-style flip)
        better_hi = any(ml[k] < cons[2][k] for k in range(2, 6))
        if better_hi and cons[2][1] < ml[1]:
            flips += 1
    print(f"\nShapes that beat consec at some MISS_k (k>=2) while losing MISS_1: "
          f"{flips}/{len(data)-1}")
    print("  (tournament-style high-order flip: AP sacrifices low level, wins high.)")


if __name__ == "__main__":
    main()
