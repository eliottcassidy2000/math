#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 1 (kps-S24-wf7): adversarial cap-check + error-vs-margin verdict on WIDE configs.

Findings to confirm/refute:
 - The signed error p0(E)-p0_decorr can EXCEED the fixed Q-margin cap-Q(k-1) (refutes the prompt's
   "<=0.012, 10x slack" premise). Locate the peak error and its sign.
 - BUT the ACTUAL p0(E) stays < cap_k on all WIDE configs (the real LRC bound). Push hard for a
   cap-violation; report the closest approach (sup p0, margin to cap).

Strategy: target the error/p0 PEAK structure found in v2 -- a small bounded base (so the decorr
baseline is small) PLUS an adjacent far cluster that itself sweeps all 6 sectors well. Adversarially
maximize p0 over wide configs with: (i) structured stretched-AP far clusters not containing 0,
(ii) random wide configs, (iii) the v2 peak family. Exact rationals.
"""
from __future__ import annotations
import sys, random, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l; den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e; x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps); num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi; mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def main():
    print("=" * 80)
    print("THREAD 1: ADVERSARIAL CAP-CHECK on WIDE LRC(14) configs (kps-S24-wf7)")
    print("=" * 80)
    for k in (8, 9, 10, 11, 12):
        cap = float(CAP[k]); margin = float(MARGIN[k])
        print("-" * 80)
        print(f"k={k}: cap={cap:.6f} Q(k-1)={float(QVAL[k]):.6f} Q-margin={margin:.6f}")
        random.seed(2024 + k)
        sup_p0 = 0.0; arg = None; overcap = 0; N = 0
        nz_over_margin = 0  # count configs where (p0 - best-decorr-split) > Q-margin

        def consider(E):
            nonlocal sup_p0, arg, overcap, N
            E = tuple(sorted(set(E)))
            if len(E) != k or E[0] != 0:
                return
            if max(E) <= 14:
                return  # not wide
            if reduce(gcd, [e for e in E if e]) != 1:
                return
            pv = float(p0_fast(E)); N += 1
            if pv > sup_p0:
                sup_p0 = pv; arg = E
            if pv > cap:
                overcap += 1

        # (i) stretched APs not containing the pin densely: {0} u {s, s+g, s+2g, ...}
        for s in range(15, 40):
            for g in range(1, 8):
                E = [0] + [s + g * j for j in range(k - 1)]
                consider(E)
        # (ii) small base + adjacent far cluster (the v2 peak shape)
        for bsize in range(1, k):
            rfar = k - bsize
            for _ in range(300):
                B = sorted(set([0] + random.sample(range(1, 15), max(0, bsize - 1))))
                if len(B) != bsize:
                    continue
                start = max(B) + random.choice([1, 2, 3, 5, 15 - max(B) + random.randint(1, 10)])
                start = max(start, 15 - 0)
                far = [start + j for j in range(rfar)]  # adjacent far cluster
                consider(B + far)
                far2 = [start + 2 * j for j in range(rfar)]
                consider(B + far2)
        # (iii) fully random wide configs
        for _ in range(8000):
            top = random.randint(15, 50)
            rest = random.sample(range(1, top), k - 2)
            E = [0] + rest + [top]
            consider(E)
        # (iv) two adjacent clusters split by a big gap (resonant: cluster ratio ~1)
        for c1 in range(1, k):
            c2 = k - c1
            if c1 < 1 or c2 < 1:
                continue
            for gap in (15, 16, 20, 28):
                for _ in range(150):
                    B = sorted(set([0] + (random.sample(range(1, 14), c1 - 1) if c1 > 1 else [])))
                    if len(B) != c1:
                        continue
                    s = max(B) + gap
                    C = sorted(set([s + x for x in ([0] + random.sample(range(1, 14), c2 - 1) if c2 > 1 else [0])]))
                    consider(B + C)

        print(f"  scanned {N} wide configs")
        print(f"  sup ACTUAL p0 = {sup_p0:.6f}  margin_to_cap = {cap - sup_p0:.6f}  OVER_CAP = {overcap}")
        print(f"  argmax p0 E = {arg}")
        print(f"  => REAL bound p0<cap holds with margin {cap-sup_p0:.6f}? {overcap == 0}")
    print()
    print("=" * 80)
    print("VERDICT: report sup p0 vs cap (real bound) and whether any wide config violates it.")
    print("=" * 80)


if __name__ == "__main__":
    main()
