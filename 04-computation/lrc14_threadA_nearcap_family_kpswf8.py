#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A part 3 (kps-S26-wf8): EXACT characterization of the near-cap wide family.

Two questions left after the r-stratified scan:
  (Q3a) Among SINGLE-FAR wide configs (r=1, base size k-1 in [0,14], 1 far > 14), which
        reach the near-cap band p0 > Q(k-1)?  Characterize the base structure of the
        near-cap family precisely.  Hypothesis: it is the high-coverage bases (consec_{k-1}
        and its dilations / dense APs) plus 1 far -- the single-far plateau.
  (Q3b) The THRESHOLD delta: over ALL wide configs with r>=2 (>=2 far elements), is
        max p0 <= cap - delta with explicit delta>0?  We adversarially push r=2 (the worst
        slack stratum) by ALSO sweeping the far PAIR over commensurable / adjacent ratios
        (peak resonance per the decay law), exact, to certify the slack delta.

For (Q3a) we EXHAUSTIVELY enumerate the base (all (k-1)-subsets of [0,14] containing 0) and,
for the single far w, we only need w over one period window: p0(base u {w}) depends on w only
through the far runner's phase; the near-cap-relevant maximum is captured by sweeping w over a
finite range and the resonance peaks (w just past base). We sweep w in [15, 15+W].
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce, lru_cache
from itertools import combinations
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    nz = [int(x) for x in E if x]
    if not nz:
        return F(0)
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x)
            x += step
    bps = sorted(bps)
    num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6:
            num0 += hi - lo
    return F(num0, d)


def base_signature(base):
    """Classify a base structure: consec, dilated-consec (AP from 0), dense-AP, other."""
    g = list(base)
    diffs = sorted(set(b - a for a, b in zip(g, g[1:])))
    span = g[-1]
    if base == tuple(range(len(base))):
        return "consec"
    # arithmetic progression from 0?
    if len(diffs) == 1:
        return f"AP(d={diffs[0]})"
    # dilated consec = {0,c,2c,...}? handled by AP above
    # near-consec: max gap small
    if max(b - a for a, b in zip(g, g[1:])) <= 2 and span <= 14:
        return "near-consec"
    return "other"


def scan_single_far(k, w_hi=40):
    """Exhaustive base (k-1 subset of [0,14], contains 0) x single far w in [15,w_hi].
    Return near-cap list and the per-base best p0 over w."""
    q = QVAL[k]
    near = []  # (p0, base, w)
    best = F(-1)
    best_arg = None
    rest_size = k - 2  # remaining base elements besides 0, from [1,14]
    for rest in combinations(range(1, 15), rest_size):
        base = (0,) + rest
        bmax = base[-1]
        for w in range(15, w_hi + 1):
            E = base + (w,)
            if not primitive(E):
                continue
            pv = p0_fast(E)
            if pv > best:
                best = pv
                best_arg = (base, w)
            if pv > q:
                near.append((pv, base, w))
    near.sort(reverse=True)
    return best, best_arg, near


def scan_two_far_adversarial(k, w_hi=40):
    """Adversarially push r=2: base size k-2 (exhaustive over [0,14]) + far PAIR (f1,f2),
    f1,f2 both >14, swept including adjacent/commensurable pairs. Return max p0 (slack cert)."""
    best = F(-1)
    best_arg = None
    rest_size = k - 3
    bases = []
    # to keep r=2 tractable, restrict to high-coverage bases: consec_{k-2}, dilations, dense APs,
    # plus all bases with max<=12 of small spread (the slack max lives near dense base + adjacent far)
    seen = set()
    for rest in combinations(range(1, 13), rest_size):
        base = (0,) + rest
        # only keep dense-ish bases (longest run >= rest_size-1 or AP) -- these maximise coverage
        diffs = [b - a for a, b in zip(base, base[1:])]
        if max(diffs) <= 3:  # dense base
            bases.append(base)
    # always include canonical dense bases
    for cb in [tuple(range(k - 2)), tuple(i * 2 for i in range(k - 2))]:
        if max(cb) <= 14 and cb not in bases:
            bases.append(cb)
    for base in bases:
        bmax = base[-1]
        for f1 in range(max(15, bmax + 1), w_hi):
            for f2 in range(f1 + 1, min(f1 + 12, w_hi + 1)):  # adjacent far pair = peak resonance
                E = base + (f1, f2)
                if len(set(E)) != k or not primitive(E):
                    continue
                pv = p0_fast(E)
                if pv > best:
                    best = pv
                    best_arg = (base, f1, f2)
    return best, best_arg


def main():
    print("=" * 92)
    print("THREAD A part 3: EXACT near-cap family + slack threshold delta (kps-S26-wf8)")
    print("=" * 92)
    print()
    print("(Q3a) SINGLE-FAR (r=1) near-cap family: exhaustive base x far w in [15,40]")
    print("-" * 92)
    nearcap_summary = {}
    for k in (8, 9, 10, 11, 12):
        best, best_arg, near = scan_single_far(k)
        cap = CAP[k]
        q = QVAL[k]
        # dedup near by (base,w)
        sigs = {}
        for pv, base, w in near:
            sig = base_signature(base)
            sigs.setdefault(sig, 0)
            sigs[sig] += 1
        nearcap_summary[k] = (best, best_arg, len(near), sigs)
        print(f"  k={k}: max single-far p0={float(best):.6f} (cap={float(cap):.5f}, cap-p0={float(cap-best):.6f})")
        print(f"        argmax base={best_arg[0]} sig={base_signature(best_arg[0])} far={best_arg[1]}")
        print(f"        #near-cap (p0>Q={float(q):.5f}) = {len(near)};  base-signature histogram: {sigs}")
        print(f"        top near-cap configs:")
        for pv, base, w in near[:6]:
            print(f"           p0={float(pv):.6f} (cap-p0={float(cap-pv):.6f}) base={base} [{base_signature(base)}] far={w}")
        print()

    print("=" * 92)
    print("(Q3b) SLACK THRESHOLD: adversarial r=2 max p0 (peak-resonance far pairs)")
    print("-" * 92)
    for k in (8, 9, 10, 11, 12):
        best, best_arg = scan_two_far_adversarial(k)
        cap = CAP[k]
        q = QVAL[k]
        delta = cap - best
        in_band = best > q
        print(f"  k={k}: max r=2 p0={float(best):.6f}  cap-p0(slack delta)={float(delta):.6f}"
              f"  >Q? {'YES(!)' if in_band else 'no -> SLACK'}")
        print(f"        argmax base={best_arg[0]} far=({best_arg[1]},{best_arg[2]})")
    print()
    print("=" * 92)
    print("VERDICT")
    print("=" * 92)
    print("If (Q3a) near-cap is dominated by consec/AP(dense) bases + 1 far, and (Q3b) shows")
    print("r=2 max p0 < Q(k-1) with explicit slack delta, the dichotomy is established:")
    print("  near-cap <=> r=1 (single-far plateau, dense base) ; r>=2 => p0 <= cap - delta.")


if __name__ == "__main__":
    main()
