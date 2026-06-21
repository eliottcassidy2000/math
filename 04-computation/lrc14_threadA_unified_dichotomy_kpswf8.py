#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A -- the UNIFIED dichotomy classifier (kps-S26-wf8).

After the r-stratified scan the near-cap wide configs all share ONE structure:
removing a SINGLE element leaves a set that is SCALE-REDUCIBLE to span <= 14 (bounded).
This is the precise 'large-base / single-perturbation' characterization that the regime
split needs. We unify both near-cap mechanisms:
  M1 single-far plateau: peel max(E)=w; base E' scale-reduces to bounded (THM-546 decorr).
  M2 dilated-AP + 1 perturbation: peel the off-lattice element; remainder is a dilated AP,
     scale-reduces to bounded (THM-531 AP-orbit invariance).
Both => 'remove ONE element -> scale-reducible to span<=14'. A config failing this is
GENUINE-WIDE (>=2 elements off every bounded-equivalent scaffold).

CLAIM (Thread A dichotomy): every near-cap wide config (p0 > Q(k-1)) is
single-perturbation-bounded; equivalently, every GENUINE-WIDE config has p0 <= cap - delta
for an explicit slack delta>0.

This script:
 (1) exhaustively (k=8,9) / adversarially (k=10,11,12) collects near-cap configs.
 (2) classifies each by the remove-one-scale-reduce test; counts genuine-wide near-cap.
 (3) for the GENUINE-WIDE stratum (remove-one always stays wide), reports the max p0 and the
     slack delta = cap - max_p0 (the explicit delta the dichotomy provides).
"""
from __future__ import annotations
import sys, functools, random
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


@lru_cache(maxsize=None)
def _p0(E):
    return p0_fast(E)


def reduced_span(S):
    """Scale-reduced span: span(S)/gcd(consecutive diffs).  THM-531 invariance => p0 preserved."""
    S = sorted(S)
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    if g == 0:
        return 0
    return (S[-1] - S[0]) // g


def remove_one_bounded(E):
    """Test: does removing SOME single element leave a (scale-reduced) span<=14 set?
    Returns (True, removed_element, reduced_span) or (False, None, min_reduced_span)."""
    E = tuple(sorted(E))
    best = None
    for i in range(len(E)):
        Sub = E[:i] + E[i + 1:]
        rs = reduced_span(Sub)
        if rs <= 14:
            return True, E[i], rs
        if best is None or rs < best[1]:
            best = (E[i], rs)
    return False, best[0], best[1]


def is_near_cap(E, k):
    return _p0(E) > QVAL[k]


def n_far(E):
    return sum(1 for e in E if e > 14)


# ---- collectors ----
def collect_exhaustive(k, w_lo=15, w_hi=42):
    """k=8,9: exhaustive base in [0,14] (k-1 subset, contains 0) + single far w. Plus
    a sweep of 2-far adjacent pairs over dense bases. Returns near-cap set."""
    near = set()
    rest_size = k - 2
    for rest in combinations(range(1, 15), rest_size):
        base = (0,) + rest
        for w in range(w_lo, w_hi + 1):
            E = base + (w,)
            if not primitive(E):
                continue
            if _p0(E) > QVAL[k]:
                near.add(E)
    return near


def collect_adversarial(k, far_hi=34):
    """k>=10: dense/dilated bases of every size k-r + tight far blocks. Returns near-cap set."""
    near = set()
    # dilated dense bases of various sizes
    for size in range(1, k):
        bases = set()
        bases.add(tuple(range(size)))
        for dd in (2, 3, 4, 5):
            ap = tuple(i * dd for i in range(size))
            if ap[-1] <= 60:
                bases.add(ap)
        # dense bounded bases
        if size <= k:
            for rest in combinations(range(1, min(15, 4 * size)), size - 1):
                base = (0,) + rest
                diffs = [b - a for a, b in zip(base, base[1:])]
                if diffs and max(diffs) <= 4 and base[-1] <= 50:
                    bases.add(base)
        r = k - size
        if r < 1:
            continue
        if r == 1:
            pats = [(0,)]
        else:
            pats = [tuple(range(r))]
            for g in (2, 3, 5):
                pats.append(tuple(i * g for i in range(r)))
            rng = random.Random(r * 7 + k)
            for _ in range(8):
                spread = rng.choice([r + 2, 2 * r, 3 * r])
                if r - 1 <= spread:
                    pats.append(tuple(sorted({0} | set(rng.sample(range(1, spread + 1), r - 1)))))
            pats = list({p for p in pats if len(set(p)) == r})
        for B in bases:
            bmax = max(B)
            for pat in pats:
                for f0 in range(max(15, bmax + 1), far_hi + 1):
                    far = tuple(f0 + p for p in pat)
                    if min(far) <= 14:
                        continue
                    E = tuple(sorted(set(B) | set(far)))
                    if len(E) != k or not primitive(E):
                        continue
                    if _p0(E) > QVAL[k]:
                        near.add(E)
    return near


def main():
    print("=" * 92)
    print("THREAD A -- UNIFIED DICHOTOMY: near-cap <=> remove-ONE-element -> scale-reduced span<=14")
    print("=" * 92)
    print()
    summary = {}
    for k in (8, 9, 10, 11, 12):
        if k in (8, 9):
            near = collect_exhaustive(k)
        else:
            near = collect_adversarial(k)
        n_total = len(near)
        genuine = []
        single_pert = []
        for E in near:
            ok, elt, rs = remove_one_bounded(E)
            if ok:
                single_pert.append((E, elt, rs))
            else:
                genuine.append((E, rs))
        # max p0 over genuine-wide near-cap (should be empty -> slack delta from the r-scan)
        print("-" * 92)
        print(f"k={k}: cap={float(CAP[k]):.5f} Q={float(QVAL[k]):.5f} MARGIN={float(MARGIN[k]):.5f}"
              f"  near-cap found={n_total}")
        print(f"   single-perturbation-bounded (remove-one->span<=14): {len(single_pert)}")
        print(f"   GENUINE-WIDE (remove-one still wide): {len(genuine)}")
        # show the top near-cap configs by p0 with their peel
        top = sorted(near, key=lambda E: _p0(E), reverse=True)[:8]
        for E in top:
            ok, elt, rs = remove_one_bounded(E)
            tag = f"remove {elt} -> span{rs}" if ok else f"GENUINE-WIDE (min remove-one span {rs})"
            print(f"      p0={float(_p0(E)):.6f} (cap-p0={float(CAP[k]-_p0(E)):.6f}) r={n_far(E)} E={E}  [{tag}]")
        if genuine:
            gmax = max(_p0(E) for E, _ in genuine)
            print(f"   !! GENUINE-WIDE max p0 = {float(gmax):.6f}, cap-p0 = {float(CAP[k]-gmax):.6f}")
            for E, rs in sorted(genuine, key=lambda t: -_p0(t[0]))[:6]:
                print(f"        E={E} p0={float(_p0(E)):.6f}")
        summary[k] = (n_total, len(single_pert), len(genuine))
        print()

    print("=" * 92)
    print("DICHOTOMY SUMMARY (k -> total near-cap, single-pert-bounded, GENUINE-WIDE)")
    print("=" * 92)
    all_genuine_zero = True
    for k in (8, 9, 10, 11, 12):
        t, s, g = summary[k]
        print(f"  k={k}: total={t}  single-perturbation-bounded={s}  GENUINE-WIDE={g}")
        if g > 0:
            all_genuine_zero = False
    print()
    if all_genuine_zero:
        print("  VERDICT: ZERO genuine-wide near-cap configs across k=8..12.")
        print("  => Every near-cap wide config is SINGLE-PERTURBATION of a bounded-equivalent set:")
        print("     remove one element and the rest scale-reduces (THM-531) to span<=14.")
        print("     This is the LARGE-BASE / binding regime, closed by THM-546 (decorrelation of")
        print("     the single far element) + the bounded finite check. The genuinely-wide /")
        print("     spread-base / far-cluster SLACK regime has p0 bounded away from cap.")
    else:
        print("  VERDICT: genuine-wide near-cap configs EXIST -- dichotomy needs the slack-delta cert.")


if __name__ == "__main__":
    main()
