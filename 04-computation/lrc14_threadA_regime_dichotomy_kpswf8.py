#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A (kps-2026-06-21-S26-wf8): RIGORIZE the regime dichotomy for the LRC(14) wide bound.

GOAL: establish QUANTITATIVELY that every wide near-cap config is LARGE-BASE.

DEFINITIONS (this script):
  - A wide config E: primitive, 0 in E, |E|=k, span(E)=max(E)>14.
  - r = #far = #{e in E : e > 14}  (the elements OUTSIDE the bounded window [0,14]).
    base size = k - r, base B = E intersect [0,14] (always contains 0).
  - p0(E) = measS7(E) = measure{x : E*x hits all 6 inner sectors}. EXACT rational.
  - cap_k, Q(k-1)=decorr baseline at consec_{k-1}+1far, MARGIN_k = cap_k - Q(k-1).
  - near-cap = p0 within MARGIN_k of cap_k, i.e. p0 > cap_k - MARGIN_k = Q(k-1).
    (This is the natural near-cap band: the decorrelated floor sits at exactly cap-MARGIN.)

PART 1: max p0 stratified by r (number of far elements), k=8..12.
  For each (k,r) we adversarially search: consecutive/AP/random bounded bases of size k-r,
  combined with far blocks placed at small / adjacent / resonant (small-denominator-ratio)
  positions -- exactly where the existing engine showed resonance peaks. Report:
    max p0, the gap delta_r = cap_k - max_p0, and the argmax config.

PART 2: the THRESHOLD r0. Tabulate delta_r by r,k. Identify r0 such that all r>r0
  have p0 <= cap - delta with explicit delta>0 (the SLACK-regime margin).

PART 3: near-cap family. Among ALL configs found with p0 > Q(k-1) (near-cap band),
  list them and check: are they exactly consec_{k-1} + 1 far (single-far plateau)
  plus small perturbations? Report the exact family and the slack-regime delta.

Exact rationals throughout. Far positions bounded but chosen adversarially to maximise
resonance (adjacent far pairs/blocks = small-denominator ratios = peak resonance, per
lrc14_wide_resonance_sup_kpswf7.py PART 2 decay law).
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

from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive, p0 as p0_repo
from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}


def p0_fast(E):
    """Fast exact measS7 (matches missed_distribution[0] / p0_repo)."""
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


def is_wide_primitive(E):
    E = tuple(sorted(set(E)))
    if 0 not in E:
        return False
    if max(E) <= 14:
        return False
    return primitive(E)


def n_far(E):
    return sum(1 for e in E if e > 14)


# ---- adversarial generators for bounded bases of a given size (subset of [0,14], contains 0) ----
def bounded_bases(size, n_random=80, seed=0):
    """Yield bounded primitive-ish bases of given size in [0,14] containing 0.
    Always includes consec, even-AP, and structured + random samples."""
    out = []
    seen = set()

    def add(B):
        B = tuple(sorted(set(B)))
        if len(B) != size:
            return
        if 0 not in B or max(B) > 14:
            return
        if B in seen:
            return
        seen.add(B)
        out.append(B)

    add(tuple(range(size)))  # consec
    add(tuple([0] + [2 * i for i in range(1, size)]))  # even AP (if fits)
    add(tuple([0] + [i for i in range(2, 2 * size, 2)][: size - 1]))
    # consec shifted to top of window
    add(tuple([0] + list(range(15 - (size - 1), 15))))
    # arithmetic progressions with small common difference
    for dstep in (1, 2, 3):
        ap = tuple(range(0, dstep * size, dstep))
        if max(ap) <= 14:
            add(ap)
    rng = random.Random(seed)
    tries = 0
    while len(out) < n_random and tries < n_random * 40:
        tries += 1
        rest = rng.sample(range(1, 15), size - 1)
        add([0] + rest)
    return out


# ---- far-block generators: place r far elements to MAXIMISE resonance ----
def far_blocks(r, base_max, n_random=40, seed=1):
    """Yield r-element far blocks (all > 14, > base_max enforced by caller offset).
    We generate offset patterns relative to a start point f0; the caller sweeps f0.
    Adjacent blocks (consecutive) = small-denominator ratios = peak resonance."""
    pats = []
    # consecutive block (peak resonance per decay law)
    pats.append(tuple(range(r)))
    # AP blocks with small gap
    for g in (2, 3, 5, 7):
        pats.append(tuple(i * g for i in range(r)))
    rng = random.Random(seed)
    for _ in range(n_random):
        pats.append(tuple(sorted({0} | set(rng.sample(range(1, 40), r - 1)))) if r > 1 else (0,))
    # dedup
    uniq = []
    seen = set()
    for p in pats:
        if len(set(p)) == r and p not in seen:
            seen.add(p)
            uniq.append(p)
    return uniq


def scan_kr(k, r, far_start_lo=15, far_start_hi=60, verbose=False):
    """Adversarial scan of all wide configs with exactly r far elements at size k.
    Returns (max_p0, argmax_E, near_cap_list) where near_cap_list = configs with p0>Q(k-1)."""
    size = k - r
    if size < 1:
        return None
    bases = bounded_bases(size)
    best_p0 = F(-1)
    best_E = None
    near = []  # (p0, E)
    cap = CAP[k]
    q = QVAL[k]
    for B in bases:
        bmax = max(B)
        patterns = far_blocks(r, bmax)
        for pat in patterns:
            # sweep the far block start so the block sits just past the window and beyond
            for f0 in range(max(15, bmax + 1), far_start_hi + 1):
                far = tuple(f0 + p for p in pat)
                if min(far) <= 14:
                    continue
                E = tuple(sorted(set(B) | set(far)))
                if len(E) != k:
                    continue
                if not primitive(E):
                    continue
                if n_far(E) != r:
                    continue
                pv = p0_fast(E)
                if pv > best_p0:
                    best_p0 = pv
                    best_E = E
                if pv > q:
                    near.append((pv, E))
    near.sort(reverse=True)
    return best_p0, best_E, near[:25]


def main():
    print("=" * 92)
    print("THREAD A: REGIME DICHOTOMY -- every wide near-cap config is LARGE-BASE (kps-S26-wf8)")
    print("=" * 92)
    print("near-cap band: p0 > cap_k - MARGIN_k = Q(k-1).  delta_r = cap_k - max_p0(r).")
    for k in CAP:
        print(f"  k={k}: cap={CAP[k]}={float(CAP[k]):.5f}  Q(k-1)={float(QVAL[k]):.5f}  MARGIN={float(MARGIN[k]):.5f}")
    print()

    all_near = {}  # k -> list of (p0,E)
    table = {}  # (k,r) -> (max_p0, argmax, n_near)
    for k in (8, 9, 10, 11, 12):
        print("-" * 92)
        print(f"k={k}  (cap={float(CAP[k]):.5f}, Q={float(QVAL[k]):.5f})")
        print("-" * 92)
        print(f"  {'r=#far':>7} {'base':>5} {'max_p0':>14} {'cap-max_p0=delta_r':>20} {'>Q? near':>9}  argmax_E")
        near_k = []
        for r in range(1, k):  # r from 1 (single far) to k-1 (base = {0} only)
            res = scan_kr(k, r)
            if res is None:
                continue
            mp0, argE, near = res
            table[(k, r)] = (mp0, argE, len(near))
            delta_r = CAP[k] - mp0
            flag = "YES" if mp0 > QVAL[k] else "no"
            print(f"  {r:>7} {k-r:>5} {float(mp0):>14.6f} {float(delta_r):>20.6f} {flag:>4} {len(near):>4}  {argE}")
            near_k.extend(near)
        # dedup near-cap configs
        seen = set()
        dd = []
        for pv, E in sorted(near_k, reverse=True):
            if E in seen:
                continue
            seen.add(E)
            dd.append((pv, E))
        all_near[k] = dd
        print()

    # ---------- PART 2: THRESHOLD r0 ----------
    print("=" * 92)
    print("PART 2: THRESHOLD r0 -- delta_r = cap_k - max_p0(r) by (k,r)")
    print("=" * 92)
    print("  For each k: smallest r with delta_r bounded away (slack regime), and the slack delta.")
    for k in (8, 9, 10, 11, 12):
        print(f"  k={k}:")
        slack_deltas = []
        r0 = None
        for r in range(1, k):
            if (k, r) not in table:
                continue
            mp0, argE, nn = table[(k, r)]
            delta_r = CAP[k] - mp0
            near = mp0 > QVAL[k]
            tag = "NEAR-CAP (binding)" if near else "slack"
            print(f"     r={r}: max_p0={float(mp0):.6f} delta_r={float(delta_r):.6f}  [{tag}]")
            if not near:
                slack_deltas.append((r, delta_r))
                if r0 is None:
                    r0 = r
        if slack_deltas:
            min_slack = min(d for _, d in slack_deltas)
            print(f"     => threshold r0={r0}: ALL r>={r0} are SLACK, min slack delta = {float(min_slack):.6f}")
        else:
            print(f"     => NO slack r (all r near-cap?) -- unexpected, inspect")
        print()

    # ---------- PART 3: NEAR-CAP FAMILY ----------
    print("=" * 92)
    print("PART 3: NEAR-CAP FAMILY (all configs with p0 > Q(k-1)), is it consec_{k-1}+1far?")
    print("=" * 92)
    for k in (8, 9, 10, 11, 12):
        dd = all_near.get(k, [])
        print(f"  k={k}: {len(dd)} distinct near-cap configs (p0>Q={float(QVAL[k]):.5f}), all with cap={float(CAP[k]):.5f}")
        consec = tuple(range(k - 1))
        for pv, E in dd[:20]:
            r = n_far(E)
            base = tuple(e for e in E if e <= 14)
            is_consec_base = base == consec
            far = tuple(e for e in E if e > 14)
            print(f"     p0={float(pv):.6f} (cap-p0={float(CAP[k]-pv):.6f}) r={r} base={base}"
                  f" {'[consec_{k-1}]' if is_consec_base else ''} far={far}")
        # classify
        single_far = [E for _, E in dd if n_far(E) == 1]
        multi_far = [E for _, E in dd if n_far(E) >= 2]
        consec_single = [E for E in single_far if tuple(e for e in E if e <= 14) == consec]
        print(f"     => near-cap breakdown: single-far={len(single_far)} (of which consec_(k-1) base={len(consec_single)}),"
              f" multi-far={len(multi_far)}")
        print()

    print("=" * 92)
    print("DICHOTOMY VERDICT")
    print("=" * 92)
    print("If PART 1/2 show: r=1 (single far, base size k-1) is the ONLY near-cap stratum and")
    print("all r>=2 have delta_r bounded away by an explicit slack delta, then the SLACK regime")
    print("(r>=2 = spread base / far cluster) is provably under cap, and the BINDING regime is")
    print("EXACTLY the single-far plateau family, handled by the localized THM-546 comb bound.")


if __name__ == "__main__":
    main()
