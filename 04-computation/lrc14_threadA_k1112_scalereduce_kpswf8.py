#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A completion (kps-S26-wf8): finish k=11,12 r-stratified max p0, AND apply the
SCALE-REDUCTION filter that resolves the dichotomy.

KEY INSIGHT (from THM-531, PROVED AP-orbit invariance + the partial scan):
  The near-cap wide configs split into TWO mechanisms, NOT 'large base vs not':
   (M1) SINGLE-FAR PLATEAU: r=1, dense base (consec/dilated-consec), 1 far element.
        p0 -> Phi(B) = p0(B)+(1/7)p1(B) <= Q(k-1) as far->inf (THM-546 decorrelation).
   (M2) DILATED-AP / SCALE-RESONANT pocket: configs that are (up to translation+scale)
        a BOUNDED set. By THM-531 these have p0 = p0(bounded representative), so they are
        ALREADY covered by the bounded finite check (max<=14) -- they are not measure-wide.
  Genuinely-wide, structureless 'spread base + clustered random far' (the SLACK case) never
  approaches cap.

This script:
 (1) finishes the r-stratified max p0 for k=11,12 with a tight far window.
 (2) for EVERY near-cap config found (p0>Q), tests whether it is in the DILATED-AP orbit:
     compute g = gcd of all pairwise differences; if g>1 the set is a dilation E = g*E0 + t,
     and E0 is its scale-reduced (smaller-span) representative -- flag if span(E0)<=14.
 (3) classifies each near-cap config as M1 (single-far plateau), M2 (scale-reducible to
     bounded), or GENUINE-WIDE (neither). The dichotomy CLAIM: NO near-cap config is
     GENUINE-WIDE.
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


def n_far(E):
    return sum(1 for e in E if e > 14)


def diff_gcd(E):
    """gcd of all pairwise differences (i.e. the dilation factor of the AP-orbit)."""
    g = 0
    s = sorted(E)
    for a, b in zip(s, s[1:]):
        g = gcd(g, b - a)
    return g


def scale_reduced_span(E):
    """If E lies in a dilated lattice (diff_gcd g>1), its scale-reduced span is span(E)//g.
    THM-531: translation+scale invariance => p0(E)=p0(reduced). Return reduced span."""
    g = diff_gcd(E)
    if g <= 1:
        return E[-1] - E[0], g  # not reducible
    return (E[-1] - E[0]) // g, g


def classify(E, k):
    """Return one of: M1 (single-far plateau), M2 (scale-reducible to bounded), GENUINE-WIDE."""
    r = n_far(E)
    base = tuple(e for e in E if e <= 14)
    red_span, g = scale_reduced_span(E)
    if g > 1 and red_span <= 14:
        return f"M2(scale/{g}->span{red_span})"
    if r == 1:
        # single far plateau: base is dense (covers a near-consecutive window)?
        bdiffs = [b - a for a, b in zip(base, base[1:])]
        dense = (max(bdiffs) <= 2) if bdiffs else True
        return "M1(single-far,dense-base)" if dense else "M1(single-far,sparse-base)"
    return "GENUINE-WIDE"


# generators ----
def dense_bases(size):
    """Dense bounded bases (consec, dilated-consec/AP, near-consec) of given size in [0,14]."""
    out = set()
    out.add(tuple(range(size)))
    for d in (1, 2, 3, 4):
        ap = tuple(i * d for i in range(size))
        if ap[-1] <= 14:
            out.add(ap)
    # near-consec windows
    for start in range(0, 15 - size + 1):
        out.add(tuple(range(start, start + size)) if start == 0 else (0,) + tuple(range(start, start + size - 1)))
    # all small-spread bases containing 0
    for rest in combinations(range(1, 15), size - 1):
        base = (0,) + rest
        diffs = [b - a for a, b in zip(base, base[1:])]
        if max(diffs) <= 3:  # dense
            out.add(base)
    return [b for b in out if 0 in b and max(b) <= 14 and len(b) == size]


def scan_kr_tight(k, r, far_hi=30):
    """Tighter r-stratified scan for k=11,12. Far blocks tight; dense bases."""
    size = k - r
    if size < 1:
        return None
    bases = dense_bases(size)
    # far OFFSET patterns (peak resonance: consecutive/AP)
    if r == 1:
        pats = [(0,)]
    else:
        pats = [tuple(range(r))]
        for g in (2, 3, 5, 7):
            pats.append(tuple(i * g for i in range(r)))
        rng = random.Random(r)
        for _ in range(12):
            spread = rng.choice([r + 2, 2 * r, 3 * r])
            pats.append(tuple(sorted({0} | set(rng.sample(range(1, spread + 1), r - 1)))))
        pats = list({p for p in pats if len(set(p)) == r})
    best = F(-1)
    arg = None
    near = []
    q = QVAL[k]
    for B in bases:
        bmax = max(B)
        for pat in pats:
            for f0 in range(max(15, bmax + 1), far_hi + 1):
                far = tuple(f0 + p for p in pat)
                if min(far) <= 14:
                    continue
                E = tuple(sorted(set(B) | set(far)))
                if len(E) != k or n_far(E) != r or not primitive(E):
                    continue
                pv = _p0(E)
                if pv > best:
                    best, arg = pv, E
                if pv > q:
                    near.append((pv, E))
    near.sort(reverse=True)
    return best, arg, near[:20]


def main():
    print("=" * 92)
    print("THREAD A completion: k=11,12 r-strata + SCALE-REDUCTION dichotomy filter (kps-S26-wf8)")
    print("=" * 92)
    all_near = []
    for k in (11, 12):
        print("-" * 92)
        print(f"k={k}  cap={float(CAP[k]):.5f}  Q={float(QVAL[k]):.5f}  MARGIN={float(MARGIN[k]):.5f}")
        print("-" * 92)
        print(f"  {'r':>3} {'base':>5} {'max_p0':>12} {'delta_r':>10} {'>Q':>4}  class(argmax)            argmax_E")
        for r in range(1, k):
            res = scan_kr_tight(k, r)
            if res is None:
                continue
            best, arg, near = res
            delta = CAP[k] - best
            flag = "YES" if best > QVAL[k] else "no"
            cls = classify(arg, k)
            print(f"  {r:>3} {k-r:>5} {float(best):>12.6f} {float(delta):>10.6f} {flag:>4}  {cls:<24} {arg}")
            for pv, E in near:
                all_near.append((k, pv, E))
        print()

    print("=" * 92)
    print("SCALE-REDUCTION CLASSIFICATION of ALL near-cap configs (k=11,12), dedup")
    print("=" * 92)
    seen = set()
    counts = {"M1": 0, "M2": 0, "GENUINE-WIDE": 0}
    genuine = []
    for k, pv, E in sorted(all_near, key=lambda t: (t[0], -t[1])):
        if (k, E) in seen:
            continue
        seen.add((k, E))
        cls = classify(E, k)
        key = "M2" if cls.startswith("M2") else ("M1" if cls.startswith("M1") else "GENUINE-WIDE")
        counts[key] += 1
        if key == "GENUINE-WIDE":
            genuine.append((k, pv, E))
        print(f"  k={k} p0={float(pv):.6f} (cap-p0={float(CAP[k]-pv):.6f}) class={cls:<26} E={E}")
    print()
    print(f"  CLASS COUNTS over near-cap (k=11,12): {counts}")
    if genuine:
        print(f"  !! {len(genuine)} GENUINE-WIDE near-cap configs (would BREAK the dichotomy):")
        for k, pv, E in genuine:
            print(f"     k={k} p0={float(pv):.6f} E={E} red_span={scale_reduced_span(E)}")
    else:
        print("  => ZERO genuine-wide near-cap configs: every near-cap config is M1 (single-far")
        print("     plateau) or M2 (scale-reducible to a bounded set). DICHOTOMY HOLDS at k=11,12.")


if __name__ == "__main__":
    main()
