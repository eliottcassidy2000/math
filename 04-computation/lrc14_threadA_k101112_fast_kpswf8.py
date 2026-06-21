#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A -- FAST k=10,11,12 near-cap collection + unified dichotomy (kps-S26-wf8).

The full collect_adversarial is too slow (LCM blowup at large far/spread). We use the
THM-546 decay law to JUSTIFY a tight far window: resonance peaks for far elements ADJACENT
to the base/each other (small-denominator ratios); large/spread far DECORRELATES -> lower p0.
So the near-cap sup lives in: far values just past the window, tight far blocks.

We collect near-cap (p0 > Q) over:
  - all base SIZES s = k-r, r=1..k-1
  - dense/dilated bounded bases (consec, 2..6 dilations, all max-gap<=3 bases with span<=24)
  - tight far blocks (consecutive / small-gap AP), far start in [15, 30]
Then classify each by remove-one-scale-reduce (the unified dichotomy test).
Value cap keeps LCM small => fast & exact. This is a SUP search (we want the worst near-cap),
not a full enumeration, and the decay law guarantees the sup is in this window.
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
    S = sorted(S)
    g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def remove_one_bounded(E):
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


def n_far(E):
    return sum(1 for e in E if e > 14)


def dense_bases(size, span_cap=24):
    out = set()
    out.add(tuple(range(size)))
    for dd in (2, 3, 4, 5, 6):
        ap = tuple(i * dd for i in range(size))
        if ap[-1] <= span_cap:
            out.add(ap)
    # all max-gap<=3 bases with span<=span_cap (dense / dilated-dense)
    if size >= 2:
        for rest in combinations(range(1, span_cap + 1), size - 1):
            base = (0,) + rest
            diffs = [b - a for a, b in zip(base, base[1:])]
            if max(diffs) <= 3:
                out.add(base)
    else:
        out.add((0,))
    return [b for b in out if 0 in b and len(b) == size]


def far_patterns(r):
    if r == 1:
        return [(0,)]
    pats = [tuple(range(r))]
    for g in (2, 3, 4, 5):
        pats.append(tuple(i * g for i in range(r)))
    rng = random.Random(r)
    for _ in range(6):
        spread = rng.choice([r + 1, r + 3, 2 * r])
        if r - 1 <= spread:
            pats.append(tuple(sorted({0} | set(rng.sample(range(1, spread + 1), r - 1)))))
    return list({p for p in pats if len(set(p)) == r})


def collect_near(k, far_hi=30, value_cap=70):
    near = set()
    for r in range(1, k):
        size = k - r
        bases = dense_bases(size)
        pats = far_patterns(r)
        for B in bases:
            bmax = max(B)
            for pat in pats:
                for f0 in range(max(15, bmax + 1), far_hi + 1):
                    far = tuple(f0 + p for p in pat)
                    if min(far) <= 14 or max(far) > value_cap:
                        continue
                    E = tuple(sorted(set(B) | set(far)))
                    if len(E) != k or not primitive(E):
                        continue
                    if _p0(E) > QVAL[k]:
                        near.add(E)
    return near


def main():
    print("=" * 92)
    print("THREAD A FAST k=10,11,12: near-cap collection + unified remove-one dichotomy (kps-S26-wf8)")
    print("=" * 92)
    print("(far window tight by THM-546 decay; SUP search for the worst near-cap config)")
    print()
    summary = {}
    for k in (10, 11, 12):
        near = collect_near(k)
        genuine, single = [], []
        for E in near:
            ok, elt, rs = remove_one_bounded(E)
            (single if ok else genuine).append((E, elt, rs))
        print("-" * 92)
        print(f"k={k}: cap={float(CAP[k]):.5f} Q={float(QVAL[k]):.5f} MARGIN={float(MARGIN[k]):.5f}  near-cap={len(near)}")
        print(f"   single-perturbation-bounded={len(single)}  GENUINE-WIDE={len(genuine)}")
        for E in sorted(near, key=lambda E: _p0(E), reverse=True)[:8]:
            ok, elt, rs = remove_one_bounded(E)
            tag = f"remove {elt}->span{rs}" if ok else f"GENUINE-WIDE(min span {rs})"
            print(f"      p0={float(_p0(E)):.6f} (cap-p0={float(CAP[k]-_p0(E)):.6f}) r={n_far(E)} E={E} [{tag}]")
        if genuine:
            gmax = max(_p0(E) for E, _, _ in genuine)
            print(f"   !! GENUINE-WIDE max p0={float(gmax):.6f} cap-p0={float(CAP[k]-gmax):.6f}")
            for E, _, rs in sorted(genuine, key=lambda t: -_p0(t[0]))[:6]:
                print(f"        E={E} p0={float(_p0(E)):.6f} min-remove-one-span={rs}")
        summary[k] = (len(near), len(single), len(genuine))
        print()
    print("=" * 92)
    print("SUMMARY k -> (near-cap, single-pert-bounded, GENUINE-WIDE)")
    for k in (10, 11, 12):
        print(f"  k={k}: {summary[k]}")
    allzero = all(summary[k][2] == 0 for k in summary)
    print()
    print("  VERDICT:", "ZERO genuine-wide near-cap at k=10,11,12 -- dichotomy holds." if allzero
          else "genuine-wide near-cap EXIST -- inspect.")


if __name__ == "__main__":
    main()
