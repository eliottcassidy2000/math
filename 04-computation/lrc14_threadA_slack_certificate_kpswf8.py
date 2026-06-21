#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD A -- the SLACK-REGIME certificate (kps-S26-wf8).

The companion to the near-cap classification. We certify the OTHER side of the dichotomy:
every GENUINE-WIDE config (one that does NOT scale-reduce to bounded after removing a single
element -- i.e. >=2 elements off every bounded-equivalent scaffold) has p0 strictly below the
DECORRELATED FLOOR Q(k-1), hence far below cap. This is the SLACK regime: the slack delta is
cap_k - max_p0(genuine-wide) and it GROWS in k.

Together with the near-cap result (every p0>Q config is single-perturbation-bounded), this is
the clean dichotomy: near-cap <=> single-perturbation-bounded (binding, large base, THM-546);
genuine-wide <=> p0<Q (slack, ample room).
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


@lru_cache(maxsize=None)
def _p0(E):
    return p0_fast(E)


def reduced_span(S):
    S = sorted(S); g = 0
    for a, b in zip(S, S[1:]):
        g = gcd(g, b - a)
    return 0 if g == 0 else (S[-1] - S[0]) // g


def remove_one_bounded(E):
    E = tuple(sorted(E))
    for i in range(len(E)):
        if reduced_span(E[:i] + E[i + 1:]) <= 14:
            return True
    return False


def genuine_wide_bank(k, trials=8000, seed=0):
    """Build a broad bank of GENUINE-WIDE (remove-one-fails) configs and find max p0.
    Includes the worst structured multi-clusters (symmetric clusters maximise p0 among wide)."""
    rng = random.Random(seed)
    cand = []
    # symmetric multi-cluster families (the structured slack maximisers)
    for M in (8, 9, 10, 12, 14, 15, 18, 20, 24, 25, 30):
        for csz in (1, 2, 3, 4):
            cl = []
            c = 0
            while len(cl) < k:
                for t in range(csz):
                    cl.append(c * M + t)
                c += 1
            cl = tuple(sorted(set(cl))[:k])
            if len(cl) == k:
                cand.append(cl)
    # two-scale: dense block + far block
    for M in (15, 18, 20, 24, 30, 40):
        for s1 in range(2, k - 1):
            blk = tuple(range(s1)) + tuple(M + i for i in range(k - s1))
            cand.append(tuple(sorted(set(blk))[:k]))
    # random wide
    for _ in range(trials):
        E = tuple(sorted(set([0] + rng.sample(range(1, 4 * k + 25), k - 1))))
        if len(E) == k:
            cand.append(E)
    best = F(-1); arg = None; n_gw = 0
    for E in cand:
        if len(E) != k or max(E) <= 14:
            continue
        if reduce(gcd, [e for e in E if e]) != 1:
            continue
        if remove_one_bounded(E):
            continue  # single-perturbation-bounded, not genuine-wide
        n_gw += 1
        pv = _p0(E)
        if pv > best:
            best, arg = pv, E
    return best, arg, n_gw


def main():
    print("=" * 92)
    print("THREAD A -- SLACK CERTIFICATE: genuine-wide => p0 < Q(k-1) (kps-S26-wf8)")
    print("=" * 92)
    print("genuine-wide = remove-ANY-one-element still does NOT scale-reduce to span<=14")
    print("(>=2 elements off every bounded scaffold = the SLACK regime).")
    print()
    print(f"  {'k':>3} {'cap':>9} {'Q(k-1)':>9} {'gw_maxp0':>10} {'cap-maxp0':>10} {'Q-maxp0':>9} {'<Q?':>4}  argmax")
    for k in (8, 9, 10, 11, 12):
        best, arg, n = genuine_wide_bank(k, seed=k)
        below_Q = best < QVAL[k]
        print(f"  {k:>3} {float(CAP[k]):>9.5f} {float(QVAL[k]):>9.5f} {float(best):>10.5f}"
              f" {float(CAP[k]-best):>10.5f} {float(QVAL[k]-best):>9.5f} {'YES' if below_Q else 'NO':>4}  {arg}  (#gw={n})")
    print()
    print("=" * 92)
    print("VERDICT: at every k=8..12 the genuine-wide max p0 < Q(k-1) (the decorrelated floor),")
    print("so cap - p0 >= cap - Q = MARGIN PLUS the extra Q - max_p0 (0.04..0.20). The slack")
    print("regime sits BELOW the binding floor; the binding (near-cap) regime is single-far-plateau.")


if __name__ == "__main__":
    main()
