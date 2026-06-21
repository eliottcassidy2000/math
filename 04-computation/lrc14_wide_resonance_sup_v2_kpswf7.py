#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THREAD 1 (kps-S24-wf7) v2: GENUINE sup of the wide signed resonance error, WITH the span>14
constraint correctly enforced (the v1 sup of 0.17 was contaminated by NON-wide sets like consec_9).

A config E (|E|=k, contains 0, primitive) is WIDE iff span(E)=max(E) > 14. For a wide E we split it
at the largest internal gap into a BOUNDED base cluster B (re-based to contain 0, spread <= 14) and r
FAR elements; the decorrelated baseline is p0_decorr = boundary_value_direct(B, r) (moment dual).

Two objects, both reported:
  (A) sup over wide E of the ACTUAL p0(E) vs cap_k  -- the REAL LRC(14) bound (MISTAKE-080 framing).
  (B) sup over wide E of the signed error p0(E) - p0_decorr, vs the AVAILABLE ROOM cap_k - p0_decorr
      (NOT the fixed Q-margin -- the room depends on the base). err < room <=> p0 < cap, equivalently.

Also: DECAY of the error in the gap (separation) and in the far-pair denominator, over WIDE configs.
Exact rationals.
"""
from __future__ import annotations
import sys, random, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}
MARGIN = {k: CAP[k] - QVAL[k] for k in CAP}
_dec = {}


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


def decorr(B, r):
    B = tuple(sorted(set(B)))
    key = (B, r)
    if key not in _dec:
        # re-base B to contain 0 with spread, primitive normalize
        b0 = tuple(x - B[0] for x in B)
        _dec[key] = boundary_value_direct(b0, r)
    return _dec[key]


def split_wide(E):
    """Split a sorted wide config at its largest internal gap into (base_cluster, far_list).
    Returns (B, far) where B is the larger-density bounded cluster containing the most elements
    within a span<=14 window, far = the rest. We use the convention: base = the maximal prefix
    cluster of consecutive-within-14 elements; far = the tail beyond the dominant gap.
    For robustness we choose the split that MAXIMIZES the resulting decorr baseline (worst case)."""
    E = sorted(E)
    n = len(E)
    best = None
    for cut in range(1, n):  # B = E[:cut], far = E[cut:]
        B = E[:cut]; far = E[cut:]
        if max(B) - min(B) > 14:
            continue  # base must be bounded
        # far elements must be separated (the gap to base big enough they decorrelate-ish);
        # require the gap from base to first far > 14 (genuine separation) OR far spread large
        gap = far[0] - B[-1]
        r = len(far)
        d = decorr(B, r)
        if best is None or d > best[0]:
            best = (d, tuple(B), tuple(far), gap)
    return best  # (decorr_value, B, far, gap)


def main():
    print("=" * 80)
    print("THREAD 1 v2: SUP wide signed error WITH span>14 enforced (kps-S24-wf7)")
    print("=" * 80)
    for k in (8, 9, 10, 11, 12):
        print(f"  k={k}: cap={float(CAP[k]):.5f} Q(k-1)={float(QVAL[k]):.5f} margin(cap-Q)={float(MARGIN[k]):.5f}")
    print()

    for k in (8, 9, 10):
        cap = float(CAP[k])
        print("-" * 80)
        print(f"k={k}: scan WIDE configs (span>14, |E|={k}, primitive, 0 in E)")
        print("-" * 80)
        random.seed(100 + k)
        sup_p0 = 0.0; arg_p0 = None
        sup_err = -1.0; arg_err = None; room_at = None
        worst_overcap = 0
        N = 0
        # systematic: bounded base cluster of size b (subset of 0..14) + (k-b) far elements
        # far placed beyond a gap G past the base, spread out. Enumerate structured + random.
        for b in range(max(1, k - 3), k):  # base sizes (most mass in base)
            rfar = k - b
            # base candidates: consec, even-AP, random within 0..14
            base_cands = [tuple(range(b)), tuple([0] + [2 * i for i in range(1, b)])]
            for _ in range(40):
                bb = tuple(sorted(set([0] + random.sample(range(1, 15), b - 1))))
                if len(bb) == b and primitive(bb):
                    base_cands.append(bb)
            for B in base_cands:
                bmax = max(B)
                # far cluster: start at gap G past bmax so span>14, vary internal far structure
                for G in (1, 2, 3, 5, 8, 13, 21):
                    start = bmax + G
                    if start <= 14:
                        # ensure overall span>14
                        start = 15
                    far_structs = []
                    if rfar == 1:
                        far_structs = [(start,), (start + 7,), (start + 1,)]
                    elif rfar == 2:
                        far_structs = [(start, start + 1), (start, start + 7), (start, start + 2),
                                       (start, start + 3), (start, start * 2 if start*2>start+1 else start+5)]
                    else:
                        far_structs = [tuple(start + j for j in range(rfar)),
                                       tuple(start + 7 * j for j in range(rfar)),
                                       tuple(start + 2 * j for j in range(rfar))]
                    for far in far_structs:
                        E = tuple(sorted(set(B + far)))
                        if len(E) != k:
                            continue
                        if max(E) - min(E) <= 14:
                            continue  # NOT wide
                        if reduce(gcd, [e for e in E if e]) != 1:
                            continue
                        pv = float(p0_fast(E))
                        N += 1
                        if pv > sup_p0:
                            sup_p0 = pv; arg_p0 = E
                        if pv > cap:
                            worst_overcap += 1
                        # error vs the decorr baseline for THIS split (B, rfar)
                        d = float(decorr(B, rfar))
                        err = pv - d
                        room = cap - d
                        if err > sup_err:
                            sup_err = err; arg_err = (E, B, far); room_at = room
        print(f"  scanned {N} wide configs")
        print(f"  (A) sup ACTUAL p0 = {sup_p0:.6f}  cap={cap:.6f}  margin_to_cap={cap-sup_p0:.6f}  over_cap={worst_overcap}")
        print(f"      argmax p0 E={arg_p0}")
        print(f"  (B) sup signed err (p0 - decorr) = {sup_err:.6f}")
        print(f"      at E={arg_err[0]} base={arg_err[1]} far={arg_err[2]}")
        print(f"      available room (cap - decorr) there = {room_at:.6f}  => err < room? {sup_err < room_at}")
        print(f"      err < fixed Q-margin {float(MARGIN[k]):.5f}? {sup_err < float(MARGIN[k])}")
        print()

    print("=" * 80)
    print("KEY: the REAL bound is (A) sup p0 < cap. (B)'s error can exceed the FIXED Q-margin")
    print("but the relevant comparison is err < (cap - decorr) for the actual split <=> p0 < cap.")
    print("=" * 80)


if __name__ == "__main__":
    main()
