#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: does the Delsarte dual gK8 close the WIDE region directly?

gK8 = (10,0,0,1,0,0,10) weights the MISS-DISTRIBUTION q=(q0,...,q6), q_t = meas{exactly t of the
6 inner sectors missed}; q0 = p0 = measS7. L_yK8 = 10*q0 + q3 + 10*q6. The pointwise majorant
10*p0 <= L_yK8 is trivial (q3,q6>=0); the CONTENT is the moment bound max_E L_yK8 <= 10*cap_k,
which gives p0 <= cap. mac-mini: this MAY close the wide region (clears E* + 1600 wide, 30% slack).

This script computes the EXACT miss-distribution q and L_yK8 for the BINDING WIDE configs:
  - genuine-wide maximizers (HYP-2807): k=10,11,12 generalized doublets
  - mac-mini's k=12 breaker E*
  - single-far near-cap plateau (consec_{k-1}+far) -- the THM-563 binding family
  - dilated even-AP configs
and checks L_yK8 <= 10*cap_k with margin. If it holds with healthy margin across ALL binding wide
families, gK8 closes the wide region directly (superseding dichotomy+doublet+R-tail). Exact.
"""
from __future__ import annotations
import sys, functools
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import CAP
ALL_INNER = 0b1111110


def miss_dist(E):
    """Exact miss-distribution q=(q0,...,q6): q_t = meas{phi: exactly t inner sectors missed}."""
    nz = [int(x) for x in E if x]
    if not nz:
        return [F(0)] * 7
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
    q = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nz:
            mask |= 1 << ((e * midnum // den2) % 7)
        occ = bin(mask & ALL_INNER).count("1")
        missed = 6 - occ
        q[missed] += F(hi - lo, d)
    return q


def LyK8(E):
    q = miss_dist(E)
    return 10 * q[0] + q[3] + 10 * q[6], q[0]


def check(E, k, label):
    ly, p0 = LyK8(E)
    cap = CAP[k]
    bound = 10 * cap
    ok = ly <= bound
    print(f"   {label:28s} {E}")
    print(f"      p0={float(p0):.5f}  L_yK8={float(ly):.5f}  10cap={float(bound):.5f}  "
          f"margin={float(bound-ly):+.5f}  {'OK' if ok else '*** VIOLATES ***'}")
    return ok


def main():
    print("=" * 80)
    print("Does gK8 (10,0,0,1,0,0,10) close the WIDE region? L_yK8 <= 10*cap on binding wide configs")
    print("claude-opus 2026-06-22")
    print("=" * 80)
    allok = True
    # genuine-wide maximizers (HYP-2807)
    gw = {10: (0, 1, 3, 5, 7, 9, 11, 13, 15, 17),
          11: (0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22),
          12: (0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 18, 20)}
    print("\n[A] genuine-wide maximizers (HYP-2807):")
    for k, E in gw.items():
        allok &= check(E, k, f"gw-max k={k}")
    # E* breaker
    print("\n[B] mac-mini k=12 breaker E*:")
    allok &= check((0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18), 12, "E* (k=12)")
    # single-far near-cap plateau (binding family, consec_{k-1}+far)
    print("\n[C] single-far near-cap plateau (consec_{k-1} + far, THM-563 binding):")
    for k in (10, 11, 12):
        E = tuple(range(k - 1)) + (3 * k,)  # consec_{k-1} + one far
        allok &= check(E, k, f"single-far k={k}")
    # dilated even-AP + far doublet (the dilated binding case)
    print("\n[D] dilated even-AP families:")
    for k in (10, 12):
        E = tuple(range(0, 2 * (k - 1), 2)) + (2 * (k - 1) + 4,)
        allok &= check(E, k, f"even-AP+far k={k}")
    print("\n" + "=" * 80)
    print(f"VERDICT: gK8 L_yK8<=10cap on all tested binding wide configs: {'HOLDS' if allok else 'FAILS somewhere'}")
    print("If HOLDS with margin: gK8 may close the WHOLE wide region via the miss-distribution moment")
    print("bound -- a cleaner completion than dichotomy+doublet+R-tail. Needs: max_E L_yK8<=10cap over")
    print("ALL wide E (the Delsarte LP feasibility for wide q-moments).")


if __name__ == "__main__":
    main()
