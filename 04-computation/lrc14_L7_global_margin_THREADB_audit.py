#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_L7_global_margin_THREADB_audit.py  (THREAD B audit, 2026-06-21)

Adversarial global check of the L7 atlas: over ALL bounded bases (span<=B) AND the
small-denominator resonance atlas (p/q in (1,2.15], q<=QMAX), is

    p0_inf(B, p/q)  <  cap_k

with a STRICTLY POSITIVE margin?  Reports the GLOBAL sup p0_inf per k and the worst
(base, p/q).  This is the binding finite check behind L7 (the balanced two-far case).

Resonance R(p/q) decays ~1/q (proved D<=14/p), so the sup over all q is attained at
small q; QMAX=6 captures p/q in {2/1,3/2,4/3,5/3,5/4,7/4,...} which dominate.
EXACT rationals throughout.
"""
import os, importlib.util, sys, time
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

_d = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("atl", os.path.join(_d, "lrc_q108_L7_resonance_atlas_kps.py"))
atl = importlib.util.module_from_spec(spec); spec.loader.exec_module(atl)

CAP = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91)}
QMAX = 6

def prim(E):
    return reduce(gcd, [x for x in E if x != 0], 0) == 1

def atlas_pairs():
    out = []
    for q in range(1, QMAX + 1):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1:
                continue
            r = Fr(p, q)
            if Fr(1) < r <= Fr(43, 20):
                out.append((p, q))
    return out

def main():
    pairs = atlas_pairs()
    print(f"L7 GLOBAL MARGIN AUDIT: bounded bases x atlas {pairs}")
    # (k, nbase, base span bound) : nbase = k-2 (two far elements)
    for k, nbase, Bspan in [(8, 6, 12), (9, 7, 13), (10, 8, 14)]:
        cap = CAP[k]
        t0 = time.time()
        glob = Fr(0); arg = None; nbases = 0
        for tail in combinations(range(1, Bspan + 1), nbase - 1):
            B = (0,) + tail
            if not prim(B):
                continue
            nbases += 1
            for (p, q) in pairs:
                val = atl.p0_inf(list(B), p, q)
                if val > glob:
                    glob = val; arg = (B, p, q)
        dt = time.time() - t0
        print(f"k={k}: scanned {nbases} bases x {len(pairs)} ratios ({dt:.0f}s)")
        print(f"   GLOBAL sup p0_inf = {glob} = {float(glob):.5f}  @ base={arg[0]} p/q={arg[1]}/{arg[2]}")
        print(f"   cap_{k} = {cap} = {float(cap):.5f}   MIN margin = {float(cap-glob):.5f}   strictly under cap: {glob < cap}")

if __name__ == "__main__":
    main()
