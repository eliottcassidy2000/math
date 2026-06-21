#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_base_domination_audit_THREADC.py  (THREAD C audit, 2026-06-21)

ADVERSARIAL: the L7 atlas-exhaustive script only checks ONE base B per k (the dense
even AP [0,2,...,2(k-3)]). The closure DEPENDS on the claim that this base is the
WORST (maximizes p0_inf) over ALL primitive bounded bases B subset {0,...,14} of the
right size. If some other base B' gives a LARGER p0_inf at some atlas ratio p/q, the
single-base atlas is INSUFFICIENT.

This script EXHAUSTIVELY enumerates bounded bases B subset {0,...,14} with 0 in B
(WLOG translate so min=0), of size k-2 (so that |E|=k with 2 far elements), and for
each computes sup_{p/q in small atlas} p0_inf(B,p,q). It reports whether the even-AP
base attains the global sup, and the actual global-max base. EXACT rational.

For tractability we sweep the binding atlas ratios {2/1,3/2,4/3,5/3,5/4,7/5,7/4,5/2-cap}
and a few more, and restrict the base search appropriately per k. We focus on k=9,10
(binding) and k=8.
"""
import itertools, os, importlib.util
from fractions import Fraction as Fr
from math import gcd

_d = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location("atl", os.path.join(_d, "lrc_q108_L7_resonance_atlas_kps.py"))
atl = importlib.util.module_from_spec(spec); spec.loader.exec_module(atl)
p0_inf = atl.p0_inf
P2_decorrelated = atl.P2_decorrelated

CAP = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91), 11: Fr(66, 91), 12: Fr(6, 7)}

# binding small-q atlas ratios in (1, 43/20]
ATLAS = []
for q in range(1, 9):
    for p in range(q + 1, int(Fr(43, 20) * q) + 1):
        if gcd(p, q) == 1 and Fr(1) < Fr(p, q) <= Fr(43, 20):
            ATLAS.append((p, q))

def base_sup(B):
    """sup over atlas ratios of p0_inf(B,p,q), exact."""
    best = (Fr(0), 0, 0)
    for (p, q) in ATLAS:
        v = p0_inf(B, p, q)
        if v > best[0]:
            best = (v, p, q)
    return best

def enumerate_bases(k, max_elt=14, sample_cap=None):
    """bases B subset {0..max_elt}, 0 in B, |B|=k-2 (so E has k elements w/ 2 far)."""
    size = k - 2
    rest = list(range(1, max_elt + 1))
    out = []
    for combo in itertools.combinations(rest, size - 1):
        B = (0,) + combo
        # primitivity of B itself: gcd of all = 1 (0 ignored)
        g = 0
        for b in B:
            g = gcd(g, b)
        if g != 1:
            continue
        out.append(B)
        if sample_cap and len(out) >= sample_cap:
            break
    return out

def main():
    print("=" * 84)
    print("BASE DOMINATION AUDIT: is the even-AP base the WORST (max p0_inf) over all bounded bases?")
    print("=" * 84)
    for k, evenAP in [(8, (0, 2, 4, 6, 8, 10)), (9, (0, 2, 4, 6, 8, 10, 12))]:
        size = k - 2
        nbases = sum(1 for _ in itertools.combinations(range(1, 15), size - 1))
        print(f"\n--- k={k}: |B|={size}, bases B subset {{0..14}} with 0 in B (~{nbases} combos) ---")
        ap_sup = base_sup(evenAP)
        print(f"  even-AP base {evenAP}: sup p0_inf = {float(ap_sup[0]):.5f} ({ap_sup[0]}) at {ap_sup[1]}/{ap_sup[2]}")
        bases = enumerate_bases(k)
        global_best = (Fr(0), None, 0, 0)
        nbeats = 0
        beaters = []
        for B in bases:
            bs = base_sup(B)
            if bs[0] > global_best[0]:
                global_best = (bs[0], B, bs[1], bs[2])
            if bs[0] > ap_sup[0]:
                nbeats += 1
                if len(beaters) < 8:
                    beaters.append((B, bs))
        print(f"  enumerated {len(bases)} primitive bases.")
        print(f"  GLOBAL MAX p0_inf = {float(global_best[0]):.5f} ({global_best[0]}) at base {global_best[1]}, ratio {global_best[2]}/{global_best[3]}")
        print(f"  cap_{k} = {float(CAP[k]):.5f}; global-max margin = {float(CAP[k]-global_best[0]):.5f}")
        if nbeats == 0:
            print(f"  => even-AP IS the worst base (0 bases beat it). Single-base atlas SUFFICIENT.")
        else:
            print(f"  => {nbeats} bases BEAT the even-AP! Single-base atlas INSUFFICIENT. Examples:")
            for B, bs in beaters:
                print(f"       base {B}: sup p0_inf={float(bs[0]):.5f} at {bs[1]}/{bs[2]} (even-AP was {float(ap_sup[0]):.5f})")
        # Critical: does the GLOBAL max stay below cap?
        if global_best[0] >= CAP[k]:
            print(f"  *** CAP VIOLATION at global max base {global_best[1]} ***")
        else:
            print(f"  Global-max base still < cap_{k}: SAFE (margin {float(CAP[k]-global_best[0]):.5f})")

if __name__ == "__main__":
    main()
