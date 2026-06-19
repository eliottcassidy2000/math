#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) — does the gcd-3 stratum carry the DANGEROUS near-AP corrections,
and the unit stratum the harmless ones?  kind-pasteur-2026-06-19-S13.

Method (per task ANGLE): the correction L_y(E) - L_y^inf is driven by orbit
COLLISIONS, i.e. by additive relations among the speeds e_i.  We probe which
collisions matter by:
  (A) For the AP consec, perturb ONE element e -> e+27 (a pure +C=+27 shift =
      same residue mod 27, but breaks the small integer relations while
      PRESERVING the residue/stratum).  Measure the L_y drop.
  (B) For the AP consec, replace one element by a different residue in the SAME
      stratum vs a DIFFERENT stratum, and measure L_y.
  (C) Directly: classify each element of consec by stratum and measure the L_y
      sensitivity dL = L_y(consec) - L_y(consec minus e plus e') when e' is chosen
      to land in unit vs gcd-3 vs gcd-9 stratum.

The diagnostic: a "dangerous" perturbation keeps L_y HIGH (near cap); a
"harmless" one drops it.  If gcd-3 perturbations keep L_y high (the AP-like
collisions survive) and unit perturbations drop it (collisions break), then the
gcd-3 stratum is where the near-AP danger lives — matching V* sitting in the
gcd-3 blind spot (HYP-2083).
"""
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, "04-computation")
from lrc14_summand_strata_relation_lattice_kps_S13 import L_y, stratum, C, g_poly

def report(k, base, perturbs, title):
    Lbase, _ = L_y(base, k)
    print(f"\n--- {title}  (k={k}) ---")
    print(f"  base E = {base}   L_y = {float(Lbase):.5f}")
    for label, E in perturbs:
        E = sorted(set(E))
        if len(E) != k or 0 not in E:
            print(f"    [{label}] SKIP |E|={len(E)} E={E}")
            continue
        L, _ = L_y(E, k)
        dL = float(L - Lbase)
        print(f"    [{label:28s}] E={E}  L_y={float(L):.5f}  dL={dL:+.5f}")

if __name__ == "__main__":
    print("="*78)
    print("STRATUM CORRECTION SENSITIVITY — which stratum carries the danger?")
    print("="*78)
    print(f"C={C}=3^3. unit=gcd1, gcd3, gcd9. AP consec elements 0..k-1 sit in low residues.")

    # ===== k=8 =====
    k = 8
    base = list(range(8))  # 0..7
    # (A) +27 shift of one element: preserves residue & stratum, breaks integer relations.
    A = []
    for e in range(1, 8):
        E2 = [x for x in base if x != e] + [e + 27]
        A.append((f"shift {e}->+27 (resid {e%27}, {stratum(e)})", E2))
    report(k, base, A, "(A) +27 shift (residue/stratum PRESERVED, integer relations broken)")

    # (B) replace top element 7 by a higher value in chosen stratum.
    #     Candidate replacement residues mod 27: unit {10,11}, gcd3 {12,15}, gcd9 {9,18}.
    B = []
    for cand in [10, 11, 13, 9, 18, 12, 15, 21, 24]:
        E2 = [x for x in base if x != 7] + [cand]
        B.append((f"7->{cand} ({stratum(cand)})", E2))
    report(k, base, B, "(B) replace top elt 7 -> candidate residue (vary stratum)")

    # (C) replace a MIDDLE element by same-residue+27 vs cross-stratum jump
    # ===== k=9 (the TIGHTEST: margin to cap only 0.0014) =====
    k = 9
    base9 = list(range(9))
    A9 = []
    for e in [3, 6, 1, 2, 9 if 9 in base9 else 8]:
        if e not in base9: continue
        E2 = [x for x in base9 if x != e] + [e + 27]
        A9.append((f"shift {e}->+27 ({stratum(e)})", E2))
    report(k, base9, A9, "(A) k=9 +27 shift on gcd-3 elements {3,6} vs unit elements {1,2}")

    # ===== Direct: per-stratum L_y of consec with each element pushed far (e+27*3) =====
    print("\n" + "="*78)
    print("DIRECT STRATUM TEST: push element e by +81 (=3*27, preserves residue),")
    print("compare L_y drop for gcd-3 elements {3,6} vs unit elements {1,2,4,5,7}.")
    print("="*78)
    for k in [8, 9]:
        base = list(range(k))
        Lbase, _ = L_y(base, k)
        print(f"\n k={k}  consec L_y={float(Lbase):.5f}")
        rows = []
        for e in range(1, k):
            E2 = sorted(set([x for x in base if x != e] + [e + 81]))
            L, _ = L_y(E2, k)
            rows.append((e, stratum(e), float(Lbase - L)))
        # aggregate by stratum
        from collections import defaultdict
        agg = defaultdict(list)
        for e, st, drop in rows:
            agg[st].append(drop)
            print(f"   push elt {e} (+81, resid {e%27}, {st:5s}) -> L_y drop = {drop:+.5f}")
        print(f"   --- mean L_y drop by stratum:")
        for st, ds in sorted(agg.items()):
            print(f"       {st:6s}: mean drop {sum(ds)/len(ds):+.5f} over {len(ds)} elts")
