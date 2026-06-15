#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TWO DIMENSIONS of the OCF non-spectral content (monad-explorer-2026-06-15-S? deep)

The growth law (THM-505) reports dim_nonspec(H)(n) = #partitions(odd>=3,<=n) - 3
= q(n)-3 = 3,5,7 at n=8,9,10, by RANK of the individual OCF carriers
{c7, D33, D35, ...}.  But H = I(Omega,2) = 1 + 2*a1 + 4*a2 + 8*a3 depends ONLY on
the LEVEL-SUMS a_j (the number of disjoint odd-cycle j-packings).  Within a
cospectral class, dH = 2 da1 + 4 da2 + 8 da3 is an IDENTITY.  So the honest
"number of independent non-spectral invariants needed to determine H" is

    dim_func(H)(n) = rank{ da1, da2, da3, ... }   (level-sums)

NOT rank of the individual carriers.  These differ whenever a level contains >1
independent carrier (e.g. a2 = D33+D35 at n>=8).  This script measures BOTH and
the relation H in span(...) for each, to settle which quantity the "growth law"
is really counting.

Reuses the exact carrier machinery of ocf_nonspectral_n10_monad.py.
"""
import sys
from collections import defaultdict
from fractions import Fraction
import numpy as np

from ocf_nonspectral_n10_monad import (
    sample_sigs, code_to_adj, analyze, matrix_rank_Q, within_class_deltas,
    partitions_odd_ge3_upto,
)


def run(n, samples, p2cap, seed=12345):
    print(f"\n{'='*64}\n n = {n}   samples={samples}  cap={p2cap}\n{'='*64}", flush=True)
    # phase 1: bucket codes by trace signature
    sig_codes = defaultdict(list)
    for sig, code in sample_sigs(n, samples, batch=20000, seed=seed):
        sig_codes[sig].append(code)
    colliding = {s: cs for s, cs in sig_codes.items() if len(cs) >= 2}
    print(f"  cospectral classes with >=2 members: {len(colliding)}", flush=True)

    # phase 2: full analysis on colliding members (capped per class)
    classes = {}
    nclass = 0
    for sig, codes in colliding.items():
        recs = []
        for code in codes[:p2cap]:
            recs.append(analyze(code_to_adj(code, n), n, sig))
        classes[sig] = recs
        nclass += 1
    # add level-sum fields
    for sig, recs in classes.items():
        for r in recs:
            # a1 non-spectral part within class = c7 + c9 + c11 (c3,c5 spectral)
            r["a1ns"] = r.get("c7", 0) + r.get("c9", 0) + r.get("c11", 0)
            r["a2sum"] = (r["D33"] + r["D35"] + r["D37"] + r["D55"]
                          + r["D57"] + r["D39"])
            r["a3sum"] = r["T333"] + r["T335"]
    nmembers = sum(len(v) for v in classes.values())
    ok = all(r["ok_ocf"] for v in classes.values() for r in v)
    print(f"  analyzed members: {nmembers}   OCF holds on all: {ok}", flush=True)

    # OCF individual non-spectral carriers at this n (same set as growth-law script)
    ocf = ["c7"]
    if n >= 9:  ocf.append("c9")
    if n >= 11: ocf.append("c11")
    if n >= 6:  ocf.append("D33")
    if n >= 8:  ocf.append("D35")
    if n >= 10: ocf += ["D37", "D55"]
    if n >= 12: ocf += ["D57", "D39"]
    if n >= 9:  ocf.append("T333")
    if n >= 11: ocf.append("T335")

    # level-sum carriers present at this n
    levels = ["a1ns"]            # level 1, non-spectral part
    levels.append("a2sum")       # level 2
    if n >= 9: levels.append("a3sum")   # level 3

    drows_c = within_class_deltas(classes, ocf)
    drows_l = within_class_deltas(classes, levels)

    rk_c  = matrix_rank_Q(drows_c, ocf)
    rk_cH = matrix_rank_Q(drows_c, ocf + ["H"])
    rk_l  = matrix_rank_Q(drows_l, levels)
    rk_lH = matrix_rank_Q(drows_l, levels + ["H"])

    _, qn = partitions_odd_ge3_upto(n)
    growth = qn - 3

    print(f"\n  --- CARRIER-SPACE basis {ocf} ---")
    print(f"    rank(carriers)        = {rk_c}   (growth law predicts {growth})")
    print(f"    rank(carriers + H)    = {rk_cH}  "
          f"({'H in carrier span' if rk_cH==rk_c else 'H NEEDS more'})")

    print(f"\n  --- LEVEL-SUM basis {levels} ---")
    print(f"    rank(level-sums)      = {rk_l}")
    print(f"    rank(level-sums + H)  = {rk_lH}  "
          f"({'H in level-sum span' if rk_lH==rk_l else 'H NEEDS more'})")

    # are the individual carriers within a level actually independent?
    if n >= 8:
        for pair in (["D33", "D35"], ["D33", "D35", "D37", "D55"]):
            pair = [c for c in pair if c in ocf]
            if len(pair) >= 2:
                rp = matrix_rank_Q(within_class_deltas(classes, pair), pair)
                print(f"    rank within level-2 {pair} = {rp} "
                      f"({'all independent' if rp==len(pair) else 'dependent'})")

    print(f"\n  >>> CARRIER-SPACE dim = {rk_c}  vs  H-FUNCTIONAL dim = {rk_l}")
    if rk_l < rk_c:
        print(f"  >>> MISMATCH: H depends on only {rk_l} level-sums, but the growth "
              f"law counts {rk_c} individual carriers.")
    return rk_c, rk_l, growth


if __name__ == "__main__":
    # n: (samples, cap)
    plan = {8: (200000, 200), 9: (300000, 60), 10: (1500000, 12)}
    if len(sys.argv) > 1:
        ns = [int(x) for x in sys.argv[1:]]
    else:
        ns = [8, 9, 10]
    summary = []
    for n in ns:
        s, c = plan.get(n, (200000, 50))
        summary.append((n, *run(n, s, c)))
    print(f"\n{'='*64}\n SUMMARY\n{'='*64}")
    print(f"  {'n':>3} {'carrier-dim':>12} {'H-func-dim':>11} {'growth q(n)-3':>14}")
    for n, rc, rl, g in summary:
        print(f"  {n:>3} {rc:>12} {rl:>11} {g:>14}")
