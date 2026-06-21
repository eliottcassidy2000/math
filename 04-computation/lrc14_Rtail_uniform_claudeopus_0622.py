#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claude-opus 2026-06-22: is the GENERAL bounded-base R-tail UNIFORMLY BOUNDED over (base, gap)?

THM-564 closes the adjacent doublet {M,M+1} via M*(p0-Phi)=P(M)+R(M), R=M*(d2-d_inf)=O(1/M),
sup|R|~0.4-0.7 (consec base). HYP-2807: the genuine-wide max is a GENERALIZED doublet {M,M+g}
(any base, any gap). To complete the proof we need the R-tail UNIFORM over (base B, gap g):
   R_g(M) = M*(d2_g(M) - d_inf_g),   d2_g(M)=p0(B+{M,M+g})-p0(B+{M})-p0(B+{M+g})+p0(B).

This script computes sup_{M in [15,F]} |R_g(M)| for MANY bounded bases B (consec, even-AP,
top-cluster, random) of size k-2 and gaps g=1,2,3, and reports the GLOBAL sup. If it stays
bounded (~<2) UNIFORMLY, the R-tail closes with a uniform cutoff f0 = ceil(G/min-margin).
d_inf_g estimated from a late-M window (exact value = THM-564 Phi_frozen; estimate suffices for
the sup over [15,F]). Exact rationals.
"""
from __future__ import annotations
import sys, functools, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from lrc14_threadA_regime_dichotomy_kpswf8 import p0_fast, CAP


def base_family(k, rng, n_random=12):
    size = k - 2
    fams = {}
    fams["consec"] = tuple(range(size))
    fams["even-AP"] = tuple([0] + [2 * i for i in range(1, size)]) if 2 * (size - 1) <= 14 else None
    fams["top-cluster"] = tuple([0] + list(range(15 - (size - 1), 15)))
    # a couple of structured non-consec
    fams["even+oddbridge"] = tuple(sorted(set([0] + [2 * i for i in range(1, size - 1)] + [ (2*(size-2)) -1 ])))[:size] if size>=3 else None
    out = {}
    for name, B in fams.items():
        if B is None:
            continue
        B = tuple(sorted(set(B)))
        if len(B) == size and 0 in B and max(B) <= 14:
            out[name] = B
    cnt = 0
    while cnt < n_random:
        S = tuple(sorted(rng.sample(range(1, 15), size - 1)))
        B = (0,) + S
        if len(B) == size:
            out[f"rand{cnt}"] = B
            cnt += 1
    return out


def Rtail_sup(B, g, F_win=200, late=150):
    p0B = p0_fast(B)
    A = {}
    for f in range(15, F_win + g + 1):
        A[f] = p0_fast(tuple(sorted(B + (f,))))
    d2 = {}
    for M in range(15, F_win + 1):
        d2[M] = (p0_fast(tuple(sorted(B + (M, M + g)))) - A[M] - A[M + g] + p0B)
    d_inf = sum(d2[M] for M in range(late, F_win + 1)) / (F_win + 1 - late)
    supR = F(0)
    argM = None
    for M in range(15, F_win + 1):
        r = abs(M * (d2[M] - d_inf))
        if r > supR:
            supR, argM = r, M
    return supR, argM, d_inf


def main():
    rng = random.Random(99)
    print("=" * 78)
    print("GENERAL bounded-base R-tail: sup_M |R_g(M)| over (base, gap)  claude-opus 2026-06-22")
    print("=" * 78)
    global_sup = F(0)
    global_at = None
    for k in (9, 10, 11):
        fams = base_family(k, rng)
        print(f"\nk={k}  (base size {k-2}):")
        for name, B in fams.items():
            row = []
            for g in (1, 2, 3):
                supR, argM, d_inf = Rtail_sup(B, g)
                row.append(f"g={g}:{float(supR):.3f}@{argM}")
                if supR > global_sup:
                    global_sup, global_at = supR, (k, name, B, g, argM)
            print(f"   {name:14s} {B}:  " + "  ".join(row))
    print("\n" + "=" * 78)
    print(f"GLOBAL sup_M |R_g(M)| over all tested (base, gap, k) = {float(global_sup):.4f}")
    print(f"   attained at k={global_at[0]} base={global_at[1]}={global_at[2]} g={global_at[3]} M={global_at[4]}")
    print(f"=> general bounded-base R-tail UNIFORMLY BOUNDED by ~{float(global_sup):.2f}"
          f" (if this holds over ALL bases, cutoff f0 = ceil(G/margin) is uniform).")


if __name__ == "__main__":
    main()
