#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TASK 2: test the organizing principle -- is the R-ODD (epsilon=-1) eigenspace
of the cap the OBSTRUCTION?  (mac-mini-2026-06-29-S8)

PRINCIPLE (user): the LRC 'signed certificate' = the R-odd (eps=-1) eigenspace of
the cap; the R-symmetric half data on (0,1/2) is the SOS/Brouwer-provable piece;
the even/Brouwer (+) odd/Borsuk-Ulam split IS the +-1 eigenspace split of the
reversal R.  FALSIFIABLE: compute the R-odd component of the cap for the binding
LRC config and check it carries the obstruction (the deviation / the hard part),
while the R-even component is the SOS-provable bulk.

We use kps's pairwise sector co-emptiness matrix M (6x6 on inner sectors 1..6),
M[a][b] = meas{t: inner sectors a AND b both EMPTY}, and the inclusion-exclusion
   cap = 1 - S1 + S2 - S3 + S4 - S5 + S6,  S_j = sum_{|T|=j} meas(T all empty).
R is the time reflection a->6-a on sectors (the complement t->-t).  On inner
sectors {1..6} it acts as (1 5)(2 4) fixing 3,6 (since sector 0's partner is the
danger sector).  R-even eigenspace dim 4, R-odd dim 2 = span(e1-e5, e2-e4).
We decompose M = M_even (+) M_odd and the cap's S2 = (1/2)(1^T M 1 - tr M) into
R-even/odd parts, and the FULL cap via the R-orbit decomposition of subsets.
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)
from fractions import Fraction as F


def sector_empty_measure(E, subset):
    """meas{ t in [0,1): for every a in subset, NO nonzero e in E lands in inner sector a }.
    sector of e at t = floor(7*frac(e t)); inner sectors 1..6; speed 0 ignored."""
    nz = [int(x) for x in E if x]
    if not nz:
        return F(1) if not subset else F(0)
    from functools import reduce
    from math import gcd
    l = reduce(lambda a, b: a // gcd(a, b) * b, nz)
    d = 7 * l
    den2 = 2 * l
    bps = {0, d}
    for e in nz:
        step = l // e
        x = 0
        for _ in range(7 * e + 1):
            bps.add(x); x += step
    bps = sorted(bps)
    target = set(subset)
    num = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        occ = set((e * midnum // den2) % 7 for e in nz)
        if not (target & occ):     # none of the subset's sectors are occupied
            num += hi - lo
    return F(num, d)


def coemptiness_matrix(E):
    M = [[F(0)] * 6 for _ in range(6)]
    for a in range(1, 7):
        for b in range(a, 7):
            m = sector_empty_measure(E, {a, b}) if a != b else sector_empty_measure(E, {a})
            M[a - 1][b - 1] = m; M[b - 1][a - 1] = m
    return M


def cap_incl_excl(E):
    """cap = sum_{T subset {1..6}} (-1)^|T| meas(T all empty), split by R-orbit parity.
    Returns (cap, cap_from_R_fixed_subsets, cap_from_R_paired_subsets)."""
    # R on inner sectors: a -> 6-a, i.e. 1<->5, 2<->4, 3->3, 6->6 (partner 0 gone)
    Rmap = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1, 6: 6}
    cap = F(0); cap_fixed = F(0); cap_paired = F(0)
    for r in range(0, 7):
        for T in itertools.combinations(range(1, 7), r):
            Tset = frozenset(T)
            m = sector_empty_measure(E, Tset)
            term = ((-1) ** r) * m
            cap += term
            RT = frozenset(Rmap[a] for a in T)
            if RT == Tset:
                cap_fixed += term
            else:
                cap_paired += term
    return cap, cap_fixed, cap_paired


def R_decompose(M):
    """R = (1 5)(2 4) fix 3,6 on indices 0..5 (sectors 1..6).
    Even projector P+ = (I+R)/2, odd P- = (I-R)/2.  Return M_even=P+ M P+, M_odd=P- M P-."""
    import numpy as np
    R = np.zeros((6, 6))
    perm = {0: 4, 1: 3, 2: 2, 3: 1, 4: 0, 5: 5}  # 0-indexed: s1<->s5,s2<->s4,s3,s6 fix
    for i in range(6):
        R[perm[i]][i] = 1
    Mn = np.array([[float(x) for x in row] for row in M])
    Pp = (np.eye(6) + R) / 2
    Pm = (np.eye(6) - R) / 2
    Meven = Pp @ Mn @ Pp
    Modd = Pm @ Mn @ Pm
    return Mn, Meven, Modd, Pp, Pm


def main():
    import numpy as np
    print("=" * 84)
    print("TASK 2: is the R-ODD eigenspace of the cap the OBSTRUCTION? (mac-mini-S8)")
    print("=" * 84)

    configs = {
        "consec_8 (k=8 base)": (0, 1, 2, 3, 4, 5, 6, 7),
        "minimizer {1,5,7,8,9}": (0, 1, 5, 7, 8, 9),
        "consec_9 (k=9 base)": (0, 1, 2, 3, 4, 5, 6, 7, 8),
        "j=4 min {1,11,12,13}": (0, 1, 11, 12, 13),
    }
    for name, E in configs.items():
        M = coemptiness_matrix(E)
        Mn, Meven, Modd, Pp, Pm = R_decompose(M)
        S2 = sum(M[a][b] for a in range(6) for b in range(a + 1, 6))
        trM = sum(M[a][a] for a in range(6))
        cap, cap_fixed, cap_paired = cap_incl_excl(E)
        # eigenvalues of full, even, odd blocks
        evfull = sorted(np.linalg.eigvalsh(Mn), reverse=True)
        # odd block lives on span(e1-e5,e2-e4): build 2x2
        v1 = np.array([1, 0, 0, 0, -1, 0]) / np.sqrt(2)   # s1 - s5
        v2 = np.array([0, 1, 0, -1, 0, 0]) / np.sqrt(2)   # s2 - s4
        Vodd = np.column_stack([v1, v2])
        Modd2 = Vodd.T @ Mn @ Vodd
        ev_odd = sorted(np.linalg.eigvalsh(Modd2), reverse=True)
        tr_odd = float(np.trace(Modd2))
        tr_even = float(trM) - tr_odd
        print(f"\n--- {name}: E={E} ---")
        print(f"  cap (incl-excl) = {float(cap):.6f}  [R-fixed subsets: {float(cap_fixed):.6f}, "
              f"R-paired: {float(cap_paired):.6f}]")
        print(f"  S2={float(S2):.5f}  trM={float(trM):.5f}  "
              f"=> tr(M_even)={tr_even:.5f}, tr(M_odd)={tr_odd:.5f}")
        print(f"  R-odd contribution to S2 = -tr(M_odd)/2 = {-tr_odd/2:.5f}  "
              f"(S2 even part = {(float(S2)+tr_odd/2):.5f})")
        print(f"  full M eigenvalues: {[f'{x:.4f}' for x in evfull]}")
        print(f"  R-ODD 2x2 block eigenvalues: {[f'{x:.4f}' for x in ev_odd]}  (Perron={evfull[0]:.4f} is R-{'EVEN' if abs(evfull[0]-ev_odd[0])>1e-6 else 'ODD'})")

    print("\n" + "=" * 84)
    print("Test: if the Perron (bulk) mode is R-EVEN and the binding deviation /")
    print("hardness sits in the R-ODD block, the user's principle is CONFIRMED:")
    print("R-even = SOS-provable bulk, R-odd = the signed obstruction.")
    print("=" * 84)


if __name__ == "__main__":
    main()
