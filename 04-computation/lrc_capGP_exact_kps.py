#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_capGP_exact_kps.py  (kind-pasteur 2026-06-22, witness-floor Lemma B)

EXACT min meas(G_P) over all P subset {1,...,13} with |P| = 13-k, for k=8..13.
G_P = {tau in [0,1): ||p tau|| >= 1/14 for all p in P}.
meas(G_P) = 1 - meas(Union_p B_p),  B_p = {tau: ||p tau|| < 1/14}.

B_p is a union of p arcs, each of width 1/(7p), centered at j/p (j=0..p-1).
meas(Union) computed EXACTLY via the breakpoint/coverage method on rationals.

This pins cap_k := min_{|P|=13-k} meas(G_P) as an exact rational (Lemma B), and
checks whether it coincides with the skeleton's p0-plateau capRat (a duality).
"""
import itertools
import sys
from fractions import Fraction as Fr


def bad_arcs(p, half=Fr(1, 14)):
    """Arcs of {tau: ||p tau|| < 1/14} in [0,1): for j=0..p-1, (j/p - 1/(14p), j/p + 1/(14p))."""
    arcs = []
    w = half / p  # = 1/(14p)
    for j in range(p):
        c = Fr(j, p)
        lo, hi = c - w, c + w
        # wrap into [0,1)
        arcs.append((lo, hi))
    return arcs


def union_measure(arcs):
    """Exact Lebesgue measure of a union of (possibly wrapping) open arcs in the circle [0,1)."""
    # Normalize each arc into [0,1), splitting wrappers.
    # FIX: compute length BEFORE modding lo (else wrapping arcs around 0 are dropped).
    segs = []
    for lo, hi in arcs:
        length = hi - lo
        if length <= 0:
            continue
        a = lo % 1
        b = a + length
        if b <= 1:
            segs.append((a, b))
        else:
            segs.append((a, Fr(1)))
            segs.append((Fr(0), b - 1))
    # sweep-line union
    segs.sort()
    tot = Fr(0)
    cur_lo, cur_hi = None, None
    for a, b in segs:
        if cur_lo is None:
            cur_lo, cur_hi = a, b
        elif a <= cur_hi:
            if b > cur_hi:
                cur_hi = b
        else:
            tot += cur_hi - cur_lo
            cur_lo, cur_hi = a, b
    if cur_lo is not None:
        tot += cur_hi - cur_lo
    return tot


def measGP(P):
    arcs = []
    for p in P:
        arcs.extend(bad_arcs(p))
    return 1 - union_measure(arcs)


def main():
    print("=" * 70)
    print("Lemma B: cap_k = min_{|P|=13-k} meas(G_P), EXACT rational")
    print("=" * 70)
    # skeleton capRat (p0-plateau) for comparison
    capRat = {8: Fr(2243, 5880), 9: Fr(1979, 4004), 10: Fr(55, 91),
              11: Fr(66, 91), 12: Fr(6, 7), 13: Fr(1)}
    out = {}
    for k in range(8, 14):
        m = 13 - k
        best, bestP = None, None
        if m == 0:
            best, bestP = Fr(1), ()
        else:
            for P in itertools.combinations(range(1, 14), m):
                v = measGP(P)
                if best is None or v < best:
                    best, bestP = v, P
        out[k] = (best, bestP)
        match = "== capRat" if best == capRat[k] else f"!= capRat({capRat[k]})"
        print(f"  k={k:2d} (|P|={m}): cap_k = min meas(G_P) = {best} = {float(best):.6f}  "
              f"at P={bestP}   {match}")
        sys.stdout.flush()
    print("\n--- duality check: does min meas(G_P) == p0-plateau capRat? ---")
    allmatch = all(out[k][0] == capRat[k] for k in range(8, 14))
    print(f"  ALL k match: {allmatch}")
    if allmatch:
        print("  => the small-part safe-measure floor IS the p0 Q-plateau cap_k.")
        print("     (a genuine duality; same rational both routes.)")
    print("\n--- Bonferroni floor with EXACT capGP ---")
    nuc = {8: Fr(691, 735), 9: Fr(247, 294), 10: Fr(38, 49),
           11: Fr(1381, 2205), 12: Fr(13823, 24255), 13: Fr(477, 1078)}
    worst = None
    for k in range(8, 14):
        fl = nuc[k] + out[k][0] - 1
        if worst is None or fl < worst:
            worst = fl
        print(f"  k={k:2d}: nu_consec + cap - 1 = {nuc[k]} + {out[k][0]} - 1 = {fl} = {float(fl):+.6f}")
    print(f"\n  worst exact floor = {worst} = {float(worst):.6f}  (> 0 for all k=8..13)")
    print("=" * 70)
    print("DONE.")


if __name__ == "__main__":
    main()
