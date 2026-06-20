#!/usr/bin/env python3
"""
half_tiling_maxH_kpswf.py
================================================================
THREAD 1 -- HYP-2688: the GLOBAL H-maximizer is attained by a grid-symmetric
(phi-self-converse) tiling, so a max-H search only needs the 2^half(n)
fundamental domain instead of the full 2^m tiling cube.

We REUSE the canonical machinery from half_tiling_framework_kps.py
(tiles, tiling_to_tournament, count_hamiltonian_paths, rho, involution_structure).

Tasks:
 (a) n=8: exhaust all 2^12=4096 grid-symmetric tilings; report max grid-sym H;
     compare to global max H(8)=661.
 (b) n=9: exhaust all 2^16=65536 grid-symmetric tilings; report max grid-sym H;
     compare to global max H(9)=3357 (canon THM-329 / A038375).
 (c) n=5,6,7: full 2^m exhaustion; find ALL maximizing tilings; report how many
     are grid-sym vs not (does maximizer SET merely CONTAIN a grid-sym point, or
     is it ENTIRELY grid-sym?).
 (d) Proof-attempt logic: confirm maximizer set is rho-invariant (H=H^op), and
     report the structural status.
 (e) Speedup table 2^half vs 2^m for n=5..14.

KEY SUBTLETY (tested below): "grid-symmetric tiling" means the tiling is fixed by
rho for the FIXED canonical base path P0 = n->...->1. HYP-2688 as stated is about
THIS specific representation. We also test the WEAKER claim (d): whether the
global maximizer tournament admits SOME phi-self-converse tiling representative.

Output is exhaustive/verified where the cube is enumerable.
"""

import sys, os, itertools
from math import comb, floor
from collections import Counter

# import the canonical framework helpers
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from half_tiling_framework_kps import (
    tiles, rho, tiling_to_tournament, count_hamiltonian_paths,
    involution_structure,
)

# Canonical global max H (A038375 / THM-329). a(n) for n vertices.
# a(1..11) = 1,1,3,5,15,45,189,661,3357,15745,95095
GLOBAL_MAX_H = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357,
                10:15745, 11:95095}

def half_formula(n):
    return floor((n-1)**2 / 4) if n >= 1 else 0


# ----------------------------------------------------------------------
# Grid-symmetric tiling enumeration: assign a free bit to each rho-orbit
# representative, propagate the SAME bit to the orbit's other cell.
# ----------------------------------------------------------------------

def grid_sym_max_H(n, return_argmax=False):
    """Exhaust all 2^half(n) grid-symmetric tilings. Return (max_H, count_tilings,
    spectrum Counter, list_of_argmax_bitvectors_optional)."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    # orbit reps: each fixed cell is its own orbit; each pair contributes one rep.
    # We enumerate one bit per orbit and build the full tile-bit dict by copying.
    orbit_cells = []  # list of lists of tile-positions sharing one free bit
    for f in fixed:
        orbit_cells.append([f])
    for (p, q) in pairs:
        orbit_cells.append([p, q])
    assert len(orbit_cells) == half, (len(orbit_cells), half)

    spec = Counter()
    best = -1
    argmax = []
    for bits in itertools.product([0, 1], repeat=half):
        t = {}
        for i, cells in enumerate(orbit_cells):
            for c in cells:
                t[c] = bits[i]
        T = tiling_to_tournament(t, tl, n)
        H = count_hamiltonian_paths(T, n)
        spec[H] += 1
        if H > best:
            best = H
            argmax = [bits] if return_argmax else []
        elif H == best and return_argmax:
            argmax.append(bits)
    return best, 2**half, spec, argmax


# ----------------------------------------------------------------------
# Full-cube exhaustion: every tiling, classify grid-sym vs not, find all maxima.
# ----------------------------------------------------------------------

def full_cube_analysis(n):
    """Exhaust all 2^m tilings. Return dict with:
       global_max_H (over all tilings = true global max for n vertices),
       n_max_total, n_max_gridsym, n_max_nongridsym, full_spectrum."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rho_map = {p: rho(*p, n) for p in tl}
    spec = Counter()
    best = -1
    n_max_total = 0
    n_max_gridsym = 0
    n_max_nongridsym = 0
    # We do two passes to keep memory tiny: pass1 find max, pass2 count.
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        H = count_hamiltonian_paths(T, n)
        spec[H] += 1
        if H > best:
            best = H
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        H = count_hamiltonian_paths(T, n)
        if H == best:
            n_max_total += 1
            is_gs = all(t[p] == t[rho_map[p]] for p in tl)
            if is_gs:
                n_max_gridsym += 1
            else:
                n_max_nongridsym += 1
    return {
        'global_max_H': best, 'm': m, 'half': half,
        'n_max_total': n_max_total,
        'n_max_gridsym': n_max_gridsym,
        'n_max_nongridsym': n_max_nongridsym,
        'spectrum': spec,
    }


# ----------------------------------------------------------------------
# (d) weaker claim: does the global-max TOURNAMENT admit SOME phi-self-converse
# tiling rep? Equivalent (for fixed base path) to: is the global maximum H value
# also achieved by at least one grid-symmetric tiling? We answer this directly by
# comparing grid_sym_max_H to the true global max value.
# ----------------------------------------------------------------------


def main():
    print("=" * 72)
    print("HALF-TILING MAX-H SEARCH  --  HYP-2688  --  kind-pasteur-S-wf")
    print("=" * 72)
    print("\nGlobal max H (A038375 / THM-329): a(n) by #vertices n:")
    print("  ", {k: GLOBAL_MAX_H[k] for k in sorted(GLOBAL_MAX_H)})

    # ---------- (c) full-cube exhaustion for n=5,6,7 ----------
    print("\n" + "=" * 72)
    print("[c] FULL-CUBE exhaustion (n=5,6,7): ALL maximizing tilings,")
    print("    grid-sym vs non-grid-sym among the maximizer SET.")
    print("=" * 72)
    print(f"{'n':>3} {'m':>3} {'2^m':>8} {'globalH':>8} {'#max':>6} "
          f"{'#max_GS':>8} {'#max_nonGS':>11} {'GS attains?':>12}")
    fullres = {}
    for n in [5, 6, 7]:
        r = full_cube_analysis(n)
        fullres[n] = r
        gs_attains = (r['n_max_gridsym'] > 0)
        match_canon = (r['global_max_H'] == GLOBAL_MAX_H[n])
        print(f"{n:>3} {r['m']:>3} {2**r['m']:>8} {r['global_max_H']:>8} "
              f"{r['n_max_total']:>6} {r['n_max_gridsym']:>8} "
              f"{r['n_max_nongridsym']:>11} {str(gs_attains):>12} "
              f"{'(canon OK)' if match_canon else '(canon MISMATCH!)'}")
    # interpretation line
    for n in [5, 6, 7]:
        r = fullres[n]
        if r['n_max_nongridsym'] == 0 and r['n_max_gridsym'] > 0:
            verdict = "maximizer set ENTIRELY grid-sym"
        elif r['n_max_gridsym'] > 0:
            verdict = "maximizer set CONTAINS grid-sym (mixed) -- speedup OK"
        else:
            verdict = "maximizer set has NO grid-sym point -- HYP-2688 FAILS"
        print(f"    n={n}: {verdict}")

    # ---------- (a) n=8 grid-sym exhaustion ----------
    print("\n" + "=" * 72)
    print("[a] n=8 : exhaust 2^half(8)=2^12=4096 grid-symmetric tilings.")
    print("=" * 72)
    best8, ntil8, spec8, _ = grid_sym_max_H(8)
    print(f"    grid-sym tilings enumerated = {ntil8}  (= 2^{half_formula(8)})")
    print(f"    MAX grid-sym H(8)           = {best8}")
    print(f"    global max H(8)             = {GLOBAL_MAX_H[8]}")
    print(f"    grid-sym ATTAINS global max = {best8 == GLOBAL_MAX_H[8]}")
    top8 = dict(sorted(spec8.items(), reverse=True)[:8])
    print(f"    top grid-sym H values (H:count) = {top8}")

    # ---------- (b) n=9 grid-sym exhaustion ----------
    print("\n" + "=" * 72)
    print("[b] n=9 : exhaust 2^half(9)=2^16=65536 grid-symmetric tilings.")
    print("=" * 72)
    best9, ntil9, spec9, _ = grid_sym_max_H(9)
    print(f"    grid-sym tilings enumerated = {ntil9}  (= 2^{half_formula(9)})")
    print(f"    MAX grid-sym H(9)           = {best9}")
    print(f"    global max H(9) (canon)     = {GLOBAL_MAX_H[9]}")
    print(f"    grid-sym ATTAINS global max = {best9 == GLOBAL_MAX_H[9]}")
    top9 = dict(sorted(spec9.items(), reverse=True)[:8])
    print(f"    top grid-sym H values (H:count) = {top9}")

    # ---------- (d) proof-attempt structural checks ----------
    print("\n" + "=" * 72)
    print("[d] PROOF-ATTEMPT structural facts")
    print("=" * 72)
    print("  Fact 1 (H=H^op): the maximizer SET is closed under rho (converse-relabel),")
    print("    because H(T)=H(T^op) always (reverse a Ham path -> Ham path of T^op),")
    print("    and rho realizes the converse-relabel as a coordinate involution (THM-280).")
    print("    => the set of maximizing tilings is a UNION of rho-orbits.")
    print("  Fact 2 (needed): HYP-2688 needs >=1 rho-FIXED maximizer (a grid-sym tiling).")
    print("    A rho-orbit of size 2 (a non-grid-sym tiling + its mirror) does NOT by")
    print("    itself give a fixed point. So Fact 1 alone is INSUFFICIENT.")
    print("    Verdict per n (does some grid-sym tiling attain the global max?):")
    summary = {}
    for n in [5, 6, 7]:
        r = fullres[n]
        summary[n] = (r['n_max_gridsym'] > 0)
    summary[8] = (best8 == GLOBAL_MAX_H[8])
    summary[9] = (best9 == GLOBAL_MAX_H[9])
    for n in sorted(summary):
        print(f"      n={n}: grid-sym attains global max = {summary[n]}")

    # ---------- (e) speedup table ----------
    print("\n" + "=" * 72)
    print("[e] SPEEDUP TABLE  2^half(n)  vs  2^m,  m=C(n-1,2)")
    print("=" * 72)
    print(f"{'n':>3} {'m=C(n-1,2)':>11} {'half(n)':>8} {'2^m (full)':>16} "
          f"{'2^half (fund)':>16} {'speedup=2^(m-half)':>20}")
    for n in range(5, 15):
        m = comb(n-1, 2)
        h = half_formula(n)
        print(f"{n:>3} {m:>11} {h:>8} {2**m:>16} {2**h:>16} {2**(m-h):>20}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
