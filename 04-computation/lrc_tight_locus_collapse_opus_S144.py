"""
lrc_tight_locus_collapse_opus_S144.py   (opus-2026-07-07-S144, HYP-5137 part 4)

THE TIGHT-LOCUS COLLAPSE CONJECTURE -- the surviving form of the ladder question after
the mu = M refutation (93 separations at |S| >= 4, all at LOOSE sets, M ~ 0.2+):

  (TLC)  mu(S) = M(S) on the NEAR-TIGHT covering sets -- the configurations that
         actually matter for LRC(14).  If TLC holds on the tight locus, the fractional
         (Motzkin/LP) relaxation is faithful exactly where LRC(14) binds, and
         MOTZKIN-14-on-the-locus recovers the constant.

TESTS (exact, max-cycle-mean engine; window graph w = max(S) so bounded-max families only):
 (1) The bounded repo extremal 13-families: GW {1..11,13,24} (M = 1/14 EXACTLY TIGHT),
     prim-sat 2{1..12}+{13} (M = 1/13), parity record {2,4,..,22,11,13} (M = 1/12).
     Plus the k <= 8 tight prefixes AP_k (M = 1/(k+1)) as controls.
 (2) NEAR-TIGHT PERTURBATIONS of GW: single-element moves (the adversary class that
     found every past counterexample) -- does collapse survive the neighborhood?
 (3) THE SEPARATION BOUNDARY: for the 26 |S|=4 exceptions, how far above the LRC bar
     do they sit?  min M over exceptions vs 1/5 (=1/(k+1) at k=4).  Is there any
     exception with M <= 1/(k+1)+eps (i.e. near-tight)?  [data check, no new compute]
 (4) |S| = 3 extension: max <= 18 (states cap) -- does the |S|=3 rigidity persist?
"""
from fractions import Fraction as F
from math import gcd
import sys, time, itertools
import numpy as np

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_graph_interpretation_ladder_opus_S141 import M_exact, witness
from lrc_mu_eq_M_maxcycle_opus_S144 import (build_window_graph, positive_cycle_exists,
                                            extract_beating_cycle, certify_witness_indep)

def test_set(S, label, verbose=True):
    S = sorted(S)
    M = M_exact(S)
    states, t0, t1 = build_window_graph(S)
    if len(states) > 400_000:
        print(f"    {label:44s} SKIP ({len(states)} states)")
        return None
    wa, wq, wv = witness(S)
    wok = certify_witness_indep(S, wa, wq, wv)
    pos, iters = positive_cycle_exists(t0, t1, M.numerator, M.denominator)
    if pos:
        dens, L, B = extract_beating_cycle(states, t0, t1, M.numerator, M.denominator, max(S))
        if verbose:
            print(f"    {label:44s} M = {str(M):>7}  *** mu > M: periodic density {dens} ***")
        return ("SEP", M, dens)
    else:
        if verbose:
            print(f"    {label:44s} M = {str(M):>7}  mu = M exactly "
                  f"({len(states)} states, wit-cert {wok})")
        return ("EQ", M, M)

def main():
    t00 = time.time()
    print("=" * 100)
    print("(1) TIGHT-LOCUS COLLAPSE: the bounded repo extremal 13-families + AP controls")
    print("=" * 100)
    GW = list(range(1, 12)) + [13, 24]
    primsat = [2 * j for j in range(1, 13)] + [13]
    parityrec = [2 * j for j in range(1, 12)] + [11, 13]
    for S, nm in [(GW, "GW {1..11,13,24}  (TIGHT: M = 1/14)"),
                  (primsat, "prim-sat 2{1..12}+{13}  (M = 1/13)"),
                  (parityrec, "parity record {2..22 even,11,13}  (M = 1/12)")]:
        test_set(S, nm)
    for k in range(3, 9):
        test_set(list(range(1, k + 1)), f"AP_{k} (control, M = 1/{k+1})")

    print()
    print("=" * 100)
    print("(2) NEAR-TIGHT NEIGHBORHOOD of GW: single-element moves (the adversary class)")
    print("=" * 100)
    n_eq = n_sep = 0
    seps = []
    base = set(GW)
    for e in sorted(base):
        for delta in (-2, -1, 1, 2, 3):
            e2 = e + delta
            if e2 < 1 or e2 in base:
                continue
            S = sorted(base - {e} | {e2})
            g = 0
            for i in range(12):
                g = gcd(g, S[i + 1] - S[i])
            if g != 1 or max(S) > 30:
                continue
            r = test_set(S, f"GW[{e}->{e2}]", verbose=False)
            if r is None:
                continue
            if r[0] == "EQ":
                n_eq += 1
            else:
                n_sep += 1
                seps.append((S, r[1], r[2]))
    print(f"    single-move GW neighborhood: mu = M on {n_eq}; separations {n_sep}")
    for S, M, dens in seps[:10]:
        print(f"      *** {S}: M = {M}, mu >= {dens}")

    print()
    print("=" * 100)
    print("(3) WHERE THE SEPARATIONS LIVE: |S|=4 exceptions' M values vs the tight bar 1/5")
    print("=" * 100)
    exc4 = [(1,2,9,14),(1,3,4,5),(1,3,4,7),(1,3,4,12),(1,3,4,14),(1,3,11,12),(1,3,13,14),
            (1,4,5,9),(1,7,8,9),(1,9,10,14),(1,11,12,13),(2,3,5,8),(2,5,7,12),(2,6,7,9),
            (2,6,7,13),(2,6,9,11),(3,4,7,10),(3,4,7,11),(3,5,8,9),(3,5,8,11),(3,5,8,13),
            (3,8,11,14),(3,9,10,13),(4,5,9,13),(4,5,9,14),(5,8,9,14)]
    Ms = [(M_exact(list(S)), S) for S in exc4]
    Ms.sort()
    print(f"    min M over 26 exceptions: {Ms[0][0]} at {Ms[0][1]}   (tight bar 1/5 = 0.2)")
    print(f"    max M: {Ms[-1][0]} at {Ms[-1][1]}")
    below = [x for x in Ms if x[0] < F(1, 5)]
    at = [x for x in Ms if x[0] == F(1, 5)]
    print(f"    exceptions with M < 1/5: {len(below)}; with M = 1/5 exactly: {len(at)}"
          f" {[s for _, s in at]}")
    print(f"    [every 4-element set has M >= 1/5 iff LRC(5) tight; separations AT the bar")
    print(f"     would kill even the tight-locus collapse at k=4 -- check above]")

    print()
    print("=" * 100)
    print("(4) |S| = 3 rigidity extension: max <= 18")
    print("=" * 100)
    n_eq = n_sep = 0
    first_seps = []
    for S in itertools.combinations(range(1, 19), 3):
        g = gcd(gcd(S[0], S[1]), S[2])
        if g != 1:
            continue
        states, t0, t1 = build_window_graph(list(S))
        if len(states) > 400_000:
            continue
        M = M_exact(list(S))
        pos, _ = positive_cycle_exists(t0, t1, M.numerator, M.denominator)
        if pos:
            n_sep += 1
            if len(first_seps) < 5:
                dens, L, B = extract_beating_cycle(states, t0, t1, M.numerator,
                                                   M.denominator, max(S))
                first_seps.append((S, M, dens))
        else:
            n_eq += 1
    print(f"    |S|=3, max <= 18, primitive: mu = M on {n_eq}; separations {n_sep}")
    for S, M, dens in first_seps:
        print(f"      *** {S}: M = {M}, mu >= {dens}")
    if n_sep == 0:
        print("    => the |S| = 3 rigidity persists through max = 18 (conjecture: all |S|=3)")

    print(f"\n[total {time.time()-t00:.0f}s]")

if __name__ == "__main__":
    main()
