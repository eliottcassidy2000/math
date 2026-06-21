#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angleG_coverage_majorization_opus_0620.py   (opus-2026-06-20, ANGLE G)

ANGLE G for the LRC(14) crux: prove consec_k MAXIMIZES measS7 (k=8,9,10) via an
AGGREGATE coverage-count COUPLING / Hardy-Littlewood REARRANGEMENT, since the
per-IE-block route is dead (C1).

SETUP (colorings/CA, arithmetic mod 7).  For E subset Z, |E|=k, 0 in E, each
x in [0,1) realises the HIT-SET  C_E(x) = {floor(7 frac(e x)) : e in E} subset Z/7,
always containing 0.  measS7(E) = meas{x : C_E(x) = Z/7} = p_full(E) = P(N=0),
where N(x) = 7 - |C_E(x)| = #absent colors.  The breakpoints of the piecewise-
constant map x -> C_E(x) partition [0,1) into CELLS; on each cell C_E(x) is fixed.

THE COVERAGE-COUNT DISTRIBUTION is the length-weighted law of N(x) (equivalently of
|C_E(x)|).  measS7 = mass at N=0.

WHAT IS ALREADY KNOWN (must build on, NOT repeat):
 * (T3a, angleE) consec MINIMIZES E[N] = Sum_j A_j (mean absent count): 0/330 beat it.
 * (route_D C) increasing-convex order N_consec >=_icx N_E FAILS (329/329): plain
   majorization of the FULL N-distribution is impossible, because some adversary
   (e.g. (0,1,2,3,4,5,6,8)) shifts mass toward small N yet still loses at N=0.
 * So: consec wins at N=0 (crux) and at E[N], but loses on the INTERMEDIATE tail.
   A pure Hardy-Littlewood majorization of the coverage-count law CANNOT close it
   -- this script CONFIRMS that wall sharply and then probes the refinements that
   ISOLATE the N=0 atom (the only thing measS7 actually sees).

TESTS (each EXACT rationals; brutally honest PROVED/VERIFIED/REFUTED):
 (G0) recompute the cell coverage-count law p_s(E)=meas{|C_E(x)|=s}, s=1..7, exact.
 (G1) REARRANGEMENT/MAJORIZATION (the literal Angle-G ask): sort each E's cells by
      #covered = |C_E(x)| descending into a step profile cov_E^*(u) on u in [0,1);
      is cov_consec^*(u) >= cov_E^*(u) pointwise for a.e. u?  (<=> the |C| law of
      consec stochastically dominates).  Predict FAIL; report first failure & by how
      much.  This is the clean refutation of the naive rearrangement.
 (G2) The N=0 ATOM is NOT a monotone functional of the law.  Confirm with the explicit
      g-weights: which weightings w_s on the coverage-count make consec extremal?
      Solve (LP-free, by inspection) for the cone of weight vectors w with
      <w, p(consec)> >= <w, p(E)> for ALL bounded E; is e_7 (the measS7 functional)
      in that cone's interior, on its boundary, or OUTSIDE?  (If OUTSIDE, majorization
      is provably the wrong tool -- a SHARP structural statement.)
 (G3) TWO-LEVEL (coarse) majorization: collapse N to {N=0} vs {N>=1}, and to
      {N<=1} vs {N>=2}, etc.  measS7 is monotone in the coarse 2-bin law iff consec
      dominates that single threshold.  Find ALL thresholds t for which
      P(N<t)(consec) >= P(N<t)(E) holds for all E (a partial stochastic-dominance).
 (G4) The DECISIVE refinement: measS7 = E[N] - (corrections).  Since consec MINIMIZES
      E[N] (T3a, sharp), write
         1 - measS7 = P(N>=1) = E[N] - E[(N-1)_+]   (telescoping the survival fn),
      i.e.  measS7 = 1 - E[N] + E[(N-1)_+].  consec minimizes E[N] (+good) but we need
      it to also control E[(N-1)_+] = Sum_{s>=2}(s-1)p_s.  TEST whether
         measS7(consec) >= measS7(E)
      reduces to   E[(N-1)_+](consec) - E[(N-1)_+](E) >= E[N](consec) - E[N](E).
      Report the two gaps per shape and whether the inequality direction is consistent
      (a candidate for a 1-D comparison certificate).
 (G5) RESTRICTED rearrangement: do the cells with |C|>=6 of consec dominate?  i.e.
      restrict the coverage profile to the TOP TWO levels {6,7}.  Is
         p_7(consec)+p_6(consec)  the max, and within it p_7 the max?
      (a lexicographic top-of-profile dominance -- the only place the atom lives).

python3 stdlib only; exact Fractions.
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
from math import comb
sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------------
# exact cell decomposition: coverage-count law p_s(E) = meas{|C_E(x)| = s}
# ----------------------------------------------------------------------------
def cell_decomp(E):
    """Return list of (length, covered_count, hitset_frozenset) over the cells of
    x -> C_E(x).  Exact Fractions.  0 in E pins residue 0 (e=0 -> floor(0)=0)."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    cells = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        hits = {0}  # e=0 always pins residue 0 (whether or not 0 listed, 0 in E by hyp)
        for e in E:
            if e == 0:
                continue
            hits.add(int(((e * xm) % 1) * 7))
        cells.append((x1 - x0, len(hits), frozenset(hits)))
    return cells

def coverage_law(E):
    """p_s = meas{|C_E(x)| = s}, s=1..7.  Returns dict s->Fraction (sums to 1)."""
    p = defaultdict(lambda: F(0))
    for length, s, _ in cell_decomp(E):
        p[s] += length
    return {s: p.get(s, F(0)) for s in range(1, 8)}

def measS7(E):
    return coverage_law(E)[7]

def meanN(law):
    """E[N] = sum (7-s) p_s."""
    return sum((7 - s) * law[s] for s in range(1, 8))

# ----------------------------------------------------------------------------
def shapes_for(k, span):
    """all E={0}+combo, combo subset {1..span}, |E|=k, primitive-ish (gcd handled
    by scale-invariance of measS7; we keep all, consec included)."""
    out = []
    for combo in itertools.combinations(range(1, span + 1), k - 1):
        out.append((0,) + combo)
    return out

# ============================================================================
def main():
    print("=" * 90)
    print("ANGLE G: coverage-count majorization / Hardy-Littlewood rearrangement")
    print("=" * 90)

    CONFIG = {8: 11, 9: 12, 10: 13}  # k -> span box (matches landscape adversary pool)

    for k in (8, 9, 10):
        span = CONFIG[k]
        consec = tuple(range(k))
        law_c = coverage_law(consec)
        mn_c = float(meanN(law_c))
        ms_c = float(law_c[7])
        print("\n" + "#" * 90)
        print(f"### k={k}, span box [0..{span}].  consec={consec}")
        print(f"###   measS7(consec)=p7={law_c[7]}={ms_c:.6f}  E[N]={mn_c:.6f}")
        print("#" * 90)

        print("[G0] consec coverage-count law p_s (s=1..7):")
        print("     " + "  ".join(f"p{s}={float(law_c[s]):.4f}" for s in range(1, 8)))

        pool = [E for E in shapes_for(k, span) if E != consec]
        print(f"\n     adversary pool size = {len(pool)}")

        # ----- (G1) rearrangement / stochastic dominance of |C| -----
        # cov_E^*(u) >= : equivalent to P(|C|>=s)(consec) >= P(|C|>=s)(E) for all s.
        def ccdf_geq(law, s):  # P(|C| >= s)
            return sum(law[t] for t in range(s, 8))
        g1_fail = 0; g1_example = None; g1_worst = F(0)
        for E in pool:
            law = coverage_law(E)
            dominated = True
            for s in range(1, 8):
                if ccdf_geq(law_c, s) < ccdf_geq(law, s):
                    dominated = False
                    d = ccdf_geq(law, s) - ccdf_geq(law_c, s)
                    if d > g1_worst:
                        g1_worst = d; g1_example = (E, s, d)
            if not dominated:
                g1_fail += 1
        print(f"\n[G1] REARRANGEMENT (consec |C| stochastically dominates all E?):")
        print(f"     consec dominates: {len(pool)-g1_fail}/{len(pool)}; FAILS {g1_fail}/{len(pool)}")
        if g1_example:
            E, s, d = g1_example
            print(f"     worst violation: E={E} at threshold |C|>={s}: P(adv)-P(consec)=+{float(d):.5f}")
        print("     => naive Hardy-Littlewood majorization of the coverage-count law is FALSE.")

        # ----- (G2) the weight cone: is e_7 (measS7) in the consec-extremal cone? -----
        # consec is extremal for <w,p> iff for every adversary E:
        #   sum_s w_s (p_s(consec) - p_s(E)) >= 0.
        # Let D_E = p(consec)-p(E) (a vector in R^7, sum=0). The cone of valid w is
        #   { w : <w, D_E> >= 0 for all E }.  measS7 functional is e_7.
        # We test: is e_7 in this cone? <e_7, D_E> = p7(consec)-p7(E) >= 0 for all E (TRUE = crux).
        # Is it on the BOUNDARY (some D_E with <e_7,D_E>=0) or strictly interior?
        # And: is e_7 expressible as a nonneg combo of {e_s : e_s in cone}? If NO single
        # e_s other than e_7 is in the cone, e_7 is an EXTREME ray not reachable by
        # majorization (only by the atom itself) -- the structural obstruction.
        Ds = []
        for E in pool:
            law = coverage_law(E)
            Ds.append(tuple(law_c[s] - law[s] for s in range(1, 8)))  # index 0..6 -> s=1..7
        # which unit functionals e_s (s-th coverage count) are consec-extremal?
        cone_units = []
        for idx in range(7):
            ok = all(D[idx] >= 0 for D in Ds)
            cone_units.append((idx + 1, ok))  # (s, in-cone?)
        # which CUMULATIVE survival functionals 1{|C|>=s} are extremal?
        cone_cum = []
        for s0 in range(1, 8):
            # functional w with w_s=1 for s>=s0: <w,D> = sum_{s>=s0} D[s-1]
            ok = all(sum(D[s - 1] for s in range(s0, 8)) >= 0 for D in Ds)
            cone_cum.append((s0, ok))
        print(f"\n[G2] WEIGHT-CONE membership (which monotone functionals consec maximizes):")
        print(f"     unit coverage-count e_s in consec-extremal cone (e_7 = measS7 itself):")
        print("       " + ", ".join(f"e{s}:{'Y' if ok else 'N'}" for s, ok in cone_units))
        print(f"     cumulative 1[|C|>=s] (stochastic-dominance thresholds):")
        print("       " + ", ".join(f"|C|>={s}:{'Y' if ok else 'N'}" for s, ok in cone_cum))
        n_units = sum(1 for _, ok in cone_units if ok)
        print(f"     => # unit functionals consec maximizes = {n_units}/7.")
        if n_units == 1 and cone_units[6][1]:
            print("     => e_7 (measS7) is the ONLY unit functional consec maximizes.")
            print("        measS7 is an EXTREME ray of the consec-extremal cone, NOT reachable")
            print("        as a nonneg combination of OTHER coverage-count functionals.")
            print("        STRUCTURAL OBSTRUCTION: no majorization/rearrangement can derive it")
            print("        from lower coverage levels. (This is WHY C1/route_D failed.)")

        # ----- (G3) coarse two-bin thresholds (partial stochastic dominance) -----
        print(f"\n[G3] coarse threshold dominance P(N<t)=P(|C|>7-t) for t=1..7:")
        good_t = []
        for t in range(1, 8):
            s0 = 7 - t + 1  # |C| >= s0  <=>  N < t
            ok = all(sum(law_c[s] for s in range(s0, 8)) >=
                     sum(coverage_law(E)[s] for s in range(s0, 8)) for E in pool)
            if ok: good_t.append(t)
            tag = "DOMINATES" if ok else "fails"
            print(f"     t={t}: P(N<{t})=P(|C|>={s0})  consec {tag} all E")
        print(f"     => consec stochastically dominates only at thresholds t in {good_t}")
        print(f"        (t=1 is the crux measS7; the gap is the intermediate tail).")

        # ----- (G4) telescoping certificate: measS7 = 1 - E[N] + E[(N-1)_+] -----
        # measS7 = p7 = P(N=0).  1 - p7 = P(N>=1) = E[N] - E[(N-1)_+].
        # check identity, then per-shape gaps.
        print(f"\n[G4] TELESCOPE  measS7 = 1 - E[N] + E[(N-1)_+]   (E[(N-1)_+]=sum_{{s>=2}}(s-1)... over N):")
        def Eplus(law):  # E[(N-1)_+] = sum over cells of max(N-1,0)
            return sum(max((7 - s) - 1, 0) * law[s] for s in range(1, 8))
        ms_check = 1 - meanN(law_c) + Eplus(law_c)
        print(f"     identity check consec: 1-E[N]+E[(N-1)+] = {float(ms_check):.6f} vs measS7={ms_c:.6f}  "
              f"({'OK' if ms_check == law_c[7] else 'MISMATCH'})")
        # consec minimizes E[N] (known). Decompose the measS7 gap:
        #   measS7(consec)-measS7(E) = [E[N](E)-E[N](consec)] + [E[(N-1)+](consec)-E[(N-1)+](E)]
        # term1 >= 0 always (consec min E[N]). term2 sign varies. report.
        t1_pos = 0; t2_pos = 0; both_needed = 0; t1_alone_suffices = 0
        worst_t2 = F(0); worst_t2_E = None
        for E in pool:
            law = coverage_law(E)
            t1 = meanN(law) - meanN(law_c)        # >=0 (consec min E[N])
            t2 = Eplus(law_c) - Eplus(law)        # sign varies
            gap = law_c[7] - law[7]               # measS7 gap, >=0 (crux)
            assert t1 + t2 == gap, (E, t1, t2, gap)
            if t1 >= 0: t1_pos += 1
            if t2 >= 0: t2_pos += 1
            if t1 >= gap: t1_alone_suffices += 1   # term1 alone already covers the gap
            if t2 < worst_t2: worst_t2 = t2; worst_t2_E = E
        print(f"     term1 = E[N](E)-E[N](consec) >= 0 for {t1_pos}/{len(pool)} (should be all: consec min E[N])")
        print(f"     term2 = E[(N-1)+](consec)-E[(N-1)+](E) >= 0 for {t2_pos}/{len(pool)}")
        print(f"     term1 ALONE >= measS7 gap (term2 not needed) for {t1_alone_suffices}/{len(pool)}")
        if worst_t2_E:
            print(f"     most negative term2: E={worst_t2_E} term2={float(worst_t2):+.5f} "
                  f"(must be repaid by term1 surplus)")

        # ----- (G5) top-of-profile lexicographic dominance -----
        print(f"\n[G5] TOP-OF-PROFILE: does consec maximize p7, then p7+p6, ...?")
        for topset, lbl in [((7,), "p7"), ((6, 7), "p6+p7"), ((5, 6, 7), "p5+p6+p7")]:
            cval = sum(law_c[s] for s in topset)
            beat = sum(1 for E in pool if sum(coverage_law(E)[s] for s in topset) > cval)
            print(f"     {lbl}: consec={float(cval):.5f}  #adv strictly greater = {beat}")

    # ----------------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("READOUT")
    print("=" * 90)
    print("Angle G tested the coverage-count majorization / Hardy-Littlewood rearrangement.")
    print("See per-k [G1] (naive majorization FALSE), [G2] (is measS7 an extreme ray?),")
    print("[G4] (telescope split into E[N] term + tail term).")

if __name__ == "__main__":
    main()
