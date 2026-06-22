#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_thread2_verdict_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 2 VERDICT)

THE THREAD-2 DELIVERABLE.  Three findings, all EXACT, on the TRUE GOOD set
   GOOD = {x: maxgap{frac(e_i x): e_i in E} > 1/7}   (genuinely x->-x symmetric, ghat REAL):

FINDING 0 (correction to the kpswf12 spectral setup).  The kpswf12 engine used GOOD_proxy =
   cover_set(E)^c = {>=1 of 7 half-open sectors MISSED}.  This proxy is NOT x->-x symmetric
   (sym-diff with its reflection = 31/245 at E={0..7}), so its ghat is NOT real and its R' is a
   DIFFERENT, larger-deviation quantity.  On the TRUE maxgap GOOD, ghat IS real and the floor is
   much tighter: R'_true in [0.85, 1.01] vs the proxy's [0.59, 1.21].  The complement-even premise
   of the prompt is EXACT only for GOOD_true.

FINDING 1 (cluster complement = x->-x; the even part is the WHOLE cluster up to reflection).
   GOOD_true(-E) = GOOD_true(E) and GOOD_true(N-E) = GOOD_true(E) EXACTLY (N=max+min), because
   reflecting all phases ph->1-ph preserves every cyclic gap.  So GOOD (hence ghat, hence SPEC)
   depends on the cluster ONLY through its complement-symmetric (even) orbit.  The complement-ODD
   part of the cluster is invisible to GOOD: it is mean-zero in the strongest sense -- it does not
   enter GOOD at all.  (Consec clusters {0..k-1} are self-complement-even.)

FINDING 2 (the frequency grading: there is NO universal mean-zero Z_2 channel).
   SPEC = S_even(n) + S_odd(n) (parity) = S[7|n] + S[7nmid] (apex).  EXACT via shift-projectors:
     S_odd = (meas(GP cap GOOD) - meas(GP cap (GOOD+1/2)))/2.
   S_odd is EXACTLY 0 for some (P,E) (e.g. P={1..5} or {2..6}, E={0..7}) but NOT others
   (P={1,2,3,12,13}: S_odd=1/546; coprime P: S_odd dominates).  G_P is not 1/2-balanced in general
   (fails on a generic test set), so the parity channel is NOT a uniform controlled/mean-zero split.
   CONCLUSION: the controllable handle is FINDING 1 (cluster-level complement-even), NOT a frequency
   Z_2 grading.

This file: (i) reprints the EXACT bank table on GOOD_true; (ii) the parity + 7 grading EXACT;
(iii) a RANDOM-sample floor scan of R'_true over many (P, cluster) to estimate the true uniform
floor c (the crux quantity), separating complement-even (consec/symmetric) vs general clusters.
"""
import sys, random
from fractions import Fraction as F
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import safe_set, meas, intersect
from lrc14_true_maxgap_good_kpswf13 import good_true
from lrc14_spec_grading_kpswf13 import shift, parity_graded_spec_exact, sevens_graded_spec_exact

def rprime_true(P, E):
    gp = safe_set(P); good = good_true(sorted(set(E)))
    base = meas(gp) * meas(good)
    if base == 0:
        return None
    return meas(intersect(gp, good)) / base

def main():
    print("#" * 96)
    print("# THREAD 2 VERDICT (kind-pasteur kpswf13): complement-even reframe on the TRUE GOOD")
    print("#" * 96)

    # admissible LRC(14): |P| small part, E cluster; the floor family is consec clusters k=8..13.
    bank = [
        ([1, 2, 3, 12, 13], list(range(8)),  "k=8 consec"),
        ([1, 2, 3, 4, 5],   list(range(8)),  "k=8 consec smallP"),
        ([1, 2, 3, 12],     list(range(9)),  "k=9 consec"),
        ([1, 3, 4, 5],      list(range(9)),  "k=9 floor-ish"),
        ([1, 2, 3],         list(range(10)), "k=10 consec"),
        ([2, 3, 4, 5, 6],   list(range(8)),  "killer P"),
        ([5, 7, 11],        list(range(10)), "coprime P"),
        ([1, 2, 6], [0, 4, 6, 8, 10, 12, 14, 15, 16, 17], "wide"),
    ]

    print("\n" + "=" * 96)
    print("BANK on GOOD_true: R'_true, exact parity grading (S_even,S_odd), 7-grading")
    print("=" * 96)
    print(f"  {'case':<22}{'R_true':>9}{'SPEC':>10}{'S_even':>10}{'S_odd':>10}{'S[7|n]':>10}{'S[7nmid]':>10}")
    for P, E, lab in bank:
        pg = parity_graded_spec_exact(P, sorted(set(E)))
        sg = sevens_graded_spec_exact(P, sorted(set(E)))
        if not pg:
            continue
        print(f"  {lab:<22}{float(pg['Rp']):>9.5f}{float(pg['SPEC']):>10.5f}"
              f"{float(pg['S_even']):>10.5f}{float(pg['S_odd']):>10.5f}"
              f"{float(sg['S_sev']):>10.5f}{float(sg['S_nonsev']):>10.5f}")

    # ----- the EXACT consec floor on GOOD_true: min R'_true over admissible P, k=8..13 -----
    print("\n" + "=" * 96)
    print("EXACT CONSEC FLOOR on GOOD_true:  min R'_true over consec E (k=8..13) x admissible P")
    print("  (admissible: |P| = 13-k, P subset of {1..13}; the LRC(14) small-part budget)")
    print("=" * 96)
    import itertools
    gmin = (F(10), None, None)
    for k in range(8, 14):
        psz = 13 - k
        E = list(range(k)); good = good_true(E); mG = meas(good)
        kmin = (F(10), None)
        # to keep runtime sane, scan P combinations of {1..13}
        combos = list(itertools.combinations(range(1, 14), psz)) if psz >= 1 else [()]
        for P in combos:
            if not P:
                continue
            gp = safe_set(list(P)); base = meas(gp) * mG
            if base == 0:
                continue
            R = meas(intersect(gp, good)) / base
            if R < kmin[0]:
                kmin = (R, P)
        print(f"   k={k:2d} (|P|={psz}):  min R'_true = {float(kmin[0]):.5f}  at P={kmin[1]}")
        if kmin[0] < gmin[0]:
            gmin = (kmin[0], kmin[1], k)
    print(f"\n   *** GLOBAL consec floor on GOOD_true: R'_true >= {gmin[0]} = {float(gmin[0]):.5f}"
          f"  (k={gmin[2]}, P={gmin[1]}) ***")
    mP = F(14249, 252252)
    print(f"   m_P = 14249/252252 = {float(mP):.6f} (witness-floor target)")
    print(f"   ratio (R'_true-floor)/m_P = {float(gmin[0] / mP):.1f}x")

    # ----- random scan over general (non-consec) clusters to test uniform floor -----
    print("\n" + "=" * 96)
    print("RANDOM SCAN of R'_true over GENERAL clusters (k in 8..11, random E subset of small span,")
    print("   random admissible P).  Estimating the uniform floor c and S_odd statistics.")
    print("=" * 96)
    random.seed(20260622)
    worst = (F(10), None)
    sodd_zero = 0; sodd_nonzero = 0; n_neg = 0; tot = 0
    sodd_abs_max = 0.0
    for trial in range(220):
        k = random.choice([8, 9, 10, 11])
        span = random.choice([k - 1, k, k + 2, k + 5, 2 * k])
        E = sorted(random.sample(range(0, span + 1), min(k, span + 1)))
        if len(E) < 4:
            continue
        psz = random.choice([3, 4, 5])
        P = sorted(random.sample(range(1, 14), psz))
        try:
            pg = parity_graded_spec_exact(P, E)
        except Exception:
            continue
        if not pg:
            continue
        tot += 1
        R = pg['Rp']
        if R < worst[0]:
            worst = (R, (tuple(P), tuple(E)))
        if pg['S_odd'] == 0:
            sodd_zero += 1
        else:
            sodd_nonzero += 1
        sodd_abs_max = max(sodd_abs_max, abs(float(pg['S_odd'])))
        if pg['SPEC'] < 0:
            n_neg += 1
    print(f"   trials with valid baseline: {tot}")
    print(f"   worst R'_true in random scan: {float(worst[0]):.5f}  at P={worst[1][0]} E={worst[1][1]}")
    print(f"   fraction with SPEC<0 (anti-correlated): {n_neg}/{tot} = {n_neg/max(tot,1):.2%}")
    print(f"   S_odd EXACTLY 0: {sodd_zero}/{tot};  S_odd != 0: {sodd_nonzero}/{tot}")
    print(f"   max |S_odd| over scan: {sodd_abs_max:.5f}  (vs typical |SPEC|~0.02-0.07)")

    print("\n" + "=" * 96)
    print("THREAD-2 VERDICT (summary)")
    print("=" * 96)
    print("""  (0) CORRECTION: the prompt's 'ghat REAL' premise is EXACT only for the TRUE maxgap GOOD,
      not the kpswf12 cover^c proxy (which is x->-x ASYMMETRIC, sineE~1.6e-2).  On GOOD_true the
      floor is MUCH tighter (R'_true in ~[0.85,1.01]); the consec floor on GOOD_true is far above
      m_P -- the Node-3 floor is even safer than the proxy indicated.

  (1) The complement-EVEN part of the cluster determines GOOD ENTIRELY (the strong form of the
      conjecture): GOOD_true(E)=GOOD_true(N-E) EXACTLY, so ghat and SPEC depend only on the
      complement-symmetric orbit of E.  The complement-ODD part of the cluster is invisible to GOOD
      -- mean-zero in the strongest sense.  This is the genuine uniform handle: it is the SAME
      x->-x = T^op symmetry that makes ghat real, read on the cluster side.

  (2) There is NO universal mean-zero FREQUENCY Z_2 grading.  S_odd (parity) and S[7|n] are EXACTLY
      0 for some (P,E) but order-|SPEC| for others; G_P is not 1/2-balanced in general.  So the
      finer-grading conjecture (split SPEC into controlled + mean-zero by n-parity / n mod m) is
      REFUTED as a universal statement.  The control lives at the CLUSTER level (finding 1), and in
      the L2/Cauchy-Schwarz tail already established, not in a frequency parity.""")
    print("DONE.")

if __name__ == "__main__":
    main()
