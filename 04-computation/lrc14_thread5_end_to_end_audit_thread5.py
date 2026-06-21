#!/usr/bin/env python3
"""THREAD 5 -- COMPLETENESS CRITIC end-to-end audit of the LRC(14) sector/cap/dichotomy route.

This is a DOCUMENTATION script: it records the logical chain, what is PROVED vs VERIFIED vs
ASSUMED at each link, and the single highest-priority remaining gap, with the supporting
computational findings from this session. It does NOT recompute (the checks live in the
companion thread5 scripts and their .out files).

CHAIN (current route; NOT the older covering/THM-527 route):

 [L0] LRC(14) <=> meas(S7(E)) <= cap_k for every primitive co-offset set E (0 in E, |E|=k).
      STATUS: ASSUMED (HYP-2603 + the slow/fast glue THM-527-A). This is a HYPOTHESIS, not a
      theorem. It is the upstream reduction that turns the runner problem into the p0/cap
      problem. cap_8..13 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7, 1. caps verified >= (k-6)/7.

 [L1] k<=7: pigeonhole (THM-530). STATUS: PROVED unconditional.

 [L2] k=13: cap=1 trivial; k<=7 pigeonhole; so the live rows are k=8..12.

 [L3] BOUNDED span<=14: max p0(E) <= cap_k over ALL primitive E with max(E)<=14.
      STATUS: PROVED (finite, exhaustive, exact). Independently reconfirmed this session
      (lrc_bounded_full_recheck_thread5.out): k=8..12, 0 violations, consec the binding
      maximizer (k=8 binding margin 319/5880=0.0543; positive at every k). SOLID LINK.

 [L4] WIDE span>14: the dichotomy (HYP-2788):
      (a) near-cap (p0 > Q(k-1)) ==> single-perturbation-bounded (remove ONE element ->
          reduced span<=14), hence reduces to SINGLE-FAR over a bounded (possibly DILATED) base.
      (b) genuine-wide (remove-one-fails) ==> p0 < Q(k-1) < cap (slack floor).
      STATUS: VERIFIED by ADVERSARIAL HEURISTIC BANK, NOT proved. The bank truncates far
      positions (far_start_hi=34) and bounded-base shapes (consec/AP/45-random). Stress-tested
      this session over a WIDER bank (far up to 120, 20k genuine-wide configs): 0 genuine-wide
      with p0>Q at k=8 (lrc_dichotomy_stress2_thread5.out). Slack floor PROVED only k=8,9
      (multi-cluster error aggregation p0 <= p0_inf + sum (6/49)V_i/g_i is verified-numeric,
      not closed-form for k>=10). Both halves rest on finite generated banks.

 [L5] SINGLE-FAR closed by THM-563: w*Delta_w periodic in w (period 7*lcm(B)) via the exact
      Dedekind/sawtooth identity; sup_w Delta_w*w = finite period-max; if period-max(B) <
      15*margin(B) then Delta_w < margin for all w>=15.
      STATUS: identity+periodicity PROVED. The closure condition period-max(B) < 15*margin(B)
      is VERIFIED ONLY for: consec_{k-1} k=8,9,10 (exact), and the SMALL-period bounded bases
      (the committed general script scans ONLY periods <= 30000 and SKIPS the rest).
      THIS SESSION'S KEY FINDINGS:
       (i)  the committed general script (lrc_periodmax_general_macmini_0621s6.py) only loops
            k=8,9 top-20; the "broad 133-base, worst 10.16, ALL CLOSE" output line was produced
            by an UNCOMMITTED broad variant and SKIPS 267/400 bases (period>30000).
       (ii) skipped-base audit (lrc_periodmax_skipped_audit_thread5.out): for k=8,9,10 the
            GLOBAL smallest-margin base is CONSECUTIVE and lives in the CHECKED (small-period)
            region; every SKIPPED huge-period base has margin >= 0.17 (>> consec's 0.12-0.18).
            So "huge-period => safe" is TRUE in the correlation sense, but never proved.
       (iii) the crude a-priori bound period-max <= sum_arcs range(S_j) = (6/49)*#arcs is
            PROVABLY INSUFFICIENT: it FAILS for 687/3003 bases at k=8 (worst ratio 33.3 >> 15)
            and 2570/3432 at k=9 (lrc_periodmax_apriori_skipped_thread5.out). So the skipped
            region CANNOT be closed by the margin alone -- signed cancellation among arcs is
            ESSENTIAL, and that cancellation has never been verified for the skipped bases.
       (iv) DIRECT exact period-max for 5 feasible skipped k=8 bases (P up to 60060): ALL CLOSE,
            worst ratio 6.56 < 15 (lrc_periodmax_spotcheck2_thread5.out). Reassuring but partial
            (the P=180k/360k/2.5M bases remain unscanned).

 [L6] DILATED base (continuous period-max): for a single-perturbation reduction that yields a
      dilated base d*X, the deviation f(s)=s*Delta_s(X) is periodic in REAL s, sup = contmax,
      and contmax < 14*margin closes it for real far speed s>=14.
      STATUS: VERIFIED ONLY for consec_{k-1}, k=8,9,10 (lrc_continuous_periodmax script). The
      same coverage gap as L5: not checked for general bounded X.

 [L-ledger] THM-527/530/531/546/547/557 + THM-534 Delsarte L_y<=cap: the reductions and the
      per-E dual bound. THM-534's per-E bound meas(S7)<=L_y is PROVED; "consec maximizes L_y"
      (HYP-2607) is OPEN but only needed for the (already separately-closed-by-L3) bounded rows.
      THM-546 (|Delta_w|<=(6/49)V/w) PROVED but LOSSY (O(1) in the ungapped regime); THM-563
      replaces it for single-far. THM-531 scale-invariance PROVED.

CIRCULARITY: none found. L3 (bounded) is self-contained finite; L4 reduces wide->single-far;
L5 closes single-far; L6 closes the dilated variant. No link uses LRC(14) to prove itself. The
ONLY non-self-contained dependency is L0 (HYP-2603) which is upstream, not circular.

UNBOUNDED-SPAN ESCAPE: the partition span<=14 (L3) / span>14 (L4) is exhaustive. Within span>14,
the dichotomy partition near-cap / genuine-wide is exhaustive (a config is in exactly one). The
ESCAPE RISK is not a missed region but a missed CONFIG inside L4/L5: the dichotomy banks and the
THM-563 period scan are both TRUNCATED finite searches, so a hypothetical adversarial wide config
with large far position or huge-lcm base is not provably covered.

SINGLE HIGHEST-PRIORITY REMAINING GAP (the load-bearing unproven step):
  L4 (the HYP-2788 dichotomy) is the most load-bearing UNPROVEN link, because:
   - everything downstream (L5 single-far, L6 dilated) is CONDITIONAL on the dichotomy having
     reduced the wide region to single-far over a bounded base;
   - the dichotomy is established ONLY by truncated adversarial banks (far<=34 originally,
     <=120 in my stress test) with NO closed-form proof of either half;
   - its genuine-wide slack-floor closed form (multi-cluster error aggregation) is proved only
     k=8,9; the near-cap=>single-perturbation half has NO proof at all.
  CLOSELY SECOND: L5's general bounded-base period-max(B) < 15*margin(B) is unproven for the
  huge-period (skipped) bases, AND the a-priori margin bound provably cannot close them
  (signed cancellation required, never verified there). This is a self-contained FINITE-but-
  -uncomputed gap; L4 is a genuinely OPEN structural lemma.
  ABOVE BOTH: L0 (HYP-2603) is an unproven HYPOTHESIS -- if one counts it, it is the true
  foundation gap, but the prompt treats the reduction as given, so within the dichotomy
  architecture the binding gap is L4.

VERDICT: the chain does NOT yet prove LRC(14). It is a complete ARCHITECTURE with one PROVED
finite half (L3), one PROVED mechanism (L5 periodicity), and the wide-region closure resting on
an UNPROVEN dichotomy (L4) plus an UNCOMPUTED finite tail (L5/L6 general-base period-max). No
counterexample exists anywhere; LRC(14) is almost certainly TRUE.
"""
print(__doc__)
