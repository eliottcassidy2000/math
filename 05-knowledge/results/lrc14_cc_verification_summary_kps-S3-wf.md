# Adversarial verification: cluster-collapse-threedistance advance on LRC(14) S3

Session: kps-S3-wf (2026-06-18). Verifier: kind-pasteur.

## VERDICT
The claimed result is a TRUE PARTIAL advance, exactly as self-declared. It does NOT close S3.
No counterexample to the PROVED pieces was found. The self-declared gap (limit infimum = 1,
no proof that finite covering+primitive S3 sets stay strictly above 1) is REAL and correctly
characterized. One sub-claim (v=max determinism) is correctly self-flagged REFUTED and I
independently reproduced its failure. One reported empirical CONSTANT is wrong (see below).

## PROVED pieces — independently re-derived exactly, 0 failures
1. CLUSTER-COLLAPSE LEMMA: J0=(1/(14 V0), 13/(14 Vmax)) safe for all u in [V0,Vmax];
   width(J0)*7*Vmax = (13 - Vmax/V0)/2; >1 iff Vmax<11*V0. Confirmed exactly; band contained in
   safe_components(L) and midpoint safe for 0/4000 random L. HOLDS.
2. ARC-WIDTH IMPLICATION C(S) => M(S)>=1/14: logically SOUND. v-danger teeth are isolated open
   arcs of width 1/(7v) separated by safe gaps 6/(7v)>0, so any arc of width >1/(7v) contains a
   v-safe point. Strict > is conservative (>= would already suffice). No violations found.

## Empirical criterion check on real covering+primitive S3 sets
- Across all valid S3 sets generated (hundreds confirmed under heavy machine load; thousands in
  prior sessions), C(S) (best_margin = max_v W(S\{v})*7*v > 1) held with ZERO failures.
- ZERO sets with M(S) < 1/14. LRC(14) holds on every S3 set tested.

## v=max determinism sub-claim: REFUTED (confirmed)
- Found S3 sets where removing v=max gives margin <= 1 yet C(S) holds via a DIFFERENT v.
- Concrete: S=[1,2,3,8,9,10,11,12,13,14,17,19,20]. Removing v=max=20 -> margin 0.636 (<1).
  v=19 -> 0.972 (<1). Witness is the INTERIOR large speed v=14 -> margin 49/40=1.225. C holds.
  M(S)=2/21~0.0952 >= 1/14. The witness need not be max(S).

## CORRECTION to a reported empirical constant
- Claim package reports "realized global min best-margin ~1.336".
- I found covering+primitive S3 sets with best_margin well below 1.336:
    * 49/40 = 1.225  at S=[1,2,3,8,9,10,11,12,13,14,17,19,20] (witness v=14), M=2/21.
    * 196/171 = 1.1461... at S=[1,2,5,8,9,11,12,13,14,17,19,20,28] (witness v=28), M=3/29.
  (21 of 2719 boundary-regime sets fell below 1.336.) So the stated 1.336 floor is NOT a valid
  lower bound; the true realized floor is <= 196/171 = 1.146. This narrows the gap to the limit
  value 1 but does NOT breach it (C never failed; not a counterexample to anything PROVED).
- The floor kept dropping the harder I searched (1.336 -> 1.225 -> 1.146), all strictly > 1 — the
  signature of an infimum at 1 that real covering+primitive sets approach but (empirically) never
  reach. This is exactly the unproven crux.

## Mechanism of the limit-1 gap (independent derivation)
- In the carry-phase limit the cluster danger teeth become c open arcs of width 1/7. If they could
  be EVENLY spread, widest-gap margin = 7/c - 1: c=3 -> 4/3 (= the reported test-shape value),
  c>=4 -> <1. The abstract infimum reaches 1 precisely because free (even) tooth placement is
  attainable only as V0->inf with unconstrained carry phase; actual covering+primitive arithmetic
  (and the max-over-all-v in C, which can drop a tooth) keeps the realized margin strictly above 1.
  This is exactly the unproven lifting (limit 1 -> realized >1) the package self-declares.

## Bottom line
- C(S) => M(S)>=1/14: PROVED (sound implication).
- Cluster-collapse lemma: PROVED.
- "C(S) holds for every covering S3 set" (which would CLOSE S3): NOT proved. Still the open crux.
  Empirically true (0 failures, broad search) but the limit reduction's infimum is exactly 1 with
  no rigorous lift to the realized floor. The advance is PARTIAL, as claimed.
- One empirical constant (1.336 floor) is too high; true realized floor <= 49/40 = 1.225.
