# Unknotting number one: exact status + in-repo decider design (klein-S691, 2026-07-28)

> **CORRECTION (root, 2026-07-28; MISTAKE-318).**  The first landed
> `TRUE_CERTIFIED` rule proved only that one crossing change produces an
> unknot, hence `u(K)<=1`; it did not prove that the input was nontrivial.
> An explicit six-crossing diagram of the unknot defeats greedy input
> reduction but becomes greedily reducible after changing crossing four, so
> the old code falsely returned `TRUE_CERTIFIED`.  The repaired rule requires
> an independent input-nontriviality certificate as well as the
> changed-diagram unknot certificate.  Without the former it returns
> `UNKNOWN`, correctly retaining the `u=0` versus `u=1` ambiguity.  The
> regression is now built into the 16-check ordinary/optimized suite.

**Prompt (owner puzzle bundle / Epoch FrontierMath open problem):** an
algorithm taking a PD diagram (≤ 100 crossings), returning True iff the
represented knot has unknotting number 1, under an hour on a typical
machine. (https://epoch.ai/frontiermath/open-problems/unknotting-number)

## Status ledger (checked 2026-07-28)

- **Decidability of u(K)=1 is OPEN.** Not known decidable, per the problem
  page itself and the literature; a provably-correct always-terminating
  decider would be a genuine theorem, not an engineering artifact.
- **Known theory:** u=0 (unknot recognition) decidable, in NP ∩ co-NP
  (Lackenby for co-NP); *diagrammatic* unknotting variants are NP-hard
  (de Mesmay–Rieck–Sedgwick–Tancer, "the unbearable hardness of
  unknotting"); **genus-one knots: u=1 IS decidable** (Coward–Lackenby,
  "Unknotting genus one knots"); Lackenby et al. 2025 study practical
  unknotting via RL (Experimental Mathematics).
- **Montesinos:** u(K)=1 ⟹ the double branched cover Σ(K) is ±d/2-surgery
  on a knot (d = det K) ⟹ H₁(Σ) cyclic; **Lickorish linking-form
  obstruction:** the linking form must represent ±2/d; **Murasugi:**
  |σ(K)| ≤ 2u(K); Ozsváth–Szabó d-invariant refinements exist (heavy).
- **Consequence for the benchmark prompt:** any honest submission is a
  three-valued engine (True-certified / False-certified / UNKNOWN) plus an
  argument that UNKNOWN is empty on the benchmark's test distribution —
  or a new theorem. State this loudly; never present the pipeline as a
  proven decider.

## The in-repo pipeline (04-computation/unknot1_decider.py)

1. **PD parsing + faces** (rotation system) → checkerboard graph.
2. **Exact invariants:** Goeritz matrix; signature via Gordon–Litherland
   (integer congruence diagonalization, no floats); det; Alexander
   polynomial (Seifert matrix route; evaluation-interpolation with exact
   rationals).
3. **False-certificates:** |σ| ≥ 4 (Murasugi); Lickorish linking form
   ±2/d test via Smith normal form of the Goeritz matrix (convention
   fixed against the classic u(7₄)=2 example before it is allowed to
   decide); optional Alexander/Nakanishi module checks.
4. **True-certificates:** first certify that the input is nontrivial
   (`det!=1` or `sigma!=0` in the current engine).  Then, for each crossing
   change in the given diagram, use the det/signature screens and greedy
   R1/R2 reduction to certify the changed diagram as the unknot.  Together
   these give `1<=u<=1`.  A changed-diagram certificate without input
   nontriviality gives only `u<=1` and therefore returns `UNKNOWN`.
5. **Escalation lane (not implemented day one):** bounded Reidemeister
   exploration before crossing changes (the u=1 change need not be visible
   in the input diagram — the classic completeness gap); d-invariant
   obstructions; SnapPy/Regina outsourcing where available.

## Results (core landed, klein-S691 close-out)

`04-computation/unknot1_decider.py` (stdlib-only, subsecond suite + example):
- Sanity 16/16: trefoil det 3, σ = −2 (Knot Atlas convention; mirror +2),
  TRUE_CERTIFIED with explicit certificate (flip crossing 1, R2, R1 → 0
  crossings); figure-8 det 5, σ = 0, TRUE_CERTIFIED; 7₄ det 15, |σ| = 2,
  linking form misses ±2/15 ⟹ **Lickorish fires decisively:
  FALSE_CERTIFIED** (the calibration suite passed, so the linking-form
  test earns decision rights, not report-only).
- Hostile regression:
  `[[1,11,2,10],[6,10,7,9],[3,8,4,9],[11,5,12,4],[7,2,8,3],[5,1,6,12]]`
  is an unknot produced from the one-crossing unknot by reverse
  `R1,R2,R2,R3,R3`.  Greedy input reduction stalls, while crossing four
  changes to an R1/R2-reducible unknot.  The old engine returned
  TRUE_CERTIFIED; the repaired engine returns UNKNOWN because no input
  nontriviality certificate exists.
- The owner's 11-crossing example: **det 43, σ = −2, verdict UNKNOWN** —
  all 11 in-diagram crossing changes screen out (flipped dets 3..59,
  none 1), no obstruction fires (H₁ cyclic, form contains ±2, |σ| < 4).
  det 43 ≠ 1 ⟹ NOT the Conway/Kinoshita–Terasaka knot (a first-recall
  guess corrected by computation).
- Brief-data corrections found by the build (hostile-control culture
  working): the fig-8 PD in the brief was invalid (an arc thrice) —
  rebuilt from the braid closure (σ₁σ₂⁻¹)²; the 7₄ PD was non-planar —
  rebuilt as pretzel P(3,1,3), pinned by det 15.
- Conventions: Goeritz + Gordon–Litherland, exact Fraction congruence
  diagonalization, both checkerboard colorings cross-checked;
  ETA_SIGN/TYPE_T0 calibrated on the trefoil.
- Limitations (unchanged from the design): no Alexander route yet; TRUE
  needs a separately certified nontrivial input and the unknotting crossing
  visible in the given diagram with monotone R1/R2 finishing; no
  d-invariants; UNKNOWN is an honest possible outcome — as it must be while
  u=1 decidability is open.

## Repo-native remark

THM-2176/2191 (Gordian continuation cocycle, catalytic group length,
closure of 9_10) already live one street over: the continuation-kernel
machinery is the repo's own language for crossing-change moves. If a
session ever pushes the u=1 lane seriously, start from THM-2191's stable
envelope rather than raw diagram search, and respect MISTAKE-230–235
(no syntax-only bridges between knot moves and other carriers).
