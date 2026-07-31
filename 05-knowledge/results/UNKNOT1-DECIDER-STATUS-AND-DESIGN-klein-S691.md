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
- **Known theory:** u=0 (unknot recognition) is decidable, in NP and
  conditionally in co-NP; *diagrammatic* unknotting variants are NP-hard
  (de Mesmay–Rieck–Sedgwick–Tancer, "the unbearable hardness of
  unknotting").  Coward--Lackenby give strong structure and recognition
  results for genus-one knots, but their paper is not cited here as a
  general practical decider for `u=1`.  Applebaum et al. study practical
  Reidemeister/crossing-change search by reinforcement learning; its output
  is an upper bound, not an exact unknotting-number oracle.
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
   (integer congruence diagonalization, no floats); determinant.  An
   Alexander-polynomial route is a proposed extension, not implemented.
3. **False-certificates:** |σ| ≥ 4 (Murasugi); Lickorish linking form
   ±2/d test.  The implementation derives cyclicity and a linking
   generator from determinant divisors and adjugate entries; it does not
   compute a full Smith normal form.  Optional Alexander/Nakanishi module
   checks remain unimplemented.
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
  FALSE_CERTIFIED**.  The decision is justified by the obstruction theorem
  and exact implementation; the calibration suite is a regression control,
  not the source of validity.
- Hostile regression:
  `[[1,11,2,10],[6,10,7,9],[3,8,4,9],[11,5,12,4],[7,2,8,3],[5,1,6,12]]`
  is an unknot produced from the one-crossing unknot by reverse
  `R1,R2,R2,R3,R3`.  Greedy input reduction stalls, while crossing four
  changes to an R1/R2-reducible unknot.  The old engine returned
  TRUE_CERTIFIED; the repaired engine returns UNKNOWN because no input
  nontriviality certificate exists.
- The owner's 11-crossing example is exactly **K11n3**, with signed DT code
  `[4,8,10,-14,2,-16,-20,-6,-22,-12,-18]`, determinant `43`, signature
  `-2`, and verdict `UNKNOWN`.  All eleven crossing changes visible in this
  diagram have determinant different from one, so this diagram has
  diagrammatic unknotting number greater than one.  That does not decide
  the knot invariant: the [Knot Atlas K11n3 page](https://katlas.org/wiki/K11n3)
  records `u(K11n3)` in `{1,2}`.
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

## Exact algorithmic boundary

Visible-crossing search cannot be complete in general: Taniyama constructs,
for every nontrivial knot and every `n`, a diagram whose diagrammatic
unknotting number is at least `n`
([arXiv:0805.3174](https://arxiv.org/abs/0805.3174)).  The strongest simple
general result is a positive semidecision: first certify the input is
nontrivial with an exact unknot recognizer, fairly enumerate Reidemeister
paths, flip every crossing of every reached diagram, and apply the exact
unknot recognizer again.  It halts with TRUE whenever `u(K)=1`, but need not
halt when `u(K)>=2`.

There is a sharp promised-class exception: for an alternating knot with
`u(K)=1`, every alternating diagram contains an unknotting crossing
([McCoy, arXiv:1312.1278](https://arxiv.org/abs/1312.1278)).  Lackenby's
2026 hierarchy algorithm strengthens the exact unknot-recognition inner
oracle ([arXiv:2607.23350](https://arxiv.org/abs/2607.23350)); it does not
supply a finite outer list of all possible unknotting crossing arcs or a
one-hour worst-case bound for the benchmark.

## Repo-native remark

THM-2176/2191 (Gordian continuation cocycle, catalytic group length,
closure of 9_10) already live one street over: the continuation-kernel
machinery is the repo's own language for crossing-change moves. If a
session ever pushes the u=1 lane seriously, start from THM-2191's stable
envelope rather than raw diagram search, and respect MISTAKE-230–235
(no syntax-only bridges between knot moves and other carriers).
