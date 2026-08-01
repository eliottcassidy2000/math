---
source: opus-2026-08-01-S4 (FC and planar JC in tandem; the isolated-coincidence framing; OpenAI ten-advances lesson)
status: >
  SYNTHESIS + one structural link + honest calibration. Both the Factorial Conjecture (FC(3), bare) and the
  planar Jacobian Conjecture (JC(2)) reduce, after all STRUCTURED cases are closed, to the same shape: "does an
  ISOLATED period/reducibility COINCIDENCE exist that every structured method misses?" For FC(3) the structured
  cases -- composition (Pakovich, kps: saddle weight R-independent) and isolated-max-modulus (opus: saddle cap,
  all degrees) and transversal families (kps: capacity=P through D=6) -- are closed; the residual is an
  isolated, non-composition, non-transversal psi (a Kontsevich-Zagier period-rigidity coincidence). For JC(2)
  the structured (tame/reducible) strata are closed (the repo rigidity cage); the residual is the Abhyankar-Moh
  leading-form non-descent. NEW structural link: the FC(3) residual lives EXACTLY on the fold {J=0} (opus:
  max-ridge = fold), while a KELLER map has J = const != 0 -- NO fold. So FC-counterexamples need a fold and
  Keller/JC maps forbid one: a fold/no-fold DUALITY on the Jacobian locus. CALIBRATION from the verified
  landmark (Alpoge/Fable JC(3) counterexample, constructive): isolated coincidences DO exist at n>=3; the
  honest expectation is set by dimension/symmetry ({e}<S_3<SO(3) for FC; n=2|n>=3 for JC). HONEST: FC(2)-
  homogeneous is Liu-Sun 2020 Thm 2.6 (known, cited); the tandem is a framing + one link, not a proof of either.
tags: [factorial-conjecture, jacobian-conjecture, JC2, keller, fold, isolated-coincidence, kontsevich-zagier, pakovich, abhyankar-moh, tandem, openai-lesson, liu-sun]
related: [fc3-decisive-max-ridge-is-the-fold-equals-kps-rank-drop-locus, FC2-FC3-tandem-and-pakovich-saddle-weight-R-independent-kps-S160, factorial-conjecture-exponential-periods-and-the-repo-state]
external: ["Alpoge/Claude Fable 5 JC(3) counterexample (Anthropic, 2026)", "OpenAI 'Ten advances in mathematics' 2026-07-28 (validation not established)", "Liu-Sun 2020 Thm 2.6 (FC homogeneous)"]
---

# FC and planar JC in tandem: the isolated coincidence, and the fold / no-fold duality

## 1. The methodological frame (the OpenAI ten-advances lesson)

The July-2026 landscape: the JACOBIAN CONJECTURE was DISPROVED by a **constructive** counterexample (Alpoge
with Claude Fable 5, Anthropic) -- an explicit `C^3` map, `det = -2`, weighted grading `(-1,1,2)`, non-
injective; **checkable by hand** (= repo THM-1300). By contrast the OpenAI "Ten advances in mathematics"
(2026-07-28) is a curated list whose "research status, AI contribution and independent validation ... have not
been established," and the Weil "10 Erdos problems" episode collapsed on one check (the problems were merely
"open on one person's site"). **Lesson, and the discipline of this lane: a constructive, hand-checkable object
(a counterexample; an exact identity) is worth more than a claimed proof; fetch the primary source; do not
carry a capstone stronger than it supports.** (This hour I nearly wrote "FC(3) closed" and stopped on reading
kps-S160's honest caveats.) FC(2)-homogeneous is Liu-Sun 2020 Thm 2.6 -- cited, not claimed.

## 2. FC(3): all STRUCTURED cases closed, one residual

`FC(3)-homogeneous <=> int_{Delta_2} phi^m dsigma = 0 forall m => phi=0` on the triangle (Liu-Sun form). In
the `C_3`-covariant family (only the `3|m` "leaks" survive):

- **Isolated-max (Morse) -- opus saddle cap, ALL degrees.** `int_T phi^{3k} ~ (3/k) M^{3k} sum_j p_j
  e^{3ik theta_j} != 0` (the `C_3`-orbit of the max-modulus point). Non-vanishing leak -> not a counterexample.
- **Composition `psi = R o W` -- kps-S160, ALL R.** The saddle weight `C = exp(psi_{D-1}/(D psi_D))` is
  R-INDEPENDENT (depends only on `W`), so `R` can only coarsen `W`'s roots-of-unity cancellation, never refine
  it -- no composition closes the `3|k` leak.
- **Transversal families -- kps-S159/S160, D<=6.** leak-Jacobian rank `= P` (no rank drop, no family).

The residual is exactly what is left uncovered: an **isolated, non-composition, non-transversal** `phi` -- a
genuine Kontsevich-Zagier period-rigidity coincidence. And (opus decisive step) it is pinned geometrically: its
max-modulus ridge lies on the fold `{J=0}` = kps's rank-drop discriminant locus, so it is a NON-Morse (critical
curve) `phi` on `{det J=0}`. Three independent methods (analytic saddle, algebraic transversality, arithmetic
ideal) converge on "no counterexample," leaving only this coincidence -- **strong evidence bare FC(3) is TRUE,
not a proof.**

## 3. JC(2): the same shape

`JC(2)`: every planar Keller pair `(P,Q)` (`{P,Q}=const != 0`) is a coordinate system. The repo's rigidity
cage closes every STRUCTURED stratum -- equivariant/graded (THM-1345/1370), quasi-homogeneous single-face
(THM-2113), one-fiber-linear (THM-2063), power-free top face (THM-2102), and the entire THM-2696-2722 slice
atlas -- each survivor is a TAME triangular automorphism. The residual is the **Abhyankar-Moh leading-form
non-descent**: propagate `{P_top,Q_top}=0` down the weight filtration. Same shape as FC(3): structured
(tame/composition) cases closed; an isolated non-tame coincidence is the only opening.

## 4. NEW link: the fold / no-fold duality on the Jacobian locus

The two residuals meet on the Jacobian:

```
   FC(3) counterexample : lives on the FOLD  {J = 0}     (opus: the max-modulus ridge is a fold).
   Keller / JC(2) map   : has   J = const != 0  -- NO FOLD anywhere.
```

So **an FC-counterexample REQUIRES a fold, and a Keller map FORBIDS one** -- they occupy complementary parts
of the Jacobian locus. (Honesty: `J` here is the REAL Jacobian of the single map `phi: R^2 -> R^2`, while the
Keller `det` is the COMPLEX Jacobian of a PAIR `C^2 -> C^2`; so this is a structural ANALOGY on `{det J=0}`,
not a literal duality -- but the shared "degenerate vs non-degenerate Jacobian" axis is real.) FC probes
whether a degenerate (folded)
polynomial can hide a period coincidence; JC probes whether a non-degenerate (fold-free) polynomial map can
fail to invert. The shared object is `{det J = 0}` -- FC's counterexample locus and the boundary Keller maps
avoid. (It also says the FC period-rigidity coincidence and the JC(2) Abhyankar-Moh obstruction are NOT the
same coincidence -- they sit on opposite sides of `{det J=0}` -- so a solution of one does not mechanically
give the other, matching the known "GM => JC needs GM(4), which is false" non-transfer.)

## 5. Calibration and honest status

Alpoge's JC(3) counterexample proves isolated period/reducibility coincidences DO exist -- at `n>=3`, on a
weighted grading the lower dimensions forbid (THM-1370: positive-graded Keller is invertible in every
dimension; the counterexample's forced grading `(-1,1,2)` is indefinite, impossible to realize in a way that
survives the FC/JC-2 symmetry). The honest expectation is therefore dimension/symmetry-calibrated:

```
   FC:  arc (n=2) / discrete S_3 (n=3) caps  -> FC(2),FC(3) lean TRUE; only the isolated KZ coincidence opens.
   JC:  n=2 open (rigidity cage, no counterexample) ; n>=3 FALSE (Alpoge, constructive).
```

Both conjectures are now "structured cases closed, one isolated coincidence remains," the coincidence is a
Kontsevich-Zagier period-rigidity phenomenon, and it lives on `{det J=0}` for FC and is dual to the fold-free
Keller locus for JC. **This is a unifying framing and one structural link, not a proof of either.** The
concrete joint frontier: decide whether the fold `{det J=0}` can carry the FC coincidence (kps's discriminant +
opus's constant-curvature/uniform-winding obstruction), and whether the fold-free Keller locus can carry the
JC(2) non-descent (Abhyankar-Moh). The constructive standard set by Alpoge is the bar: an explicit object, or a
real proof -- not a capstone.
