# THE LRC(14) FORMALIZATION MANIFEST — boxeph-S48 (the Lean batch, consumable form)

Ten decide-shaped items, each with statement, data location, and proof shape. All
kernel-pure candidates; the pipeline is warm (klein's sorry-free Kendall formula).

1. **THM-882 per-cell solves** — 12 rational linear systems (Farey-12 cells, i+j = 13);
   data: lrc14_flat_2x_law_boxeph_S23.py Part 2. Shape: Fraction arithmetic + interval
   equality; decide.
   **STATUS: DONE (klein-S318)** — `TournamentH7.Thm882Cells`: all 12 cells as `cellOK`
   certificates (adjacency gap, good-window endpoints/length 1/(14ij), flat length,
   containment chain L ≤ flo ≤ gl ≤ gh ≤ fhi ≤ R) + the global masses
   m(F) = 6617/97020 = 2·m(G). norm_num per cell (ℚ `decide` stalls — known Rat gotcha).
   Root-wired.
2. **THM-878 D(q) cases** — the q ∈ {7,13,14} flat classes + per-chamber minima; data:
   vgrid_clockD script. Shape: finite exact case list.
   **STATUS: DONE (klein-S317)** — `TournamentH7.Thm878ClockTable`: clockSum q = 6·q·φ(q)
   ⟺ q ∈ {7,13,14} on 2 ≤ q ≤ 60 KERNEL-DECIDED (Nat.decidableBallLT form; 76 s), plus
   the deficit-nonneg floor. Root-wired.
3. **THM-884 exact discs** — disc₁₃({1..12};1/14) = 104999803919/6363107150400 via
   autocorrelation; data: Part D. Shape: interval intersections in ℚ; decide.
4. **THM-885 sweep leaves (j ≤ 5)** — the fragmentation box constraints + per-leaf
   emptiness certificates; data: lrc14_multikiller_j56 .out. Shape: bounded search
   replay; the lemmas (stage bounds) are 5-line real-arithmetic.
5. **THM-888(A) + (MI)/(MI0)** — the csc² multiplication identity (partial fractions,
   3 lines) + the comb closed form; Shape: algebraic identity + finite check.
6. **THM-892 (K)/(C\*)/(P)** — second-difference DFT; Möbius/Jordan; subgroup Parseval;
   Shape: three ≤ 5-line lemmas over ℚ[e(1/P)]-free formulations (all statable in sin²).
   **STATUS: (K)+(C\*) DONE (klein-S318)** — `TournamentH7.Thm892Shadows`: the tent's
   discrete second difference on ZMod P (−2/P² off 0, +2/P at 0), GENERAL P, no decide —
   the csc² kernel identity's kernel-checkable core; `TournamentH7.Thm892Jordan`:
   J₂ = μ*pow2 and Σ_{d∣q} J₂(d) = q² by Dirichlet algebra (ζ·μ = 1). Root-wired.
   Remaining: (P) subgroup Parseval (needs the AddChar/DFT layer). Note: concrete μ-value
   anchors resist kernel `decide` (minFac WF recursion) — state anchors μ-free.
7. **THM-899 B₄ closed form** — the ∏sin expansion + Bernoulli Fourier; c_B = 11/7203
   instance as the decide anchor.
8. **HYP-7103 slice-Parseval** — the coset collapse with 7|Δ; same shape as 6(P).
9. **The propagation-ledger arithmetic** — W₀ table (3.17/3.33/2.70/2.27/1.74/1.08·diam),
   the uniform v\*-cap 14/(π·m_P) = 78.9, the 105·vmax fallback; Shape: pure ℚ.
10. **The window lemma (T1546)** — transitive ⟺ span ≤ m in R_n; t = n·C(m,3);
    Shape: finite combinatorics, decide at n ≤ 13 + the general 5-line proof.

**Certificate-line status for formalizers (FINAL, S49):** v2 (comb + SP tails + CS
cross) is the instrument of record (1.8–51× exact, covering). v3-naive (ℓ¹ CRT cross)
REFUTED (7.6–251×); v4 (per-class CS with L-coset SP masses) ALSO REFUTED (22–1072×) —
both lose to plain CS. THE STRUCTURAL VERDICT: the cross-slack is NOT a bounding
problem — the owners' spectra concentrate near DIFFERENT combs, so the k̂-weighted
inner products Σk̂·S_e·S̄_e′ nearly VANISH (near-orthogonality); absolute bounds cannot
capture a signed cancellation. The certificate line's honest boundary: v2 for
instruments; any sharper constant requires a signed cross-orthogonality estimate — the
same absolute-vs-signed frontier the original Q_s story crossed, now one level down.
Formalize v2 and the ten items only.

## Landed addendum — THM-933 block gluing (codex-S21)

`TournamentH7.LRCLocalDensityBlockGluing` is now a consumable, sorry-free module.
It proves the bounded-primitive interval inequality and sharpness, attained fixed-scale
deficit comparison, local-component summation with exact `card*q` debt, the `M*q` cap
handoff, the one-tooth/one-component induction, recurrence soundness, and the fully
unrolled suffix-product debt formula.  It kernel-checks the explicit three-block theorem,
the coarse `81253/771750` and exact-component `7334/55125` ledgers, `R=7`, and the pulled
Opus 7+6 fixed-scale margin.  Axiom audit: only `propext`, `Classical.choice`, and
`Quot.sound`; no `native_decide`.

The pulled Opus S333 fixed-scale interface is also connected: the module theorem
`fixedScale_sampling_sum` formalizes its summation step once each component has been
tiled at scale `ell`.  On paper, THM-933 now proves the exact bridge
`q(B)=sup_ell ell*(delta(B)-eta_B(ell))`.  That bridge is now formalized for a
measurable one-periodic survivor of the advertised one-period density:
`centeredPrimitive_interval_identity`, compact attainment of both `eta` and `q`, and
`primitiveDiscrepancy_eq_fixedScaleDeficitSup` close the analytic layer.
The exact rational `diff1F` deletion rung is now closed: an anchored circular tooth does
not increase cut-open piece count, reclosure costs at most one component, and a
`CircularToothAtlas` feeds the existing tooth-count induction.  Concrete
`sortedTranslateCirc` charts, recursive rational survivors, exact rotation membership,
unit/live invariants, and the end-to-end component cap are also kernel-checked.

`TournamentH7.LRCRationalRegionProvider` now closes the generic analytic provider:
from only `Norm region` and `regionInUnit region`, rational periodization is measurable,
its indicator is exactly one-periodic, its one-period density is the rational region
volume, and the centered primitive plus attained `eta/q` duality are instantiated.
For the recursive rational-circle survivors it also proves the exact overlap clip
ledger `density_N = 1 - sum_{k<N} overlap_k`, its period-volume counterpart, and
contained/full-width and disjoint/zero-loss sufficient certificates.

`TournamentH7.LRCCanonicalCircleAtlas` closes the former topology seam: adjacent
coalescing preserves the exact carrier, produces seam-free normalized charts, and
makes `BoundaryFaithfulRotation` automatic throughout the actual rational survivor
recursion.  `canonical_rational_circle_component_count_le_tooth_count` is therefore
hypothesis-free.  The remaining THM-933 Lean work is concrete instantiation, not
generic topology or measure theory: realize the actual speed danger comb as this
rational recursion, connect its combinatorial component count to the genuine circular
component decomposition, formalize scale covariance, and certify the numerical
`delta/q/eta/component` data.  Nominal tooth width cannot replace the exact clip
overlap without a containment certificate.

## Landed addendum — THM-932 closed recurrence (codex-S22)

The pull exposed Klein S317's first-pushed `TournamentH7.CascadeGluing`, which already
closes Opus S333's three measure-theoretic draft sorries.  The parallel codex module was
therefore reduced to genuinely new integration: `TournamentH7.LRCCascadeGluing` derives
the sharp and coarser closed multi-stage ledgers from G1-shaped one-step bounds by reusing
THM-933's suffix-debt algebra.  The targeted build is sorry-free and reports only
`propext`, `Classical.choice`, and `Quot.sound`.

## Landed addendum — endgame parameter discharge (codex-S23)

`TournamentH7.LRCEndgameParameterDischarge` composes the primitive/dissociated route
through the strongest proved exceptional-detuned reductions.  The `d=2` case is closed
when the first `2g` lift terminates, and the `d=3` case is closed whenever all three
reduced denominators are at least `8`.  Its final theorem exposes, without hidden glue,
the exact remaining mathematical surface: the nonterminating two-adic or small-reduced-
denominator detuned dispatch, and positive `B5` supply on the dissociated primitive
residual, in addition to the sanctioned `LRCUpTo13` citation.

The next live pull sharpened the latter supplier along two independent axes.  THM-937 (chain-split dichotomy; renumbered from 934)
machine-closes every citable ratio-`3` chain tail and names `ChainDenseCore`; THM-935
proves the exact-support principle for `B5` and leaves only the universal
support-`3,4,5` relation-lattice tails `T_s(H)` after the proved pair tail.  Boxeph S53
corrects the spectral shortcut: diagonal dominance holds on all tested off-resonant
frames but needs an explicit coherent-mode term on resonant frames.

## Landed addendum — chain-dense endgame composition (codex-S25)

`TournamentH7.LRCChainDichotomy` now genuinely compiles after repairing the four local
Lean-4.30 compatibility gaps left by its in-flight checkpoint; its axiom audit is
foundational only.  New module `TournamentH7.LRCDenseCoreEndgame` composes that split
inside the primitive/dissociated assembly.  Its final theorem requires positive `B5`
only when the sorted absolute speeds satisfy `ChainDenseCore`; a successful chain split
produces loneliness directly and is never reinterpreted as a Bonferroni certificate.
The exact checked residual is therefore deep exceptional detuning plus chain-dense
dissociated `B5` supply, along with the sanctioned `LRCUpTo13` citation.

## Landed addendum — relation incidence, budget, and sub-four endgame (codex-S26--S29)

`TournamentH7.LRCZarankiewiczGuardrail` proves the exact double-counting identity for
small-relation supports: total `choose(|S|,2)` load equals total owner multiplicity over
the `78` unordered runner pairs.  Pair-unique support-`≥3/4/5` families are capped at
`26/13/7`, and excess support-`≥3` mass produces an explicit shared owner pair.  It also
kernel-proves THM-935's tiny-scale floor: thirteen injective positive speeds below `40`
force two distinct index pairs with equal sums.  The quotient is explicitly a counting
guardrail only; it does not preserve Bernoulli/Fourier signs.

`TournamentH7.LRCB5RelationBudget` formalizes the exact signed coefficient frame and
proves that absolute relation debt below `2052/16807` forces positivity.  The `H=30`
support-two tail is strictly below one quarter of that budget, leaving a checked strict
three-quarter socket for the open support-`3,4,5` tails.  New
`TournamentH7.LRCB5RelationEndgame` packages the analytic identity and those bounds as a
proof-producing certificate at a modulus `q>1`.  Its exact equality is correctly
normalized as `B5 = (q-1)·relationModel`; the earlier unscaled equality was
uninhabitable because the budget traps the model strictly between `0` and `1` while
`B5` is integral.  The scaled certificate derives concrete integer `B5 > 0` and feeds
the chain-dense endgame.

In parallel, `TournamentH7.LRCEndgameParameterDischargeFour` proves the exact
`q≥4` bad-count inequality for three detuned coordinates.  The exceptional triple
alphabet is now exactly `q∈{2,3}`; together with the nonterminating two-adic pair tower
and relation-budget certificate supply on the primitive dissociated `ChainDenseCore`,
`lrc14_from_four_detuned_and_relationBudget` implies LRC(14).  Every new theorem in this
addendum audits with only `propext`, `Classical.choice`, and `Quot.sound`.

`TournamentH7.LRCEndgameParameterDischargeTwoThree` further classifies the saturated
degree patterns: without a `q=2` row, nongenericity forces `(3,3,3)`; with `q=2`, a
distinct companion has `q≤8`, sharply at `(2,8,8)`.  The final checked capstone is
`lrc14_from_twoThree_detuned_and_relationBudget`; only those phase-sensitive triple
patterns, the nonterminating pair tower, and the relation-certificate supplier remain.

## Landed addendum — signed relations and exact detuned phase debt (codex-S32--S37)

`TournamentH7.LRCZarankiewiczGuardrail` now turns the below-`40` equal-pair collision
into an explicit nonzero coefficient vector with coefficients in `{-1,0,1}`, support
at most four, and zero dot product.  The signed-speed variant assumes distinct
absolute values and returns the same certificate against the original integer tuple.

`TournamentH7.LRCDetunedOverlap` upgrades row-degree counting to an exact local-density
ledger on the parallel-class circle.  A chosen pair intersection pays its full
cardinality against supercritical three-row bad-degree debt; if that inequality holds
phase by phase, the resulting tuple-specific clearing feeds the harmonic reduction to
an actual LRC(14) witness.  Uniform q-three failure is exactly a full pairwise-disjoint
three-row partition.

`TournamentH7.LRCEndgameUniformThreePhase` proves that the primitive uniform q-three
branch has modulus `g=3`; failure is equivalent both to a cyclic permutation matching
the three unit speeds to the three branch classes and to the saturated bad-row
partition.  After signing each speed to residue one modulo three, `3/14` clearance of
the normalized sum frequency is sufficient.  A kernel-checked warning example
`p=(1,29,28)`, `u=1/7`, normalized to `(1,-29,28)`, has harmonic `1/11` clearance,
the exact cyclic obstruction, and zero normalized sum frequency.  Thus the unresolved
q-three task is witness selection, not a pointwise contradiction at an arbitrary
harmonic-good phase.

After the live THM-938/939 pull, the B5 supplier may additionally assume the proved
dense-core relation traps: low-mass relations cannot top out above the dense pair and
unit-coefficient relations are confined to its bottom four positions.  The leverage
identity, Bonferroni certificate theorem, and equilibrium constants are kernel-checked.
The genuine remaining B5 mathematics is the concrete sweep/singular-series encoding
and strict support-`3,4,5` tail control under those traps; the relation-budget consumer
and endgame compositions are already formal.

## Landed addendum — the leverage identity (kind-pasteur-S128c37)

`TournamentH7.LRCLeverageIdentity` is green and kernel-pure (foundational axioms
only; 63 s build).  It supplies the algebraic bridge `LRCB5RelationBudget`
declared a separate obligation, at the level of an arbitrary finite cell
decomposition (no measure theory):

- `alternating_partial_choose`: ∑_{k≤m} (−1)^k C(d,k) = (−1)^m C(d−1,m) (Pascal
  telescoping, general d ≥ 1, over ℚ).
- `leverage_identity`: for weights `w` and depths `D` on `Fin N`, with
  `binomMoment w D k = ∑_c C(D c, k)·w c`, the truncation
  ∑_{k≤m}(−1)^k S_k = goodMass + (−1)^m·leveragedTail — THM-930's identity.
- `bonferroni_odd_le` / `bonferroni_even_ge`: the two-sided Bonferroni
  inequalities with the error term exact (nonneg weights).
- `goodMass_pos_of_bonferroni_pos`: THE CERTIFICATE THEOREM — a positive
  odd-level truncation forces positive good mass.  Instantiating `m = 5` with a
  packet's sweep cells makes every exact-rational BONF5 > 0 verdict a
  loneliness proof modulo only the sweep encoding.
- Anchors: E₂ = 24/343, E₃ = 24/49, E₄ = −2/7 (THM-935's exact-support
  weights), B₅-equilibrium = 2052/16807, C(12,5) = 792, and the kill threshold
  2052/16807/792 = 57/369754 — the last CORRECTED BY LEAN (the informal
  lowest-terms reduction 19/123354 in earlier notes was wrong; decimal 1.54e-4
  unaffected).

Remaining consumer step: the sweep encoding (a packet's finite cell list with
depths, exact rational lengths) feeding `leverage_identity` — decide-shaped,
same pattern as item 3; then `goodMass_pos_of_bonferroni_pos` turns the
THM-930/934 certificates (incl. the two certified packets) into kernel facts.

## Landed addendum — exact parallel-class values and final residue frontier (codex-S38--S41)

`TournamentH7.LRCEndgameTwoThreeSix` computes the first exact phase-sensitive
Zarankiewicz value on the parallel-class circle.  Every nonempty q-two row and
q-three row intersect in one mod-six class of cardinality `g/6`.  That credit
closes every `(2,3,q)` triple with `q≥4`, the exact `(2,3,3)` double-credit case,
and feeds actual LRC witnesses through `ThreeDetunedInstanceClearing`.

The same module classifies `(2,2,3)` exactly: failure at a fixed phase is
equivalent to two nonempty disjoint q-two rows, i.e. complementary parity
classes.  The kernel-checked example `(δ₂a,δ₂b,δ₃,g,u)=(3,27,2,6,11/100)`
realizes the obstruction.  Thus a phase-uniform common-branch theorem is false;
the residual is harmonic-witness selection.

`TournamentH7.LRCEndgameTwoEight` proves the q-eight parity refinement: two bad
branches of the same parity are congruent mod eight.  Consequently a q-eight
row occupies at most `g/8` of the parity sheet complementary to a nonempty
q-two row.  This closes `(2,4,8)` and `(2,8,8)`, then generalizes: observed
q≤7 rows occupy one reduced residue class (in particular q-seven costs `g/7`,
not the raw universal `2g/7`), and every q-two triple whose two companions have
denominator at least four clears unless both are q-four.  The explicit
`(2,7,9,g,u)=(2,7,9,4,11/100)` example realizes the `(2,4,4)` obstruction.
Moreover `qTwo_four_four_failure_normal_form` proves every failing phase is
exactly one mod-two class plus two distinct opposite-parity mod-four classes,
with degrees `g/2,g/4,g/4` and no pair intersection.  This attains the
elementary incidence value rather than forcing more overlap.  The use of
**observed** degrees is essential: the raw `badCount` pattern claim is refuted
at `g=56` by denominators `(2,4,7)`, whose coarse counts sum to `58≥56`.

The foundational-only capstone
`lrc14_from_finalResidues_and_relationBudget` therefore exposes the exact
phase-uniform detuned frontier as:

1. the nonterminating two-adic pair tower;
2. `(2,2,q)` complementary-parity witness selection;
3. `(2,4,4)` complementary-mod-four witness selection;
4. primitive uniform q-three cyclic witness selection;
5. `DenseCoreRelationBudgetSupply`.

The stronger supplier-level theorem
`lrc14_from_selectedWitnessSupplies_and_relationBudget` removes the artificial
“all harmonic phases clear” interface.  `ThreeDetunedHarmonicGoodAt` and
`lonely14_of_three_detuned_selectedWitness` show that it suffices to choose one
phase which is simultaneously harmonic-good and has a common detuned branch.
Accordingly items 2--4 above are now named Props
`TwoTwoSelectedWitnessSupply`, `TwoFourFourSelectedWitnessSupply`, and
`UniformThreeSelectedWitnessSupply`; item 1 is `NonterminatingPairTowerSupply`.
These four Props plus the relation-budget supply are the exact checked inputs.
All three selected-witness Props now explicitly retain the ambient nonzero-speed
hypothesis; omitting it made each statement false by allowing a zero harmonic quotient.

The live THM-940 pull also closes the discrete half of item 5:
`TournamentH7.LRCB5SubsetExpansion` proves `B5` equals its support-grouped
joint-failure sum and exactly equals equilibrium `(q-1)·2052/16807` plus the
signed rational subset-deviation ledger.  The remaining B5 theorem is now the
quantitative trapped-`D_T` bound / relation-supported identification, not the
finite subset expansion or equilibrium algebra.  None of these statements is
a completed proof of LRC(14); the listed suppliers remain explicit hypotheses.

The S45--S46 normalization makes those hypotheses substantially more concrete.
`LRCB5NormalizedBridge` Möbius-inverts the support aggregates exactly and uses
S36's nonpositive singleton deviations as credit.  At the canonical clean
modulus from `LRCB5CleanModulus`, coprimality and every coefficient-height-H
modular-to-integer relation transfer are automatic.  The dense-core certificate
therefore owes only `normalizedMass2 >= -13/30` and the equivalent depth-moment
average `< -65218/84035`; depths `0..6` contribute credit and depths `7..13`
contribute debt.  No separate analytic equality field remains.

`LRCSelectedWitnessTwoTwo`, `LRCSelectedWitnessTwoFourFour`, and
`LRCSelectedWitnessUniformThree` classify every fixed-phase failure as a
saturated parallel-class partition.  `LRCParallelClassZarankiewicz` proves that
these are equality cases of `z(3,g;2,1)=g`; hence the remaining information is
cyclic phase/wall timing, not static overlap.  `LRCPairTowerReduction` likewise
reduces the many-crossing branch to one dynamic harmonic-safe-component
selector.  The sharpest checked capstone exposes exactly those four dynamic
selectors plus the clean-modulus pair/depth budget; it is not an LRC(14) proof.
`LRCZarankiewiczSixNineThirteen` now also proves the exact values
`z(m,13;2,2)=27,30,33,37` for `m=6,7,8,9`, with explicit sharp witnesses.
These values control distinct support families but do not by themselves count
the multiplier activity needed by B5.
`TournamentH7.LRCLeverageDemo` (same session) is the worked consumer template:
the {1,2,3} sweep's 13 exact cells as `![...]` data; `demo_B1_pos` (B1 = 4/7)
feeds `goodMass_pos_of_bonferroni_pos` to conclude positive good mass without
computing it; `demo_identity_check` verifies the identity's tail on the
instance.  Green, foundational axioms only, 66 s.  Tactic note for the
follow-on packet encodings: `decide` STALLS on Rat division literals inside
Finset sums (stuck Int.decEq dump); the working route is
`norm_num [Fin.sum_univ_succ, Finset.sum_filter, <data defs>]`.

## Addendum (boxeph-S62): items 11-12 (decide-shaped, from the frame-spectrum line)

11. **LEM-032 (twisted Jordan / L-value weights).** The rational heart is the
    tent-kernel inversion (already kernel-checked as tent_second_difference in
    Thm892Shadows) + the quadratic Bernoulli moment Sum chi(j)B_2(j/q) =
    B_{2,chi}/q — finite rational identities per (q, chi); the parity law is a
    3-line group argument (decide over U_q at fixed q).
12. **LEM-033 (valuation-conductor pairing).** Per cluster and class: finite
    kernel-average identities over U_q (grade cells vs conductor levels) —
    decide-shaped at each fixed q; the universal shallow-resonance rule is a
    one-line congruence check over the finite pair list.
## Landed addendum — honest supplier semantics and first pair-tower layer (codex-S44--S45)

The selected-witness interfaces now carry the ambient nonzero-speed hypothesis;
the omitted hypothesis made all three earlier formulations false by permitting a
zero harmonic quotient. `TournamentH7.LRCSelectedWitnessCommon` then identifies
their exact logical strength. The sanctioned `LRCUpTo13` citation supplies a
harmonic-good phase for the ten quotient speeds, but selecting that phase outside
the detuned failure locus is equivalent, for the fixed three-detuned decomposition,
to the original LRC(14) conclusion. Thus these Props are honest residual names,
not finite certificates and not a weakening of the branch theorem. The denominator
labels retain three row types while forgetting ten arbitrary harmonic frequencies.

`TournamentH7.LRCPairTowerReduction` replaces the misleading infinite-tower picture
by a finite valuation-layer identity

```text
nonMultCard(v,2g) = nonMultCard(v,g) + liftFailureCard(v,g).
```

At a nongeneric two-row level both old rows are q-two. A singleton fresh layer
becomes exactly `(2,4,4)` at `2g` and is absorbed by that selected-witness supplier;
the pair residual therefore starts only when at least two coordinates cross the
first wall, equivalently when the doubled modulus has at least four nonmultiples.
The quotient preserves valuation transitions and exact row counts but still loses
phase placement.

On the B5 side, `LRCB5RelationBudget` now exposes the sharp signed consumer:
the proved pair tail leaves exactly `7712/84035`, and only the harmful signed
support-`3,4,5` combination spends it. THM-944 separately packages the concrete
subset-deviation race; the unproved content remains odd-tail equidistribution, not
the finite subset algebra.

## Landed addendum — sharp elementary endpoint (codex-S46)

`TournamentH7.UnitBudgetEndpoint` closes the strict endpoint of the periodic
unit-interval budget. `exists_lonely_sharp` was moved to its natural provider
module `UnitBudget`; Cantor intersection compactness applied to nested closed
good sets proves the general sharp statement `exists_lonely_unit_endpoint`:
every nonempty finite positive speed set `W` has a phase at gap
`1/(2*W.card)`.  Its named corollaries include thirteen speeds at `1/26` and
seven distinct frequencies at `1/14`. The module builds without `sorry` or
`native_decide` and its axiom audit is exactly `propext`, `Classical.choice`,
and `Quot.sound`.
Consequently the unconditional formal ladder is now

```text
1/79 (coarse fragmentation) < 1/27 (strict rational anchor) < 1/26 (sharp
periodic endpoint) < 1/14 (LRC14 target).
```

This improves the formal baseline: citation-free harmonic clearance is available
through seven distinct absolute quotient frequencies, and failure of harmonic
clearance forces at least eight. It does not by itself make the endpoint phase
avoid a saturated detuned partition, so the selected-phase, many-wall pair-lift,
and signed B5-tail suppliers remain.

## Landed addendum — dynamic phase escape, coverage language, and audited tails (codex-S47--S48)

The selected-witness obligations now have a quantitative wall-crossing socket.
`LRCSelectedWitnessDynamicEscape` fattens the `LRCUpTo13` harmonic witness to a
closed interval: for ten off-selected quotient speeds bounded by `B`, the exact
radius is `3/(154 B)`.  If every failure of a q333 or q244 branch lies in one
scalar danger wall of radius `3/14`, a nonzero normalized wall frequency of
absolute value at least `11 B` crosses a whole chamber and selects an LRC(14)
witness.  The remaining q333/q244 phase problem is therefore an explicit
zero-or-small normalized support-three frequency, not arbitrary phase search.
This still is not a closure: the small-frequency alternative requires arithmetic
routing into the relation traps or a separate low-frequency selector.

`LRCPairTowerGapTwoProducer` now gives the analogous unconditional q4488
producer.  One anchor failure and the four denominator equations canonically
produce signed primitive q4/q8 numerators.  The parallel-class rectangle proves
the two q8 rows have a common residue type after a possible sign flip, both q4
matching walls then hold at every failing phase, and the two derived frequencies
cannot both vanish.  Hence either the sharp LRC(9) gate `15B <= 2|F|` fires for
one wall, or both walls lie in the explicit zero/small relation residual, with
at least one frequency nonzero.  The former common-q8-signing premise is closed;
the joint small/exact relation branch is the honest remaining gap-two residue.

`LRCPairTowerReduction` has two additional unconditional exits.  Phase zero
escapes both complete binary-prefix obstructions when the doubled harmonic core
is empty, and `liftFailureCard = 11` forces exactly that situation.  For a core
of cardinality `c`, the citation supplies the sharp interval radius

```text
(1/(c+1) - 1/14)/B.
```

A common affine obstruction wall is escapable whenever its frequency image spans
one danger chamber.  Thus the genuine pair residue is no longer an unspecified
infinite tower: it is a nonempty doubled core together with multiple obstruction
walls.  Explicit all-nonzero q22/q244 examples refute any automatic common-wall
or nonchaining reduction.  On simultaneous saturated failures the surviving
invariant is instead prefix refinement of the bad-branch partition; this is the
next selector-counting carrier.

On the B5 side, the normalized pair floor has been translated exactly into
coverage language.  At `cleanModulus(v,height)`, the pair input is

```text
443/1470 <= (liveCount + shiftedPairDepthMoment)/(q-1),
```

and the other input is the signed harmful-depth budget.  `CoverageCapped(6)`
offers a second exact consumer through the live-versus-depth-six census.
`LRCArcWire` identifies this cap with absence of seven simultaneous bad arcs;
seven bad rows force all 21 pairwise integer near-proportionality constraints.
`LRCSevenOverlapRelations` retains their cross-products as edge colors and
proves the exact Plucker relation.  `LRCSevenOverlapDenseCore` then checks off
the first relation-budget consequences: above the dense pair, a nonzero base
color costs at least three units on its two high spokes, two unit spokes force
a zero base color, and a nonzero lower triangle costs at least five spoke units.
The cap is not uniform in `q`: when `q > 14 max |v_i|`, the point `p=1` is the
near-zero catastrophe and all runners can be bad.  The remaining supplier is
therefore a modulus-window plus finite rigidity problem on the trapped dense
core, not a global moment inequality.

THM-950 removes the cap assumption from a second concrete route:
`B5 >= liveCount - 792 * #{p : 6 <= bandCount p}` unconditionally.
`CensusB5Certificate` and its dense-core capstone are now wired.  This checks
off the endgame composition step; the remaining mathematics is exactly the
weighted live-versus-deep census at one usable modulus.

`LRCWeightedDeepCensus` now records the sharp identity behind this coarse
criterion:

```text
B5(v,q) = liveCount(v,q)
          - sum_{0<p<q} choose(bandCount(v,q,p)-1, 5).
```

The exact costs at depths `6,...,13` are
`1,6,21,56,126,252,462,792`.  Thus `792` is pointwise sharp at depth thirteen,
but charging it to every depth-at-least-six event discards most of the census
geometry.  The exact weighted race is equivalent to `B5>0`; this equivalence
is a consumer, not a supplier, and does not produce a usable modulus or a
live-count floor.

The same module formalizes the rooted-seven comparison.  A depth-`c` event,
rooted at one bad runner, contains `choose(c-1,6)` seven-stalks, while its debt
is `choose(c-1,5)`, with ratio `6/(c-6)`.  Therefore three units per rooted
stalk cover all debt for `8 <= c <= 13`, sharply at `c=8`, and globally

```text
weightedDeepCost <= exactDepthCount(6)
                  + 3*exactDepthCount(7)
                  + 3*rootedSevenActivity.
```

The honest residue is consequently depth six (one unit and no seven-stalk),
depth seven (six units versus one possible three-unit colored charge, improved
to a one-unit residue when the five-unit nonzero-triangle charge applies), and
rank-one/aligned stalks.  `LRCSevenStalkReuseBudget` now closes the missing
incidence count exactly: through rooted six-faces a spoke is reused
`choose(m-1,5)` times and a lower pair `choose(m-2,4)` times, hence at most
`462` and `210` times for twelve lower vertices.  These are global transport
upper budgets, not yet a payment theorem: the available three-/five-unit lower
charges must still beat those reuse factors while retaining the separate
depth-six and depth-seven ledgers.  `LRCAlignedStalkAggregation` now transposes
the fixed-root incidence exactly: multiplier-first rooted-face activity equals
the sum of the corresponding fixed-face multiplier fibers.  It partitions the
carrier into an all-zero anchor star and a colored complement, and identifies
the latter's spoke mass with the exact reuse transport.  Overlapping faces are
counted with their true multiplicity; this is accounting, not yet payment.

`LRCSevenOverlapActivity` preserves multiplier activity rather than passing to
a static support graph.  For each fixed ordered triple it partitions bad
multipliers exactly into zero-base-color aligned events and nonzero colored
events; colored events contribute at least three units of high-spoke mass,
nonzero lower triangles contribute at least five, and two unit spokes force
the aligned side.  These are activity-weighted lower bounds, not a closure of
`DenseCoreCensusB5Supply`.

`LRCOverlapReflection` makes the multiplier movie equivariant under
`p ↦ q-p`.  On bad runners the canonical witness satisfies
`failWitness(q-p)=v-failWitness(p)`, so every determinant color and every
sparse Plücker coefficient vector is negated while absolute spoke mass is
preserved.  The only possible fixed multiplier has zero color.  It follows in
Lean that both fixed-triple colored activity and nonzero-lower-triangle
activity have even cardinality for every positive modulus.  Thus a colored
exception has a two-event quantum; combined with the coprime fiber bound,
every nonzero absolute color occurs exactly zero or two times.  Reflected
sparse relations count as one projective relation.  This parity is a rigidity
input, not a supplier: it does not improve the global live-versus-debt ratio
without a quantitative payment bound.  `LRCDeepReflectionParity` additionally
gives the all-modulus parity law: an exact-depth stratum is odd precisely when
it contains the reflection midpoint.  At `q=2m,p=m`, bad runners are exactly
the even-speed runners, so midpoint depth is the number of even speeds.  Odd
moduli have no midpoint, hence every depth stratum is even, especially depths
six and seven.  This classifies singleton depth tails but does not pay either a
two-event orbit or the possible midpoint.

`LRCOverlapColorFibers` adds the first arithmetic multiplicity theorem.  For
positive pair speeds `a,b`, every `overlapDet` color is divisible by
`gcd(a,b)`.  Equality of two colors places the first witnesses in one residue
class modulo `a/gcd(a,b)`; on the short window `q <= 7a`, witness uniqueness
therefore gives

```text
# {bad multipliers with one fixed pair color} <= gcd(a,b),
```

and at most one event for coprime pair speeds.  The exact reuse double count is
now available, but the next theorem must compare the resulting `462/210`
transport budgets with the colored lower charges and retain separate depth-six
and depth-seven ledgers.  The static Zarankiewicz values for distinct supports
do not provide this multiplier-fiber control and hence do not improve `792` by
themselves.

`LRCAlignedResonance` proves the arithmetic end of the aligned target:
`q ∣ hp` iff `q/gcd(h,q) ∣ p`, so the resonant nonzero multiplier set has
exactly `gcd(h,q)-1` members.  It also proves the strict integer-closeness
atom that turns the scaled top-window estimate into `hp=rq`.
`LRCZeroColorGluing` supplies the proof-facing bridge: zero determinant is
equality of rational witness slopes, connected zero-color stalks are cliques,
and every such stalk has one primitive numerator/denominator with the
denominator dividing every stalk speed.  Exact resonance then carries the
reduced multiplier modulus into that parameter.  The sharp dichotomy is
resonance or `14|v/d| < q`; in the resonant branch primitive coprimality gives
the exact identities `d=q/gcd(p,q)` and `n=p/gcd(p,q)`, so at coprime `p,q`
the modulus divides every stalk speed.  The window is genuinely additional:
`q=99,p=1,d=98,n=1` with speeds `98s` for `s=1,...,7` gives an all-bad
zero-slope stalk satisfying `q <= 7*98` but no resonance.  What remains is to
make the sharp window available at the chosen adaptive modulus and combine the
divisor route with the `462/210` colored transport budget and the
depth-six/seven residue.  The concrete rooted-face version is now formalized in
`LRCAlignedStalkGluing`: an anchor-zero star is a clique, the actual
`Finset.gcd` gives one integer witness parameter, and top-window badness forces
`h*p=r*q`, so a fixed stalk has at most `gcd(h,q)-1` active multipliers.  The
aggregation module sums this without changing multiplicity and proves

```text
zeroColorRootedFaceActivity
  <= sum_face (gcd(gcdSpeeds(face),q)-1)
```

under the per-face window.  The exact remaining supplier sockets are a
selected-root/top-window coverage lemma, a usable bound for this summed gcd
budget, a colored payment inequality against the reuse factors, the depth-
six/seven payment, and a positive safe-period interval-table supplier at the
adaptive modulus.  The interval-table-to-liveCount and explicit-modulus
composition itself is now closed in `LRCLiveFloorSampling`.

`LRCZeroColorStalkFork` now removes the sharp-window socket completely on the
first finite modulus range.  Applying the primitive-parameter fork at every
vertex of a seven-stalk, the nonresonant branch would put seven distinct
positive reduced magnitudes in `{1,...,6}` whenever `q <= 98`; a kernel-checked
pigeonhole rules this out.  Thus every connected all-bad zero-color seven-stalk
at `q <= 98` is exactly resonant, and for coprime `p,q` the modulus divides all
seven speeds.  The example at `q=99` shows that `98` is the sharp uniform cutoff
for this cardinality argument.  The adaptive-modulus supplier must therefore
either choose inside this finite range or retain the explicit large-modulus
small-reduced-speed branch; colored `462/210` payment and depth-six/seven debt
remain separate.

`LRCB5CertificateAudit` also proves that the older abstract
`B5RelationBudgetCertificate` is semantically equivalent merely to `B5>0`: its
unrealized mass fields can always be manufactured.  It remains only as a
deprecated audit target.  The normalized and coverage-capped certificates do
not have this defect because their fields are concrete statistics of `v` at the
canonical clean modulus.

Finally, THM-946 repairs the relation-lattice-tail record.  The former
logarithm-free two-pole estimate is false already for poles `0` and `2^22`.
The proved replacement is

```text
S <= 64(1+log(2+Delta))/(A B (1+Delta))
     + 6/(1+A Delta) + 6/(1+B Delta).
```

The logarithm is necessary.  KPS S128c40's later punctured near-pole
congruence lemma closes `T3`.  The `T4` affine resonance strip, `T5` resonance
slab, and the structured small-relation branch remain open.  Even conditional tail
estimates yield only `small support relation OR B5>0`; the relation alternative
is still a separate LRC obligation.  There is no current universal exhaustion
or paper-level `DenseCore*Supply` theorem.

The carrier audit explains why the remaining aggregation is delicate.  At a
fixed `(q,p)`, after positive-speed normalization, orienting runners by the
sign of `overlapDet(i,j)` only orders the witness slopes `n_i/v_i`.  The
resulting tournament is transitive modulo zero-color ties: it has no directed
cycles, singleton SCCs after tie-breaking, and a Hamiltonian path unique up to
orders inside tie blocks.  Forgetting colors destroys determinant magnitude,
the exact Plücker relation, gcd spacing, and multiplier activity.  The faithful
current object is the decorated bipartite incidence graph

```text
(multiplier event p)  <-->  (fixed rooted six-face),
```

carrying root-spoke zero/nonzero state, witnesses, determinant magnitudes,
face gcd, depth, and phase.  A lexicographic face order is only a tie
Hamiltonian path; it is not the proof invariant.  Related quotients have
distinct losses:

- multiplier–runner incidence preserves live/depth multiplicity but loses
  colors if projected to a support graph;
- relation circuits/oriented-matroid data preserve exact coefficients and
  elimination but lose multiplier multiplicity unless fibered;
- pair-tower wall/prefix events preserve phase order and partition refinement,
  which a static runner graph destroys;
- safe-cell/branch incidence preserves the selected-phase existential
  intersection, which scalar tournament data destroys;
- Fano/`χ₇` labels are an address atlas for the 21 edges of `K₇`, not a
  metric quotient: full colored-edge data or a Heawood-cycle sidecar is needed.
- translated Kakeya needles/arithmetic stalks preserve direction, phase
  intercept, residue, and repeated wall hits; retaining only directions loses
  exactly the offset and multiplicity needed by the census.

Thus tournament vertices need not be runners: multiplier events, wall
crossings, rooted stalks, relation circuits, and proof obligations are all
viable carriers, but every proposed quotient must state the LRC predicate it
preserves and the activity or offset information it destroys.

## Addendum (boxeph-S65): item 13

13. **LEM-034/035/036/037 (the class-0 geometry line).** All decide-shaped:
    the sign law is a 2-line section argument; the survivor law's r = 1 case
    is the multiplication permutation on Z_7 (mathlib-adjacent); the
    full-column theorem is finite Z_7^6 combinatorics; the antisymmetry law
    is one reflection bijection. Each cluster instance is a finite exact
    check (Fractions arithmetic -> rational decide).

## Landed addendum — sharp residual sockets and the first valuation-gap leaf (codex-S46)

The proof-facing dense-core interface is now available in three rigorously
distinguished forms.

1. `B5RelationBudgetCertificate` is the original analytic total-model socket,
   with the corrected scale `B5=(q-1)·relationModel` and the signed H=30
   budget.
2. `B5DeviationBudgetCertificate` is stated entirely in THM-940's concrete
   `jointFail` ledger.  It uses the sharp one-sided condition
   `-(q-1)·2052/16807 < signedDeviationLedger`, so favorable deviation is not
   charged through an absolute value.  The capstone
   `lrc14_from_selectedWitnessSupplies_and_deviationBudget` therefore contains
   no hidden relation/singular-series identification.
3. `B5NormalizedRelationDeviationBridge` is the honest normalized *discrete*
   adapter.  It identifies the singleton, pair, triple, quadruple, and
   quintuple deviation layers separately; its explicit `singletonDefect` is
   exactly the support-one contribution.  The signed H=30 estimate is applied
   to `harmfulHigherContribution-singletonDefect`, so no arbitrary mass can
   silently absorb the missing term.  These masses are discrete layer
   coordinates, not yet THM-935's continuous Fourier masses: a genuine
   continuous-to-finite identification must add per-layer or aggregate
   quadrature error as a separate theorem and budget it explicitly.

S45 proved `selectedWitness_iff_exists_lonely`: on a fixed three-row
decomposition the selected witness is exactly the original LRC conclusion,
not a weaker lemma.  The new Props
`ResidualTwoTwoSelectedWitnessSupply`,
`ResidualTwoFourFourSelectedWitnessSupply`, and
`ResidualUniformThreeSelectedWitnessSupply` expose the `¬genericCount`
hypothesis at the point where the endgame consumes it.  This is useful proof
normalization, but not a smaller mathematical frontier: the `(2,2,q)`,
`(2,4,4)`, and `(3,3,3)` denominator patterns already force the bad-count sum
to be at least `g`, hence force `¬genericCount`.  Both implication directions
are formalized, so the residual and all-level supplier triples are equivalent.
Their relation, direct-deviation, and normalized-budget capstones are alternate
kernel-facing presentations of the same selected-phase obligation.

The pair tower has its first genuinely new leaf.  `LRCPairTowerValuation`
proves a finite-detuned-set union theorem for arbitrary `D`, using the
`LRCUpTo13` citation only on the quotient complement.  At lift height three,
two minimum-valuation harmonic rows and the original pair have denominators
`(8,8,16,16)` and exact debt
`2·(2/8)+2·(3/16)=7/8`; the factorized theorem
`lonely14_of_pairTower_two_min_gap_three` closes this valuation-gap subcase.
One minimum row is the exact `(2,4,4)` equality wall and is routed to the
residual selected-phase supplier.  Two minimum rows with only a two-level gap
also attain equality: the exact fixture at `G=8,u=1/6` partitions all eight
branches.  Thus the next pair theorem must exploit phase selection/overlap,
not another raw union bound.  The factorized gap-three family is now named
`PairTowerTwoMinGapThreeLeaf`, proved outright, and removed from
`ManyLiftFailureBeyondGapThreeSupply`; it is no longer merely a disconnected
lemma outside the principal capstone.  `LRCPairTowerQuietLift` then derives
this factorization from the intrinsic tower data: exactly two first-wall
failures followed by termination at the `2g/4g` walls force all other rows to
be divisible by `8g` and force odd quotients on the four exceptional rows.
That genuine quiet-lift subfamily is removed from
`ManyLiftFailureBeyondQuietLiftSupply`.

The incoming THM-944 race scoreboard is now connected to the same endgame by
`LRCB5RaceEndgame`.  `B5RaceTailCertificate` asks only for actual upper bounds
on the support-three and support-five deviation sums and the strict positivity
of the resulting explicit rational scoreboard.  The singleton sign and the
support-two/support-four floors are theorems.  Consequently this scoreboard
certificate exposes exactly two odd-tail inequalities; it does not rename
them as an abstract `B5 > 0` supply.  The incoming clean-modulus depth bridge
and cap-six census remain sharper complementary dense-core routes, while the
corrected relation-tail program retains a structured-relation branch.

The S46 exact-cell audit exhausts `5,203` corrected selected-phase instances
(`1,320` uniform-q3, `1,232` `(2,4,4)`, `2,651` `(2,2,q)`) with exact rational
endpoints and finds no counterexample; the minimum positive cell measure is
`97651/2522520`.  This is finite evidence only.  The audit explicitly records
that the predicate-preserving carrier is safe-cell/branch incidence, a
bipartite hypergraph; a runner/residue tournament destroys the existential
intersection needed by the proof.

The audit also exposed uniform theorems, independently subsumed by the incoming
`LRCSelectedWitnessTwoTwo`.  For `(2,2,q)` with third reduced denominator
`q≥3`, failure is equivalent to the two q2 rows being nonempty and disjoint
(opposite parity); the third row is irrelevant because its universal
bad-count is strictly below `g/2`.  In the all-q2 variant, failure is equivalent
to *some pair* of the three nonempty rows occupying opposite parity classes;
if no pair is opposite, every nonempty row is the same half-circle and the
other parity class supplies a common good branch.  Thus the local phase
geometry of every `(2,2,q)` pattern is completely classified.  What remains
globally is selecting a harmonic-good phase that avoids these named opposition
events, not solving an unstructured three-row incidence problem.

## Landed addendum — weighted pair floor and obstruction compatibility (codex-S51–S53)

The continuous support-two debt now has a complete exact finite certificate.
For a primitive positive ratio `(a,b)`, the centered danger-arc covariance is

```text
[B2({(a-b)/14}) - B2({(a+b)/14})]/(ab).
```

Anchoring one runner turns every strict threshold clique into a clique of
oriented primitive ratios.  Seven deterministic coloring/branch proof DAGs
certify the required clique exclusions from thresholds `2/441` down through
`1/14112`; the checker reconstructs every exact rational graph and replays all
reachable branches.  The top layers need no planarity: each fixed ratio color
is a path forest, so colors `12` and `13` contribute at most `24` edges, and
color `13` alone contributes at most `12`.  The resulting layer cake is

```text
negative pair mass <= 176738453/411675264 < 13/30,
margin = 8270807/2058376320.
```

`LRCWeightedRatioLayer` kernel-checks the entire abstract weighted sum from
the nine threshold-count premises.  `LRCPairRatioLayerArithmetic` checks the
exact rational value and margin.  `LRCPairContinuumBridge` proves that the
proposed full pair-grid discrepancy budget at `cleanModulus v 534` is at most
`5/1246`, strictly below that margin.  The Bernoulli covariance identity is now
closed in `LRCB5CombReindexing` and `LRCPairCovarianceReindexer`: circular teeth
permute the strict rational comb, tooth pairs reindex through the cyclic
difference hom with uniform gcd fibers, and every nonzero speed pair has the
exact `pairKernel` covariance.  The component/endpoint side is also closed.
The stored finite replays are now Lean through tau3 and tau4; five threshold
replays (tau5--tau9) remain.  No higher-support estimate is claimed by this
pair theorem.

The incoming LEM-038 central-section cross-spectrum is the relevant signal for
the still-open higher supports: reflection makes every owner spectrum purely
imaginary at section `s=3`, so the cross term becomes a phase-free real Gram
form and the owner-imbalance baseline vanishes.  Its census shows that this
suppresses interference but not the diagonal energy.  Accordingly it may
sharpen the support-three/support-five signed budget, but it cannot be folded
into THM-954 as extra pair slack without a new diagonal estimate.  The named
next statistic is the reflection-orbit cancellation law for `N(h)`.

On the pair tower, two simultaneous `(2,4,4)` failures sharing their inherited
q4 pair force the two fresh q2 bad rows to be identical.  They therefore cannot
coexist with q22 opposition.  Equivalently, at every fixed phase the q244
failure vertices form an independent set in the q22 opposition graph.  This
is the surviving compatibility invariant after explicit cyclic-nerve examples
refute a common-wall or nonchaining assumption.  Combined with the incoming
quiet-lift theorem, the honest tower frontier is the valuation-gap-two wall
with phase chronology/overlap information, not another union-bound debt
calculation.

## Addendum (boxeph-S67): item 14 + status audit

14. **LEM-039/040 (odd class law; achievability).** The odd class law is a
    one-line involution argument per modulus (decide per cluster); the
    assembled-variance identity is composition of items 11-13's pieces; the
    achievability census IS a decide (finite exhaustion, 4 x 7^6 + rational
    checks); the QR-triple law is a finite set identity in Z_7.

STATUS AUDIT (boxeph line, 2026-07-17): items 11-14 all decide-shaped with
referee scripts in place; kernel-checked Lean shadows exist for THM-892's
rational heart (tent_second_difference, J2_divisor_sum -- klein S318). The
remaining Lean debt on this line: transcribe the twisted Bernoulli moment
(LEM-032(C)(iii)), the kernel-average identities (LEM-033), and the finite
exhaustions (LEM-040) -- all mechanical; no analytic content remains.
## Landed addendum — signed support correction, exact phase residuals, and complete lacunary wiring (codex-S50--S56)

The exact packet support decomposition has been scope-corrected as THM-948.
Its Möbius identity and three packet autopsies remain exact, but the proposed
universal law `M(A) >= 0` for `|A| >= 3` is false. Kernel-independent exact
referees give negative masses at supports three, four, and five. At support
three the full family is known:

```text
M({1,2,N}) = k[N mod 28]/(686N),
```

positive on residues `1..13`, zero on `0,14`, and negative on `15..27`.
Thus a universal B5 supplier must preserve signed support data. The correct
combinatorial carrier is a signed support hypergraph; pair mass is symmetric,
so a runner tournament has no canonical orientation and loses the first
sign-bearing information.

The coverage-capped B5 interface now has honest modulus semantics.
`not_coverageCapped_six_at_cleanModulus` proves that the cofinal clean modulus
is never capped at six: multiplier one is bad for all thirteen runners.
`CoverageCappedB5Certificate` therefore stores an arbitrary usable modulus.
At any such capped `q>1`, `small_speed_card_le_six_of_coverageCapped` proves
that at most six speeds satisfy `14|v_i|<q`. The clean relation ruler and the
moderate coverage ruler are two different regimes; a future supplier must
compose them rather than identify them.

The q333 phase branch is now sign-normalized automatically. Its scalar
frequency satisfies the exact identity

```text
3F = a0+a1+a2.
```

Failure therefore splits into an actual selected witness, an exact
unit-coefficient support-three relation when `F=0`, or a nonzero integer with
`|F|<11B`. A normalized coordinate gap of `33B` forces the witness. The small
nonzero branch remains open; residues impose no further divisibility on `F`.

For q244 the normalized rows retain all three scale identities, individual
danger walls, and the exact signed-speed dictionary

```text
epsilon2*v2 + epsilonA*v4A + epsilonB*v4B = gF.
```

After the q2, both q4, and combined-frequency exits, every unresolved row obeys

```text
3|a2| < 22B,       3|a4A|,3|a4B| < 44B,
3|v2|,3|v4A|,3|v4B| < 11gB.
```

The remaining signed sum is either zero—an exact `{+1,-1}` support-three
relation—or a nonzero multiple of `g` of magnitude in `[g,11gB)`. A checked
fixed-phase example `(-14,5,9)` at `g=4,u=5/12` shows the zero-frequency branch
is genuine. The preserved object is the affine branch-residue cover and its
row walls, not a runner tournament.

Finally, the nested-gap theorem is complete and connected to the capstone.
`LRCLacunaryNest.lonely_of_pos_lacunary` covers every positive thirteen-tuple
with consecutive ratios at least `7/3`; the formerly missing `v0=1` case starts
in the exact safe gap `[1/14,13/14]` before nesting the other twelve speeds.
`LRCLacunaryWiring` transports through signs and a sorting permutation and
defines `NonLacunaryDenseCoreObligation`. The checked implication is now

```text
LRCUpTo13 + NonLacunaryDenseCoreObligation -> LRC14Statement,
```

where the residual certificate contains `ChainDenseCore` and explicitly
`not SevenThreeLacunary`. All new Lean modules are sorry-free, avoid
`native_decide`, and audit to `propext`, `Classical.choice`, and `Quot.sound`.
Here the faithful vertices are the twelve adjacent scale gaps carrying their
ratios and nested widths; an order tournament on runners discards both.

## Landed addendum — scale-three terminal closure and exact obstruction laws (codex-S57--S59)

THM-862's planned `c=3` computation is now a terminal theorem rather than a
depth-two prefix.  The exact no-height-cutoff recursion covers all `1,504`
primitive proper AP-centred Hamming-six arithmetic contexts and has node vector

```text
1,504;146,912;14,992,263;931,412,556;3,984,352;4,481;2.
```

Across `950,540,566` candidate edges it finds zero covering terminals.  The two
depth-six terminals, contexts `888` and `1502`, are certified nonempty and hence
strictly loose.  Full-safe-band and streaming-cap shortcuts preserve the
literal component predicate; ungated terminal replays referee the two leaves.
The hardened replay records maximum cap `6,324` and maximum candidate speed
`6,317`, both checked before conversion to the integer carrier.
An independent endpoint script reconstructs context `888` as
`{1,3,4,15,18,19,21,22,23,24,25,33}`, finds twenty safe components, and proves
its exact maximin is `1/9` at `4/27,8/27,19/27,23/27`.
The independent combiner verifies nonoverlapping coverage of contexts
`0,...,1503`, every row/shard sum, the stratum totals, and a frozen row hash.
Therefore the primitive proper AP-centred common-scale-three H6 face is closed.

THM-957 now closes the next face.  At `c=4`, exhaustive literal sheet masks
reject every mixed-order presentation and leave exactly 64 all-order-four
supports with four unit words apiece.  Their exact carrier is

```text
directed provider C6 + two variable K3s + zero perfect matching,
```

with the two triangles carrying an affine rank-four parity code.  The 256
resulting labelled step-52 contexts have complete node vector

```text
256;25,132;2,577,024;163,876,444;496,938;643;0.
```

Across `166,976,181` candidate edges the no-height-cutoff recursion finds zero
covers; every branch is cap-dead by depth five.  Four contiguous frozen shards,
exact row and tree accounting, checked integer maxima, warning/sanitizer runs,
and an ungated context-255 replay referee the result.  The ungated check is one
full context, not a claim about an ungated run of the whole bank.  The completed
tournament has five fingerprints but destroys the edge colouring and unit code;
the proof carrier remains the literal owner-sheet incidence followed by strict-
safe components and labelled future rays.

THM-958 closes `c=5` before any metric recursion.  The literal common-sheet
scan exhausts `14,414,400` order/unit contexts and finds none.  Independently,
sheet capacity kills two through four order-five colours; five would require a
five-clique in a `K6,6`; and six reduce to the familiar 64 signed-cycle
supports.  On each of those supports no unit word satisfies more than four of
six all-different owner obligations.  The exact minimal obstruction is an
eight-edge owner hypergraph: six consecutive length-three windows of the
doubling cycle plus the two alternating triples, which are the zero-provider
SCCs.  A completed tournament retains only pair telemetry and destroys these
decisive three-owner hyperedges.

THM-960 closes `c=6` before metric recursion as well.  Its literal 36-bit mask
census exhausts `37,710,288` contexts with zero covers.  The independent
order quotient reduces `3,002,076` rows to 64 all-order-six signed-doubling
supports.  Each owner obligation is an affine four-flat in `F_2^6`; the pair
nerve is octahedral, opposite flats are disjoint, every other pair meets in
four words, and no triple meets.  Hence 16 unit words satisfy zero owners and
48 satisfy exactly two.  The same eight cycle triples from `c=5` reappear,
strengthened by three antipodal pair obstructions.  The completed tournament
again sees the signed cycle but loses the decisive affine equations and nerve.

THM-962 closes `c=7`.  Capacity forces all six effective orders to seven; a
row-product character eliminates 922 of 924 supports, and a square-sum
contradiction over `F_7` kills the quadratic-residue support and its nonsquare
coset.  The terminal opposite-pair square obstruction is production-green in
`LRCScaleSevenSquareSum.lean`; the literal-mask-to-character reduction is the
separate frozen exact certificate.

THM-963 closes `c=8`.  The hereditary-lcm bank has `215,728,128` contexts;
capacity and owner-local feasibility reduce it to the recurring 64 signed
cycles.  Their owner-obligation intersection graph is `K3,3`, so a
distance-two incompatible pair occurs in every owner triple.  A separately
implemented Python referee reproduces the divisor counts, the `66=64+2`
capacity split, the `3040/960/96` unit profile, and the pair spectrum
`6x0,6x8,3x16`.  `LRCScaleEightOwnerNerve.lean` now formalizes the terminal
quotient without `native_decide`: a 16,384-row ordinary-`decide` core proves
distance-two owner obligations disjoint, and a six-cycle pigeonhole lemma
kills every triple and the sixfold intersection.  Its public theorems audit
to the standard trio.  The native certificate remains responsible for
reducing all 215,728,128 literal contexts to that signed-cycle quotient.

THM-969 closes `c=9`.  The `482,294,736` literal contexts reduce to `1,186`
scalar rows and then `76` owner-local contexts: 64 all-order-nine signed
cycles and one 12-context mixed orbit.  Exact replay of `3,048,192` unit
fibres has zero survivors.  The faithful owner-obligation nerves are `3K2`
and, on the mixed orbit's order-nine owners, `2K2`.  Optimized, unoptimized,
and sanitizer C++ runs agree with an independent Python reconstruction.
`LRCScaleNineOwnerNerve.lean` kernel-checks both finite nerve tables and their
abstract empty-total-intersection consumers without `native_decide`; the raw
bank-to-nerve completeness reduction remains computational.

THM-970 closes `c=10`.  The two-prime hereditary grammar has `3,002,076`
divisor contexts and `821,620,800` literal contexts.  Scalar capacity leaves
`1,200` rows; owner-local feasibility leaves exactly the 64 all-order-ten
sign transversals of `F_13^*/{+-1}`; and the remaining `64*4^6=262,144`
global unit words have no survivor.  On the six projective owner classes the
zero pairs form a multiplication-by-two `C6`, so the nonempty-intersection
nerve is `K6` minus `C6`, the triangular prism, with no triple faces.  Raw
runners and completed tournaments lose the disjoint-obligation meaning of
those cycle edges.

`LRCScaleTenProjectivePrism.lean` is the in-progress kernel transcription of
this terminal quotient.  It encodes the exact order-ten mask table, proves
that full capacity forces a partition, and reduces the projective classes to
the multiplication-by-two cycle.  The native `821,620,800 -> 64` reduction
remains external; this module is not counted green here until its
resource-guarded direct check finishes.

THM-974 closes `c=11`.  The prime-order hereditary grammar has 57 divisor
words and `1,636,866,000` labelled literal contexts.  Scalar capacity and
owner-local feasibility leave 66 all-order-eleven supports.  Literal replay
of all `66,000,000` remaining unit words gives zero survivors: 64 supports
have pairwise-disjoint owner obligations, while the quadratic-residue and
nonresidue supports have intersection graph `3K2`.  An independent Python
referee proves the exact provider/unit/owner/sheet covariance before reducing
to seven multiplication orbits.

Thus the AP-centred common-scale H6 faces `c=3,...,11` are closed.  The live
pull extended this: THM-976 independently certifies complete owner
orthogonality at `c=12`; THM-860 makes `c=13` primitive-impossible; and THM-977
closes `c=14` already at the owner-local gate, where every scalar row misses at
least two of fourteen sheets.  THM-978 now closes `c=15`: its 2,184 scalar rows
have respectively 0, 1, 2, or 4 feasible owners, never all six.  The faithful
terminal datum is the labelled feasibility subset and maximum-union vector;
the induced owner tournament is always transitive and loses the absolute
coverage threshold.  THM-980 closes `c=16` at the same pre-nerve layer: its
2,540 scalar rows have 0, 1, or 2 feasible owners, so at least four owner
projections are empty.  THM-981 closes `c=17`: its complete
`22,303,024,128`-context bank reduces to the two all-order-seventeen
quadratic-residue/nonresidue rows, and every one of their twelve owner fibres
has exact maximum union `16/17`.  Independent frozen C++ and Python
reconstructions agree.  THM-982 now fully closes `c=18`: independent
literal-search C++ and algebraic-CRT Python certificates traverse the exact
`27,490,799,952`-context bank and leave no row feasible at more than four
owners.  THM-983 uniformly excludes every prime common scale `p>=19` by an
exact residue-capacity recurrence plus structural arguments at the only
threshold exceptions `p=23,29`.  THM-986 closes composite `c=20`, and THM-988
closes `c=21`; on its two all-order rows a cubic character on the induced
`F_7^*` marks proves symbolically that the scalar-tight masks cannot partition
the twenty-one sheets.  THM-989 now closes `c=22`: independent algebraic-CRT
Python and literal-search C++ certificates agree on the exact owner-local
deficit and, under common serialization, every reachable-mask bank.  Since
`c=23` is prime, THM-990 independently closes `c=24`: its two implementations
reduce 154,461,339,648 state contexts to 66,984 scalar rows, make every row
miss at least two owners, and hash all 101,961,528 reachable masks.  THM-992's
structural prime-square argument reduces `c=25` to 36
rows and proves, via `Z/25 -> Z/5` fibre incidences, that every owner misses at
least three sheets.  Its frozen primary awaits an independent replay.  Scale
26 is multiple-of-thirteen excluded; the next untreated scale is `c=27`.

`LRCScaleTwelveOwnerOrthogonality.lean` now kernel-checks THM-976's terminal
quotient: every realized mask has size two, full coverage forces a partition,
and one exact two-provider core proves owner obligations zero and four
disjoint for every sign transversal.  It does not replay the native
`2,413,458,432 -> 64` reduction, the exact owner size 48, or all fifteen
pairwise intersections.  The generic `LRCPreNerveProjection` module
formalizes the logic shared by the later deficits: a global word projects to
a local witness at every owner, so one empty projection is terminal; a
two-owner Boolean counterexample proves the converse false.

The remaining composite scales, the finite ramified H5 bank, and non-AP,
deep, and higher-sheet branches remain parts of the global `n=12` problem.
None of these facewise results proves uniform sporadic emptiness.

`LRCSporadicDiscreteCap.lean` separately kernel-checks a conditional terminal
arithmetic sharpening.  THM-668 ruler data `q<=2b`, strict `mu>1/n`, and the
THM-759 completion inequality force

```text
mu >= 1/n + 1/(2nb),
w <= 2nb^2/(2b+n+1),
```

hence `w<=24b^2/(2b+13)<12b` at `n=12`.  The module is sorry-free, uses no
`native_decide`, and audits to the standard foundational trio.  It is
deliberately conditional on the two cited bridges and is not a proof that the
global sporadic branch is empty.

The false `q<=25` proposal now has a much stronger exact guardrail.  The row

```text
(43,55,61,70,73,79,83,99,103,104,109,113,156)
```

is primitive, covering, zero-free at every `q=15,...,25`, and has a complete
signed deck at every such modulus.  It has first rational witness `2/27`, exact
`M=43/199`, diameter `113`, ratio below four, maximum common-prime packet three,
all deletion gcds one, and every deletion strictly above `1/13`; reattaching
the deleted runner privately blocks its displayed witness.  Translation by
`k*lcm(2,...,25)` preserves the entire bounded-period blocker code,
primitivity, diameter, and packet cap while its ratio tends to one and an
explicit clearance tends to `1/2`.  The strongest exact local replacement is

```text
C_q=N_q-|B_q(S)| > N_q-phi(q)/2,
```

equivalent to a witness for a zero-free covering row at `15<=q<=28`.  The
faithful object is an irredundant runner--signed-card--modulus incidence
hypergraph; two transitive runner gauges flip 47 of 78 edges and retain none of
the cross-modulus private-owner circuit.

`LRCQ25Obstruction.lean` now kernel-checks the full explicit statement: the
base row and every common translate are covering and primitive, zero-free on
`15<=q<=25`, and have no inclusive-band witness for any `2<=q<=25`; multiplier
`2` at `q=27` gives the honest LRC14 witness and `Mreach` floor.  The module is
sorry-free, avoids `native_decide`, and audits to the standard foundational
trio.

Finally, THM-948's support-three structure now includes the complete doubling
locus.  For `d=gcd(a,b)`, `A=a/d`, `B=b/d`,

```text
M({a,2a,b})=C(A mod 14,B mod 28)/(686AB),
C=0  iff  7|A or 14|B,
sign C=centeredSign_14(A)*centeredSign_28(B) when C!=0.
```

Fourier coefficients of `g(x)=1_{||x||<1/14}-1/7` and
`H(x)=g(x)g(2x)`, followed by the periodic Bernoulli identity, prove the
formula; all `392` residue cells, `5,028` endpoint triples, `267` midpoint
all-face triples, and `54` dilations pass exact referees.  This removes two
entire divisor hyperplanes from the B5 support-three tail and gives every other
doubling triple an exact sign.  The carrier is a signed support hyperedge with
its `(A mod14,B mod28)` Fourier cell, not a runner tournament.

`LRCExactDoublingTriple.lean` now proves the denominator-cleared `/48`
Bernoulli coefficient identity, audits all `392` cells once, and lifts
integrality, zero/sign, reflection, and the infinite negative family to
arbitrary residues.  It is sorry-free and avoids `native_decide`.  The analytic
Fourier/integral equality connecting concrete support mass to this closed
Bernoulli interface remains an explicit paper-level boundary.

## Concurrent mainline integration — witness rigidity, deep census, and owner survivors

THM-949 and THM-950 materially sharpen the residual while this manifest is
being assembled.  Inside `q<=14v`, a bad runner has an integer witness in
`[1,v]`; simultaneous bad runners with speed ratio at least three inherit a
factor-three witness ladder, so seven consecutive overlaps force a top witness
at least `729` times the bottom one.  For `q<=7v`, a fixed witness supports at
most one bad multiplier.  Consequently a windowed bottom runner injects the
deep-multiplier set into a finite witness range.  The unconditional Lean theorem

```text
B5 >= liveCount - 792 * #{p : bandCount(p) >= 6}
```

replaces the former need for a global pointwise cap by an exact `792:1` census
target.  The remaining composition problem is to combine this injection and
ladder pinning with the nonlacunary dense-core strata strongly enough that live
multipliers beat deep ones uniformly.

THM-951 separately closes the proof plumbing for explicit capped rows:
`CoverageCapped(6)` plus `deepSixCount<liveCount` now yields a lonely time
end-to-end, and a finite `decide` example passes through the B5 funnel without
direct witness search.  The dense-core relation endgame also consumes the
unconditional weighted census certificate above, so the unresolved content is
the live-versus-deep supply theorem, not another assembly interface.  On the
scale side, `window_tail_glue` and `norm_glue` now compose any certified base
window with a `7/3`-nested tail.  THM-955 now proves the cluster-gap width under
its necessary positive-numerator hypothesis; its original unconditional wording
and raw-tooth constant derivation were false and have been corrected.  THM-959
then proves a direct nested induction for every supplied partition into blocks
of size at most six.  Using THM-955 directly gives
`ell=(3/4,5/4,2,3,6,13)`,
`mu=(2/11,17/154,1/14,5/119,1/36,1/85)`, and worst junction
`G0(6,6)=1105`.  The earlier `2198` was valid but needlessly weakened both
THM-955 constants; the original `2030` was invalid because it reversed the
monotonicity in `M/m`.  The old claim that the complement is exactly one
seven-comparable block has been withdrawn.  What remains is the partition
trichotomy and THM-955's L2 periodic-comb enumeration/discrepancy plus L3
`arcSafe` assembly; L1 is production-green in `LRCClusterGapBrick`.

LEM-034 and LEM-035 independently refine the boundary side of the same state.
Endpoint classes adjacent to a section have forced signs, while class-zero
crossings of a `7m` owner have a co-lattice of simultaneous owners.  Genuine
survivors can therefore be reassigned by the minimum-owner convention.  On
`{1,...,6,7M}`, survivor columns reduce to multiplication permutations of
`{1,...,6}` with explicit carry and boundary corrections.  The formalization
interface must retain both the intrinsic surviving boundary and its attributed
owner; either quotient alone loses information consumed by the signed comb
ledger.

## S54 addendum — exact pair ledgers and the two finite THM-954 sockets

The strict-open pair-grid bridge is no longer an open formalization item.
`LRCRationalOpenComb` constructs the exact rational danger teeth and normalized
pair merge, identifies their real carrier, and proves the circle volume.
`LRCOpenPairLedger` then proves the exact complete-grid count with the mandatory
zero atom, the sharp component/endpoint budget, the all-pair discrepancy, and
the clean-534 consumer.  This is the same strict endpoint technology used by
the corrected positive-window cluster work, but it does not imply THM-959's
still-open partition trichotomy.

The continuous pair formula is no longer an open analytic bridge.
`LRCB5MergeLength` replaces the linear merge by the quadratic all-clips sum;
`LRCB5PairOverlapSum` evaluates the cyclic ledger by Bernoulli/Raabe;
`LRCB5DifferenceFibers` proves uniform gcd fibers; and
`LRCB5CombReindexing` closes the split-boundary tooth identity and the full
`PairOverlapReindexing` statement.  `LRCPairCovarianceReindexer` consequently
exports the unconditional pair-correlation and family `pairKernel` identities.

For the finite layer certificates, `LRCPairTopClassification` proves the top
ratio colors `{12,13}` and `{13}` and their actual `24/12` path caps on thirteen
distinct magnitudes.  `LRCAnchoredCliqueTransfer` and
`LRCPairRatioQuotient` prove the concrete absolute-speed anchor injection and
reduce every middle threshold to a finite allowed-ratio clique exclusion.  A
finite-cover theorem separates envelope enumeration from the stored graph
replays.  Tau3 is closed by its edgeless quotient.  Tau4 is closed by an exact
60-ratio replay whose active component maps to the eight-class Wagner circle;
the graph contains 18 four-cycles, so a raw `K_{2,2}` Zarankiewicz cap is false
and Wagner triangle-freeness is the faithful certificate.  The exact remaining
THM-954 producer socket is the five finite replays tau5--tau9.

This completion is deliberately local.  The `(4,4,8,8)` valuation-gap-two wall
now has an unconditional phase-independent producer and a large-frequency
escape.  `LRCSelectedWitnessGapTwoResidual` identifies the exact affine pencil

```text
4(F_a-F_b) = a4b-a4a.
```

Thus a normalized q4 gap at least `60B` fires one sharp selected-witness gate;
the honest q4488 residue is the close nonzero q4 pencil/joint zero-small
relation branch, not a missing matching-wall producer.  The q22 opposition
selector and the q244/q333 zero-or-small branches also remain, as does the
dense-core signed support-3/4/5 payment.  The incoming scale closures and frame
bounds are compatible structural inputs but do not replace those obligations.

The sparse rational-pair ledger has also crossed its last interval-to-phase
boundary.  `LRCSparseBranchLattice` proves the residue bijection on every
positive Bezout branch, injectivity from width `<q`, disjointness of the
`k=0,±1` multiplier families, and the complete finite-`q` joint-failure count.
The positive contribution is the `Int.toNat` of the floor difference.  This
truncation is essential: at the integral zero-width boundary `(i',j')=(1,13),
q=14`, the raw difference is `-1` while the actual branch is empty.  The older
raw-floor wording in S49 is superseded by this kernel-checked formula.

`LRCRealRelationLock` adds the continuous coefficient-weight lock used by all
of these scalar branches: at an arbitrary real phase, every integer speed
relation with total absolute coefficient weight at most fourteen transfers
exactly to the selected integer witnesses.  It closes the real-error-to-exact-
relation step, but it does not eliminate a nonzero small frequency by itself.

## Landed addendum — THM-985, the canonical two-circle classifier

`LRCTwoCircle` proves the forward inclusion of the integer and half-resonance
circles in the canonical `{1,...,13}` depth-six set.  The kernel-pure
`LRCTwoCircleII` now proves the full converse and exposes the exact statement

```text
bandCount canonical q p >= 6
  iff (84p < q or 84(q-p) < q) or 84|2p-q| < q,
```

for `0<p<q`.  Its middle-case theorem `compat_card_le` checks in one finite
kernel decision that every normalized coprime witness for `k₀=3,...,8` has at
most four compatible later partners, so those cases have fewer than six
failures; the `k₀=1,2` cases give the two circles and `k₀>=9` is cardinally
impossible.  Exact reconstruction is now a referee rather than the logical
basis.  THM-987 closes the remaining card computation in Lean for every
`q>=2`:

```text
#deep(q) = 2B + (B + 1 - (q+B) % 2),    B = (q-1)/84.
```

The integer and half circles are proved disjoint.  Its `q=14` corollary has
one deep and six live multipliers and kernel-checks loneliness of the tight
family through the exact census.  `LRCTwoCircleCount` now packages the same
geometry as reusable low/high/half Finsets, with membership equivalences,
all pairwise disjointness, parity and compact card formulas, the full set
identity, and `canonicalDeepMultipliers 14 = {7}`.  THM-985/987 remain canonical-family
theorems, not the uniform weighted/deep or trapped-core supplier required by
LRC(14).

## Addendum (boxeph-S69): item 15

15. **LEM-042 (pair-overlap law).** (A) exact trapezoid-sum formula: finite
    rational arithmetic per pair (decide); (B) the 1/49 integral identity:
    one ring computation; (C) consecutive-floor: finite-form induction on
    n = floor((2a'+1)/14) (decide-adjacent).  `LRCPairOverlapArcs` turns these
    into restricted-volume credits.  `LRCC8Consecutive` uses the residue-
    `3 mod 7` strict edge in every seven-edge Hunter path to prove positive
    safe measure for every consecutive block `v,...,v+7`.  This is the live
    generic shifted eight-comb theorem.  It does not add the remaining five
    combs, handle the inherited admissible-window problem, or prove LRC(14).

## Landed addendum — congruence-averaging core (kind-pasteur-S128c44)

`TournamentH7.LRCCongruenceAveraging` (58 s, foundational axioms only): the
arithmetic core of THM-952's averaging lemma — `lavSum_mul_unit` (unit rotation
preserves the reciprocal least-abs-value sum; `Units.mulLeft_bijective`),
`leastAbsVal_neg` (the fold pairing), and `lavSum_89 = 1 + 2·H₄₄` EXACTLY,
kernel-decided at the adversarial modulus of the cont.43 referee.  The exact
harmonic form is stronger than the paper's logged bound for any instance; the
general folded identity (odd n) is the named next rung.

`LRCCongruenceAveraging` extended (kind-pasteur-S128c45): **the general folded
identity `lavSum_odd : lavSum (2m+1) = 1 + 2·H_m` is kernel-green** (foundational
axioms only) — the averaging lemma's exact harmonic law at every odd modulus, via
the val-transfer (`sum_nbij'`, klein's Thm892 template) and `sum_range_reflect`.
The 89-anchor is now a corollary of the general law.  Tactic notes: `sum_nbij'`
takes PLAIN maps (no membership binder); numeral rewrites under `NeZero`
instances are motive-incorrect — reindex the range argument with `show` instead.

## Landed addendum — eight-consecutive Hunter assembly (codex-S68)

`LRCC8Consecutive` closes the previously separated `c=8` analytic interfaces.
On the length-one window, each clipped danger comb has measure at most `1/7`.
Each consecutive overlap has the exact closed-form lower credit

```text
1/49 + r(6-r)/(49k(k+1)),  r = k mod 7.
```

Across seven consecutive edges the baseline is exactly `1/7`, and at least one
edge has positive excess.  Path Hunter therefore bounds the eight-event union
strictly below the unit-window measure.  The parallel S77 module proves this
measure theorem compactly; `LRCC8ConsecutiveWitness` is only the proof-facing
consumer that extracts a point and packages a literal `Lonely 14` witness for
every positive block `base,...,base+7`.  The ordered danger-comb path is the
faithful carrier; a tournament or unordered Zarankiewicz support graph forgets
the numeric credits and Hamiltonian charging order.  This closes one concrete
block exit but does not pay the general depth-six/seven or `462/210` census
debts.

The same module now contains `c7_consecutive_good_pos`.  For seven combs the
union budget is exactly one, so any single adjacent overlap credit crosses the
wall; `consecutive_credit_closed` supplies at least `1/49`.  This closes the
consecutive c7 measure theorem, and `LRCC8ConsecutiveWitness` now extracts its
literal `Fin 7` witness alongside the existing `Fin 8` theorem.  Neither adds
the remaining runners of a thirteen-speed family.

THM-979's explicit sampling modulus is a distinct next bridge.  It can turn a
continuous positive `B5` floor into a finite grid, but its proposed modulus is
far larger than the sharp `q<=98` zero-stalk resonance range.  The two routes
must be joined through large-modulus payment or a formal sampling/T_s theorem,
not by silently applying the finite-window fork.

## Landed addendum — live-count exit and exact-zero relation routing (codex-S69)

`LRCLiveCountLonely` removes a needless consumer-side hypothesis stack.  By
definition, `liveCount(v,q)>0` supplies `0<p<q` with every residue in the
inclusive safe band.  The module derives `q>0` and nonzero speeds from that
same witness, feeds `mreach_ge_of_pairsum_band`, and obtains
`∃ t, Lonely 14 v t`.  No `CoverageCapped` or depth-six comparison is required
after liveness is known.  The unresolved theorem is therefore the live floor,
not its conversion to a witness.

`LRCSelectedWitnessRelationRouter` closes the exact-zero residues of all three
selected support-three circles.  The q333, q244, and either q4 side of q4488
yield signed unit relations on the original selected speeds.  After sorting
absolute speeds, a factor-three dense ladder forces the top relation position
to at most `lastDensePair+1`, strengthening the generic relation localization.
The module also records two guardrails: real witness locking requires a single
common phase, while the matching-wall phases are shifted; and an explicit
q4488 residue packet has one genuinely nonzero small frequency.  Thus exact
zero is routed, but small nonzero frequency and signed B5 payment remain.

THM-984 was audited at this integration point.  Its intended endpoint bound
must use `E=2*sum|v_i|`.  `LRCGridSampling` now proves the abstract interval
and separated-family grid counts in Lean.  `LRCLiveFloorSampling` now proves
the missing application/consumer wire from any separated safe-period interval
table to `liveCount >= q*mu0-error`, defines the strict explicit
`q0=ceil((error+1)/mu0)`, and composes positivity with both direct loneliness
and the capped-five census.  `LRCGridCount` supplies an independent rational
one-interval kernel.  The remaining producer is exactly the finite separated
interval-table decomposition of `safePeriod`, with its total-length and
component-count bounds.  The original statement that the measure-zero family
`(1,...,13)` has zero live points at every modulus is false—its q14 histogram
has six live multipliers and `B5=5`.  The two positive-measure computations
remain useful evidence, not a formal trapped-core closure.  Moreover thirteen
distinct positive magnitudes force `E>=182`, so `ceil(2E/mu0)>=364`; this route
cannot be inserted into the sharp `q<=98` stalk fork by modulus size alone.

`LRCPairOverlapArcs.good_interval_of_margin` and `good_measure_of_margin`
also close the Lipschitz conversion in-kernel: a phase with strict margin
`M>1/14` contains an explicit safe interval and supplies a quantitative
positive measure floor.  Consequently the hard trapped-core supplier is now
the equality/rigidity split—classify every residual with `M=1/14` and prove a
uniform or packet-wise strict margin on the rest—not another measure-to-grid
adapter.
