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
walls that have not yet been reduced to one common affine wall.

On the B5 side, the normalized pair floor has been translated exactly into
coverage language.  At `cleanModulus(v,height)`, the pair input is

```text
443/1470 <= (liveCount + shiftedPairDepthMoment)/(q-1),
```

and the other input is the signed harmful-depth budget.  `CoverageCapped(6)`
offers a second exact consumer through the live-versus-depth-six census.
`LRCArcWire` identifies this cap with absence of seven simultaneous bad arcs;
seven bad rows force all 21 pairwise integer near-proportionality constraints.
The cap is not uniform in `q`: when `q > 14 max |v_i|`, the point `p=1` is the
near-zero catastrophe and all runners can be bad.  The remaining supplier is
therefore a modulus-window plus finite rigidity problem on the trapped dense
core, not a global moment inequality.

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

The logarithm is necessary.  This leaves `T3` open at a punctured near-pole
congruence sum and leaves `T4/T5` open at affine resonance strips and slabs.
Even conditional tail estimates yield only `small support relation OR B5>0`;
the structured relation branch is still a separate LRC obligation.  The natural
quotient here has slice cosets and pole windows as vertices and retains affine
offsets; a runner tournament destroys those offsets, so the honest higher object
is an oriented-matroid resonance arrangement unless a preservation theorem for
an elimination-order tournament is found.

## Addendum (boxeph-S65): item 13

13. **LEM-034/035/036/037 (the class-0 geometry line).** All decide-shaped:
    the sign law is a 2-line section argument; the survivor law's r = 1 case
    is the multiplication permutation on Z_7 (mathlib-adjacent); the
    full-column theorem is finite Z_7^6 combinatorics; the antisymmetry law
    is one reflection bijection. Each cluster instance is a finite exact
    check (Fractions arithmetic -> rational decide).
