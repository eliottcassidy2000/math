# HYP-7161 — The endgame parameter discharge: maximal unconditional composition (death-star-2026-07-16-S33)

**Status:** CLAIMED — in progress (S33)

**Claim (verify-first: read the build, not this stub, for the final state):**
The two parameters of `lrc14_grand_assembly` are `cite : LRCUpTo13` (the sanctioned
citation node — stays named per owner directive) and `hresidual : ResidualObligation`.
S33 works the second: compose the existing sorry-free reducers
(`residualObligation_of_ledger` / `_of_B5` / the dissociated chain
`LRCDissociatedAssembly` → `LRCDissociatedRigidity`, the detuned-dispatch reduce,
kind-pasteur's HYP-7145 stratification: generic bulk + opus THM-928(A) lacunary ratio ≥ 19
+ structured sliver) into the sharpest machine-checked surface, prove new glue where a
finite gap remains, and report the exact residual honestly.

**Known at claim time:**
- census leg CLOSED: `hdistinct20_from_data` (native_decide, 6084 rows) + `hwindowW_closed 22`
  kernel-pure (LEM-024 route).
- peel/far leg reduced: `CoveringFarLonely` ⟸ top-ratio-91 (cite) + dominant peel (cite)
  + spread13 + detuned + coarse + common-residue — remaining class = ResidualObligation
  (covering ∧ gap>13 ∧ compressed ∧ distinct |speeds| ∧ max ≥ 23 ∧ no detuned ∧ no coarse
  ∧ no common residue).
- Discrete reducers exist: HasLiveRuler (pair-sum ledger), B5 > 0 (discrete Bonferroni).

**Needs evidence:** whether the proved strata (lacunary ≥ 19, rigidity, block-gluing
THM-933) jointly cover the residual class, or what exact sub-class remains.

**S33 mid-session state (checkpoint):** frontier mapped — the fleet's minimal analytic
target is `SafeMeasureFloorPrimitive` (opus S207; "prove the floor, not the rigidity",
kps S127); the 7-wall crossing has Hunter (klein) + telescoping engine (mac-mini)
proved with the pair-floor per-cell obligation open; the paper program's remaining
items are codex's K6 formula + the manifest's ten Lean items (boxeph S49/S50).
NEW THIS SESSION (→ THM-937): `LRCChainDichotomy.lean` — the split-search dichotomy
wiring kps-S19's `cite_chain_lonely` into `ResidualObligation` (never wired before):
LRC14Statement ⟸ cite + DenseCoreObligation, where the dense core adds to the
residual the explicit certificate {dense pair ratio < 3 at j; ratios ≥ 3 above j;
entry fee 2(12−k)w(k+1) < 21(k+2)w(k) failing ∀ k ≥ j}. Referee: 200k-tuple
exhaustiveness PASS, 50k planted citable tails all closed
(`chain_dichotomy_referee_deathstar_S33.out`). Lean build in progress at checkpoint
time — VERIFY the build report in the session log before consuming.

## Collaborative machine-checked discharge — codex-S23

`TournamentH7.LRCEndgameParameterDischarge` now gives the sharpest composition
available from the proved detuned reducers.  It closes two pieces of the exceptional
dispatch without adding hypotheses:

- if `nonMultCard v g = 2` and every `g`-multiple is already a `2*g`-multiple, the
  existing two-detuned lift theorem supplies the lonely time;
- if `nonMultCard v g = 3` and every exceptional coordinate has reduced denominator
  at least `8`, the existing three-detuned coarse theorem supplies the lonely time.

Consequently `lrc14_from_deep_detuned_and_dissociated_B5` reduces the complete sharp
assembly to the sanctioned `cite : LRCUpTo13` plus exactly two named suppliers:

1. `DeepExceptionalDetunedDispatch`: the nonterminating two-adic tower and the
   three-detuned branch with a reduced denominator below `8`;
2. `DissociatedB5Supply`: positive depth-five Bonferroni mass on the primitive,
   covering, gap-family, compressed, distinct-absolute-speed, nonstructured,
   nondetuned residual.

No generic branch, primitive peel, witness-attainment step, or assembly implication
remains implicit.  The subsequent live pull adds two orthogonal refinements.  THM-937
closes every residual family with a citable ratio-`3` tail, leaving its explicit
`ChainDenseCore`.  THM-935 proves that `B5` only sees exact-support `2`--`5` relations
and reduces universal algorithmic exhaustion to the support-`3,4,5` lattice tails
`T_s(H)`.  Boxeph S53 also refutes the unqualified S52 diagonal-dominance conjecture:
cross terms are negative off resonance but positive on six resonant frames, where a
mode-coherence correction is required.  These three facts sharpen the second supplier,
not the now-closed detuned dispatch glue.

**Assumption challenge.**  The quotient vertices here are divisibility levels and
detuned coordinates, not runners, arcs, or fixed circle sections.  This quotient
preserves `nonMultCard`, reduced denominators, and immediate lift termination, while it
forgets phase chronology and circle-component incidence; those data must remain in the
actual tuple `v` and in the downstream lonely-time consumer.

## Dense-core composition and verified THM-937 repair — codex-S25

The clean build found that the THM-937 checkpoint had not completed elaboration; its
four local compatibility failures are now repaired, and `lonely_or_denseCore`,
`residualObligation_of_denseCore`, and `lrc14_of_denseCore` audit with foundational
axioms only.

`TournamentH7.LRCDenseCoreEndgame` then composes the two reductions at the correct
logical level.  Inside `ResidualObligationPrimitive`, it first branches on a `d=2,3`
detuning and consumes `MultiDetunedDispatch`; on the complementary dissociated branch,
it sorts absolute speeds and consumes `lonely_or_denseCore`.  A successful chain split
is transported directly back to the original tuple.  Only the other branch calls the
new `DenseCoreDissociatedB5Supply`.  Thus no invalid implication “lonely implies
positive B5” is used.

The current narrowest checked theorem is
`lrc14_from_deep_detuned_and_denseCore_dissociated_B5`.  Besides the sanctioned
`LRCUpTo13` citation, its exact suppliers are:

1. `DeepExceptionalDetunedDispatch` (nonterminating two-adic lift or a three-detuned
   small reduced denominator);
2. positive `B5` only on the primitive, dissociated residual whose sorted absolute
   speeds satisfy `ChainDenseCore`.

THM-935 identifies the likely mathematical content of supplier 2—universal
support-`3,4,5` relation-lattice tails after the proved support-2 tail—but that
relation-mass theorem is not yet the Lean supplier.

## Sub-four detuned discharge and relation-budget capstone — S27/S29

`TournamentH7.LRCEndgameParameterDischargeFour` sharpens the exceptional triple
alphabet by exact floor arithmetic.  For every reduced detuning denominator `q ≥ 4`,
it proves `7 * badCount ≤ 2 * g`; consequently an exactly-three-detuned tuple with all
three reduced denominators at least `4` is generic and therefore lonely.  It also proves
that a nonmultiple's reduced denominator is at least `2`, so the complementary
`HasSubFourDetuningDenominator` condition is exactly the existence of a `q=2` or `q=3`
coordinate.  All statements audit with foundational axioms only.

The strengthened theorem
`lrc14_from_four_detuned_and_denseCore_dissociated_B5` composes this with the chain split.
Besides `LRCUpTo13`, its exact two suppliers are now:

1. a nonterminating two-adic pair tower or a non-generic triple containing reduced
   denominator `2` or `3`;
2. positive `B5` on the primitive dissociated `ChainDenseCore`.

`TournamentH7.LRCB5RelationBudget` and `TournamentH7.LRCB5RelationEndgame` replace the
second raw supplier by THM-935-shaped proof data.  The exact scaled equality
`B5=(q-1)·relationModel`,
the proved support-two quarter budget at `H=30`, and a strict three-quarter bound on the
support-`3,4,5` mass imply integer `B5 > 0`; the capstone
`lrc14_from_four_detuned_and_relationBudget` then gives LRC(14).  The universal
`T_s(H)` bounds and the analytic equality with concrete discrete `B5` remain genuine
mathematics, not hidden axioms.

## Exact degree-pattern normalization — S31

`TournamentH7.LRCEndgameParameterDischargeTwoThree` exhausts the remaining triple
degree arithmetic.  For `q ≥ 3`, one-coordinate bad degree is at most `g/3`, with
equality only at `q=3`; hence a non-generic triple with no `q=2` row is exactly the
uniform `(3,3,3)` pattern.  If a `q=2` row is present, nongenericity forces a distinct
companion row with `q≤8`; `(2,8,8)` shows this cutoff is sharp for degree-only bounds.

The sharpest checked theorem is now
`lrc14_from_twoThree_detuned_and_relationBudget`.  After `LRCUpTo13`, its suppliers are:

1. the nonterminating two-adic pair tower, or a triple of pattern `q=2` with a distinct
   `q≤8` companion, or `(3,3,3)`;
2. THM-935 relation-budget certificates on the primitive dissociated `ChainDenseCore`.

The denominator quotient preserves bad-neighborhood row degrees but destroys cyclic
phase intersections.  Closing these saturated patterns therefore needs a parity/two-
adic or explicit intersection argument; no unsupported Zarankiewicz `K_{2,t}`-free
assumption is introduced.

## Exact phase, topology, and relation sharpening — S32--S37

The next formalization pass separates three previously conflated residuals.

First, THM-933's generic rational-tooth machinery is now closed at the topology,
analytic-provider, and exact-density-recursion levels.  Canonical adjacent coalescing
makes every rational survivor seam-free and `BoundaryFaithfulRotation` automatic;
the exact density is `1 - sum overlap_k`, with each overlap an explicit rational clip
sum.  This does **not** yet identify an arbitrary actual speed-block safe set with the
rational recursion or instantiate its component-measure and numerical certificates.
Those concrete realization steps remain formalization work, not an open seam lemma.

Second, `TournamentH7.LRCDetunedOverlap` now proves the local-density debt principle
on the parallel-class circle: pair-overlap cardinality pays the same amount of excess
three-row bad degree, phase by phase, and the certificate produces an actual lonely
time through `ThreeDetunedInstanceClearing`.  In the primitive uniform q-three case,
`TournamentH7.LRCEndgameUniformThreePhase` proves `g=3` and the equivalences

```text
no common good branch
  <-> cyclic permutation obstruction
  <-> three full pairwise-disjoint bad rows partition the branch circle.
```

After signing to residue one modulo three, `3/14` clearance of the normalized sum
frequency suffices.  The kernel-checked tuple `p=(1,29,28)` at `u=1/7`, normalized to
`(1,-29,28)`, simultaneously has harmonic `1/11` clearance, the cyclic obstruction,
and zero sum frequency.  Hence existing structural hypotheses do not contradict the
obstruction at an arbitrary cited witness; the remaining q-three problem is to select
a different harmonic-good phase.

Third, the below-`40` Sidon collision now returns an explicit signed support-at-most-
four coefficient certificate.  THM-939 proves complementary exclusion results on the
chain-dense core: below-mass-two relations cannot top out above the dense pair, and
unit relations cannot top out four positions above it.  The newly pulled leverage
identity closes the alternating-binomial and Bonferroni certificate algebra, including
the exact equilibrium and kill threshold.  THM-940 additionally proves the concrete
finite subset expansion of `B5` and its exact rational subset-deviation identity.
What is still missing is the quantitative trapped-deviation / relation-supported
identification and strict support-`3,4,5` tail budget under those relation traps.

The sharp honest residual after these pulls is therefore:

1. the nonterminating q-two two-adic pair tower;
2. `(2,2,q)` witness selection: at a bad phase the two q-two rows are exactly
   complementary parity classes, and `(3,27,2; g=6,u=11/100)` realizes this;
3. `(2,4,4)` witness selection: every other q-two triple is closed by exact CRT,
   one-class q≤7 rigidity, or q-eight parity slicing, while
   `(2,7,9; g=4,u=11/100)` realizes the surviving obstruction; every failure is
   forced into the exact disjoint mod-two/mod-four partition of degrees
   `g/2,g/4,g/4`;
4. primitive `g=3` cyclic q-three witness selection (the scalar `3/14` gate is only a
   sufficient selector, not an automatic consequence);
5. trapped relation-budget supply: THM-940 supplies the concrete subset ledger; the
   support-`3,4,5` trapped-deviation strict three-quarter bound remains.

`lrc14_from_finalResidues_and_relationBudget` machine-checks exactly this capstone.
This is a sharpening of the interface, not a proof of LRC(14).

The stronger theorem `lrc14_from_selectedWitnessSupplies_and_relationBudget`
eliminates the phase-uniform formulation entirely.  For each surviving triple family,
`TwoTwoSelectedWitnessSupply`, `TwoFourFourSelectedWitnessSupply`, or
`UniformThreeSelectedWitnessSupply` chooses one phase satisfying both
`ThreeDetunedHarmonicGoodAt` and `HasThreeDetunedGoodBranch`; the generic theorem
`lonely14_of_three_detuned_selectedWitness` consumes that single witness.  Thus the
five precise remaining inputs are those three selection Props,
`NonterminatingPairTowerSupply`, and `DenseCoreRelationBudgetSupply`.
The three witness Props explicitly assume all ambient speeds are nonzero; the previous
omission admitted a zero quotient and made the harmonic-good conclusion false.
