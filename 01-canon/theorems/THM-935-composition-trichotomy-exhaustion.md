---
id: THM-935
title: THE COMPOSITION TRICHOTOMY — exhaustion of the middle strata at the algorithmic level. (I) THE RELATION-MASS IDENTITY (PROVED + refereed to 10⁻⁸): for any 13-packet, B_m = B_m^eq + Σ_{s=2}^{m} (−1)^s E_s^{(m)} M_s, where M_s = Σ over relations h·v = 0 of exact support s of ∏_{i∈supp} c(h_i) (c(h) = sin(2πhλ)/(πh)), and E_s^{(5)} = 24/343, 24/49, −2/7, 1 for s = 2,3,4,5 — with the SUPPORT PRINCIPLE: relations among ≥ 6 speeds are INVISIBLE to BONF5 (they enter only S₆⁺), so the level-5 certificate is decided ENTIRELY by pair ratios and 3-, 4-, 5-term linear relations; (II) THE TRICHOTOMY ALGORITHM: compute B₅ exactly (depth-spectrum sweep); if > 0 → CERTIFIED; else enumerate small relations (reduced ratios ≤ 30; linear forms with coefficients ≤ 8) → STRUCTURED with explicit witness, routed to the classified strata (support 2 → dilate/lacunary = opus THM-928(A)/codex THM-933 territory; support 3–5 → linear-forms/near-AP = the covering program); else → MIDDLE. BATTERY (26 packets: tight, deep well, geometric ×2/×3, opus-30Z, exact APs, near-AP perturbations, 12 small-scale randoms at [13,182]–[80,1120], 6 large randoms): 9 CERTIFIED, 17 STRUCTURED-WITH-WITNESS, **ZERO MIDDLE** — exhaustion holds; every BONF5 failure surrenders its witness (opus-30Z surrenders exactly the 15/14 horizon escape; the one ratio-clean failure surrenders the support-3 form −5·134 − 2·690 + 2·1025 = 0); (III) PROVED FLOORS: (a) the s = 2 TAIL LEMMA: per pair, Σ_m |c(ma)c(mb)| ≤ 1/(6ab), so ratio-dissociation beyond H bounds the s=2 deviation by (24/343)·13/H — at H = 30 that is 0.0303, only 25% of the budget B₅^eq = 2052/16807 ≈ 0.12209; (b) the PIGEONHOLE FLOOR: max(V) < 40 forces a Sidon collision (78 pairwise sums in [3,77]) = a ±1-coefficient support-≤4 witness — tiny scales can NEVER be middle; (c) ERRATUM FIX to THM-930: B₅^eq = 2052/16807 (not 0.0821, which was the certified packet's value); the depth-13 kill threshold is 2052/(16807·792) ≈ 1.54×10⁻⁴ and the deep well is 71× over (not 106×); (IV) OBSERVATIONS: the certification TRANSITION ZONE localizes at scale ≈ [80, 1120] (barely-positive/negative certificates ±0.0003–0.02 there — the natural boundary for the block-gluing composition); geometric ratio-3 chains certify at plain level 5 while ratio-2 fails (the pure-BONF5 lacunary boundary sits in (2,3), complementing codex's cascade threshold 7)
status: (I) PROVED (grouping the singular series by exact support; E_s from the binomial partial-sum identity) + refereed exactly at n = 3 (errors 7.5e-9, 9.4e-9 = box truncation only); (II) battery-complete, zero middle; (III)(a),(b) PROVED. UPDATED by THM-2051: the coarse UNIVERSAL alternative is now proved at the explicit horizon H=2^20 by whole-product Fejer--BV approximation, so the T_s(H), s=3,4,5, absolute lattice tails are no longer required for that dichotomy. They remain open if one wants termwise mass bounds or a substantially smaller horizon. The structured small-relation branch remains open.
source: kind-pasteur-2026-07-16-S128 (cont.36; owner: work the middle strata composition and prove exhaustion)
depends_on:
  - THM-930 (leverage identity/depth spectrum — the sweep engine + the erratum fixed here)
  - THM-934 (generic certification — the bulk this composes with)
related:
  - THM-932/933 (opus/codex block-gluing: the SCALE-GAP composition axis; this file is the RELATION-CONTENT axis — together the two coordinates of the exhaustion program)
  - THM-928(A) (lacunary cascade — where support-2 witnesses route)
  - THM-897/926/929 (the wall, the blockers, the anti-Lee-Yang extremizer — the structured sliver's certificates)
  - THM-2051 (proved whole-product universal exhaustion at height 2^20)
script: 04-computation/composition_trichotomy_kps_S128c36.py -> 05-knowledge/results/composition_trichotomy_kps_S128c36.out
---

# THM-935 — the composition trichotomy

> **2026-07-21 closure of the coarse middle (THM-2051).** If no exact
> support-two-through-five relation has all coefficient magnitudes at most
> `2^20`, then continuous `BONF5>0`. Thus the universal
> `small relation OR positive BONF5` alternative sought here is proved, with a
> deliberately large explicit horizon. The absolute `T_s(H)` program below is
> still mathematically open and would give finer information; the remaining
> LRC obstruction is classification/descent on the relation-rich branch.

## (I) The identity

Group the classical singular-series expansion of μ(∩_A B_v) by the exact support S of the
relation h (h_i ≠ 0 exactly on S): summing the alternating binomial series over the 13−s
free slots gives

> **B_m = B_m^eq + Σ_{s=2}^{m} (−1)^s E_s^{(m)} M_s,  E_s^{(m)} = Σ_{j≤m−s} (−1)^j C(13−s,j)(2λ)^j.**

At m = 5: E = (24/343, 24/49, −2/7, 1). Corollary (support principle): B₅ never sees
relations of support ≥ 6. Corollary (budget): B₅ ≤ 0 requires Σ_s |E_s||M_s| ≥ 2052/16807.

## (II) The algorithm and the battery

CERTIFY (exact sweep) | SURRENDER (small-relation witness, box H₂ = 30 / H₃₄ = 8) | MIDDLE.
26-packet battery: 9 | 17 | **0**. Witnesses match the known strata exactly (the 30Z packet
surrenders (420,450) = 15/14 — THM-926's horizon escape, found automatically).

## (III) Proved floors

- s=2 tail: Σ_m |c(ma)c(mb)| ≤ Σ_m 1/(π²m²ab) = 1/(6ab); 78 pairs, ratio-dissociated
  beyond H ⟹ s=2 deviation ≤ (24/343)(13/H). H = 30 spends 25% of budget.
- Pigeonhole: max(V) < 40 ⟹ Sidon collision ⟹ ±1-coefficient support-≤4 witness.
- B₅^eq = 2052/16807 exactly (THM-930 erratum: threshold 1.54×10⁻⁴; deep well 71× over).

## Named next
- T_s(H) for s = 3,4,5: uniform tails for rank-(s−1) relation-lattice sums with first
  minimum > H — closes the gap from algorithmic to universal exhaustion (THM-920-style).
- The transition zone [80,1120]: map it; feed the block-gluing (THM-932/933) composition
  through it (their scale-gap axis × this relation axis = the full 2D program).
- Ratio-2 vs ratio-3 pure-BONF5 boundary (the plain-certificate lacunary constant ∈ (2,3)).

## Formal incidence and budget guardrails — codex-S26/S28

`TournamentH7.LRCZarankiewiczGuardrail` formalizes the honest combinatorial role of
the parallel-class/Zarankiewicz quotient.  For a family of relation supports, the sum
of `choose(|S|,2)` is exactly the sum of owner multiplicities over the `78` unordered
runner pairs.  Hence an owner-multiplicity cap `m` gives total pair load at most
`78m`; in the pair-unique case, support-`≥3`, support-`≥4`, and support-`≥5` families
have cardinality at most `26`, `13`, and `7`.  More than `26` support-`≥3` relations
produces an explicit pair shared by two distinct supports.  The same module proves the
tiny-scale floor: thirteen injective positive speeds below `40` give two distinct
two-index sets with equal speed sum.  The proof uses the sharp range `[3,77]` (the
former `[3,79]` wording was a typo).  This collision is now converted to an explicit
nonzero coefficient vector in `{-1,0,1}`, with support at most four and zero dot
product.  A signed version applies to integer tuples with distinct absolute values,
so the formal result carries the relation coefficients rather than only the two
equal-pair supports.

This quotient preserves exact-support size, pair-owner incidence, multiplicity load,
and the collision graph.  It destroys relation coefficients, Fourier/Bernoulli signs,
coefficient heights, and phase chronology.  This loss is essential: the S20 divisor-
sheet bank contains identical Zarankiewicz/parity data with opposite signed
contributions, so crossing values cannot by themselves prove `B5 > 0`.

`TournamentH7.LRCB5RelationBudget` therefore keeps the signed masses separate.  It
checks the exact algebraic consumer

```text
B5model = 2052/16807 + (24/343)M2 - (24/49)M3 - (2/7)M4 - M5
```

and proves that absolute relation debt below `2052/16807` forces positivity.  At
`H=30`, the proved pair tail `(24/343)(13/30)` is strictly less than one quarter of
the equilibrium budget, leaving a machine-checked strict three-quarter socket for the
still-open `T_s(H)`, `s=3,4,5`, bounds.  The module does not claim that the analytic
singular-series relation-mass identity is already formalized.

THM-940's `TournamentH7.LRCB5SubsetExpansion` now closes the exact **discrete**
identity that was formerly missing: concrete integer `B5` is the alternating
support-grouped sum of joint-failure counts and, over `ℚ`, is exactly
`(q-1)·2052/16807` plus the signed subset-deviation ledger.  Thus the remaining
bridge is narrower: identify/bound the trapped relation-supported contribution to
those concrete deviations (especially support `3,4,5`), rather than derive the
finite subset expansion or equilibrium constant.

`TournamentH7.LRCLeverageIdentity` independently closes the finite depth-spectrum
algebra: the exact alternating-binomial leverage identity, two-sided Bonferroni
inequalities, the positive odd-truncation certificate theorem, the coefficients
`E_2,E_3,E_4`, equilibrium `2052/16807`, leverage `792`, and exact kill threshold
`57/369754` are kernel-checked.  This removes the abstract Bonferroni algebra from the
residual, but not the packet sweep encoding or the singular-series relation-mass
identification consumed by `LRCB5RelationBudget`.

THM-939's `TournamentH7.LRCDenseCoreRelationTrap` further narrows where the remaining
tails must be proved: above the dense pair there are no top-supported below-mass-two
relations, and unit-coefficient relations are trapped in the bottom four positions.
`TrappedDenseCoreB5Supply` exports exactly these theorem hypotheses to the B5 supplier.

`TournamentH7.LRCB5RelationEndgame` is the proof-producing consumer of that missing
bridge.  Its older `B5RelationBudgetCertificate` contains a modulus, four nominal
real mass coordinates, the scaled equality `B5=(q-1)·relationModel`, the sharp
one-sided pair bound `mass2 ≥ -13/30`, and the signed harmful-tail bound
`harmfulHigherContribution(mass3,mass4,mass5) < 7712/84035`.  Those coordinates
are not constrained to equal concrete relation masses: the formal semantic audit
`nonempty_b5RelationBudgetCertificate_iff` proves that inhabiting the structure is
equivalent to already having some positive concrete `B5`.  It is therefore an
obsolete abstract compatibility interface, not a realization of the missing bridge.
The `(q-1)` factor is mandatory: the former unscaled interface could not contain an
integer `B5` under its strict subunit debt budget.  The corrected certificate proves
integer `B5 > 0`, turns certificate supply on the primitive
dissociated `ChainDenseCore` into `DenseCoreDissociatedB5Supply`, and machine-checks
`lrc14_from_finalResidues_and_relationBudget`.  Exact **observed-row** arithmetic
(not the coarser raw `badCount` ledger) has reduced the triple dispatch to
`(2,2,q)`, `(2,4,4)`, or the uniform `(3,3,3)` pattern.  The q-two/q-three CRT
intersection is exactly `g/6`; q-eight costs only `g/8` inside the complementary
parity sheet; hence every other q-two triple has a common branch.  Explicit checked
phases show that both surviving q-two patterns can genuinely cover the branch circle,
so their remaining task is harmonic-witness selection, not another phase-uniform
cardinality bound.  A failing `(2,4,4)` phase is further forced to be the exact disjoint
mod-two/mod-four partition of degrees `g/2,g/4,g/4`.  The qualifier “observed-row” is
necessary: raw `badCount` arithmetic falsely retains `(2,4,7)` (at `g=56`, coarse count
`28+14+16=58`).  The q-three branch is primitive modulus `3`, and failure is exactly a
cyclic matching, equivalently a saturated pairwise-disjoint partition of the three
parallel-class rows.  A normalized sum-frequency clearance of `3/14` is sufficient,
but the checked `(1,29,28)` example at `u=1/7` shows that an arbitrary harmonic-good
witness can realize the obstruction with zero normalized sum frequency.  Thus, after
`LRCUpTo13`, the remaining mathematics is exposed exactly as the nonterminating
two-adic pair tower, witness selection for `(2,2,q)`, `(2,4,4)`, and uniform q-three,
and the trapped subset-deviation bound supplying relation-budget certificates; no
sign, cast, finite subset-expansion, or endgame composition step is left informal.

The sharper capstone `lrc14_from_selectedWitnessSupplies_and_relationBudget` removes
the artificial phase-uniform clearing interface.  Its three triple suppliers select
one phase that is simultaneously harmonic-good and branch-good, which
`lonely14_of_three_detuned_selectedWitness` turns directly into an LRC witness.
Together with `NonterminatingPairTowerSupply` and `DenseCoreRelationBudgetSupply`,
the named Props `TwoTwoSelectedWitnessSupply`, `TwoFourFourSelectedWitnessSupply`,
and `UniformThreeSelectedWitnessSupply` are now the exact unchecked inputs.
Each selected-witness Prop includes the ambient nonzero-speed hypothesis; without it a
zero harmonic quotient makes `ThreeDetunedHarmonicGoodAt` impossible.

## S45--S46 exact normalization

The former absolute quarter/three-quarter certificate is no longer the sharpest
consumer.  `TournamentH7.LRCB5NormalizedBridge` starts directly from THM-940,
normalizes each support aggregate by `q-1`, and performs the triangular Möbius
inversion in the kernel.  Singleton deviations need only be nonpositive and hence
help.  The exact residual inequalities are

```text
normalizedMass2 >= -13/30,
(24/49)M3 + (2/7)M4 + M5 < 7712/84035.
```

The second inequality is equivalent to a depth-spectrum statement.  For
`P(d)=C(d,5)-C(d,4)+C(d,3)-(319/343)C(d,2)`, the required average is
`(q-1)^(-1) sum_p P(bandCount(p)) < -65218/84035`.  Lean checks `P(d)<=0` for
`d<=6` and `P(d)>0` for `7<=d<=13`; the mathematical crux is therefore exactly
control of depth-seven-and-higher multipliers with credit retained from lower
depths.  `LRCB5CleanModulus` supplies explicit `q_H=1 mod 14`, coprime to every
speed and larger than `H sum|v|`, so all height-H modular relations are genuine
integer relations subject to THM-939.  `NormalizedB5RelationBudgetCertificate`
now contains only H and the two inequalities above.

On the detuned side, the first two-adic sheet satisfies the exact wall-count
identity.  One crossing is `(2,4,4)`; two or more crossings reduce to a dynamic
harmonic-safe-component selector avoiding the complete dyadic prefix codes.
The three selected-witness failures are exactly parity opposition for `(2,2,q)`
and complete disjoint partitions for `(2,4,4)` and `(3,3,3)`.  They attain the
Zarankiewicz equality `z(3,g;2,1)=g`, so no static incidence improvement can
select the phase.  The final obstruction-normalized capstone retains these
dynamic selectors and the depth budget as explicit hypotheses; LRC(14) is not
yet proved.
