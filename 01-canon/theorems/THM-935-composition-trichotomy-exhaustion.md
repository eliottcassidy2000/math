---
id: THM-935
title: THE COMPOSITION TRICHOTOMY — exhaustion of the middle strata at the algorithmic level. (I) THE RELATION-MASS IDENTITY (PROVED + refereed to 10⁻⁸): for any 13-packet, B_m = B_m^eq + Σ_{s=2}^{m} (−1)^s E_s^{(m)} M_s, where M_s = Σ over relations h·v = 0 of exact support s of ∏_{i∈supp} c(h_i) (c(h) = sin(2πhλ)/(πh)), and E_s^{(5)} = 24/343, 24/49, −2/7, 1 for s = 2,3,4,5 — with the SUPPORT PRINCIPLE: relations among ≥ 6 speeds are INVISIBLE to BONF5 (they enter only S₆⁺), so the level-5 certificate is decided ENTIRELY by pair ratios and 3-, 4-, 5-term linear relations; (II) THE TRICHOTOMY ALGORITHM: compute B₅ exactly (depth-spectrum sweep); if > 0 → CERTIFIED; else enumerate small relations (reduced ratios ≤ 30; linear forms with coefficients ≤ 8) → STRUCTURED with explicit witness, routed to the classified strata (support 2 → dilate/lacunary = opus THM-928(A)/codex THM-933 territory; support 3–5 → linear-forms/near-AP = the covering program); else → MIDDLE. BATTERY (26 packets: tight, deep well, geometric ×2/×3, opus-30Z, exact APs, near-AP perturbations, 12 small-scale randoms at [13,182]–[80,1120], 6 large randoms): 9 CERTIFIED, 17 STRUCTURED-WITH-WITNESS, **ZERO MIDDLE** — exhaustion holds; every BONF5 failure surrenders its witness (opus-30Z surrenders exactly the 15/14 horizon escape; the one ratio-clean failure surrenders the support-3 form −5·134 − 2·690 + 2·1025 = 0); (III) PROVED FLOORS: (a) the s = 2 TAIL LEMMA: per pair, Σ_m |c(ma)c(mb)| ≤ 1/(6ab), so ratio-dissociation beyond H bounds the s=2 deviation by (24/343)·13/H — at H = 30 that is 0.0303, only 25% of the budget B₅^eq = 2052/16807 ≈ 0.12209; (b) the PIGEONHOLE FLOOR: max(V) < 40 forces a Sidon collision (78 pairwise sums in [3,77]) = a ±1-coefficient support-≤4 witness — tiny scales can NEVER be middle; (c) ERRATUM FIX to THM-930: B₅^eq = 2052/16807 (not 0.0821, which was the certified packet's value); the depth-13 kill threshold is 2052/(16807·792) ≈ 1.54×10⁻⁴ and the deep well is 71× over (not 106×); (IV) OBSERVATIONS: the certification TRANSITION ZONE localizes at scale ≈ [80, 1120] (barely-positive/negative certificates ±0.0003–0.02 there — the natural boundary for the block-gluing composition); geometric ratio-3 chains certify at plain level 5 while ratio-2 fails (the pure-BONF5 lacunary boundary sits in (2,3), complementing codex's cascade threshold 7)
status: (I) PROVED (grouping the singular series by exact support; E_s from the binomial partial-sum identity) + refereed exactly at n = 3 (errors 7.5e-9, 9.4e-9 = box truncation only); (II) battery-complete, zero middle; (III)(a),(b) PROVED (three lines each); the UNIVERSAL exhaustion theorem is reduced to the named T_s(H) lattice-sum tail lemma for s = 3,4,5 (the fleet's THM-920-style certified-sum machinery — one focused session)
source: kind-pasteur-2026-07-16-S128 (cont.36; owner: work the middle strata composition and prove exhaustion)
depends_on:
  - THM-930 (leverage identity/depth spectrum — the sweep engine + the erratum fixed here)
  - THM-934 (generic certification — the bulk this composes with)
related:
  - THM-932/933 (opus/codex block-gluing: the SCALE-GAP composition axis; this file is the RELATION-CONTENT axis — together the two coordinates of the exhaustion program)
  - THM-928(A) (lacunary cascade — where support-2 witnesses route)
  - THM-897/926/929 (the wall, the blockers, the anti-Lee-Yang extremizer — the structured sliver's certificates)
script: 04-computation/composition_trichotomy_kps_S128c36.py -> 05-knowledge/results/composition_trichotomy_kps_S128c36.out
---

# THM-935 — the composition trichotomy

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
singular-series identity with concrete discrete `B5` is already formalized.

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
bridge.  A `B5RelationBudgetCertificate` contains a modulus, four exact-support masses,
the equality of the signed model with the concrete integer `B5`, and the quarter / three-
quarter bounds.  It proves integer `B5 > 0`, turns certificate supply on the primitive
dissociated `ChainDenseCore` into `DenseCoreDissociatedB5Supply`, and machine-checks
`lrc14_from_twoThree_detuned_and_relationBudget`.  Exact degree arithmetic has reduced
the triple dispatch to a `q=2` row with a distinct `q≤8` companion or the uniform
`(3,3,3)` pattern.  The latter is now primitive modulus `3`, and failure is exactly a
cyclic matching, equivalently a saturated pairwise-disjoint partition of the three
parallel-class rows.  A normalized sum-frequency clearance of `3/14` is sufficient,
but the checked `(1,29,28)` example at `u=1/7` shows that an arbitrary harmonic-good
witness can realize the obstruction with zero normalized sum frequency.  Thus, after
`LRCUpTo13`, the remaining mathematics is exposed exactly as q-two phase/overlap or
two-adic dispatch, q-three witness selection, and construction of the trapped
relation-budget certificates; no sign, cast, or endgame composition step is left
informal.
