# TournamentH7 — Architecture Document

A reference for the Lean 4 + Mathlib4 formalisation of tournament
theory developed by the "Parity in Tournaments" project (Eliott Cassidy
et al.).

## Module dependency graph

```
                    ┌─────────────┐
                    │ Basic.lean  │  Tournament structure
                    └──────┬──────┘
              ┌────────────┼──────────────────┐
              │            │                  │
        ┌─────┴────┐ ┌─────┴────┐    ┌────────┴─────┐
        │ Cycles   │ │   SCC    │    │   Tilings    │
        │  .lean   │ │  .lean   │    │   .lean      │
        └─────┬────┘ └─┬────────┘    └──────┬───────┘
              │        │                    │
        ┌─────┴────────┴─────┐    ┌─────────┴─────┐
        │      OCF.lean      │    │ GridReflection │
        │  (Grinberg-Stanley)│    │  StaircaseModel│
        └──────────┬─────────┘    │      .lean     │
                   │              └──┬─────────────┘
       ┌───────────┼───────────┐     │
       │           │           │     │
   ┌───┴───┐  ┌────┴────┐  ┌───┴─────┴──┐
   │ H7    │  │  Redei  │  │ Forbidden  │  arithmetic
   │ .lean │  │  .lean  │  │   .lean    │  enumeration
   └───┬───┘  └────┬────┘  └────┬───────┘
       │           │            │
       └───────────┼────────────┘
                   │
              ┌────┴─────┐
              │   H21    │  ──┐
              │   H63    │    │  (n≤7 only)
              │HSpectrum │    │
              └──────────┘    │
                              │
       ┌──────────────────────┘
       │
   ┌───┴────────────────────────┐
   │ Universal pair (H ≠ 7,21)
   │ Finite n≤7 trio (H ≠ 7,21,63)
   │ Rédei existence + parity   │
   └────────────────────────────┘

  Project-novel extensions (oracle-S2, S3, S3-bonus):

   Tilings ─→ SelfComplementary ─→ Iso ─→ IsoProperties
                                    │       │
                                    ↓       ↓
                              IsomorphismClasses, AntiAutomorphism

   SCCounts, SCFraction, SmallTournaments        — concrete counts
   HPIPIdentity                                  — THM-326 restatement
```

## Modules summary

| Module | Content | Status |
|---|---|---|
| `Basic.lean` | `Tournament n` structure | foundation |
| `RootSigns.lean` | type-A root lattice atoms; finite-walk telescoping and 3-cycle relation | **PROVED** |
| `RootPackets.lean` | open-walk endpoint boundaries; closed packets; `DirectedCycle` to zero-root packet | **PROVED** |
| `NaturalOperationDigraphs.lean` | natural-number operation shadows; product-sum divisor layer and seed-defect padding | **PROVED** |
| `Cycles.lean` | `DirectedCycle`, `isOdd` | foundation |
| `SCC.lean` | `Reaches`, `IsSCC`, `H`, Hamilton paths | foundation |
| `OCF.lean` | OCF axiom, Moon-Moser, Moon-Camion | external classical axioms |
| `Redei.lean` | Rédei 1934 (existence + parity) | axioms + corollaries |
| `RedeiFromOCF.lean` | Rédei existence/parity and `H % 2 = 1` derived from OCF | proved from OCF |
| `H7.lean`, `H21.lean` | Universal forbidden-H proofs | proved (modulo cited axioms) |
| `H63.lean` | Finite n≤7 absence of H=63 | finite verification axiom; universal claim false at n=8 |
| `Forbidden.lean` | Generic α-enumeration + binomial bounds | mixed |
| `HSpectrum.lean` | Universal pair + finite trio bundles | bundle |
| `Tilings.lean` | `HasBasePath`, `tilde`, score formula | axioms + corollaries |
| `GridReflection.lean` | `op`, `relabel`, vertex reversal | proved infrastructure |
| `StaircaseTileModel.lean` | concrete tile coordinates, THM-280 arc statement | axiom |
| `StaircaseModel.lean` | THM-330 (SC cut theorem) — **FULLY PROVED** | proved |
| `GoodCuts.lean` | good-cut buckets, exact spectrum `{0} ∪ {2,...,n-1}`, reflection invariance | **PROVED** |
| `StaircaseConnectivity.lean` | concrete tilings induce tournaments; top good-cut bucket iff SC | **PROVED** |
| `BucketBalance.lean` | abstract finite bucket half-line and unordered balance | **PROVED** |
| `StaircaseBucketTransport.lean` | concrete staircase mask families and good-cut transport checksums | **PROVED** |
| `SelfComplementary.lean` | `IsSelfFlip`, `PaleyLike`, regular ⟹ ¬SF | axiom + proved corollaries |
| `Iso.lean` | `TournamentIso`, `≅` | proved |
| `IsoProperties.lean` | iso-invariants — **PROVED IN LEAN** | proved |
| `AntiAutomorphism.lean` | THM-316 abstract framework | **PROVED** |
| `HPIPIdentity.lean` | THM-326 restatement | proved from OCF |
| `SCCounts.lean` | SC tile counts + THM-342 diagonals | axioms |
| `SCFraction.lean` | SC tiling fractions (THM-330 cor) | axioms |
| `IsomorphismClasses.lean` | `IsoClass n`, A000568 | axioms |
| `SmallTournaments.lean` | `transitiveTournament`, `threeCycle` | **PROVED** |
| `LRCDeathChain.lean` | finite LRC14 death-chain/live-depth quotient for direct cover and survival currencies | pending local build; decidable finite target |
| `Verify.lean` | Audit ledger | reports |

## Axiom hierarchy

Sorted by Lean foundation → external classical → project-novel.

### Lean foundational axioms (always present)
`propext`, `Classical.choice`, `Quot.sound`.

### External classical axioms (cited literature)

- `ocf` — Grinberg-Stanley, arXiv:2412.10572 (2024).
- `ocf_extended` — same, truncated to α₆.
- `moonMoser` — Moon, 1962.
- `moonCamion_oddSize` — Camion, 1959; Moon, 1968.
- `redei_existence`, `redei_parity` — Rédei, 1934.
- `omegaTriangleLocalises`, `oddCyclesIn_size3`, `oddCyclesIn_size4` —
  folklore on SCC partition + small-SCC enumeration.

### Project-novel axioms

#### Independence-polynomial bounds (from Ω structure)
- `alpha_subset_bound` — α_k ≠ 0 ⟹ α_1 ≥ k.
- `alpha_chain_step` — α_k ≠ 0 ⟹ α_{k-1} ≥ k.
- `alpha_binomial_bound` — α_k ≤ C(α_1, k).
- `alpha_descent` (NEW, oracle-S2) — α_k ≥ 1 ⟹ α_j ≥ C(k, j).
  *Recommended canonical replacement.*
- `alpha_pair_bound`, `alpha_triple_subset`, etc. — specialisations.

#### H = 21 case-by-case; H = 63 finite audit
- `no_alpha_10_0`, `no_alpha_8_1`, `no_alpha_6_2`, `no_alpha_4_3`.
- `H_ne_sixtythree_le_seven_axiom` records exhaustive n≤7 absence only.
  The universal H≠63 theorem is false at n=8; see
  `05-knowledge/results/h63_counterexample_audit_s8.out`.

#### Tiling/staircase axioms
- `tilde_score_sink`, `tilde_score_source`, `tilde_score_interior`
  (oracle-2026-05-11-S1).
- `thm280_arc` (THM-280, opus-2026-04-03-S27), now stated in the
  concrete `StaircaseTileModel.lean` coordinates.

#### Counting constants
- `SCcount_2..7`, `SCtilings_3..8`, `NonSCtilings_3..8`.
- `numIsoClasses_1..7` (A000568), `numSC_1..7`.

#### THM-342 diagonal formulas
- `thm342_diag0..3`.

#### H-invariance
- `H_iso_invariant` — **NOW PROVED**.

## Theorems fully proved in Lean (no project axioms)

These depend ONLY on `propext`, `Classical.choice`, `Quot.sound` plus
external/project axioms that are already in their proof OR — best case
— on nothing at all.

### Pure Lean proofs (no project axioms)

- `H_iso_invariant` — Hamiltonian path count is iso-invariant.
- `TypeA.RootWalk.rootTotal_eq_boundary` — open root walks telescope to endpoint boundary.
- `TypeA.RootPacket.rootTotal_eq_zero` — closed root packets carry zero total root.
- `Tournament.DirectedCycle.toRootPacket_rootTotal` — directed cycles induce zero-root packets.
- `NatOperation.addShadow_iff_lt` — additive one-input shadow is strict order.
- `NatOperation.mulUnitShadow_iff_dvd` — unit-allowing multiplicative shadow is divisibility on positive targets.
- `NatOperation.mulShadow_iff_dvd_and_lt` — nonunit multiplicative shadow is proper divisibility on positive sources.
- `NatOperation.twoFactor_productSum_iff` — the two-nonunit product-sum layer is the divisor equation `a*b=r+1`.
- `NatOperation.trivial_twoFactor_productSum` — `r` ones plus factors `2` and `r+2` always solve `sum=product`.
- `NatOperation.seedDefect_add_sum` — the product-minus-sum defect is exactly the unit padding needed to repair the additive fold.
- `abstract_anti_palindrome` — anti-automorphism reverses endpoint counts.
- `epStart_sum_eq_H` / `epEnd_sum_eq_H` — endpoint fibers partition H.
- `outDegree_iso` — out-degree is iso-invariant (modulo relabel).
- `isRegular_iso` — regularity is iso-invariant.
- `SC_implies_all_cuts_crossing` (THM-330 hard direction).
- `thm330_SC_iff_all_cuts_crossing` (THM-330 full iff).
- `crossesUpward_all_implies_SC` (THM-330 easy direction).
- `reaches_zero` — every vertex reaches 0 via base path.
- `reaches_descent` — descent via base path.
- `not_SC_implies_no_crossing` — derived from easy direction.
- `apex_implies_SC` — derived from THM-330.
- `goodCuts_empty_iff_all_down` — bucket 0 is exactly the all-down tiling.
- `goodCuts_nonempty_iff_exists_upward_tile` — nonempty support iff some tile is upward.
- `goodCutCount_eq_zero_iff_all_down` — cardinality form of bucket 0.
- `goodCutCount_pos_iff_exists_upward_tile` / `goodCutCount_pos_iff_not_all_down`
  — positive bucket forms.
- `two_le_goodCutCount_of_upward_tile` — one upward tile forces at least two good cuts.
- `goodCutCount_ne_one` — THM-336 Lean core: no good-cut bucket 1.
- `goodCutCount_eq_zero_or_two_le` — strengthened THM-336 dichotomy: bucket 0 or bucket ≥ 2.
- `goodCuts_empty_or_two_le_card` — set-cardinality form of the dichotomy.
- `goodCutCount_reflect` — grid reflection preserves good-cut bucket size.
- `isGoodCut_iff_exists_upward_tile_interval` — good cuts are unions of tile intervals.
- `goodCutCount_mono` — turning more tiles upward can only add good cuts.
- `goodCutCount_bucket_bounds` — the only bucket possibilities are 0 or 2..n-1.
- `goodCutCount_eq_top_iff_all_cuts_good` — top bucket iff every legal cut is good.
- `exists_crossesCut_of_mem_cutSet` — for n≥3 every legal cut is crossed by some tile.
- `goodCuts_singleUp_eq_cutInterval` — a one-tile tiling has exactly that tile's cut interval.
- `exists_goodCutCount_eq_of_allowed` — every bucket size 2..n-1 is realised.
- `goodCutCount_spectrum` — exact good-cut bucket spectrum: precisely `{0} ∪ {2,...,n-1}`.
- `goodCutCount_allUp` / `goodCutCount_complement_of_all_down` — all-up and complement-of-bottom hit the top bucket.
- `StTiling.toTournament_hasBasePath` — the tournament induced by a concrete staircase tiling has the base path.
- `StTiling.isGoodCut_iff_crossesUpward_toTournament` — concrete good cuts are exactly THM-330 crossing-upward cuts.
- `StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected` — the top good-cut bucket is exactly strong connectivity.
- `StTiling.allUp_toTournament_stronglyConnected` / `StTiling.allDown_toTournament_not_stronglyConnected`
  — explicit SC and non-SC witnesses in the concrete tiling family.
- `BucketBalance.halfLine_balance` — oriented incident half-lines split into internal and escaping half-lines.
- `BucketBalance.crossHalf_card_eq_zero_iff` — no escaping half-lines iff every move from the bucket remains in it.
- `BucketBalance.pairHalf_mem_selfHalf` — involutive moves partner internal half-lines with internal half-lines.
- `BucketBalance.pairHalf_ne_of_fixedPointFree` — fixed-point-free moves give no self-partnered internal half-lines.
- `BucketBalance.even_card_of_fixedPointFree_involutiveOn` — finite fixed-point-free involutions have even cardinality.
- `BucketBalance.selfHalf_card_even_of_involutive_fixedPointFree` — internal half-lines are even for fixed-point-free involutive move systems.
- `BucketBalance.unordered_balance_of_even_selfHalf` — even internal half-line cardinality yields the unordered bucket balance.
- `BucketBalance.unordered_balance_of_involutive_fixedPointFree` — fixed-point-free involutive move systems satisfy unordered bucket balance.
- `BucketBalance.xorMask_involutive` — Boolean cube xor by a mask is an involution.
- `BucketBalance.xorMask_fixedPointFree_of_nonzero` — Boolean cube xor by a nonzero mask has no fixed point.
- `BucketBalance.unordered_balance_boolCube_masks` — finite Boolean cube quotients and nonzero mask families satisfy unordered bucket balance.
- `BucketBalance.fiber_eq_empty_iff` — a finite quotient bucket is empty exactly when no point maps to it.
- `BucketBalance.transportHalf_eq_empty_of_source_fiber_eq_empty` /
  `BucketBalance.transportHalf_eq_empty_of_target_fiber_eq_empty` — empty quotient fibers give zero transport rows and columns.
- `StTile.equivGapPair` — concrete staircase tiles are equivalent to the finite subtype of legal coordinate pairs.
- `StTiling.nonzeroMasks` / `StTiling.singleTileMasks` / `StTiling.complementMask` — concrete staircase mask families.
- `StTiling.transport_row_balance_allNonzeroMasks` / `StTiling.transport_row_balance_singleTileMasks` /
  `StTiling.transport_row_balance_complementMask` — THM-352 specialized to concrete staircase tiling masks.
- `StTiling.goodCutBucket_eq_zero_iff_all_down` — the finite good-cut quotient bottom bucket is exactly all-down.
- `StTiling.goodCutBucket_eq_top_iff_toTournament_SC` — the finite good-cut quotient top bucket is exactly strong connectivity.
- `StTiling.goodCutBucket_image_iff` — for `n >= 3`, the finite good-cut quotient image is exactly `{0} ∪ {2,...,n-1}`.
- `StTiling.goodCutBucket_fiber_one_eq_empty` /
  `StTiling.goodCutBucket_fiber_overTop_eq_empty` — buckets `1` and `n` are certified finite quotient gaps.
- `StTiling.transport_row_balance_goodCutBucket_allNonzeroMasks` /
  `StTiling.transport_row_balance_goodCutBucket_singleTileMasks` /
  `StTiling.transport_row_balance_goodCutBucket_complementMask` — concrete good-cut quotient row checksums.
- `H_mod_two_eq_one_from_ocf` — OCF leaves the explicit Hamiltonian-path residue `H % 2 = 1`.
- `threeCycle_isRegular` — 3-cycle is regular.
- `transitive_hasBasePath` — transitive tournament has base path.
- `transitive_not_regular` — transitive (n ≥ 2) is not regular.
- `Reaches.trans` — reachability is transitive.
- `tildeArc_exactly_one` — tilde produces a valid tournament.
- `isSelfComplementary_iff_iso_op` — clean characterisation.
- `TypeA.root_self` / `TypeA.root_swap` / `TypeA.root_add_root` /
  `TypeA.walkRootSum_append_single` / `TypeA.walkRootSum_closed` /
  `TypeA.root_cycle_sum` — type-A root-sign atoms, finite-walk
  telescoping, and the directed-triangle relation as the first closed-walk
  case.

### Lean proofs with project axioms (single-axiom dependence)

| Theorem | Axiom |
|---|---|
| `regular_not_SF` | `tilde_score_sink` |
| `regular_not_SF_id` | `tilde_score_sink` |
| `paleyLike_not_SF_id` | `tilde_score_sink` |
| `H_eq_independence_poly_at_two_truncated` | `ocf` |
| `alpha_solution_H7` | `alpha_descent`, `ocf` |
| `alpha_candidates_H21` | `alpha2_pair_bound`, `alpha_descent`, `ocf` |
| `H_pos` | `redei_existence` |
| `H_ne_two`, `H_ne_even` | `redei_parity` |

## De-axiomatisation roadmap

### Easy targets (≤ 100 LOC of Lean)
1. **`tilde_score_sink/source/interior`** — counting argument
   `outNbrs = consec ⊔ nonconsec` + `Finset.card_disjoint_union`.
2. **Concrete SC tiling counts** (`SCcount_n`, `SCtilings_n`) by
   `native_decide` once a computable `IsSC` predicate exists.

### Medium targets (200-500 LOC)
3. **`thm280_arc`** (THM-280) — case split on consecutive vs
   non-consecutive pairs in `StaircaseTileModel.lean`.
4. **`redei_existence`** — classical insertion argument (~150 lines).

### Hard targets (significant Mathlib development)
5. **`ocf`** (Grinberg-Stanley) — requires deep combinatorial development.
6. **`moonMoser` / `moonCamion_oddSize`** — classical structural theorems.

## Build & verify

```bash
# One-time setup
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
export PATH="$HOME/.elan/bin:$PATH"

# Build the project (cached oleans, ~30 s)
cd 04-computation/lean/TournamentH7
lake exe cache get
lake build TournamentH7

# Audits are printed by Verify.lean during build.
```

## Status snapshot (oracle/kind-pasteur-2026-05-29 family)

- **2982 build targets** clean for `TournamentH7` after adding
  `StaircaseBucketTransport.lean` (codex-2026-05-30).
- **20+ project-novel theorems** formalised.
- **25+ audited theorems**, with a growing axiom-free core (iso-invariants,
  THM-330, THM-316 endpoint facts, concrete examples).
- **Good-cut bucket gap completed in Lean**: nonempty good-cut support is
  equivalent to an upward tile, every tiling has `goodCutCount = 0` or
  `2 ≤ goodCutCount`, and for n≥3 every allowed size `0` or `2..n-1`
  is realised.
- **Concrete staircase connectivity bridge completed**: `StTiling.toTournament`
  is now a valid base-path tournament, good cuts translate to crossing-upward
  cuts, and `goodCutCount = n - 1` iff the induced tournament is strongly
  connected.
- **Abstract bucket balance layers formalised**: THM-348 proves any finite
  bucket map and finite move set satisfies
  `|selfHalf| + |crossHalf| = |fiber| * |moves|`; THM-350 proves the
  partner-map, finite orbit-parity, and fixed-point-free involutive unordered
  layer; THM-351 specializes it to finite Boolean cubes with nonzero xor masks.
- **THM-330 FULLY PROVED** (both directions, no project axiom).
- **THM-316 abstract anti-palindrome PROVED** by the endpoint-reversal
  bijection `σ ↦ φ * σ * vertexReversal n`.
- Single SUBMISSION.md, this ARCHITECTURE.md, and supporting reflections in `07-reflections/`.
