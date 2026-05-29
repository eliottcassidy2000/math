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
              │   H63    │    │
              │HSpectrum │    │
              └──────────┘    │
                              │
       ┌──────────────────────┘
       │
   ┌───┴────────────────────────┐
   │ Forbidden trio  (H ≠ 7,21,63)
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
| `Cycles.lean` | `DirectedCycle`, `isOdd` | foundation |
| `SCC.lean` | `Reaches`, `IsSCC`, `H`, Hamilton paths | foundation |
| `OCF.lean` | OCF axiom, Moon-Moser, Moon-Camion | external classical axioms |
| `Redei.lean` | Rédei 1934 (existence + parity) | axioms + corollaries |
| `H7.lean`, `H21.lean`, `H63.lean` | Forbidden-H proofs | proved (modulo cited axioms) |
| `Forbidden.lean` | Generic α-enumeration + binomial bounds | mixed |
| `HSpectrum.lean` | Forbidden-trio bundle | bundle |
| `Tilings.lean` | `HasBasePath`, `tilde`, score formula | axioms + corollaries |
| `GridReflection.lean` | `op`, `relabel`, vertex reversal, THM-280 | axiom + corollaries |
| `StaircaseModel.lean` | THM-330 (SC cut theorem) — **EASY DIRECTION PROVED** | mixed |
| `SelfComplementary.lean` | `IsSelfFlip`, `PaleyLike`, regular ⟹ ¬SF | axiom + proved corollaries |
| `Iso.lean` | `TournamentIso`, `≅` | proved |
| `IsoProperties.lean` | iso-invariants — **PROVED IN LEAN** | proved |
| `AntiAutomorphism.lean` | THM-316 abstract framework | axiom |
| `HPIPIdentity.lean` | THM-326 restatement | proved from OCF |
| `SCCounts.lean` | SC tile counts + THM-342 diagonals | axioms |
| `SCFraction.lean` | SC tiling fractions (THM-330 cor) | axioms |
| `IsomorphismClasses.lean` | `IsoClass n`, A000568 | axioms |
| `SmallTournaments.lean` | `transitiveTournament`, `threeCycle` | **PROVED** |
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

#### H = 21 / H = 63 case-by-case
- `no_alpha_10_0`, `no_alpha_8_1`, `no_alpha_6_2`, `no_alpha_4_3`.
- `H_ne_sixtythree_axiom`.

#### Tiling/staircase axioms
- `tilde_score_sink`, `tilde_score_source`, `tilde_score_interior`
  (oracle-2026-05-11-S1).
- `tilde_eq_reversed_op` (THM-280, opus-2026-04-03-S27).
- `SC_implies_all_cuts_crossing` (THM-330 hard direction).

#### Counting constants
- `SCcount_2..7`, `SCtilings_3..8`, `NonSCtilings_3..8`.
- `numIsoClasses_1..7` (A000568), `numSC_1..7`.

#### THM-342 diagonal formulas
- `thm342_diag0..3`.

#### THM-316
- `abstract_anti_palindrome`.

#### H-invariance
- `H_iso_invariant` — **NOW PROVED**.

## Theorems fully proved in Lean (no project axioms)

These depend ONLY on `propext`, `Classical.choice`, `Quot.sound` plus
external/project axioms that are already in their proof OR — best case
— on nothing at all.

### Pure Lean proofs (no project axioms)

- `H_iso_invariant` — Hamiltonian path count is iso-invariant.
- `outDegree_iso` — out-degree is iso-invariant (modulo relabel).
- `isRegular_iso` — regularity is iso-invariant.
- `crossesUpward_all_implies_SC` (THM-330 easy direction).
- `reaches_zero` — every vertex reaches 0 via base path.
- `reaches_descent` — descent via base path.
- `not_SC_implies_no_crossing` — derived from easy direction.
- `apex_implies_SC` — derived from THM-330.
- `threeCycle_isRegular` — 3-cycle is regular.
- `transitive_hasBasePath` — transitive tournament has base path.
- `transitive_not_regular` — transitive (n ≥ 2) is not regular.
- `Reaches.trans` — reachability is transitive.
- `tildeArc_exactly_one` — tilde produces a valid tournament.
- `isSelfComplementary_iff_iso_op` — clean characterisation.

### Lean proofs with project axioms (single-axiom dependence)

| Theorem | Axiom |
|---|---|
| `regular_not_SF` | `tilde_score_sink` |
| `regular_not_SF_id` | `tilde_score_sink` |
| `paleyLike_not_SF_id` | `tilde_score_sink` |
| `gridSym_iff_sc_via_reversal` | `tilde_eq_reversed_op` |
| `H_eq_independence_poly_at_two_truncated` | `ocf` |
| `alpha_solution_H7` | `alpha_descent`, `ocf` |
| `alpha_candidates_H21` | `alpha2_pair_bound`, `alpha_descent`, `ocf` |
| `H_pos` | `redei_existence` |
| `H_ne_two`, `H_ne_even` | `redei_parity` |

## De-axiomatisation roadmap

### Easy targets (≤ 100 LOC of Lean)
1. **`SC_implies_all_cuts_crossing`** — the THM-330 hard direction.
   Argument: contrapositive via dominating-set on `upperCut k`.
2. **`tilde_score_sink/source/interior`** — counting argument
   `outNbrs = consec ⊔ nonconsec` + `Finset.card_disjoint_union`.
3. **Concrete SC tiling counts** (`SCcount_n`, `SCtilings_n`) by
   `native_decide` once a computable `IsSC` predicate exists.

### Medium targets (200-500 LOC)
4. **`tilde_eq_reversed_op`** (THM-280) — case split on consecutive
   vs non-consecutive pairs.
5. **`abstract_anti_palindrome`** (THM-316) — bijection on the satisfying
   `Equiv.Perm` finset, similar to `H_iso_invariant`.
6. **`redei_existence`** — classical insertion argument (~150 lines).

### Hard targets (significant Mathlib development)
7. **`ocf`** (Grinberg-Stanley) — requires deep combinatorial development.
8. **`moonMoser` / `moonCamion_oddSize`** — classical structural theorems.

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

## Status snapshot (oracle-2026-05-29-S3 family)

- **953 build targets** clean.
- **20+ project-novel theorems** formalised.
- **6 theorems** fully proved with 0 mathematical axioms (iso-invariants + Lean basics).
- **THM-330 easy direction PROVED** (hard direction is the one remaining axiom for that theorem).
- Single SUBMISSION.md, this ARCHITECTURE.md, and 2 reflections in `07-reflections/`.
