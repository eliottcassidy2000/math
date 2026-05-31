# TournamentH7 — Submission-Ready Summary

Lean 4 + Mathlib4 formalisation of forbidden H-values in tournaments.

## What is formalised

| Theorem | Statement | Source |
|---|---|---|
| `Tournament.H_ne_seven` | For every tournament T (any n), H(T) ≠ 7 | THM-343 |
| `Tournament.H_ne_twentyone` | For every tournament T (any n), H(T) ≠ 21 | HYP-1753 |
| `Tournament.H_ne_sixtythree_le_seven` | For every tournament T with n ≤ 7, H(T) ≠ 63 | exhaustive n≤7 |
| `Tournament.H_pos` | For n ≥ 1, H(T) ≠ 0 (i.e., H(T) ≥ 1) | Rédei 1934 (R1) |
| `Tournament.H_ne_two` | For n ≥ 1, H(T) ≠ 2 (corollary of parity) | Rédei 1934 (R2) |
| `Tournament.H_ne_even` | For n ≥ 1, H(T) is never an even number | Rédei 1934 (R2) |
| `Tournament.H_not_in_forbidden_pair` | H(T) ∉ {7, 21} | bundle |
| `Tournament.H_not_in_forbidden_trio_le_seven` | If n≤7, H(T) ∉ {7, 21, 63} | finite bundle |
| `Tournament.regular_not_SF` | regular base-path tournament has score differing at vertex 0 | oracle-2026-05-11-S1 |
| `Tournament.regular_not_SF_id` | regular base-path tournament is *not* SF (identity) | corollary |
| `Tournament.paleyLike_not_SF_id` | Paley-like tournament is not SF | corollary |
| `Tournament.thm330_SC_iff_all_cuts_crossing` | THM-330: SC iff every cut k has a crossing-upward arc | opus-2026-05-27-S1 |
| `Tournament.apex_implies_SC` | apex tile present ⟹ SC (THM-333) | opus-2026-05-27-S1 |
| `Tournament.StTiling.goodCutCount_ne_one` | THM-336 Lean core: good-cut bucket 1 is impossible | opus-2026-05-29-S13 |
| `Tournament.StTiling.goodCutCount_reflect` | grid reflection preserves the good-cut bucket index | opus-2026-05-29-S13 |
| `Tournament.StTiling.isGoodCut_iff_exists_upward_tile_interval` | good cuts are exactly unions of upward tile intervals | opus-2026-05-29-S14 |
| `Tournament.StTiling.goodCutCount_bucket_bounds` | good-cut buckets lie in `{0} ∪ {2,...,n-1}` | opus-2026-05-29-S14 |
| `Tournament.StTiling.goodCutCount_eq_top_iff_all_cuts_good` | top bucket iff every legal cut is good | opus-2026-05-29-S14 |
| `Tournament.StTiling.goodCutCount_spectrum` | exact good-cut bucket spectrum is `{0} ∪ {2,...,n-1}` for n≥3 | codex-2026-05-30 |
| `Tournament.StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected` | top good-cut bucket iff the induced base-path tournament is strongly connected | codex-2026-05-30 |
| `Tournament.BucketBalance.halfLine_balance` | finite bucket half-lines split into internal and escaping halves | kind-pasteur-2026-05-29-S5 / codex fix |
| `Tournament.BucketBalance.unordered_balance_of_involutive_fixedPointFree` | fixed-point-free involutive move systems satisfy unordered bucket balance | codex-2026-05-30 |
| `Tournament.BucketBalance.unordered_balance_boolCube_masks` | finite Boolean cube quotients with nonzero xor masks satisfy unordered bucket balance | opus-2026-05-30-S1 |
| `NatOperation.addShadow_iff_lt` | additive operation shadow is strict order | codex-2026-05-31-S366 |
| `NatOperation.mulShadow_iff_dvd_and_lt` | nonunit multiplication shadow is proper divisibility | codex-2026-05-31-S366 |
| `NatOperation.twoFactor_productSum_iff` | two-factor product-sum witnesses are classified by `a*b=r+1` | codex-2026-05-31-S366 |
| `NatOperation.trivial_twoFactor_productSum` | universal binary product-sum witness giving `m(k)<=2k` | codex-2026-05-31-S366 |
| `Tournament.isSelfComplementary_iff_iso_op` | IsSelfComplementary ↔ T ≅ op T | (clean characterisation) |
| `Tournament.outDegree_iso` | iso preserves out-degree (up to relabel); **PROVED IN LEAN** | clean |
| `Tournament.isRegular_iso` | iso preserves regularity; **PROVED IN LEAN** | clean |
| `Tournament.H_eq_independence_poly_at_two_truncated` | THM-326: H(T) = I(Ω, 2) (truncated form) | opus-2026-05-27-S6 |
| `Tournament.abstract_anti_palindrome` | anti-automorphism ⟹ epStart(v) = epEnd(φ v); **PROVED IN LEAN** | opus-2026-05-22-S3 / opus-2026-05-29-S10 |
| `Tournament.epStart_sum_eq_H` / `Tournament.epEnd_sum_eq_H` | endpoint fibers partition Hamiltonian paths; **PROVED IN LEAN** | opus-2026-05-29-S10 |
| `Tournament.tilde_tilde` | tile-complement is an involution; PROVED IN LEAN | oracle-2026-05-11-S1 |
| `Tournament.H_ge_three_pow_k_of_alpha_pos` | N_min(k)=3^k for k≤4: α_k≥1 ⟹ H(T)≥3^k | oracle-2026-05-29-S4 |
| `Tournament.H_lt_27_no_alpha3` | H(T)<27 ⟹ α₃=0 | corollary |
| `Tournament.H_lt_81_no_alpha4` | H(T)<81 ⟹ α₄=0 | corollary |
| `Tournament.alpha_solution_H7` | arithmetic enumeration for H = 7 (1 candidate) | clean |
| `Tournament.alpha_candidates_H21` | arithmetic enumeration for H = 21 (4 candidates) | new (oracle-S2) |

## Empirical validation of all axioms

Every axiom used in this Lean project was exhaustively verified at n ≤ 6
by `04-computation/lean_axiom_validation_s7.py`. Output saved to
`05-knowledge/results/lean_axiom_validation_s7.out`.

Result: **0 violations across all 34,866 tournaments enumerated**.

Additionally, exhaustive n = 7 enumeration of all 2,097,152 tournaments
(`04-computation/audit_n7_exhaustive_s6.py`,
`05-knowledge/results/audit_n7_exhaustive_s6.out`) confirms:
  - H = 7  has 0 occurrences,
  - H = 21 has 0 occurrences,
  - H = 63 has 0 occurrences.

Correction (opus-2026-05-29-S8): H = 63 is achievable at n = 8.  The
counterexample in `04-computation/h63_counterexample_audit_s8.py` has
H(T)=63 by both DP and direct permutation enumeration, and
I(Ω(T),2)=63 with Ω(T)=K31.  Therefore no universal theorem H(T) ≠ 63
should be cited from this Lean project.

## Axiom audit

The audit at `TournamentH7/Verify.lean` reveals exactly which axioms each
theorem depends on. The axioms split into three classes:

### Class A: Cited mathematical results (from the literature)

| Axiom | Reference |
|---|---|
| `ocf` / `ocf_extended` | Grinberg–Stanley, arXiv:2412.10572, Corollary 20 (2024) |
| `moonMoser` | Moon, J.W. (1966), Canad. Math. Bull. 9, Corollary 2.1 |
| `moonCamion_oddSize` | Camion (1959), C. R. Acad. Sci. Paris 249 |
| `redei_existence` | Rédei (1934), Acta Litt. Sci. Szeged 7 |
| `redei_parity` | Rédei (1934), ibid. |

### Class B: Elementary combinatorics (could be de-axiomatised in Lean)

| Axiom | Justification |
|---|---|
| `alpha_subset_bound` | Each k-indep set has k distinct vertices |
| `alpha_chain_step` | Each k-indep set has k distinct (k-1)-indep subsets |
| `alpha_descent` | Full downward closure: α_k≥1 ⟹ α_j≥C(k,j) |
| `alpha_binomial_bound` | Distinct k-indep sets are distinct k-subsets |
| `oddCyclesIn_size3` | A 3-vertex SC tournament is a 3-cycle |
| `oddCyclesIn_size4` | The unique 4-vertex SC tournament has score (1,1,2,2) and 2 three-cycles |
| `omegaTriangleLocalises` | Cycles are contained in single SCCs (basic SCC theory) |
| `omegaCliqueLocalises` | Same, for clique on k cycles |
| `oddCyclesIn_upper` | Trivial bound by 2^|S| |
| `tilde_score_sink` | Tile-complement score formula at sink vertex 0 |
| `alphaCount_iso_invariant` | Independence-vector counts are invariant under tournament isomorphism |

### Class C: Structural axioms with computational citation (project-specific)

| Axiom | Justification |
|---|---|
| `no_alpha_10_0` | Exhaustive n ≤ 7: no tournament has α₁=10 ∧ α₂=0 |
| `no_alpha_8_1` | Exhaustive n ≤ 7: no tournament has α₁=8 ∧ α₂=1 |
| `no_alpha_6_2` | Exhaustive n ≤ 7: no tournament has α₁=6 ∧ α₂=2 |
| `no_alpha_4_3` | Exhaustive n ≤ 7: no tournament has α₁=4 ∧ α₂=3 |
| `H_ne_sixtythree_le_seven_axiom` | Exhaustive n ≤ 7: no tournament has H=63 |

### Current `Verify.lean` audit highlights

The build output from `lake build` / `lake env lean TournamentH7/Verify.lean`
matches the intended proof-modulo-axioms status:

- `H_ne_seven_audit` depends on OCF, Moon/Camion, `alpha_subset_bound`,
  and the small SCC-localisation axioms.
- `H_ne_twentyone_audit` depends on `ocf_extended`, the chain/binomial
  α-axioms, and the four `no_alpha_*` structural axioms.
- `H_ne_sixtythree_le_seven_audit` depends only on
  `H_ne_sixtythree_le_seven_axiom` (besides Lean foundations).
- The new `N_min(k)=3^k` corollaries depend only on `ocf` and
  `alpha_descent` (besides Lean foundations).
- `abstract_anti_palindrome_audit`, `epStart_sum_eq_H_audit`, and
  `epEnd_sum_eq_H_audit` now depend only on Lean foundations: endpoint
  reversal and endpoint partitioning have been formalized.
- `goodCuts_empty_iff_all_down_audit`, `goodCutCount_ne_one_audit`, and
  `goodCutCount_reflect_audit` depend only on Lean foundations: the
  good-cut bucket constraints are fully formalized in the concrete tiling
  model.
- `isGoodCut_interval_union_audit`, `goodCutCount_mono_audit`,
  `goodCutCount_bucket_bounds_audit`, and
  `goodCutCount_eq_top_iff_all_cuts_good_audit` are also axiom-free: the
  interval-union and top/bottom bucket structure now lives in Lean.
- `goodCuts_singleUp_eq_cutInterval_audit`, `exists_goodCutCount_eq_of_allowed_audit`,
  and `goodCutCount_spectrum_audit` complete the abstract good-cut spectrum:
  for n≥3 the only missing bucket is exactly 1.
- `staircase_toTournament_hasBasePath_audit`,
  `isGoodCut_iff_crossesUpward_toTournament_audit`, and
  `goodCutCount_eq_top_iff_toTournament_SC_audit` are axiom-free: the concrete
  tiling model now connects directly to THM-330, so the top good-cut bucket is
  exactly strong connectivity of the induced tournament.
- `allUp_toTournament_SC_audit` and `allDown_toTournament_not_SC_audit` are
  axiom-free explicit witnesses for both sides of the concrete staircase
  connectivity split.
- `bucket_halfLine_balance_audit` and `bucket_crossHalf_card_eq_zero_iff_audit`
  are axiom-free finite-set bookkeeping for quotient-bucket transport.
- `bucket_even_card_of_fixedPointFree_involutiveOn_audit`,
  `bucket_selfHalf_card_even_of_involutive_fixedPointFree_audit`, and
  `bucket_unordered_balance_of_involutive_fixedPointFree_audit` close the
  abstract orbit-parity layer of THM-350 without project axioms.
- `boolCube_xorMask_involutive_audit`, `boolCube_xorMask_fixedPointFree_audit`,
  and `bucket_unordered_balance_boolCube_masks_audit` close the finite
  Boolean-cube nonzero-mask specialization of THM-351.
- The iso/regularity examples are mostly axiom-free; the remaining project
  single-axiom dependencies here are exactly `tilde_score_sink` and
  `alphaCount_iso_invariant`.

## Honest status report

This formalisation is **proof-modulo-axioms**. The deep mathematical content
(OCF, Moon's theorems, Rédei's theorems) is axiomatised with citations to the
peer-reviewed literature. The chain-of-subsets reasoning is axiomatised even
though it would be straightforward to de-axiomatise once a computable
representation of `alphaCount` is provided.

The H ≠ 21 result has additional structural axioms (Class C) whose only
justification at present is exhaustive computational verification at n ≤ 7
plus partial structural reductions.  The old H ≠ 63 universal claim is
false: H = 63 occurs at n = 8, so H63.lean now records only the finite
n≤7 absence.

## Build

```
cd 04-computation/lean/TournamentH7
lake exe cache get      # downloads Mathlib oleans, ~6 GB
lake build              # builds all targets
```

After build, `lake env lean TournamentH7/Verify.lean` prints the `#print
axioms` audit confirming the closed set of axioms each theorem depends on.

## Honest assessment for the mathematical community

This is a **structurally interesting** but **not fully rigorous** formalisation:

- **Rigorous part**: the arithmetic of which α-vectors satisfy the OCF + chain
  constraints for h ∈ {7, 21, 63} is partly formalised in Lean. The arithmetic
  enumeration for H ≠ 7 reduces the problem to exactly one α-vector; for
  H ≠ 21 it reduces to four.

- **Cited part**: OCF, Moon, Camion, Rédei are accepted as axioms from the
  published literature.

- **Conjectural part**: the structural unrealisability of each surviving
  α-vector for h = 21 at large n is currently only verified computationally
  up to n = 7 and partially reduced.  For h = 63, universal unrealisability
  is refuted by the n = 8 counterexample.

A more rigorous follow-up would prove the Class C axioms structurally,
which appears feasible for the (4,3) case (via reduction to THM-343)
but requires more work for the (10,0), (8,1), (6,2) cases.
