/-
  TournamentH7.Verify — Axiom dependency audit

  Inspect after `lake build` — the `#print axioms` output should list
  exactly the axioms documented in OCF.lean, Forbidden.lean, H21.lean,
  H63.lean, and Redei.lean (plus Lean's foundational `propext`,
  `Classical.choice`, `Quot.sound`).
-/

import TournamentH7.HSpectrum
import TournamentH7.Tilings
import TournamentH7.GridReflection
import TournamentH7.StaircaseModel
import TournamentH7.SelfComplementary
import TournamentH7.AntiAutomorphism
import TournamentH7.HPIPIdentity
import TournamentH7.Iso
import TournamentH7.IsoProperties
import TournamentH7.SCCounts
import TournamentH7.SmallTournaments
import TournamentH7.ForbiddenHCounting
import TournamentH7.GoodCuts
import TournamentH7.BucketBalance

open Tournament

/-! ### Audit each individual theorem -/

theorem H_ne_seven_audit {n : ℕ} (T : Tournament n) : H T ≠ 7 := H_ne_seven T
#print axioms H_ne_seven_audit

theorem H_ne_twentyone_audit {n : ℕ} (T : Tournament n) : H T ≠ 21 := H_ne_twentyone T
#print axioms H_ne_twentyone_audit

theorem H_ne_sixtythree_le_seven_audit {n : ℕ} (hn : n ≤ 7) (T : Tournament n) :
    H T ≠ 63 := H_ne_sixtythree_le_seven hn T
#print axioms H_ne_sixtythree_le_seven_audit

theorem H_pos_audit {n : ℕ} (hn : 1 ≤ n) (T : Tournament n) : H T ≠ 0 := H_pos hn T
#print axioms H_pos_audit

theorem forbidden_pair_audit {n : ℕ} (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 := H_not_in_forbidden_pair T
#print axioms forbidden_pair_audit

theorem forbidden_trio_le_seven_audit {n : ℕ} (hn : n ≤ 7) (T : Tournament n) :
    H T ≠ 7 ∧ H T ≠ 21 ∧ H T ≠ 63 := H_not_in_forbidden_trio_le_seven hn T
#print axioms forbidden_trio_le_seven_audit

/-! ### Project-novel results — audit -/

-- (REMOVED) gridSym_iff_audit was based on the wrong tilde_eq_reversed_op
-- axiom; tile-complement and grid-reflection are different tiling
-- involutions; correct THM-280 formalisation requires a concrete
-- tile-coordinate model and is deferred.

/-- Score-formula corollary: regular tournaments are not self-flip
    (project-novel, oracle-2026-05-11-S1).  Used to prove
    Paley(p) ∉ SF for p ≡ 3 (mod 4). -/
theorem regular_not_SF_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hreg : IsRegular T)
    (v0 : Fin n) (hv0 : v0.val = 0) :
    (tilde T).outDegree v0 ≠ T.outDegree v0 :=
  regular_not_SF T hbp hn hreg v0 hv0
#print axioms regular_not_SF_audit

/-! ### THM-330 (SC Cut Theorem) - project-novel, opus-2026-05-27-S1 -/

/-- THM-330: SC iff every cut k ∈ {1, …, n-1} has a crossing-upward arc. -/
theorem thm330_audit {n : ℕ} (T : Tournament n) (hbp : HasBasePath T) :
    IsStronglyConnected T ↔ ∀ k, 1 ≤ k → k < n → CrossesUpward T k :=
  thm330_SC_iff_all_cuts_crossing T hbp
#print axioms thm330_audit

/-- THM-330 EASY direction (now FULLY PROVED in Lean — no axioms beyond foundations). -/
theorem thm330_easy_audit {n : ℕ} (T : Tournament n) (hbp : HasBasePath T)
    (h : ∀ k, 1 ≤ k → k < n → CrossesUpward T k) : IsStronglyConnected T :=
  crossesUpward_all_implies_SC T hbp h
#print axioms thm330_easy_audit

/-- Base-path descent: any vertex u reaches any v with v.val ≤ u.val. PROVED. -/
theorem reaches_descent_audit {n : ℕ}
    (T : Tournament n) (hbp : HasBasePath T) (u v : Fin n) (h : v.val ≤ u.val) :
    Reaches T u v :=
  reaches_descent T hbp u v h
#print axioms reaches_descent_audit

/-- Every vertex reaches 0 (via base path). PROVED. -/
theorem reaches_zero_audit {n : ℕ}
    (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n) (u : Fin n) :
    Reaches T u ⟨0, hn⟩ :=
  reaches_zero T hbp hn u
#print axioms reaches_zero_audit

/-- THM-333 (apex tile is SC): if the apex arc 0 → (n-1) is present, T is SC. -/
theorem apex_implies_SC_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hv0 : 0 < n) (hvn : n - 1 < n)
    (h_apex : T.arc ⟨0, hv0⟩ ⟨n - 1, hvn⟩ = true) :
    IsStronglyConnected T :=
  apex_implies_SC T hbp hn hv0 hvn h_apex
#print axioms apex_implies_SC_audit

/-! ### Self-complementary (clean isomorphism formulation) -/

/-- IsSelfComplementary ↔ T ≅ op T (via the new TournamentIso structure). -/
theorem isSelfComplementary_iff_iso_op_audit {n : ℕ} (T : Tournament n) :
    IsSelfComplementary T ↔ T ≅ op T :=
  isSelfComplementary_iff_iso_op T
#print axioms isSelfComplementary_iff_iso_op_audit

/-! ### Regular ⟹ ¬ SF chain (THM-345 candidate) -/

/-- Any regular base-path tournament is not self-flip via identity. -/
theorem regular_not_SF_id_audit {n : ℕ} (T : Tournament n)
    (hbp : HasBasePath T) (hn : 3 ≤ n) (hreg : IsRegular T) :
    ¬ IsSelfFlip_id T :=
  regular_not_SF_id T hbp hn hreg
#print axioms regular_not_SF_id_audit

/-- Any Paley-like tournament is not self-flip via identity. -/
theorem paleyLike_not_SF_audit {n : ℕ} (P : PaleyLike n) :
    ¬ IsSelfFlip_id P.T :=
  paleyLike_not_SF_id P
#print axioms paleyLike_not_SF_audit

/-! ### THM-326 (HP = IP universal identity) -/

theorem hp_ip_truncated_audit {n : ℕ} (T : Tournament n) :
    H T = 1 + 2 * alphaCount 1 T + 4 * alphaCount 2 T
            + 8 * alphaCount 3 T + 16 * alphaCount 4 T :=
  H_eq_independence_poly_at_two_truncated T
#print axioms hp_ip_truncated_audit

/-! ### THM-316 (Abstract anti-palindrome) -/

theorem abstract_anti_palindrome_audit {n : ℕ}
    (T : Tournament n) (hn : 0 < n) (φ : Equiv.Perm (Fin n))
    (hφ : IsAntiAutomorphism T φ) (v : Fin n) :
    epStart T hn v = epEnd T hn (φ v) :=
  abstract_anti_palindrome T hn φ hφ v
#print axioms abstract_anti_palindrome_audit

/-- Endpoint-start fibers partition the Hamiltonian paths. -/
theorem epStart_sum_eq_H_audit {n : ℕ} (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epStart T hn v = H T :=
  epStart_sum_eq_H T hn
#print axioms epStart_sum_eq_H_audit

/-- Endpoint-end fibers partition the Hamiltonian paths. -/
theorem epEnd_sum_eq_H_audit {n : ℕ} (T : Tournament n) (hn : 0 < n) :
    ∑ v : Fin n, epEnd T hn v = H T :=
  epEnd_sum_eq_H T hn
#print axioms epEnd_sum_eq_H_audit

/-! ### Isomorphism invariants -/

/-- Tournament isomorphism preserves vertex out-degrees (up to relabelling).
    PROVED IN LEAN (no axiom). -/
theorem outDegree_iso_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (f : TournamentIso T₁ T₂) (v : Fin n) :
    T₁.outDegree v = T₂.outDegree (f.perm v) :=
  outDegree_iso T₁ T₂ f v
#print axioms outDegree_iso_audit

/-- Tournament isomorphism preserves the regularity property. PROVED. -/
theorem isRegular_iso_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) (hreg : IsRegular T₁) :
    IsRegular T₂ :=
  isRegular_iso T₁ T₂ h hreg
#print axioms isRegular_iso_audit

/-- H is an isomorphism invariant — PROVED IN LEAN. -/
theorem H_iso_invariant_audit {n : ℕ}
    (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    H T₁ = H T₂ :=
  H_iso_invariant T₁ T₂ h
#print axioms H_iso_invariant_audit

/-- alphaCount is an isomorphism invariant. -/
theorem alphaCount_iso_invariant_audit {n : ℕ}
    (k : ℕ) (T₁ T₂ : Tournament n) (h : T₁ ≅ T₂) :
    alphaCount k T₁ = alphaCount k T₂ :=
  alphaCount_iso_invariant k T₁ T₂ h
#print axioms alphaCount_iso_invariant_audit

/-! ### Bucket-balance half-line conservation (THM-346 core) -/

/-- Ordered half-line bucket balance for any finite quotient and finite move set.
    PROVED. -/
theorem bucketBalance_halfLine_balance_audit
    {X M B : Type*} [Fintype X] [DecidableEq B]
    (q : X → B) (step : M → X → X) (moves : Finset M) (b : B) :
    (BucketBalance.selfHalf q step moves b).card +
      (BucketBalance.crossHalf q step moves b).card =
        (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.halfLine_balance q step moves b
#print axioms bucketBalance_halfLine_balance_audit

/-! ### Good-cut buckets for staircase tilings -/

/-- Good-cut bucket 0 is exactly the all-down tiling. PROVED. -/
theorem goodCuts_empty_iff_all_down_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts = ∅ ↔ ∀ t : StTile n, b t = false :=
  StTiling.goodCuts_empty_iff_all_down b
#print axioms goodCuts_empty_iff_all_down_audit

/-- Good-cut support is nonempty iff some tile is upward. PROVED. -/
theorem goodCuts_nonempty_iff_exists_upward_tile_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts.Nonempty ↔ ∃ t : StTile n, b t = true :=
  StTiling.goodCuts_nonempty_iff_exists_upward_tile b
#print axioms goodCuts_nonempty_iff_exists_upward_tile_audit

/-- Good-cut bucket 0, cardinality form. PROVED. -/
theorem goodCutCount_eq_zero_iff_all_down_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ↔ ∀ t : StTile n, b t = false :=
  StTiling.goodCutCount_eq_zero_iff_all_down b
#print axioms goodCutCount_eq_zero_iff_all_down_audit

/-- Positive good-cut count iff some tile is upward. PROVED. -/
theorem goodCutCount_pos_iff_exists_upward_tile_audit {n : ℕ} (b : StTiling n) :
    0 < b.goodCutCount ↔ ∃ t : StTile n, b t = true :=
  StTiling.goodCutCount_pos_iff_exists_upward_tile b
#print axioms goodCutCount_pos_iff_exists_upward_tile_audit

/-- Positive good-cut count iff the tiling is not all-down. PROVED. -/
theorem goodCutCount_pos_iff_not_all_down_audit {n : ℕ} (b : StTiling n) :
    0 < b.goodCutCount ↔ ¬ ∀ t : StTile n, b t = false :=
  StTiling.goodCutCount_pos_iff_not_all_down b
#print axioms goodCutCount_pos_iff_not_all_down_audit

/-- One upward tile forces at least two good cuts. PROVED. -/
theorem two_le_goodCutCount_of_upward_tile_audit {n : ℕ}
    {b : StTiling n} {t : StTile n} (ht : b t = true) :
    2 ≤ b.goodCutCount :=
  StTiling.two_le_goodCutCount_of_upward_tile ht
#print axioms two_le_goodCutCount_of_upward_tile_audit

/-- THM-336 Lean core: no tiling has exactly one good cut. PROVED. -/
theorem goodCutCount_ne_one_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount ≠ 1 :=
  StTiling.goodCutCount_ne_one b
#print axioms goodCutCount_ne_one_audit

/-- THM-336 strengthened: bucket count is 0 or at least 2. PROVED. -/
theorem goodCutCount_eq_zero_or_two_le_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ∨ 2 ≤ b.goodCutCount :=
  StTiling.goodCutCount_eq_zero_or_two_le b
#print axioms goodCutCount_eq_zero_or_two_le_audit

/-- Good-cut set form: empty or cardinality at least two. PROVED. -/
theorem goodCuts_empty_or_two_le_card_audit {n : ℕ} (b : StTiling n) :
    b.goodCuts = ∅ ∨ 2 ≤ b.goodCuts.card :=
  StTiling.goodCuts_empty_or_two_le_card b
#print axioms goodCuts_empty_or_two_le_card_audit

/-- Grid reflection preserves the good-cut bucket. PROVED. -/
theorem goodCutCount_reflect_audit {n : ℕ} (b : StTiling n) :
    b.reflect.goodCutCount = b.goodCutCount :=
  StTiling.goodCutCount_reflect b
#print axioms goodCutCount_reflect_audit

/-- Good cuts are exactly membership in an upward tile interval. PROVED. -/
theorem isGoodCut_interval_union_audit {n : ℕ} {b : StTiling n} {k : ℕ} :
    StTiling.IsGoodCut b k ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval :=
  StTiling.isGoodCut_iff_exists_upward_tile_interval
#print axioms isGoodCut_interval_union_audit

theorem mem_goodCuts_interval_union_audit {n : ℕ} {b : StTiling n} {k : ℕ} :
    k ∈ b.goodCuts ↔
      ∃ t : StTile n, b t = true ∧ k ∈ t.cutInterval :=
  StTiling.mem_goodCuts_iff_exists_upward_tile_interval
#print axioms mem_goodCuts_interval_union_audit

theorem cutInterval_subset_goodCuts_audit {n : ℕ} {b : StTiling n} {t : StTile n}
    (ht : b t = true) :
    t.cutInterval ⊆ b.goodCuts :=
  StTiling.cutInterval_subset_goodCuts_of_upward_tile ht
#print axioms cutInterval_subset_goodCuts_audit

/-- Good-cut count is monotone under turning more tiles upward. PROVED. -/
theorem goodCutCount_mono_audit {n : ℕ} {b c : StTiling n}
    (h : ∀ t : StTile n, b t = true → c t = true) :
    b.goodCutCount ≤ c.goodCutCount :=
  StTiling.goodCutCount_mono h
#print axioms goodCutCount_mono_audit

/-- The only possible buckets are 0 or at least 2, bounded above by n-1. PROVED. -/
theorem goodCutCount_bucket_bounds_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = 0 ∨
      (2 ≤ b.goodCutCount ∧ b.goodCutCount ≤ n - 1) :=
  StTiling.goodCutCount_bucket_bounds b
#print axioms goodCutCount_bucket_bounds_audit

/-- The top bucket is equivalent to every legal cut being good. PROVED. -/
theorem goodCutCount_eq_top_iff_all_cuts_good_audit {n : ℕ} (b : StTiling n) :
    b.goodCutCount = n - 1 ↔
      ∀ k, k ∈ cutSet n → StTiling.IsGoodCut b k :=
  StTiling.goodCutCount_eq_top_iff_all_cuts_good b
#print axioms goodCutCount_eq_top_iff_all_cuts_good_audit

/-- For n >= 3, every legal cut is crossed by some staircase tile. PROVED. -/
theorem exists_crossesCut_of_mem_cutSet_audit {n : ℕ} (hn : 3 ≤ n) {k : ℕ}
    (hk : k ∈ cutSet n) :
    ∃ t : StTile n, t.crossesCut k :=
  StTile.exists_crossesCut_of_mem_cutSet hn hk
#print axioms exists_crossesCut_of_mem_cutSet_audit

/-- A single-up tiling has good cuts exactly the interval crossed by that tile. PROVED. -/
theorem goodCuts_singleUp_eq_cutInterval_audit {n : ℕ} (t : StTile n) :
    (StTiling.singleUp t).goodCuts = t.cutInterval :=
  StTiling.goodCuts_singleUp_eq_cutInterval t
#print axioms goodCuts_singleUp_eq_cutInterval_audit

/-- Every allowed nonzero good-cut bucket size is realized. PROVED. -/
theorem exists_goodCutCount_eq_of_allowed_audit {n r : ℕ}
    (hn : 3 ≤ n) (hr2 : 2 ≤ r) (hrn : r ≤ n - 1) :
    ∃ b : StTiling n, b.goodCutCount = r :=
  StTiling.exists_goodCutCount_eq_of_allowed hn hr2 hrn
#print axioms exists_goodCutCount_eq_of_allowed_audit

/-- Exact spectrum of the good-cut bucket abstraction. PROVED. -/
theorem goodCutCount_spectrum_audit {n r : ℕ} (hn : 3 ≤ n) :
    (∃ b : StTiling n, b.goodCutCount = r) ↔
      r = 0 ∨ (2 ≤ r ∧ r ≤ n - 1) :=
  StTiling.goodCutCount_spectrum hn
#print axioms goodCutCount_spectrum_audit

/-- For n >= 3, the all-up tiling is in the top good-cut bucket. PROVED. -/
theorem goodCutCount_allUp_audit {n : ℕ} (hn : 3 ≤ n) :
    (StTiling.allUp n).goodCutCount = n - 1 :=
  StTiling.goodCutCount_allUp hn
#print axioms goodCutCount_allUp_audit

/-- Complementing an all-down tiling puts it in the top bucket. PROVED. -/
theorem goodCutCount_complement_of_all_down_audit {n : ℕ}
    {b : StTiling n} (hn : 3 ≤ n) (h : ∀ t : StTile n, b t = false) :
    b.complement.goodCutCount = n - 1 :=
  StTiling.goodCutCount_complement_of_all_down hn h
#print axioms goodCutCount_complement_of_all_down_audit

/-! ### Abstract bucket balances -/

/-- THM-346 Lean core: oriented half-lines split into internal and escaping
    half-lines. PROVED. -/
theorem bucket_halfLine_balance_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.selfHalf q step moves b).card +
      (BucketBalance.crossHalf q step moves b).card =
        (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.halfLine_balance q step moves b
#print axioms bucket_halfLine_balance_audit

/-- Bucket closure criterion: zero escaping half-lines iff every chosen move
    from the fiber stays in the same bucket. PROVED. -/
theorem bucket_crossHalf_card_eq_zero_iff_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.crossHalf q step moves b).card = 0 <->
      forall x, x ∈ BucketBalance.fiber q b ->
        forall u, u ∈ moves -> q (step u x) = b :=
  BucketBalance.crossHalf_card_eq_zero_iff q step moves b
#print axioms bucket_crossHalf_card_eq_zero_iff_audit

theorem bucket_crossHalf_card_le_total_audit {alpha beta move : Type}
    [Fintype alpha] [DecidableEq beta]
    (q : alpha -> beta) (step : move -> alpha -> alpha)
    (moves : Finset move) (b : beta) :
    (BucketBalance.crossHalf q step moves b).card <=
      (BucketBalance.fiber q b).card * moves.card :=
  BucketBalance.crossHalf_card_le_total q step moves b
#print axioms bucket_crossHalf_card_le_total_audit

/-! ### THM-342 (small diagonal value) -/

example : Qcount 2 1 = 1 := by
  have := thm342_diag0 1 (by omega); simp at this; exact this

/-! ### Concrete tournament examples (no axioms needed) -/

/-- The transitive tournament on 4 vertices has the base path. -/
theorem transitive_4_hasBasePath : HasBasePath (transitiveTournament 4) :=
  transitive_hasBasePath 4
#print axioms transitive_4_hasBasePath

/-- The 3-cycle tournament is regular. -/
theorem threeCycle_regular_audit : IsRegular threeCycle :=
  threeCycle_isRegular
#print axioms threeCycle_regular_audit

/-- The transitive tournament on n ≥ 2 vertices is NOT regular. -/
theorem transitive_not_regular_audit (n : ℕ) (hn : 2 ≤ n) :
    ¬ IsRegular (transitiveTournament n) :=
  transitive_not_regular n hn
#print axioms transitive_not_regular_audit

/-! ### N_min(k) = 3^k theorem -/

/-- For any tournament T with α_k ≥ 1 (k ∈ {1, 2, 3, 4}),
    H(T) ≥ 3^k.  Project-novel (oracle-S4). -/
theorem H_ge_three_pow_k_audit {n : ℕ}
    (T : Tournament n) (k : ℕ) (hk_pos : 1 ≤ k) (hk_le : k ≤ 4)
    (h : 1 ≤ alphaCount k T) :
    3 ^ k ≤ H T :=
  H_ge_three_pow_k_of_alpha_pos T k hk_pos hk_le h
#print axioms H_ge_three_pow_k_audit

/-- H(T) < 27 ⟹ no independent triple of vertex-disjoint odd cycles. -/
theorem H_lt_27_no_alpha3_audit {n : ℕ} (T : Tournament n) (hH : H T < 27) :
    alphaCount 3 T = 0 :=
  H_lt_27_no_alpha3 T hH
#print axioms H_lt_27_no_alpha3_audit

/-- H(T) < 81 ⟹ no independent quadruple of vertex-disjoint odd cycles. -/
theorem H_lt_81_no_alpha4_audit {n : ℕ} (T : Tournament n) (hH : H T < 81) :
    alphaCount 4 T = 0 :=
  H_lt_81_no_alpha4 T hH
#print axioms H_lt_81_no_alpha4_audit
