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
