/-
  TournamentH7.StaircaseModel — THM-330 (SC Cut Theorem)

  Defines the cut-crossing predicate and proves THM-330 fully in Lean.
-/

import TournamentH7.Basic
import TournamentH7.SCC
import TournamentH7.Tilings
import Mathlib.Data.Finset.Card

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Cuts -/

/-- The upward-closed interval `{k, k+1, …, n−1} ⊂ Fin n`. -/
def upperCut (n k : ℕ) : Finset (Fin n) :=
  Finset.univ.filter (fun v => k ≤ v.val)

@[simp] lemma mem_upperCut {k : ℕ} {v : Fin n} :
    v ∈ upperCut n k ↔ k ≤ v.val := by unfold upperCut; simp

/-- The lower complementary cut `{0, 1, …, k−1}`. -/
def lowerCut (n k : ℕ) : Finset (Fin n) :=
  Finset.univ.filter (fun v => v.val < k)

@[simp] lemma mem_lowerCut {k : ℕ} {v : Fin n} :
    v ∈ lowerCut n k ↔ v.val < k := by unfold lowerCut; simp

lemma upper_union_lower (k : ℕ) :
    upperCut n k ∪ lowerCut n k = (Finset.univ : Finset (Fin n)) := by
  ext v; simp; omega

lemma upper_disjoint_lower (k : ℕ) :
    Disjoint (upperCut n k) (lowerCut n k) := by
  rw [Finset.disjoint_left]; intro v hu hl; simp at hu hl; omega

/-! ### Crossing-upward predicate -/

/-- A non-consecutive pair `(i, j)` *crosses cut k upward* if i ∈ upperCut,
    j ∈ lowerCut, the pair is non-consecutive, and the arc goes j → i. -/
def CrossesUpward (T : Tournament n) (k : ℕ) : Prop :=
  ∃ (i j : Fin n), k ≤ i.val ∧ j.val < k ∧ j.val + 2 ≤ i.val
    ∧ T.arc j i = true

/-! ### Strong connectivity -/

/-- T is strongly connected: every pair is mutually reachable. -/
def IsStronglyConnected (T : Tournament n) : Prop :=
  ∀ u v : Fin n, Reaches T u v

/-! ### Reachability composition -/

/-- Composing reaches.  Reachability is transitive. -/
lemma Reaches.trans (T : Tournament n) {a b c : Fin n}
    (r1 : Reaches T a b) (r2 : Reaches T b c) : Reaches T a c := by
  induction r1 with
  | refl _ => exact r2
  | step h _ ih => exact Reaches.step h (ih r2)

/-! ### Local lemma: no crossing-upward gives all-downward arcs -/

lemma all_arcs_downward_of_no_crossing
    (T : Tournament n) (hbp : HasBasePath T) (k : ℕ) (_hk : 1 ≤ k) (_hkn : k < n)
    (h : ¬ CrossesUpward T k) :
    ∀ j ∈ lowerCut n k, ∀ i ∈ upperCut n k, T.arc j i = false := by
  intro j hjL i hiU
  simp at hjL hiU
  by_cases hcons : j.val + 1 = i.val
  · have hi_lt : i.val < n := i.is_lt
    have hj_lt : j.val + 1 < n := by rw [hcons]; exact hi_lt
    have hbp_arc : T.arc i j = true := by
      have h_eq : i = ⟨j.val + 1, hj_lt⟩ := Fin.ext hcons.symm
      rw [h_eq]; exact hbp j hj_lt
    cases hT : T.arc j i with
    | false => rfl
    | true => exact absurd (T.asym i j ⟨hbp_arc, hT⟩) (fun x => x)
  · by_cases hT : T.arc j i = true
    · exfalso; apply h
      refine ⟨i, j, hiU, hjL, ?_, hT⟩; omega
    · cases htv : T.arc j i with
      | false => rfl
      | true => exact absurd htv hT

/-! ### Base-path descent: every vertex reaches 0 -/

/-- Auxiliary: for any `k < n`, the vertex `⟨k, _⟩` reaches `⟨0, _⟩`. -/
private lemma reaches_zero_aux (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n) :
    ∀ k (hk : k < n), Reaches T ⟨k, hk⟩ ⟨0, hn⟩ := by
  intro k
  induction k with
  | zero => intro hk; exact Reaches.refl _
  | succ k ih =>
    intro hk
    -- arc ⟨k+1, hk⟩ → ⟨k, _⟩ is the base path arc.
    have hk_lt : k < n := by omega
    have h_arc : T.arc ⟨k + 1, hk⟩ ⟨k, hk_lt⟩ = true := hbp ⟨k, hk_lt⟩ hk
    have h_rec := ih hk_lt
    exact Reaches.step h_arc h_rec

/-- Base path implies every vertex reaches 0. -/
lemma reaches_zero (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n) (u : Fin n) :
    Reaches T u ⟨0, hn⟩ := by
  have h_u : u = ⟨u.val, u.is_lt⟩ := rfl
  rw [h_u]
  exact reaches_zero_aux T hbp hn u.val u.is_lt

/-- Auxiliary: for `j ≤ k < n`, vertex `⟨k, _⟩` reaches `⟨j, _⟩`. -/
private lemma reaches_descent_aux (T : Tournament n) (hbp : HasBasePath T) :
    ∀ k (hk : k < n) (j : ℕ) (hj : j < n) (_hjk : j ≤ k),
      Reaches T ⟨k, hk⟩ ⟨j, hj⟩ := by
  intro k
  induction k with
  | zero =>
    intro hk j hj hjk
    have hj_eq : j = 0 := by omega
    subst hj_eq
    exact Reaches.refl _
  | succ k ih =>
    intro hk j hj hjk
    by_cases h_eq : j = k + 1
    · subst h_eq; exact Reaches.refl _
    · have hjk' : j ≤ k := by omega
      have hk_lt : k < n := by omega
      have h_arc : T.arc ⟨k + 1, hk⟩ ⟨k, hk_lt⟩ = true := hbp ⟨k, hk_lt⟩ hk
      have h_rec := ih hk_lt j hj hjk'
      exact Reaches.step h_arc h_rec

/-- Base path implies: any vertex `u` reaches any `v` with `v.val ≤ u.val`. -/
lemma reaches_descent (T : Tournament n) (hbp : HasBasePath T)
    (u v : Fin n) (h : v.val ≤ u.val) : Reaches T u v := by
  have h_u : u = ⟨u.val, u.is_lt⟩ := rfl
  have h_v : v = ⟨v.val, v.is_lt⟩ := rfl
  rw [h_u, h_v]
  exact reaches_descent_aux T hbp u.val u.is_lt v.val v.is_lt h

/-! ### From 0, reach every vertex (PROVED — strong induction) -/

/-- Auxiliary: vertex 0 reaches every vertex.  Strong induction on v.val. -/
private lemma zero_reaches_aux
    (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n)
    (h : ∀ k, 1 ≤ k → k < n → CrossesUpward T k) :
    ∀ k (hk : k < n), Reaches T ⟨0, hn⟩ ⟨k, hk⟩ := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro hk
    by_cases hk0 : k = 0
    · subst hk0; exact Reaches.refl _
    · -- k ≥ 1.  Use crossing-upward at cut k.
      have hk_pos : 1 ≤ k := Nat.one_le_iff_ne_zero.mpr hk0
      obtain ⟨i, j, hi_ge, hj_lt, _hgap, hT⟩ := h k hk_pos hk
      -- 0 reaches j by IH (since j.val < k).
      have h_0_to_j : Reaches T ⟨0, hn⟩ j := by
        have := ih j.val hj_lt j.is_lt
        have hj_eq : j = ⟨j.val, j.is_lt⟩ := rfl
        rw [hj_eq]; exact this
      -- j reaches i (one step via hT).
      have h_j_to_i : Reaches T j i := Reaches.step hT (Reaches.refl _)
      -- i reaches ⟨k, hk⟩ by descent (since k ≤ i.val).
      have h_i_to_k : Reaches T i ⟨k, hk⟩ :=
        reaches_descent T hbp i ⟨k, hk⟩ (by show k ≤ i.val; omega)
      exact (h_0_to_j.trans T h_j_to_i).trans T h_i_to_k

/-- From vertex 0, reach every vertex v ∈ Fin n. -/
lemma zero_reaches_all
    (T : Tournament n) (hbp : HasBasePath T) (hn : 0 < n)
    (h : ∀ k, 1 ≤ k → k < n → CrossesUpward T k)
    (v : Fin n) : Reaches T ⟨0, hn⟩ v := by
  have h_v : v = ⟨v.val, v.is_lt⟩ := rfl
  rw [h_v]
  exact zero_reaches_aux T hbp hn h v.val v.is_lt

/-! ### The easy direction (PROVED): every cut crossing ⟹ SC -/

/-- **PROVED IN LEAN.** If T has the base path and every cut has a
    crossing-upward arc, then T is strongly connected. -/
theorem crossesUpward_all_implies_SC
    (T : Tournament n) (hbp : HasBasePath T)
    (h : ∀ k, 1 ≤ k → k < n → CrossesUpward T k) :
    IsStronglyConnected T := by
  intro u v
  by_cases hn : n = 0
  · subst hn; exact absurd u.is_lt (Nat.not_lt_zero _)
  · have hn_pos : 0 < n := Nat.pos_of_ne_zero hn
    -- u → 0 → v.
    have h_u_to_0 : Reaches T u ⟨0, hn_pos⟩ := reaches_zero T hbp hn_pos u
    have h_0_to_v : Reaches T ⟨0, hn_pos⟩ v := zero_reaches_all T hbp hn_pos h v
    exact h_u_to_0.trans T h_0_to_v

/-! ### Hard direction (PROVED) -/

/-- Any reachability path from a lower-cut vertex to an upper-cut vertex
    must contain a non-consecutive crossing-upward arc. -/
private lemma exists_crossing_in_reaches
    (T : Tournament n) (hbp : HasBasePath T) (k : ℕ) (_hk : 1 ≤ k) :
    ∀ {a b : Fin n}, Reaches T a b → a.val < k → k ≤ b.val →
      CrossesUpward T k := by
  intro a b r
  induction r with
  | refl v =>
    intro h1 h2
    omega
  | @step u v w h_arc r ih =>
    intro h1 h2
    by_cases hv : v.val < k
    · exact ih hv h2
    · push Not at hv
      by_cases hcons : u.val + 1 = v.val
      · exfalso
        have hv_lt : v.val < n := v.is_lt
        have hu_succ : u.val + 1 < n := by rw [hcons]; exact hv_lt
        have h_eq : v = ⟨u.val + 1, hu_succ⟩ := Fin.ext hcons.symm
        have hbp_arc : T.arc v u = true := by
          rw [h_eq]; exact hbp u hu_succ
        exact T.asym u v ⟨h_arc, hbp_arc⟩
      · refine ⟨v, u, hv, h1, ?_, h_arc⟩
        omega

/-- **PROVED (THM-330 hard direction).**

    If T has the base path AND is strongly connected, then every cut
    `k ∈ {1, …, n − 1}` has at least one crossing-upward arc. -/
theorem SC_implies_all_cuts_crossing
    (T : Tournament n) (hbp : HasBasePath T) (h : IsStronglyConnected T) :
    ∀ k, 1 ≤ k → k < n → CrossesUpward T k := by
  intro k hk1 hkn
  -- Pick v = ⟨0, _⟩ ∈ lowerCut, u = ⟨n - 1, _⟩ ∈ upperCut.
  have h0 : 0 < n := by omega
  have hn1 : n - 1 < n := by omega
  let v : Fin n := ⟨0, h0⟩
  let u : Fin n := ⟨n - 1, hn1⟩
  have hr : Reaches T v u := h v u
  exact exists_crossing_in_reaches T hbp k hk1 hr (by simp [v]; omega) (by simp [u]; omega)

/-! ### THM-330 (full iff, PROVED) -/

/-- **THM-330**: a base-path tournament is strongly connected iff every
    cut has a crossing-upward arc. -/
theorem thm330_SC_iff_all_cuts_crossing
    (T : Tournament n) (hbp : HasBasePath T) :
    IsStronglyConnected T ↔ ∀ k, 1 ≤ k → k < n → CrossesUpward T k := by
  refine ⟨?_, crossesUpward_all_implies_SC T hbp⟩
  exact SC_implies_all_cuts_crossing T hbp

/-! ### Corollaries -/

theorem not_SC_implies_no_crossing
    (T : Tournament n) (hbp : HasBasePath T) (h : ¬ IsStronglyConnected T) :
    ∃ k, 1 ≤ k ∧ k < n ∧ ¬ CrossesUpward T k := by
  by_contra hno
  push Not at hno
  have hall : ∀ k, 1 ≤ k → k < n → CrossesUpward T k :=
    fun k hk1 hkn => hno k hk1 hkn
  exact h (crossesUpward_all_implies_SC T hbp hall)

theorem apex_implies_SC
    (T : Tournament n) (hbp : HasBasePath T) (hn : 3 ≤ n)
    (hv0 : 0 < n) (hvn : (n - 1 : ℕ) < n)
    (h_apex : T.arc ⟨0, hv0⟩ ⟨n - 1, hvn⟩ = true) :
    IsStronglyConnected T := by
  apply (thm330_SC_iff_all_cuts_crossing T hbp).mpr
  intro k hk1 hkn
  refine ⟨⟨n - 1, hvn⟩, ⟨0, hv0⟩, ?_, ?_, ?_, h_apex⟩
  · show k ≤ n - 1; omega
  · show 0 < k; omega
  · show 0 + 2 ≤ n - 1; omega

end Tournament
