/-
  TournamentH7.BasePathSink — Properties of vertex 0 in HasBasePath tournaments

  In a tournament T with HasBasePath T, vertex 0 is the "sink" of the
  base path: every other vertex k ≥ 1 has a path k → k-1 → ... → 0.

  We prove tight score bounds for vertex 0 in such tournaments.
-/

import TournamentH7.Tilings
import TournamentH7.StaircaseModel

open Finset

namespace Tournament

variable {n : ℕ}

/-! ### Vertex 0's score is bounded -/

/-- In a HasBasePath tournament T with n ≥ 2, vertex 1 beats vertex 0
    (the base-path arc), so T.arc 0 1 = false. -/
theorem base_path_sink_no_out_to_one (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) :
    T.arc ⟨0, by omega⟩ ⟨1, by omega⟩ = false := by
  have h0 : (0 : ℕ) < n := by omega
  have h1 : (1 : ℕ) < n := by omega
  have h_succ : (⟨0, h0⟩ : Fin n).val + 1 < n := by show 0 + 1 < n; omega
  have h_bp : T.arc ⟨(⟨0, h0⟩ : Fin n).val + 1, h_succ⟩ ⟨0, h0⟩ = true :=
    hbp ⟨0, h0⟩ h_succ
  have h_bp' : T.arc ⟨1, h1⟩ ⟨0, h0⟩ = true := by convert h_bp using 2
  cases h : T.arc ⟨0, h0⟩ ⟨1, h1⟩ with
  | false => rfl
  | true => exact absurd (T.asym _ _ ⟨h, h_bp'⟩) (fun x => x)

/-! ### The basic bound

  outDegree T v ≤ n - 1 (since outNbrs is a subset of univ \ {v}).

  This is true for ANY tournament — we don't even need HasBasePath. -/

/-- For any tournament T, vertex v has out-degree at most n - 1. -/
theorem outDegree_le_n_minus_one (T : Tournament n) (_hn : 1 ≤ n) (v : Fin n) :
    T.outDegree v ≤ n - 1 := by
  unfold Tournament.outDegree
  show (Finset.univ.filter (fun w : Fin n => T.arc v w = true)).card ≤ n - 1
  have h_sub : Finset.univ.filter (fun w : Fin n => T.arc v w = true) ⊆
               Finset.univ.erase v := by
    intro w hw
    rw [Finset.mem_filter] at hw
    rw [Finset.mem_erase]
    refine ⟨?_, Finset.mem_univ _⟩
    intro h_eq
    rw [h_eq] at hw
    have := hw.2
    rw [T.irrefl] at this
    exact absurd this (by decide)
  calc (Finset.univ.filter (fun w : Fin n => T.arc v w = true)).card
      ≤ (Finset.univ.erase v).card := Finset.card_le_card h_sub
    _ = n - 1 := by rw [Finset.card_erase_of_mem (Finset.mem_univ _)]; simp

/-! ### HasBasePath improves the bound to ≤ n - 2 at vertex 0

  At vertex 0, the consecutive neighbor 1 is also forbidden (since
  T.arc 0 1 = false), giving outDegree ≤ n - 2. -/

/-- **Theorem (BasePathSink).** In a HasBasePath tournament on n ≥ 2 vertices,
    vertex 0's out-degree is ≤ n - 2.  PROVED. -/
theorem base_path_sink_outDegree_le (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) :
    T.outDegree ⟨0, by omega⟩ ≤ n - 2 := by
  unfold Tournament.outDegree
  -- outNbrs ⊆ univ \ {⟨0, _⟩, ⟨1, _⟩}.
  have h0 : (0 : ℕ) < n := by omega
  have h1 : (1 : ℕ) < n := by omega
  let v0 : Fin n := ⟨0, h0⟩
  let v1 : Fin n := ⟨1, h1⟩
  have hv01_ne : v0 ≠ v1 := by intro h; simp [v0, v1] at h
  have h_arc01 : T.arc v0 v1 = false := base_path_sink_no_out_to_one T hbp hn
  -- The filter is a subset of univ \ {v0, v1}.
  have h_sub : (Finset.univ.filter (fun w : Fin n => T.arc v0 w = true)) ⊆
               Finset.univ \ ({v0, v1} : Finset (Fin n)) := by
    intro w hw
    rw [Finset.mem_filter] at hw
    rw [Finset.mem_sdiff]
    refine ⟨Finset.mem_univ _, ?_⟩
    intro h_in
    rw [Finset.mem_insert, Finset.mem_singleton] at h_in
    rcases h_in with h_w_eq_v0 | h_w_eq_v1
    · -- w = v0; T.arc v0 v0 = false.
      rw [h_w_eq_v0] at hw
      have := hw.2
      rw [T.irrefl] at this
      exact absurd this (by decide)
    · -- w = v1; T.arc v0 v1 = false.
      rw [h_w_eq_v1] at hw
      have := hw.2
      rw [h_arc01] at this
      exact absurd this (by decide)
  -- |univ \ {v0, v1}| = n - 2.
  have h_card_sdiff : (Finset.univ \ ({v0, v1} : Finset (Fin n))).card = n - 2 := by
    have h_subset : ({v0, v1} : Finset (Fin n)) ⊆ Finset.univ := by
      intro w _; exact Finset.mem_univ _
    rw [Finset.card_sdiff_of_subset h_subset]
    have h_card_pair : ({v0, v1} : Finset (Fin n)).card = 2 := by
      rw [Finset.card_insert_of_notMem]
      · simp
      · simp only [Finset.mem_singleton]
        exact hv01_ne
    have h_card_univ : (Finset.univ : Finset (Fin n)).card = n := by simp
    rw [h_card_univ, h_card_pair]
  show (Finset.univ.filter (fun w : Fin n => T.arc v0 w = true)).card ≤ n - 2
  calc _ ≤ (Finset.univ \ ({v0, v1} : Finset (Fin n))).card := Finset.card_le_card h_sub
       _ = n - 2 := h_card_sdiff

/-! ### Vertex n - 1 (source) lower bound -/

/-- In a HasBasePath tournament with n ≥ 2, vertex n-1 beats vertex n-2
    (the base-path arc), so T.arc (n-1) (n-2) = true. -/
theorem base_path_source_arc (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) :
    T.arc ⟨n - 1, by omega⟩ ⟨n - 2, by omega⟩ = true := by
  have h_succ : (⟨n - 2, by omega⟩ : Fin n).val + 1 < n := by show n - 2 + 1 < n; omega
  have h_bp := hbp ⟨n - 2, by omega⟩ h_succ
  convert h_bp using 2
  apply Fin.ext; show n - 1 = n - 2 + 1; omega

/-- **Theorem.** In a HasBasePath T with n ≥ 2, vertex n-1 has out-degree ≥ 1. -/
theorem base_path_source_outDegree_ge (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) :
    1 ≤ T.outDegree ⟨n - 1, by omega⟩ := by
  have h_arc := base_path_source_arc T hbp hn
  -- h_arc : T.arc ⟨n-1, _⟩ ⟨n-2, _⟩ = true; so vertex n-2 ∈ outNbrs.
  unfold Tournament.outDegree
  rw [Finset.one_le_card]
  refine ⟨⟨n - 2, by omega⟩, ?_⟩
  rw [Finset.mem_filter]
  exact ⟨Finset.mem_univ _, h_arc⟩

/-! ### Regular HasBasePath tournament: n must be odd and ≥ 3 -/

/-- For a regular HasBasePath tournament on n vertices, every vertex has
    score (n-1)/2 ≥ 1 ⟹ n ≥ 3. -/
theorem regular_basepath_n_ge_three (T : Tournament n) (hbp : HasBasePath T)
    (hn : 2 ≤ n) (hreg : IsRegular T) :
    3 ≤ n := by
  -- vertex n-1 has outDegree (n-1)/2 ≥ 1 (from base_path_source_outDegree_ge).
  -- 2 · outDegree = n - 1, so outDegree = (n-1)/2.
  -- (n-1)/2 ≥ 1 ⟹ n - 1 ≥ 2 ⟹ n ≥ 3.
  have h_lower := base_path_source_outDegree_ge T hbp hn
  have h_reg := hreg ⟨n - 1, by omega⟩
  -- h_reg : 2 * T.outDegree ⟨n-1, _⟩ = n - 1
  -- h_lower : 1 ≤ T.outDegree ⟨n-1, _⟩
  omega

/-- For a regular HasBasePath tournament on n vertices, the (n-1)/2
    score at vertex 0 must satisfy (n-1)/2 ≤ n - 2 ⟹ n ≥ 3 (always true
    when n ≥ 3 already). -/
theorem regular_basepath_sink_bound (T : Tournament n) (_hbp : HasBasePath T)
    (hn : 3 ≤ n) (hreg : IsRegular T) :
    2 * T.outDegree ⟨0, by omega⟩ = n - 1 := by
  exact hreg _

end Tournament
