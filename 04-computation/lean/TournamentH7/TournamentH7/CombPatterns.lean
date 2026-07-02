/-
  TournamentH7.CombPatterns  (mac-mini-2026-07-02-S10)

  Playbook modules 1–2 core (CombSets + PatternOverlap avoidance), over the
  module-0 `Region` API.  All-ℚ (T1).  Contents:

    • `length_eq_zero_of_no_mem` — the empty-region bridge (module-0 gap fill)
    • `mem_comb` — membership characterization of the danger comb
    • `pattern_overlap_zero` — THM-601/605(i) avoidance at the Region level:
      for `2r(P+Q) ≤ 1` the phase `θ = 1/(2P)` makes the pattern overlap
      `length (inter (comb P r 0) (comb Q r θ)) = 0`.  The half-open interval
      orientations OPPOSE under the subtraction, so the resonance combination
      is STRICTLY inside `(−r(P+Q), r(P+Q))`, covering the boundary `P+Q = 7`.

  Sorry-free target.
-/
import Mathlib.Tactic
import TournamentH7.RatIntervals

namespace TournamentH7.CombPatterns

open LonelyRunner.RatIntervals

/-- A region no point inhabits has zero length (all intervals degenerate). -/
theorem length_eq_zero_of_no_mem {L : Region} (h : ∀ x : ℚ, ¬ mem x L) :
    length L = 0 := by
  induction L with
  | nil => rfl
  | cons p t ih =>
      have hdeg : p.2 ≤ p.1 := by
        by_contra hlt
        push_neg at hlt
        exact h p.1 ⟨p, List.mem_cons_self .., le_refl _, hlt⟩
      have ht : ∀ x : ℚ, ¬ mem x t := by
        intro x ⟨q, hq, h1, h2⟩
        exact h x ⟨q, List.mem_cons_of_mem _ hq, h1, h2⟩
      unfold length at *
      simp only [List.map_cons, List.sum_cons]
      rw [ih ht, max_eq_left (by linarith)]
      ring

/-- Membership in the comb: some tooth `k < v` has `k + φ − r ≤ v·x < k + φ + r`. -/
theorem mem_comb {x : ℚ} {v : ℕ} (hv : 0 < v) {r φ : ℚ} :
    mem x (comb v r φ) ↔ ∃ k : ℕ, k < v ∧ (k : ℚ) + φ - r ≤ v * x ∧ v * x < (k : ℚ) + φ + r := by
  have hv' : (0 : ℚ) < (v : ℚ) := by exact_mod_cast hv
  unfold mem comb
  constructor
  · rintro ⟨p, hp, h1, h2⟩
    simp only [List.map_map, List.mem_map, List.mem_range, Function.comp_apply] at hp
    obtain ⟨k, hk, rfl⟩ := hp
    refine ⟨k, hk, ?_, ?_⟩
    · have := (div_le_iff₀ hv').mp h1
      linarith [this]
    · have := (lt_div_iff₀ hv').mp h2
      linarith [this]
  · rintro ⟨k, hk, h1, h2⟩
    refine ⟨(((k : ℚ) + φ - r) / v, ((k : ℚ) + φ + r) / v), ?_, ?_, ?_⟩
    · simp only [List.map_map, List.mem_map, List.mem_range, Function.comp_apply]
      exact ⟨k, hk, rfl⟩
    · rw [div_le_iff₀ hv']
      linarith
    · rw [lt_div_iff₀ hv']
      linarith

/-- **Pattern-overlap avoidance (THM-605(i), Region level).**  If
`0 < r`, `2r(P+Q) ≤ 1`, then with `θ = 1/(2P)` no point lies in both combs. -/
theorem no_mem_inter_comb (P Q : ℕ) (r : ℚ) (hP : 0 < P) (hQ : 0 < Q)
    (hr : 0 < r) (hsum : 2 * r * ((P : ℚ) + Q) ≤ 1) :
    ∀ x : ℚ, ¬ mem x (inter (comb P r 0) (comb Q r (1 / (2 * P)))) := by
  intro x hx
  have hP' : (0 : ℚ) < (P : ℚ) := by exact_mod_cast hP
  have hQ' : (0 : ℚ) < (Q : ℚ) := by exact_mod_cast hQ
  obtain ⟨h1, h2⟩ := mem_inter.mp hx
  obtain ⟨a, ha, ha1, ha2⟩ := (mem_comb hP).mp h1
  obtain ⟨b, hb, hb1, hb2⟩ := (mem_comb hQ).mp h2
  -- the resonance combination: Q(Px − a) − P(Qx − θ − b) = Pθ + Pb − Qa
  -- with θ = 1/(2P): Pθ = 1/2, so |1/2 + Pb − Qa| < rQ + rP ≤ 1/2: contradiction
  -- (strictness: Q·[−r, r) − P·[−r−θ…): orientations oppose)
  have e1 : -r ≤ (P : ℚ) * x - a - 0 := by linarith
  have e1' : (P : ℚ) * x - a - 0 < r := by linarith
  have e2 : -r ≤ (Q : ℚ) * x - (1 / (2 * P)) - b := by linarith
  have e2' : (Q : ℚ) * x - (1 / (2 * P)) - b < r := by linarith
  set m : ℚ := (Q : ℚ) * a - (P : ℚ) * b with hm
  have key : (P : ℚ) * (1 / (2 * P)) + ((P : ℚ) * b - (Q : ℚ) * a)
      = (P : ℚ) * ((Q : ℚ) * x - (1 / (2 * P)) - b) * (-1) + (Q : ℚ) * ((P : ℚ) * x - a - 0) := by
    field_simp
    ring
  have hub : (P : ℚ) * (1 / (2 * P)) + ((P : ℚ) * b - (Q : ℚ) * a)
      < (P : ℚ) * r + (Q : ℚ) * r := by
    rw [key]
    have t1 : (P : ℚ) * ((Q : ℚ) * x - (1 / (2 * P)) - b) * (-1) ≤ (P : ℚ) * r := by
      have := mul_le_mul_of_nonneg_left e2 (le_of_lt hP')
      nlinarith
    have t2 : (Q : ℚ) * ((P : ℚ) * x - a - 0) < (Q : ℚ) * r :=
      mul_lt_mul_of_pos_left e1' hQ'
    linarith
  have hlb : -((P : ℚ) * r + (Q : ℚ) * r)
      < (P : ℚ) * (1 / (2 * P)) + ((P : ℚ) * b - (Q : ℚ) * a) := by
    rw [key]
    have t1 : -((P : ℚ) * r) < (P : ℚ) * ((Q : ℚ) * x - (1 / (2 * P)) - b) * (-1) := by
      have := mul_lt_mul_of_pos_left e2' hP'
      nlinarith
    have t2 : -((Q : ℚ) * r) ≤ (Q : ℚ) * ((P : ℚ) * x - a - 0) := by
      have := mul_le_mul_of_nonneg_left e1 (le_of_lt hQ')
      nlinarith
    linarith
  have hPθ : (P : ℚ) * (1 / (2 * P)) = 1 / 2 := by field_simp
  rw [hPθ] at hub hlb
  -- so 1/2 + (Pb − Qa) ∈ (−r(P+Q), r(P+Q)) ⊆ (−1/2, 1/2): forces the integer
  -- Pb − Qa ∈ (−1, 0): impossible for an integer combination
  have hrng : (P : ℚ) * r + (Q : ℚ) * r ≤ 1 / 2 := by nlinarith
  -- −1/2 < 1/2 + (Pb − Qa) < 1/2  ⟹  the integer Pb − Qa lies in (−1, 0): impossible
  have hI1 : (-1 : ℚ) < (P : ℚ) * b - (Q : ℚ) * a := by linarith
  have hI2 : (P : ℚ) * b - (Q : ℚ) * a < 0 := by linarith
  obtain ⟨z, hz⟩ : ∃ z : ℤ, ((z : ℚ)) = (P : ℚ) * b - (Q : ℚ) * a :=
    ⟨(P : ℤ) * (b : ℤ) - (Q : ℤ) * (a : ℤ), by push_cast; ring⟩
  rw [← hz] at hI1 hI2
  have hz1 : (-1 : ℤ) < z := by exact_mod_cast hI1
  have hz2 : z < 0 := by exact_mod_cast hI2
  omega

/-- **The avoidance theorem packaged**: the pattern overlap vanishes. -/
theorem pattern_overlap_zero (P Q : ℕ) (r : ℚ) (hP : 0 < P) (hQ : 0 < Q)
    (hr : 0 < r) (hsum : 2 * r * ((P : ℚ) + Q) ≤ 1) :
    length (inter (comb P r 0) (comb Q r (1 / (2 * P)))) = 0 :=
  length_eq_zero_of_no_mem (no_mem_inter_comb P Q r hP hQ hr hsum)

/-! ### The forced direction (THM-605(i) converse), arithmetic level.

For `1 < 2r(P+Q)`, coprime `(P,Q)`: EVERY phase admits a double-covered point.
Construction (S11 plan): `s = Pθ + z ∈ (−r(P+Q), r(P+Q))` by ceiling choice;
the case-free split `u = s/(P+Q)`, `v = −u` solves `Qu − Pv = s` with
`|u|,|v| < r` strictly; Bézout teeth `a = −z·d`, `b = z·c` make the two comb
residues exactly `u` and `v`. -/

theorem exists_double_cover (P Q : ℕ) (r θ : ℚ) (hP : 0 < P) (hQ : 0 < Q)
    (hcop : Nat.Coprime P Q) (hr : 0 < r) (hforce : 1 < 2 * r * ((P : ℚ) + Q)) :
    ∃ x : ℚ, ∃ a b : ℤ,
      -r ≤ (P : ℚ) * x - a ∧ (P : ℚ) * x - a < r ∧
      -r ≤ (Q : ℚ) * x - θ - b ∧ (Q : ℚ) * x - θ - b < r := by
  have hP' : (0 : ℚ) < (P : ℚ) := by exact_mod_cast hP
  have hPne : (P : ℚ) ≠ 0 := ne_of_gt hP'
  have hPQ : (0 : ℚ) < (P : ℚ) + Q := by
    have : (0 : ℚ) < (Q : ℚ) := by exact_mod_cast hQ
    linarith
  obtain ⟨z, hz1, hz2⟩ : ∃ z : ℤ,
      -(r * ((P : ℚ) + Q)) < (P : ℚ) * θ + z ∧ (P : ℚ) * θ + z < r * ((P : ℚ) + Q) := by
    refine ⟨⌊-(r * ((P : ℚ) + Q)) - (P : ℚ) * θ⌋ + 1, ?_, ?_⟩
    · have h := Int.lt_floor_add_one (-(r * ((P : ℚ) + Q)) - (P : ℚ) * θ)
      push_cast
      linarith
    · have h := Int.floor_le (-(r * ((P : ℚ) + Q)) - (P : ℚ) * θ)
      push_cast
      linarith
  set s : ℚ := (P : ℚ) * θ + z with hs
  set u : ℚ := s / ((P : ℚ) + Q) with hu
  have hu1 : -r < u := by
    rw [hu, lt_div_iff₀ hPQ]
    linarith
  have hu2 : u < r := by
    rw [hu, div_lt_iff₀ hPQ]
    linarith
  have huv : (Q : ℚ) * u - (P : ℚ) * (-u) = s := by
    rw [hu]; field_simp; ring
  obtain ⟨c, d, hcd⟩ : ∃ c d : ℤ, (P : ℤ) * c + (Q : ℤ) * d = 1 := by
    refine ⟨Nat.gcdA P Q, Nat.gcdB P Q, ?_⟩
    have h1 := Nat.gcd_eq_gcd_ab P Q
    have h2 : Nat.gcd P Q = 1 := hcop
    rw [h2] at h1
    push_cast at h1 ⊢
    linarith
  have hcdQ : (P : ℚ) * c + (Q : ℚ) * d = 1 := by exact_mod_cast hcd
  have hba : (P : ℚ) * ((z : ℚ) * c) - (Q : ℚ) * (-(z : ℚ) * d) = z := by
    linear_combination (z : ℚ) * hcdQ
  set x : ℚ := ((-(z : ℚ) * d) + u) / P with hxdef
  have hx1 : (P : ℚ) * x - (-(z : ℚ) * d) = u := by
    rw [hxdef]; field_simp; ring
  have hmul : (Q : ℚ) * ((-(z : ℚ) * d) + u) = (θ + ((z : ℚ) * c) + (-u)) * P := by
    linear_combination huv + hs - hba
  have hx2 : (Q : ℚ) * x - θ - ((z : ℚ) * c) = -u := by
    rw [hxdef]
    field_simp
    linear_combination hmul
  refine ⟨x, -z * d, z * c, ?_, ?_, ?_, ?_⟩
  · rw [show ((-z * d : ℤ) : ℚ) = -(z : ℚ) * d by push_cast; ring, hx1]; linarith
  · rw [show ((-z * d : ℤ) : ℚ) = -(z : ℚ) * d by push_cast; ring, hx1]; linarith
  · rw [show ((z * c : ℤ) : ℚ) = (z : ℚ) * c by push_cast; ring, hx2]; linarith
  · rw [show ((z * c : ℤ) : ℚ) = (z : ℚ) * c by push_cast; ring, hx2]; linarith

end TournamentH7.CombPatterns
