/-
  TournamentH7.LRCLem012NearAP — LEM-012 IN LEAN: the near-AP (coherent-branch) good period
  (death-star-2026-07-09-S7, HYP-5810).

  LEM-012 (klein-S196/S197; hypothesis corrected by klein-S201/MISTAKE-129 to `V ≥ Q+1`):
  a co-offset cluster containing an arithmetic progression of length `L` with at most `m ≤ 5`
  stray points admits a good period `j ≤ Q = ⌈7(L−1)/(L−k+6)⌉` — ELEMENTARY, by Dirichlet +
  a pigeonhole gap-split.  This is the branch on which every modular certificate dies (THM-675:
  the 40→41 near-dilation has B5 = 0/260) and on which nothing more is needed: the mechanism
  below is the whole proof.

  THE FORMALIZATION SHAPE (the ∀-translate real form — no `Int.fract`, no circular sorting):
  1. **Dirichlet** (Mathlib `Real.exists_nat_abs_mul_sub_round_le`, the kps-S111 pattern):
     some `j ∈ [1, Q]` has `j·d/V = r + ε` with `r ∈ ℤ`, `|ε| ≤ 1/(Q+1)`.
  2. **The AP collapses to an interval's orbit**: tooth `i` of the AP is
     `(a + i·d)·j/V = (x₀ + i·ε) + i·r`, so its translate-orbit is the orbit of
     `y_i = x₀ + i·ε ∈ [u, w]`, an interval of width `(L−1)|ε| ≤ (L−1)/(Q+1) < (6−m)/7`
     (one-sided: the `i·ε` all share a sign).
  3. **The complement window** `C = Ioo w (u+1)` has length `> (m+1)/7`; it meets NO translate
     of `[u, w]` (an integer strictly between `0` and `1` would be needed), and each stray
     orbit meets the half-open envelope `[w, w+1)` in exactly one representative.
  4. **Pigeonhole**: the `≤ m` stray representatives leave one of the `m+1` equal open pieces
     of `C` untouched — an arc of length `> 1/7` free of EVERY tooth's EVERY translate.

  The output is the free-gap form that the reach machinery consumes directly (the `hfree`
  hypothesis shape of klein-S205's `minReach_ge_of_driftGap` and of
  `LRCComposedRealization.minReach_ge_of_composed_realization`), and that kps-S99's
  good-period dispatch takes as the near-AP branch input.

  Numeric hypothesis, subtraction-free: `7·L + m·Q + m < 6·Q + 13` ⟺ `7(L−1) < (6−m)(Q+1)`.
  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCMreachConcrete

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-! ## Step 3a: an interval's orbit misses the complement window -/

/-- **The complement window dodges the interval's whole orbit.**  If `y ∈ [u, w]` with
`w − u < 1`, no integer translate of `y` lands in `Ioo w (u + 1)`: it would force an integer
strictly between `0` and `1`. -/
theorem translate_notMem_window {u w y : ℝ} (hy : y ∈ Icc u w) (hwu : w - u < 1)
    (n : ℤ) : y + (n : ℝ) ∉ Ioo w (u + 1) := by
  rintro ⟨h1, h2⟩
  obtain ⟨hu, hw⟩ := hy
  have hn0 : (0 : ℝ) < (n : ℝ) := by linarith
  have hn1 : (n : ℝ) < 1 := by linarith
  have : (0 : ℤ) < n := by exact_mod_cast hn0
  have : n < 1 := by exact_mod_cast hn1
  omega

/-- **A length-`≤ 1` window sees each orbit at most once**: two translates of the same point in
`[w, w + 1)` coincide. -/
theorem translate_unique_in_window {w z : ℝ} {n₁ n₂ : ℤ}
    (h₁ : z + (n₁ : ℝ) ∈ Ico w (w + 1)) (h₂ : z + (n₂ : ℝ) ∈ Ico w (w + 1)) : n₁ = n₂ := by
  obtain ⟨h1a, h1b⟩ := h₁
  obtain ⟨h2a, h2b⟩ := h₂
  have hlt1 : ((n₁ - n₂ : ℤ) : ℝ) < 1 := by push_cast; linarith
  have hgt1 : (-1 : ℝ) < ((n₁ - n₂ : ℤ) : ℝ) := by push_cast; linarith
  have h1 : (n₁ - n₂ : ℤ) < 1 := by exact_mod_cast hlt1
  have h2 : (-1 : ℤ) < (n₁ - n₂ : ℤ) := by exact_mod_cast hgt1
  omega

/-! ## Step 4: the finite pigeonhole on equal pieces -/

/-- **The pigeonhole gap-split.**  At most `m` points cannot touch every one of the `m + 1`
equal open pieces of `[w, w + c)`: some piece `Ioo (w + i·ℓ) (w + (i+1)·ℓ)`, `ℓ = c/(m+1)`,
is free of all of them. -/
theorem exists_free_piece (P : Finset ℝ) (m : ℕ) (hP : P.card ≤ m) (w c : ℝ) (hc : 0 < c) :
    ∃ i : ℕ, i < m + 1 ∧ ∀ x ∈ P,
      x ∉ Ioo (w + i * (c / (m + 1))) (w + (i + 1) * (c / (m + 1))) := by
  classical
  set ℓ : ℝ := c / (m + 1) with hℓ
  have hℓpos : 0 < ℓ := by
    apply div_pos hc
    positivity
  -- index map: clamp ⌊(x − w)/ℓ⌋ into Fin (m+1)
  set idx : ℝ → Fin (m + 1) := fun x =>
    ⟨min m (⌊(x - w) / ℓ⌋.toNat), by omega⟩ with hidx
  have hcard : (P.image idx).card ≤ m := le_trans (Finset.card_image_le) hP
  have hne : ∃ i : Fin (m + 1), i ∉ P.image idx := by
    by_contra hall
    push_neg at hall
    have hsub : (Finset.univ : Finset (Fin (m + 1))) ⊆ P.image idx := fun i _ => hall i
    have := Finset.card_le_card hsub
    rw [Finset.card_univ, Fintype.card_fin] at this
    omega
  obtain ⟨i₀, hi₀⟩ := hne
  refine ⟨(i₀ : ℕ), i₀.isLt, fun x hx hmem => ?_⟩
  obtain ⟨hlo, hhi⟩ := hmem
  -- x in the open piece i₀ forces idx x = i₀ — contradicting i₀ ∉ image
  apply hi₀
  have hfloor : ⌊(x - w) / ℓ⌋ = (i₀ : ℕ) := by
    rw [Int.floor_eq_iff]
    constructor
    · rw [le_div_iff₀ hℓpos]
      push_cast
      nlinarith [hlo]
    · rw [div_lt_iff₀ hℓpos]
      push_cast
      push_cast at hhi
      nlinarith [hhi]
  have hidxx : idx x = i₀ := by
    simp only [hidx]
    apply Fin.ext
    simp only [hfloor, Int.toNat_natCast]
    exact min_eq_right (by omega : (i₀ : ℕ) ≤ m)
  exact Finset.mem_image.mpr ⟨x, hx, hidxx⟩

/-! ## The main theorem -/

set_option maxHeartbeats 800000 in
/-- **LEM-012 (klein-S196/S197, MISTAKE-129-corrected), the near-AP free gap.**

Cluster teeth at ruler `V`, period `j`: the points `e·j/V` for `e ∈ E`, taken with all their
integer translates.  Suppose `E` consists of an arithmetic progression `a, a+d, …, a+(L−1)d`
(`L ≥ 2`, `d ≠ 0`) together with at most `m ≤ 5` stray points, and

  * `hnum  : 7·L + m·Q + m < 6·Q + 13`   (⟺ `7(L−1) < (6−m)(Q+1)`, the Dirichlet budget),
  * `_hQV  : Q < V`                       (the MISTAKE-129 hypothesis `V ≥ Q+1`; the free-gap
    statement is V-free, but the consumer needs `j ≤ Q < V` for `j` to be a nontrivial period).

Then some period `j ∈ [1, Q]` leaves an open arc of length `> 1/7` that avoids EVERY integer
translate of EVERY tooth — the good period, in the exact `hfree` shape the drift-embed /
composed-realization machinery consumes.  Elementary: Dirichlet + pigeonhole. -/
theorem lem012_nearAP_free_gap
    (V : ℤ) (hV : 0 < V) (a d : ℤ) (L m Q : ℕ)
    (hL2 : 2 ≤ L) (hm5 : m ≤ 5) (hQ1 : 1 ≤ Q)
    (hnum : 7 * L + m * Q + m < 6 * Q + 13)
    (_hQV : (Q : ℤ) < V)  -- MISTAKE-129: carried so the consumer's j ≤ Q is a NONTRIVIAL period (j < V)
    (strays : Finset ℤ) (hstray : strays.card ≤ m) :
    ∃ j : ℕ, 1 ≤ j ∧ j ≤ Q ∧ ∃ α g : ℝ, 1 / 7 < g ∧
      (∀ i : ℕ, i < L → ∀ n : ℤ,
        (((a + (i : ℤ) * d : ℤ) : ℝ) * (j : ℝ) / (V : ℝ) + (n : ℝ)) ∉ Ioo α (α + g)) ∧
      (∀ e ∈ strays, ∀ n : ℤ,
        (((e : ℤ) : ℝ) * (j : ℝ) / (V : ℝ) + (n : ℝ)) ∉ Ioo α (α + g)) := by
  classical
  have hVR : (0 : ℝ) < (V : ℝ) := by exact_mod_cast hV
  -- ── Step 1: Dirichlet on ξ = d/V at depth Q
  obtain ⟨j, hj1, hjQ, hround⟩ :=
    Real.exists_nat_abs_mul_sub_round_le ((d : ℝ) / (V : ℝ)) (n := Q) hQ1
  refine ⟨j, hj1, hjQ, ?_⟩
  set ξ : ℝ := (j : ℝ) * ((d : ℝ) / (V : ℝ)) with hξ
  set r : ℤ := round ξ with hr
  set ε : ℝ := ξ - (r : ℝ) with hε
  have hεabs : |ε| ≤ 1 / ((Q : ℝ) + 1) := by
    simpa [hε, hξ, hr] using hround
  -- ── Step 2: the AP orbit lies in the orbit of the interval [u, w]
  set x₀ : ℝ := (a : ℝ) * (j : ℝ) / (V : ℝ) with hx₀
  set u : ℝ := min x₀ (x₀ + ((L : ℝ) - 1) * ε) with hu
  set w : ℝ := max x₀ (x₀ + ((L : ℝ) - 1) * ε) with hw
  have huw : u ≤ w := min_le_max
  have hQR : (0 : ℝ) < (Q : ℝ) + 1 := by positivity
  have hL1R : (1 : ℝ) ≤ (L : ℝ) - 1 := by
    have h2L : (2 : ℝ) ≤ (L : ℝ) := by exact_mod_cast hL2
    linarith
  have hwidth : w - u ≤ ((L : ℝ) - 1) * |ε| := by
    rcases le_or_gt 0 ε with hpos | hneg
    · have h0 : 0 ≤ ((L : ℝ) - 1) * ε := mul_nonneg (by linarith) hpos
      rw [hw, hu, max_eq_right (by linarith : x₀ ≤ x₀ + ((L : ℝ) - 1) * ε),
        min_eq_left (by linarith : x₀ ≤ x₀ + ((L : ℝ) - 1) * ε),
        abs_of_nonneg hpos]
      linarith
    · have h0 : ((L : ℝ) - 1) * ε ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos (by linarith) (le_of_lt hneg)
      have hring : ((L : ℝ) - 1) * (-ε) = -(((L : ℝ) - 1) * ε) := by ring
      rw [hw, hu, max_eq_left (by linarith : x₀ + ((L : ℝ) - 1) * ε ≤ x₀),
        min_eq_right (by linarith : x₀ + ((L : ℝ) - 1) * ε ≤ x₀),
        abs_of_neg hneg]
      linarith
  -- numeric: width < (6−m)/7, so the window length c > (m+1)/7 ≥ (m+1) · 1/7
  have hwidth2 : w - u < (6 - (m : ℝ)) / 7 := by
    have hnumR : 7 * ((L : ℝ) - 1) < (6 - (m : ℝ)) * ((Q : ℝ) + 1) := by
      have h7L : (7 : ℝ) * (L : ℝ) + (m : ℝ) * (Q : ℝ) + (m : ℝ)
          < 6 * (Q : ℝ) + 13 := by exact_mod_cast hnum
      nlinarith
    have hb : ((L : ℝ) - 1) * |ε| ≤ ((L : ℝ) - 1) / ((Q : ℝ) + 1) := by
      rw [div_eq_mul_inv]
      apply mul_le_mul_of_nonneg_left _ (by linarith)
      rw [← one_div]
      exact hεabs
    have hd : ((L : ℝ) - 1) / ((Q : ℝ) + 1) < (6 - (m : ℝ)) / 7 := by
      rw [div_lt_div_iff₀ hQR (by norm_num)]
      linarith
    linarith
  have hm6 : (0 : ℝ) < 6 - (m : ℝ) := by
    have : (m : ℝ) ≤ 5 := by exact_mod_cast hm5
    linarith
  set c : ℝ := (u + 1) - w with hc
  have hcpos : (m + 1 : ℝ) / 7 < c := by
    have : w - u < (6 - (m : ℝ)) / 7 := hwidth2
    rw [hc]
    push_cast
    linarith
  have hc0 : 0 < c := lt_trans (by positivity) hcpos
  have hwu1 : w - u < 1 := by
    have : (6 - (m : ℝ)) / 7 ≤ 1 := by
      rw [div_le_one (by norm_num)]
      have : (0 : ℝ) ≤ (m : ℝ) := by positivity
      linarith
    linarith
  -- ── Step 3: stray representatives inside the envelope [w, w+1)
  set posOf : ℤ → ℝ := fun e => ((e : ℤ) : ℝ) * (j : ℝ) / (V : ℝ) with hposOf
  set rep : ℤ → ℝ := fun e => w + Int.fract (posOf e - w) with hrep
  have hrep_mem : ∀ e, rep e ∈ Ico w (w + 1) := by
    intro e
    constructor
    · simp only [hrep]
      have := Int.fract_nonneg (posOf e - w)
      linarith
    · simp only [hrep]
      have := Int.fract_lt_one (posOf e - w)
      linarith
  have hrep_orbit : ∀ e, ∃ n : ℤ, rep e = posOf e + (n : ℝ) := by
    intro e
    refine ⟨-⌊posOf e - w⌋, ?_⟩
    simp only [hrep, Int.fract]
    push_cast
    ring
  set P : Finset ℝ := strays.image rep with hP
  have hPcard : P.card ≤ m := le_trans (Finset.card_image_le) hstray
  -- ── Step 4: the free piece
  obtain ⟨i₀, hi₀lt, hi₀free⟩ := exists_free_piece P m hPcard w c hc0
  set ℓ : ℝ := c / (m + 1) with hℓ
  have hℓpos : 0 < ℓ := by
    apply div_pos hc0
    positivity
  have hℓ7 : 1 / 7 < ℓ := by
    rw [hℓ, lt_div_iff₀ (by positivity : (0:ℝ) < (m:ℝ) + 1)]
    calc (1 : ℝ) / 7 * ((m : ℝ) + 1) = ((m : ℝ) + 1) / 7 := by ring
      _ < c := by
          have := hcpos
          push_cast at this ⊢
          linarith
  refine ⟨w + i₀ * ℓ, ℓ, hℓ7, ?_, ?_⟩
  · -- the AP teeth: their orbits avoid the whole window Ioo w (u+1) ⊇ the piece
    intro i hiL n hmem
    -- the tooth's position equals (x₀ + i·ε) + (i·r), an orbit point of y_i ∈ [u, w]
    have halg : ((a + (i : ℤ) * d : ℤ) : ℝ) * (j : ℝ) / (V : ℝ) + (n : ℝ)
        = (x₀ + (i : ℝ) * ε) + (((i : ℤ) * r + n : ℤ) : ℝ) := by
      simp only [hx₀, hε, hξ]
      push_cast
      field_simp
      ring
    have hi0R : (0 : ℝ) ≤ (i : ℝ) := by positivity
    have hiL' : (i : ℝ) ≤ (L : ℝ) - 1 := by
      have hn : i + 1 ≤ L := hiL
      have : ((i : ℝ)) + 1 ≤ (L : ℝ) := by exact_mod_cast hn
      linarith
    have hyi : x₀ + (i : ℝ) * ε ∈ Icc u w := by
      constructor
      · rcases le_or_gt 0 ε with hpos | hneg
        · have h1 : 0 ≤ (i : ℝ) * ε := mul_nonneg hi0R hpos
          have h2 : u ≤ x₀ := min_le_left _ _
          linarith
        · have h1 : ((L : ℝ) - 1) * ε ≤ (i : ℝ) * ε :=
            mul_le_mul_of_nonpos_right hiL' (le_of_lt hneg)
          have h2 : u ≤ x₀ + ((L : ℝ) - 1) * ε := min_le_right _ _
          linarith
      · rcases le_or_gt 0 ε with hpos | hneg
        · have h1 : (i : ℝ) * ε ≤ ((L : ℝ) - 1) * ε :=
            mul_le_mul_of_nonneg_right hiL' hpos
          have h2 : x₀ + ((L : ℝ) - 1) * ε ≤ w := le_max_right _ _
          linarith
        · have h1 : (i : ℝ) * ε ≤ 0 :=
            mul_nonpos_of_nonneg_of_nonpos hi0R (le_of_lt hneg)
          have h2 : x₀ ≤ w := le_max_left _ _
          linarith
    have hsub : Ioo (w + i₀ * ℓ) (w + i₀ * ℓ + ℓ) ⊆ Ioo w (u + 1) := by
      intro x hx
      obtain ⟨hx1, hx2⟩ := hx
      constructor
      · have : (0 : ℝ) ≤ i₀ * ℓ := by positivity
        linarith
      · have hi₀m : (i₀ : ℝ) + 1 ≤ (m : ℝ) + 1 := by
          have : (i₀ : ℝ) ≤ (m : ℝ) := by exact_mod_cast Nat.lt_succ_iff.mp hi₀lt
          linarith
        have hprod : ((i₀ : ℝ) + 1) * ℓ ≤ ((m : ℝ) + 1) * ℓ :=
          mul_le_mul_of_nonneg_right hi₀m (le_of_lt hℓpos)
        have hmc : ((m : ℝ) + 1) * ℓ = c := by
          rw [hℓ]
          field_simp
        have hexp : ((i₀ : ℝ) + 1) * ℓ = i₀ * ℓ + ℓ := by ring
        have hcu : w + c = u + 1 := by rw [hc]; ring
        linarith
    rw [halg] at hmem
    exact translate_notMem_window hyi hwu1 _ (hsub hmem)
  · -- the strays: the unique window representative is dodged by the free piece
    intro e he n hmem
    have hsub2 : Ioo (w + i₀ * ℓ) (w + i₀ * ℓ + ℓ) ⊆ Ico w (w + 1) := by
      intro x hx
      obtain ⟨hx1, hx2⟩ := hx
      constructor
      · have : (0 : ℝ) ≤ i₀ * ℓ := by positivity
        linarith
      · have hi₀m : (i₀ : ℝ) + 1 ≤ (m : ℝ) + 1 := by
          have : (i₀ : ℝ) ≤ (m : ℝ) := by exact_mod_cast Nat.lt_succ_iff.mp hi₀lt
          linarith
        have hprod : ((i₀ : ℝ) + 1) * ℓ ≤ ((m : ℝ) + 1) * ℓ :=
          mul_le_mul_of_nonneg_right hi₀m (le_of_lt hℓpos)
        have hmc : ((m : ℝ) + 1) * ℓ = c := by
          rw [hℓ]
          field_simp
        have hexp : ((i₀ : ℝ) + 1) * ℓ = i₀ * ℓ + ℓ := by ring
        have hcu : w + c = u + 1 := by rw [hc]; ring
        have hu1 : u + 1 ≤ w + 1 := by linarith
        linarith
    -- the offending translate is IN the envelope, hence equals rep e — but rep e is dodged
    obtain ⟨n', hn'⟩ := hrep_orbit e
    have hmem' : posOf e + (n : ℝ) ∈ Ico w (w + 1) := hsub2 hmem
    have hrepmem : posOf e + (n' : ℝ) ∈ Ico w (w + 1) := by
      rw [← hn']
      exact hrep_mem e
    have : n = n' := translate_unique_in_window hmem' hrepmem
    subst this
    have : rep e ∈ Ioo (w + i₀ * ℓ) (w + (i₀ + 1) * ℓ) := by
      rw [hn']
      have : w + i₀ * ℓ + ℓ = w + (i₀ + 1) * ℓ := by ring
      rw [← this]
      exact hmem
    exact hi₀free (rep e) (Finset.mem_image_of_mem rep he) this

/-! ## Axiom audit -/
#print axioms translate_notMem_window
#print axioms translate_unique_in_window
#print axioms exists_free_piece
#print axioms lem012_nearAP_free_gap

end LRC14Concrete
end LonelyRunner
