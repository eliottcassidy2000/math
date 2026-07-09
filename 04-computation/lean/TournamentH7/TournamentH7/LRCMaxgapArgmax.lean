/-
  TournamentH7.LRCMaxgapArgmax — the argmax + wrapping pieces of klein-S203's `hlink`
  (kind-pasteur-2026-07-09-S102).

  `maxCircGap E Vmax j = (zipWith (·-·) cyc cyc.tail).foldl max 0` (klein's `LRCGoodPeriodMaxgap`).
  A good period gives `maxCircGap > Vmax/7 > 0`, so the widest gap is a genuine list member (the
  EXTREMAL WITNESS — cf. the saturated-witness methodology of Przybocki–Mackey–Heule–Subercaseaux
  2604.21187).  Extracting the free residue interval from it has two branches: the max gap is an
  INTERNAL consecutive gap (non-wrapping, kps-S101 `LRCTeethGap.free_translate_of_free_subInterval`)
  or the WRAPAROUND gap `(r_max, r_min + Vmax)` (wrapping, handled here).

  This file supplies the two reusable pieces:
    * `foldl_max_pos_mem` — a positive `foldl max 0` is attained by a list element (the argmax step);
    * `free_translate_wrap` — the wrapping companion to the S101 reduction: for a gap
      `(a, a+g)` straddling `1` whose complement contains every tooth (`a+g−1 ≤ c ≤ a`), no integer
      translate of a tooth lands in `(a, a+g)`.

  Self-contained (imports Mathlib only).
-/
import Mathlib

namespace LRC14Argmax

/-- **`foldl max` accumulator dichotomy.**  `L.foldl max b` is either the seed `b` or a member of
`L`.  (Induction on `L`, threading the accumulator.) -/
theorem foldl_max_eq_or_mem : ∀ (b : ℕ) (L : List ℕ), L.foldl max b = b ∨ L.foldl max b ∈ L
  | _, [] => Or.inl rfl
  | b, a :: t => by
    rw [List.foldl_cons]
    rcases foldl_max_eq_or_mem (max b a) t with h | h
    · rcases le_total a b with hab | hab
      · exact Or.inl (by rw [h, max_eq_left hab])
      · exact Or.inr (by rw [h, max_eq_right hab]; exact List.mem_cons_self)
    · exact Or.inr (List.mem_cons_of_mem a h)

/-- **A positive `foldl max 0` is a list member** — the argmax witness.  If `L.foldl max 0 > 0`
then that maximum equals some element of `L`. -/
theorem foldl_max_pos_mem (L : List ℕ) (h : 0 < L.foldl max 0) : L.foldl max 0 ∈ L := by
  rcases foldl_max_eq_or_mem 0 L with heq | hmem
  · exact absurd (heq ▸ h) (lt_irrefl 0)
  · exact hmem

/-- **`foldl max` upper bound (accumulator form).**  If the seed and every element are `≤ B`, so is
`L.foldl max b`. -/
theorem foldl_max_le_of_all (B : ℕ) : ∀ (b : ℕ), b ≤ B → ∀ (L : List ℕ), (∀ x ∈ L, x ≤ B) →
    L.foldl max b ≤ B
  | b, hb, [], _ => hb
  | b, hb, a :: t, h => by
    rw [List.foldl_cons]
    exact foldl_max_le_of_all B (max b a) (max_le hb (h a List.mem_cons_self)) t
      (fun x hx => h x (List.mem_cons_of_mem a hx))

/-- **The complement of the argmax step (the density-floor trigger).**  If every gap is `≤ B`, the
maximum circular gap is `≤ B`.  Contrapositive: `maxCircGap > B` forces some gap `> B` (a good
period).  So `¬ HasGoodPeriod` ⟺ every gap `≤ Vmax/7` — the fragmented / near-AP / density-floor
regime. -/
theorem foldl_max_le (L : List ℕ) (B : ℕ) (h : ∀ x ∈ L, x ≤ B) : L.foldl max 0 ≤ B :=
  foldl_max_le_of_all B 0 (Nat.zero_le B) L h

/-- **Wrapping-gap translate reduction** (companion to kps-S101 `free_translate_of_free_subInterval`).
For a gap `(a, a+g)` that straddles `1` (`0 ≤ a`, `a+g ≤ 2`), if every tooth `c ∈ S ⊆ [0,1)` lies in
the COMPLEMENT arc `[a+g−1, a]` (i.e. `≥` the wrap-start and `≤` the wrap-end), then no integer
translate of any tooth lands in `(a, a+g)`.  This is the wraparound branch of `hlink`'s conclusion:
the max gap runs from `r_max/Vmax = a` up over `1` to `r_min/Vmax + 1 = a+g`, and all residues sit in
`[r_min, r_max]`. -/
theorem free_translate_wrap (S : Finset ℝ) (a g : ℝ)
    (hSIco : ∀ c ∈ S, c ∈ Set.Ico (0 : ℝ) 1)
    (ha : 0 ≤ a) (_hag2 : a + g ≤ 2)
    (hfree : ∀ c ∈ S, a + g - 1 ≤ c ∧ c ≤ a) :
    ∀ c ∈ S, ∀ n : ℤ, (c + (n : ℝ)) ∉ Set.Ioo a (a + g) := by
  intro c hc n hmem
  obtain ⟨hlo, hhi⟩ := hmem
  obtain ⟨hc0, hc1⟩ := hSIco c hc
  obtain ⟨hclo, hchi⟩ := hfree c hc
  rcases lt_trichotomy n 0 with hn | hn | hn
  · -- n ≤ -1 : c + n ≤ c - 1 < 0 ≤ a < c + n
    have hn' : n ≤ -1 := by omega
    have hnr : (n : ℝ) ≤ -1 := by exact_mod_cast hn'
    linarith
  · -- n = 0 : a < c ≤ a
    subst hn
    simp only [Int.cast_zero, add_zero] at hlo
    linarith
  · -- n ≥ 1 : c + n ≥ c + 1 ≥ a + g, contradicting c + n < a + g
    have hn' : (1 : ℤ) ≤ n := by omega
    have hnr : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn'
    linarith

end LRC14Argmax
