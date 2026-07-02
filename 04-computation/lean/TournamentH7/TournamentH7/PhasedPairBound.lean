/-
klein-2026-07-02-S103 — THM-604(b), the phased pair bound (the count argument, ℚ layer).

For coprime speeds P, Q and ANY phases φ, ψ, when `2r(P+Q) ≤ 1` the combs' intersection
has length at most `2r/Q` — one deep interval. Mechanism: a positive clip between the
k-th P-interval and l-th Q-interval forces the integer `kQ − lP` into an open interval of
length `2r(P+Q) ≤ 1` (unique), and coprimality makes `(k,l) ↦ kQ − lP` injective on
`range P × range Q` — at most ONE positive clip. The shift-uniform upper half of THM-604
(attainment: OriginNestQ.lean); same residue-uniqueness engine as NestDecidable /
LRCSevenTranslate, now with free phases.
-/
import Mathlib
import TournamentH7.RatIntervals

namespace PhasedPairBound

open LonelyRunner.RatIntervals

/-- Width-after-clip is at most the right factor's width. -/
theorem widthNN_clip_le_right (p q : ℚ × ℚ) :
    max 0 ((clip p q).2 - (clip p q).1) ≤ max 0 (q.2 - q.1) := by
  apply max_le_max (le_refl 0)
  have h1 : min p.2 q.2 ≤ q.2 := min_le_right _ _
  have h2 : q.1 ≤ max p.1 q.1 := le_max_right _ _
  unfold clip
  simp only
  linarith

/-- Exclusion sum: nonnegative entries, each ≤ c, pairwise at-most-one-nonzero ⟹ sum ≤ c. -/
theorem sum_le_of_pairwise_zero {L : List ℚ} {c : ℚ} (hc : 0 ≤ c)
    (hnn : ∀ x ∈ L, 0 ≤ x) (hle : ∀ x ∈ L, x ≤ c)
    (hpw : L.Pairwise (fun a b => a = 0 ∨ b = 0)) : L.sum ≤ c := by
  induction L with
  | nil => simpa using hc
  | cons a t ih =>
    rcases List.pairwise_cons.mp hpw with ⟨hhead, htail⟩
    rcases eq_or_lt_of_le (hnn a (by simp)) with ha | ha
    · simp only [List.sum_cons, ← ha, zero_add]
      exact ih (fun x hx => hnn x (by simp [hx])) (fun x hx => hle x (by simp [hx])) htail
    · have hz : ∀ x ∈ t, x = 0 := fun x hx =>
        (hhead x hx).resolve_left (by linarith)
      have hts : t.sum = 0 := List.sum_eq_zero hz
      simp only [List.sum_cons, hts, add_zero]
      exact hle a (by simp)

/-- Length of a flatMap of regions. -/
theorem length_flatMap {α : Type*} (L : List α) (f : α → Region) :
    length (L.flatMap f) = (L.map fun a => length (f a)).sum := by
  induction L with
  | nil => rfl
  | cons a t ih =>
    rw [List.flatMap_cons, length_append, ih, List.map_cons, List.sum_cons]

/-- The meeting inequality: the (k,l) clip is positive iff
`|(kQ − lP) + (φQ − ψP)| < r(P+Q)`. -/
theorem clip_pos_iff {P Q : ℕ} (hP : 0 < (P : ℚ)) (hQ : 0 < (Q : ℚ))
    {r φ ψ : ℚ} {k l : ℚ} :
    (max ((k + φ - r) / P) ((l + ψ - r) / Q) < min ((k + φ + r) / P) ((l + ψ + r) / Q))
    ↔ |k * Q - l * P + (φ * Q - ψ * P)| < r * ((P : ℚ) + Q) := by
  simp only [max_lt_iff, lt_min_iff, abs_lt, div_lt_div_iff₀ hP hQ,
    div_lt_div_iff₀ hQ hP, div_lt_div_iff₀ hP hP, div_lt_div_iff₀ hQ hQ]
  constructor
  · rintro ⟨⟨-, h1⟩, h2, -⟩
    constructor <;> nlinarith
  · rintro ⟨h1, h2⟩
    refine ⟨⟨?_, ?_⟩, ?_, ?_⟩ <;> nlinarith

/-- Meeting uniqueness. -/
theorem meeting_unique {P Q : ℕ} (hcop : Nat.Coprime P Q)
    {r φ ψ : ℚ} (hr2 : 2 * r * ((P : ℚ) + Q) ≤ 1)
    {k l k' l' : ℕ} (hk : k < P) (hl : l < Q) (hk' : k' < P) (hl' : l' < Q)
    (h : |(k : ℚ) * Q - l * P + (φ * Q - ψ * P)| < r * ((P : ℚ) + Q))
    (h' : |(k' : ℚ) * Q - l' * P + (φ * Q - ψ * P)| < r * ((P : ℚ) + Q)) :
    k = k' ∧ l = l' := by
  set N : ℤ := (k : ℤ) * Q - l * P with hN
  set N' : ℤ := (k' : ℤ) * Q - l' * P with hN'
  have hcast : ((N : ℚ)) - ((N' : ℚ))
      = ((k : ℚ) * Q - l * P + (φ * Q - ψ * P))
        - ((k' : ℚ) * Q - l' * P + (φ * Q - ψ * P)) := by
    simp only [hN, hN']
    push_cast
    ring
  have hdiffQ : |((N : ℚ)) - ((N' : ℚ))| < 1 := by
    rw [hcast]
    calc |((k : ℚ) * Q - l * P + (φ * Q - ψ * P))
          - ((k' : ℚ) * Q - l' * P + (φ * Q - ψ * P))|
        ≤ |(k : ℚ) * Q - l * P + (φ * Q - ψ * P)|
          + |(k' : ℚ) * Q - l' * P + (φ * Q - ψ * P)| := abs_sub _ _
      _ < r * ((P : ℚ) + Q) + r * ((P : ℚ) + Q) := add_lt_add h h'
      _ ≤ 1 := by linarith
  have hNN : N = N' := by
    have h1 : |N - N'| < 1 := by exact_mod_cast hdiffQ
    rw [abs_lt] at h1
    omega
  have hDQ : ((l : ℤ) - l') * P = (Q : ℤ) * ((k : ℤ) - k') := by
    have hexp : (k : ℤ) * Q - l * P = (k' : ℤ) * Q - l' * P := by
      simpa [hN, hN'] using hNN
    linarith
  have hcop' : IsCoprime (Q : ℤ) (P : ℤ) := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    simpa [Nat.gcd_comm] using hcop
  have hQl : (Q : ℤ) ∣ ((l : ℤ) - l') :=
    hcop'.dvd_of_dvd_mul_right ⟨(k : ℤ) - k', hDQ⟩
  obtain ⟨m, hm⟩ := hQl
  have hQpos : (0 : ℤ) < (Q : ℤ) := by
    have h0 : 0 < Q := lt_of_le_of_lt (Nat.zero_le l) hl
    exact_mod_cast h0
  have hm0 : m = 0 := by
    rcases lt_trichotomy m 0 with hmlt | hme | hmgt
    · exfalso
      have h1 : m ≤ -1 := by omega
      have h2 : (Q : ℤ) * m ≤ (Q : ℤ) * (-1) :=
        mul_le_mul_of_nonneg_left h1 (le_of_lt hQpos)
      have h3 : -(Q : ℤ) < (l : ℤ) - l' := by omega
      linarith [hm]
    · exact hme
    · exfalso
      have h1 : (1 : ℤ) ≤ m := by omega
      have h2 : (Q : ℤ) * 1 ≤ (Q : ℤ) * m :=
        mul_le_mul_of_nonneg_left h1 (le_of_lt hQpos)
      have h3 : (l : ℤ) - l' < Q := by omega
      linarith [hm]
  have hE0 : ((l : ℤ) - l') = 0 := by rw [hm, hm0, mul_zero]
  have hll : l = l' := by omega
  have hkk : k = k' := by
    have hz : (Q : ℤ) * ((k : ℤ) - k') = 0 := by rw [← hDQ, hE0, zero_mul]
    rcases mul_eq_zero.mp hz with hq0 | hd0
    · omega
    · omega
  exact ⟨hkk, hll⟩

/-- Positive-sum extraction over nonnegative lists. -/
theorem exists_pos_of_sum_pos {L : List ℚ} (hnn : ∀ x ∈ L, 0 ≤ x) (h : 0 < L.sum) :
    ∃ x ∈ L, 0 < x := by
  by_contra hcon
  push_neg at hcon
  have hz : ∀ x ∈ L, x = 0 := fun x hx => le_antisymm (hcon x hx) (hnn x hx)
  rw [List.sum_eq_zero hz] at h
  exact lt_irrefl 0 h

/-- The width-positivity bridge. -/
theorem widthNN_pos_iff {P Q : ℕ} (hP : 0 < (P : ℚ)) (hQ : 0 < (Q : ℚ))
    {r φ ψ : ℚ} (k l : ℕ) :
    (0 < max 0 ((clip (((k : ℚ) + φ - r) / P, ((k : ℚ) + φ + r) / P)
        (((l : ℚ) + ψ - r) / Q, ((l : ℚ) + ψ + r) / Q)).2
      - (clip (((k : ℚ) + φ - r) / P, ((k : ℚ) + φ + r) / P)
        (((l : ℚ) + ψ - r) / Q, ((l : ℚ) + ψ + r) / Q)).1))
    ↔ |(k : ℚ) * Q - l * P + (φ * Q - ψ * P)| < r * ((P : ℚ) + Q) := by
  rw [← clip_pos_iff hP hQ]
  unfold clip
  simp only [lt_max_iff, lt_irrefl, false_or, sub_pos]

end PhasedPairBound
