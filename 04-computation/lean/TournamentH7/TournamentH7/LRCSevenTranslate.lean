/-
  TournamentH7.LRCSevenTranslate — MODULE 3 KEYSTONE (klein HYP-4007 spec: K1 + K2 engine).

  MANUAL MERGE of two concurrent same-day builds (fleet rule: never union-merge .lean):
   * klein-2026-07-02-S102b — K1-arith: the ZMod 7 residue selector (`exists_unique_translate`),
     first-pusher, kept verbatim below (import narrowed from `Mathlib` to the cached
     `Mathlib.Data.ZMod.Basic` + project modules; content untouched).
   * kind-pasteur-2026-07-02-S8 (HYP-3965) — the interval-level layer klein spec'd as "next":
     the DIAL LEMMA, interval-level K1 (`sevenTranslate_existsUnique`), cover/disjoint
     corollaries, and the K2 chain engine (`clip_chain_sum`).

  THE ℚ RE-PROOF LAYER (opus INSIGHT 1): with half-open intervals the seven `1/7`-translates
  of a `7 ∤ P` comb at `r = 1/14` tile `[0,1)` EXACTLY — no null sets, no Haar measure.

  THE DIAL LEMMA (`mem_sevenTranslate_iff_dvd`, this file's engine): for `x ∈ [0,1)`,
  membership in the `j`-th translate is the pure INTEGER congruence

      mem x (translateCirc (j/7) (comb P (1/14) ψ))  ↔  7 ∣ ⌊7(Px − ψ) + 1/2⌋ − P·j.

  The rational `x` enters only through the dial reading `N := ⌊7(Px − ψ) + 1/2⌋`; the seven
  translates are the seven residues of `N` mod 7, twisted by `P`.  K1 (exactly one `j < 7`,
  `sevenTranslate_existsUnique`) is then Euclidean division `N = 7M + i` plus invertibility
  of `P` mod 7 — decided by `omega` over the six residue cases (the same bijection klein's
  ZMod selector expresses structurally; both cores are kept, either can feed K1).

  K2 ENGINE (`clip_chain_sum`): the telescoping cursor identity — clipping an interval
  against a uniform chain of `N` abutting cells sums to exactly the covered length
  `max 0 (b − max a c)`.  This is the length-additivity engine for the partition form
  `length L = Σ_j length (inter L T_j)` (K2) whose remaining glue is the bookkeeping from
  wrapped comb teeth to chains; K3 (commensuration-ℚ, `7∣Q, 7∤P ⟹ overlap = 1/49` at every
  phase) = K2 applied to `L = wrap (comb Q r φ)` + equality of the seven terms by
  `1/7`-periodicity.  Per-instance K3 rows are already machine-checked (opus CommensurationQ).
-/
import TournamentH7.CombPatterns
import TournamentH7.RatIntervalsWrap
import Mathlib.Data.ZMod.Basic

/-! ## K1-arith (klein-2026-07-02-S102b, verbatim): the ZMod 7 residue selector -/

namespace LRCSevenTranslate

/-- K1-arith: for 7 ∤ P, every residue class mod 7 is hit by P·j for exactly one j < 7. -/
theorem exists_unique_translate (P : ℕ) (hP : ¬ (7 ∣ P)) (t : ZMod 7) :
    ∃! j : ZMod 7, (P : ZMod 7) * j = t := by
  have hne : (P : ZMod 7) ≠ 0 := by
    rw [Ne, CharP.cast_eq_zero_iff (ZMod 7) 7 P]
    exact hP
  haveI : Fact (Nat.Prime 7) := ⟨by norm_num⟩
  refine ⟨(P : ZMod 7)⁻¹ * t, ?_, fun y hy => ?_⟩
  · field_simp
  · have h1 : (P : ZMod 7)⁻¹ * ((P : ZMod 7) * y) = (P : ZMod 7)⁻¹ * t := by rw [hy]
    rwa [← mul_assoc, inv_mul_cancel₀ hne, one_mul] at h1

end LRCSevenTranslate

/-! ## The interval level (kind-pasteur-S8): dial lemma, K1, K2 engine -/

namespace LonelyRunner
namespace RatIntervals

/-! ### Tooth widths: every translated comb tooth has width `2r/v ≤ 1` -/

theorem width_translate_comb {v : ℕ} (hv : 0 < v) {r φ t : ℚ} (hr2 : 2 * r ≤ (v : ℚ)) :
    ∀ p ∈ translate t (comb v r φ), p.2 - p.1 ≤ 1 := by
  intro p hp
  have hv' : (0 : ℚ) < (v : ℚ) := by exact_mod_cast hv
  simp only [translate, comb, List.map_map, List.mem_map, List.mem_range,
    Function.comp_apply] at hp
  obtain ⟨k, hk, rfl⟩ := hp
  simp only
  have hdiff : ((k : ℚ) + φ + r) / v + t - (((k : ℚ) + φ - r) / v + t) = 2 * r / v := by
    field_simp
    ring
  rw [hdiff, div_le_one hv']
  linarith

/-! ### Membership through a circle-translated comb, as tooth arithmetic -/

theorem mem_translateCirc_comb {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {v : ℕ} (hv : 0 < v)
    {r φ t : ℚ} (hr2 : 2 * r ≤ (v : ℚ)) :
    mem x (translateCirc t (comb v r φ)) ↔
      ∃ n : ℤ, ∃ k : ℕ, k < v ∧ (k : ℚ) + φ - r ≤ (v : ℚ) * (x + (n : ℚ) - t) ∧
        (v : ℚ) * (x + (n : ℚ) - t) < (k : ℚ) + φ + r := by
  unfold translateCirc
  rw [mem_wrap hx0 hx1 (width_translate_comb hv hr2)]
  exact exists_congr fun n => by
    rw [mem_translate, TournamentH7.CombPatterns.mem_comb hv]

/-! ### THE DIAL LEMMA: translate membership = one integer congruence (r = 1/14) -/

theorem mem_sevenTranslate_iff_dvd {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {P : ℕ} (hP : 0 < P)
    (ψ : ℚ) (j : ℕ) :
    mem x (translateCirc ((j : ℚ) / 7) (comb P (1/14) ψ)) ↔
      (7 : ℤ) ∣ ⌊7 * ((P : ℚ) * x - ψ) + 1/2⌋ - (P : ℤ) * (j : ℤ) := by
  have hP' : (0 : ℚ) < (P : ℚ) := by exact_mod_cast hP
  have hP1 : (1 : ℚ) ≤ (P : ℚ) := by exact_mod_cast hP
  have hPz : (0 : ℤ) < (P : ℤ) := by exact_mod_cast hP
  have hPzne : (P : ℤ) ≠ 0 := ne_of_gt hPz
  rw [mem_translateCirc_comb hx0 hx1 hP (by linarith : 2 * (1/14 : ℚ) ≤ (P : ℚ))]
  constructor
  · rintro ⟨n, k, hk, h1, h2⟩
    have hexp : (P : ℚ) * (x + (n : ℚ) - (j : ℚ) / 7)
        = (P : ℚ) * x + (P : ℚ) * (n : ℚ) - (P : ℚ) * (j : ℚ) / 7 := by ring
    rw [hexp] at h1 h2
    have hfl : ⌊7 * ((P : ℚ) * x - ψ) + 1/2⌋
        = (P : ℤ) * (j : ℤ) + 7 * ((k : ℤ) - (P : ℤ) * n) := by
      rw [Int.floor_eq_iff]
      constructor
      · push_cast
        linarith
      · push_cast
        linarith
    rw [hfl]
    exact ⟨(k : ℤ) - (P : ℤ) * n, by ring⟩
  · rintro ⟨m, hm⟩
    have hNle := Int.floor_le (7 * ((P : ℚ) * x - ψ) + 1/2)
    have hNlt := Int.lt_floor_add_one (7 * ((P : ℚ) * x - ψ) + 1/2)
    have hmZ : ⌊7 * ((P : ℚ) * x - ψ) + 1/2⌋ = (P : ℤ) * (j : ℤ) + 7 * m := by linarith
    have hmQ := congrArg (fun z : ℤ => (z : ℚ)) hmZ
    push_cast at hmQ
    have hem0 : 0 ≤ m % (P : ℤ) := Int.emod_nonneg m hPzne
    have hemP : m % (P : ℤ) < (P : ℤ) := Int.emod_lt_of_pos m hPz
    have htn : (((m % (P : ℤ)).toNat : ℕ) : ℚ) = ((m % (P : ℤ) : ℤ) : ℚ) := by
      exact_mod_cast Int.toNat_of_nonneg hem0
    have hme : ((m % (P : ℤ) : ℤ) : ℚ) + (P : ℚ) * ((m / (P : ℤ) : ℤ) : ℚ) = (m : ℚ) := by
      exact_mod_cast Int.emod_add_mul_ediv m (P : ℤ)
    refine ⟨-(m / (P : ℤ)), (m % (P : ℤ)).toNat, ?_, ?_, ?_⟩
    · omega
    · have hexp : (P : ℚ) * (x + ((-(m / (P : ℤ)) : ℤ) : ℚ) - (j : ℚ) / 7)
          = (P : ℚ) * x - (P : ℚ) * ((m / (P : ℤ) : ℤ) : ℚ) - (P : ℚ) * (j : ℚ) / 7 := by
        push_cast
        ring
      rw [htn, hexp]
      linarith
    · have hexp : (P : ℚ) * (x + ((-(m / (P : ℤ)) : ℤ) : ℚ) - (j : ℚ) / 7)
          = (P : ℚ) * x - (P : ℚ) * ((m / (P : ℤ) : ℤ) : ℚ) - (P : ℚ) * (j : ℚ) / 7 := by
        push_cast
        ring
      rw [htn, hexp]
      linarith

/-! ### The integer core: exactly one residue class hits the dial -/

theorem dial_existsUnique (N : ℤ) {P : ℕ} (h7 : ¬ (7 ∣ P)) :
    ∃! j : ℕ, j < 7 ∧ (7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ) := by
  have hp7 : P % 7 < 7 := Nat.mod_lt _ (by norm_num)
  have hp0 : P % 7 ≠ 0 := fun h => h7 (Nat.dvd_of_mod_eq_zero h)
  have hPd : (P : ℤ) = 7 * ((P / 7 : ℕ) : ℤ) + ((P % 7 : ℕ) : ℤ) := by
    exact_mod_cast (Nat.div_add_mod P 7).symm
  have key : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔
      ((7 : ℤ) ∣ N - ((P % 7 : ℕ) : ℤ) * (j : ℤ)) := by
    intro j
    have hexp : N - (P : ℤ) * (j : ℤ)
        = (N - ((P % 7 : ℕ) : ℤ) * (j : ℤ)) - 7 * (((P / 7 : ℕ) : ℤ) * (j : ℤ)) := by
      rw [hPd]
      ring
    rw [hexp]
    constructor
    · rintro ⟨c, hc⟩
      exact ⟨c + ((P / 7 : ℕ) : ℤ) * (j : ℤ), by linarith⟩
    · rintro ⟨c, hc⟩
      exact ⟨c - ((P / 7 : ℕ) : ℤ) * (j : ℤ), by linarith⟩
  have hcases : P % 7 = 1 ∨ P % 7 = 2 ∨ P % 7 = 3 ∨ P % 7 = 4 ∨ P % 7 = 5 ∨ P % 7 = 6 := by
    omega
  rcases hcases with h | h | h | h | h | h
  -- inverse table mod 7: 1↔1, 2↔4, 3↔5, 4↔2, 5↔3, 6↔6
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 1 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((1 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 2 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((4 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 3 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((5 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 4 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((2 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 5 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((3 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega
  · have keyl : ∀ j : ℕ, ((7 : ℤ) ∣ N - (P : ℤ) * (j : ℤ)) ↔ ((7 : ℤ) ∣ N - 6 * (j : ℤ)) := by
      intro j; rw [key j, h]; norm_num
    refine ⟨((6 * N) % 7).toNat, ⟨by omega, (keyl _).mpr (by omega)⟩, ?_⟩
    rintro y ⟨hy7, hyd⟩
    have hy := (keyl y).mp hyd
    omega

/-! ### K1: the seven-translate exact-partition membership (klein HYP-4007 [K1]) -/

/-- **K1 (`sevenTranslate_mem_unique`)**: for `7 ∤ P`, every `x ∈ [0,1)` lies in EXACTLY ONE
of the seven `1/7`-translates of the `1/14`-comb of speed `P` — at every phase `ψ`.  The
half-open seven-translates tile `[0,1)` exactly: no null sets, no measure theory. -/
theorem sevenTranslate_existsUnique {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {P : ℕ} (hP : 0 < P)
    (h7 : ¬ (7 ∣ P)) (ψ : ℚ) :
    ∃! j : ℕ, j < 7 ∧ mem x (translateCirc ((j : ℚ) / 7) (comb P (1/14) ψ)) := by
  obtain ⟨j, ⟨hj7, hjd⟩, huniq⟩ := dial_existsUnique ⌊7 * ((P : ℚ) * x - ψ) + 1/2⌋ h7
  refine ⟨j, ⟨hj7, (mem_sevenTranslate_iff_dvd hx0 hx1 hP ψ j).mpr hjd⟩, ?_⟩
  rintro y ⟨hy7, hym⟩
  exact huniq y ⟨hy7, (mem_sevenTranslate_iff_dvd hx0 hx1 hP ψ y).mp hym⟩

/-- Cover half of K1: some translate contains `x`. -/
theorem sevenTranslate_cover {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {P : ℕ} (hP : 0 < P)
    (h7 : ¬ (7 ∣ P)) (ψ : ℚ) :
    ∃ j : ℕ, j < 7 ∧ mem x (translateCirc ((j : ℚ) / 7) (comb P (1/14) ψ)) := by
  obtain ⟨j, hj, _⟩ := sevenTranslate_existsUnique hx0 hx1 hP h7 ψ
  exact ⟨j, hj⟩

/-- Disjointness half of K1: no `x` lies in two distinct translates. -/
theorem sevenTranslate_disjoint {x : ℚ} (hx0 : 0 ≤ x) (hx1 : x < 1) {P : ℕ} (hP : 0 < P)
    (h7 : ¬ (7 ∣ P)) (ψ : ℚ) {j j' : ℕ} (hj : j < 7) (hj' : j' < 7) (hne : j ≠ j') :
    ¬ (mem x (translateCirc ((j : ℚ) / 7) (comb P (1/14) ψ)) ∧
       mem x (translateCirc ((j' : ℚ) / 7) (comb P (1/14) ψ))) := by
  rintro ⟨hm, hm'⟩
  obtain ⟨i, _, huniq⟩ := sevenTranslate_existsUnique hx0 hx1 hP h7 ψ
  exact hne ((huniq j ⟨hj, hm⟩).trans (huniq j' ⟨hj', hm'⟩).symm)

/-! ### K2 ENGINE: the telescoping chain-clip identity (the cursor lemma's core) -/

/-- **The chain-clip telescope**: clipping `[a, b)` against the uniform chain of `N` abutting
cells `[c + Kw, c + (K+1)w)` sums to exactly the covered length `max 0 (b − max a c)`.
This is the length-additivity engine behind K2 (`length L = Σ_j length (inter L T_j)`):
partitions contribute their exact overlap, cell by cell, with no measure theory. -/
theorem clip_chain_sum {w : ℚ} (hw : 0 < w) :
    ∀ (N : ℕ) (c a b : ℚ), a ≤ b → b ≤ c + (N : ℚ) * w →
      ((List.range N).map fun K : ℕ =>
          max 0 (min b (c + ((K : ℚ) + 1) * w) - max a (c + (K : ℚ) * w))).sum
        = max 0 (b - max a c) := by
  intro N
  induction N with
  | zero =>
      intro c a b hab hb
      simp only [Nat.cast_zero, zero_mul, add_zero] at hb
      simp only [List.range_zero, List.map_nil, List.sum_nil]
      have hbc : b - max a c ≤ 0 := by
        have := le_max_right a c
        linarith
      exact (max_eq_left hbc).symm
  | succ N ih =>
      intro c a b hab hb
      rw [List.range_succ_eq_map, List.map_cons, List.sum_cons, List.map_map]
      have htail : (List.range N).map ((fun K : ℕ =>
            max 0 (min b (c + ((K : ℚ) + 1) * w) - max a (c + (K : ℚ) * w))) ∘ Nat.succ)
          = (List.range N).map (fun K : ℕ =>
            max 0 (min b ((c + w) + ((K : ℚ) + 1) * w) - max a ((c + w) + (K : ℚ) * w))) := by
        apply List.map_congr_left
        intro K _
        simp only [Function.comp_apply, Nat.succ_eq_add_one]
        have e1 : c + (((K + 1 : ℕ) : ℚ) + 1) * w = (c + w) + ((K : ℚ) + 1) * w := by
          push_cast
          ring
        have e2 : c + ((K + 1 : ℕ) : ℚ) * w = (c + w) + (K : ℚ) * w := by
          push_cast
          ring
        rw [e1, e2]
      have hb' : b ≤ (c + w) + (N : ℚ) * w := by
        push_cast at hb
        linarith
      rw [htail, ih (c + w) a b hab hb']
      simp only [Nat.cast_zero, zero_mul, add_zero, zero_add, one_mul]
      -- goal: max 0 (min b (c+w) − max a c) + max 0 (b − max a (c+w)) = max 0 (b − max a c)
      rcases le_or_gt b (c + w) with hbc | hbc
      · rw [min_eq_left hbc]
        have hma : b - max a (c + w) ≤ 0 := by
          have := le_max_right a (c + w)
          linarith
        rw [max_eq_left hma, add_zero]
      · rw [min_eq_right hbc.le]
        rcases le_or_gt a (c + w) with hac | hac
        · have hma : max a (c + w) = c + w := max_eq_right hac
          rw [hma]
          have hmac : max a c ≤ c + w := max_le hac (by linarith)
          have h1 : max 0 (c + w - max a c) = c + w - max a c := max_eq_right (by linarith)
          have h2 : max 0 (b - (c + w)) = b - (c + w) := max_eq_right (by linarith)
          have h3 : max 0 (b - max a c) = b - max a c := by
            apply max_eq_right
            linarith
          rw [h1, h2, h3]
          ring
        · have hma : max a (c + w) = a := max_eq_left hac.le
          have hmac : max a c = a := max_eq_left (by linarith)
          rw [hma, hmac]
          have h1 : max 0 (c + w - a) = 0 := max_eq_left (by linarith)
          rw [h1, zero_add]

/-! ### K2 GLUE: the swap lemmas (phase-free sum plumbing)

The route from K1 to K2/K3 without ever sorting phase-dependent teeth:
`Σ_j length (inter L T_j) = length (inter L (T_0 ++ ⋯ ++ T_6))` by
`length_inter_append_right`, then flip to `length (inter (T_0 ++ ⋯) L)` by
`length_inter_comm` when the tiling side must drive the cursor.  Both are pure
list-sum rearrangements — no `Norm`, no ordering, no phases. -/

theorem length_inter_nil (B : Region) : length (inter B []) = 0 := by
  induction B with
  | nil => rfl
  | cons q B ih =>
      simp only [inter, List.flatMap_cons, List.map_nil, List.nil_append] at ih ⊢
      exact ih

/-- Intersection length is additive in the right argument's concatenation. -/
theorem length_inter_append_right (A B C : Region) :
    length (inter A (B ++ C)) = length (inter A B) + length (inter A C) := by
  induction A with
  | nil => simp [inter, length]
  | cons p A ih =>
      simp only [inter, List.flatMap_cons] at ih ⊢
      rw [List.map_append, length_append, length_append, length_append, length_append, ih]
      ring

/-- Clip is symmetric as a value. -/
theorem clip_comm (p q : ℚ × ℚ) : clip p q = clip q p := by
  unfold clip
  rw [max_comm, min_comm]

private theorem length_inter_cons_right (B : Region) (p : ℚ × ℚ) (A : Region) :
    length (inter B (p :: A))
      = length (B.map fun q => clip q p) + length (inter B A) := by
  induction B with
  | nil => simp [inter, length]
  | cons q B ih =>
      simp only [inter, List.flatMap_cons, List.map_cons] at ih ⊢
      rw [length_append]
      unfold length at ih ⊢
      simp only [List.map_cons, List.sum_cons, List.map_append, List.sum_append] at ih ⊢
      linarith

/-- Intersection length is symmetric — the Fubini swap for the quadratic intersection. -/
theorem length_inter_comm (A B : Region) : length (inter A B) = length (inter B A) := by
  induction A with
  | nil =>
      rw [length_inter_nil]
      rfl
  | cons p A ih =>
      rw [length_inter_cons_right B p A, ← ih]
      simp only [inter, List.flatMap_cons]
      rw [length_append]
      congr 1
      unfold length
      congr 1
      rw [List.map_map, List.map_map]
      apply List.map_congr_left
      intro q _
      simp only [Function.comp_apply]
      rw [clip_comm]

end RatIntervals
end LonelyRunner
