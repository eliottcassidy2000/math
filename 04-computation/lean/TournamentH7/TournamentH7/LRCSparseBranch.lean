/-
  TournamentH7.LRCSparseBranch — THE SPARSE BÉZOUT-BRANCH DECOMPOSITION
  (death-star-2026-07-17-S49, HYP-7245; the 29 sparse pairs of {1,…,13}).

  For a reduced-ratio pair `(g·i', g·j')` with `14 ≤ i' + j' ≤ 27` (the sparse
  regime — every sparse pair of the canonical family has `i'+j' ≤ 25`), the
  joint-failure set splits by the Bézout residue `k = j'·w_a − i'·w_b`:

  * `witness_unique` — a failing runner's witness is UNIQUE (two witnesses
    would sit `< q/7 < q` apart in `w·q`-space), so `k` is a well-defined
    function of the jointly-failing multiplier and the three branches
    `k ∈ {−1, 0, +1}` (THM-963's branch bound) PARTITION the joint-fail set.
  * `branch_zero_iff` / `branch_zero_count` — the `k = 0` branch is EXACTLY
    the GCD-speed narrow band (the locked-ray analysis needs no `i'+j' ≤ 13`
    once `k = 0` is imposed): count `2·⌊(q−1)/(14·max(i',j'))⌋`.
  * `branch_one_iff` — the `k = 1` branch, given a Bézout pair
    `j'·X₀ − i'·Y₀ = 1`: it is EXACTLY the two-band condition on the single
    integer `Z = g·p − t·q`:
    `∃t, 14|i'·Z − X₀·q| < q ∧ 14|j'·Z − Y₀·q| < q`
    — the coprimality move rewrites `(w_a, w_b) = (X₀ + i'·t, Y₀ + j'·t)`.
  * `branch_no_qmultiple` — the branch interval NEVER contains a multiple of
    `q`: `Z = m·q` in both bands would force `j'·X₀ − i'·Y₀ = 0 ≠ 1`.  So the
    `p ↦ Z` correspondence has no `p ≡ 0` boundary correction, and the branch
    count is a PURE lattice-interval count (the floor-form closed forms of the
    recon are exactly the integers of one rational interval).

  Recon (`sparsebranch_recon_deathstar_S49.out`): partition, `N⁰` formula,
  and the mirror `N⁺ = N⁻` PASS on 174/174 (pair, q) cases over all 29 sparse
  pairs; the limit law `N/(q−1) → 1/(7·max) + 2(i'+j'−14)/(14·i'·j')` is
  EXACT (0 error) at `q ≡ 1 (mod 14·i'·j')`, and reproduces boxeph LEM-044's
  μ-values (pair (8,9): 56/3528 + 21/3528 = 77/3528 ✓; pair (7,8) → 1/49,
  their `7 ∣ k` zero-excess law ✓).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCRationalLock

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **Witness uniqueness**: a failing runner has exactly one band witness. -/
theorem witness_unique (v w w' : ℤ) (q p : ℕ) (hq : 0 < q)
    (h : 14 * |v * (p : ℤ) - w * q| < q)
    (h' : 14 * |v * (p : ℤ) - w' * q| < q) :
    w = w' := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have htri : |(w - w') * (q : ℤ)|
      ≤ |v * (p : ℤ) - w' * q| + |v * (p : ℤ) - w * q| := by
    have h1 : (w - w') * (q : ℤ)
        = (v * (p : ℤ) - w' * q) - (v * (p : ℤ) - w * q) := by ring
    rw [h1]
    exact abs_sub _ _
  have hfact : |(w - w') * (q : ℤ)| = |w - w'| * q := by
    rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))]
  rw [hfact] at htri
  by_contra hne
  have h1 : 1 ≤ |w - w'| := Int.one_le_abs (sub_ne_zero.mpr hne)
  nlinarith [htri, h, h']

/-- **The `k = 0` branch is the GCD-speed narrow band** — the locked-ray
equivalence holds for ANY reduced pair once `k = 0` is imposed; no
`i' + j' ≤ 13` needed. -/
theorem branch_zero_iff (g : ℤ) (i' j' : ℕ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hcop : Nat.Coprime i' j') :
    (∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
        14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb = 0)
      ↔ (∃ s : ℤ, 14 * (max i' j' : ℕ) * |g * (p : ℤ) - s * q| < q) := by
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (0 : ℤ) < (j' : ℤ) := by exact_mod_cast hj
  constructor
  · rintro ⟨wa, wb, ha, hb, hk⟩
    have hlock : (j' : ℤ) * wa = (i' : ℤ) * wb := by linarith
    have hdvd : (i' : ℤ) ∣ (j' : ℤ) * wa := ⟨wb, hlock⟩
    have hcopZ : IsCoprime (i' : ℤ) (j' : ℤ) := by
      rw [Int.isCoprime_iff_gcd_eq_one]
      simpa using hcop
    obtain ⟨s, hs⟩ := hcopZ.dvd_of_dvd_mul_left hdvd
    refine ⟨s, ?_⟩
    have hXa : (g * i') * (p : ℤ) - wa * q = (i' : ℤ) * (g * p - s * q) := by
      rw [hs]; ring
    rw [hXa, abs_mul, abs_of_pos hiZ] at ha
    have hwb : wb = (j' : ℤ) * s := by
      have h1 : (i' : ℤ) * wb = (i' : ℤ) * ((j' : ℤ) * s) := by
        calc (i' : ℤ) * wb = (j' : ℤ) * wa := hlock.symm
          _ = (j' : ℤ) * ((i' : ℤ) * s) := by rw [hs]
          _ = (i' : ℤ) * ((j' : ℤ) * s) := by ring
      exact mul_left_cancel₀ (by linarith) h1
    have hXb : (g * j') * (p : ℤ) - wb * q = (j' : ℤ) * (g * p - s * q) := by
      rw [hwb]; ring
    rw [hXb, abs_mul, abs_of_pos hjZ] at hb
    rcases max_cases i' j' with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm]
    · calc 14 * (i' : ℤ) * |g * (p : ℤ) - s * q|
          = 14 * ((i' : ℤ) * |g * (p : ℤ) - s * q|) := by ring
        _ < q := ha
    · calc 14 * (j' : ℤ) * |g * (p : ℤ) - s * q|
          = 14 * ((j' : ℤ) * |g * (p : ℤ) - s * q|) := by ring
        _ < q := hb
  · rintro ⟨s, hs⟩
    have habs : (0 : ℤ) ≤ |g * (p : ℤ) - s * q| := abs_nonneg _
    have hmax_i : (i' : ℤ) ≤ (max i' j' : ℕ) := by
      have := le_max_left i' j'
      exact_mod_cast this
    have hmax_j : (j' : ℤ) ≤ (max i' j' : ℕ) := by
      have := le_max_right i' j'
      exact_mod_cast this
    refine ⟨(i' : ℤ) * s, (j' : ℤ) * s, ?_, ?_, by ring⟩
    · have hXa : (g * i') * (p : ℤ) - ((i' : ℤ) * s) * q
          = (i' : ℤ) * (g * p - s * q) := by ring
      rw [hXa, abs_mul, abs_of_pos hiZ]
      nlinarith [hs, habs]
    · have hXb : (g * j') * (p : ℤ) - ((j' : ℤ) * s) * q
          = (j' : ℤ) * (g * p - s * q) := by ring
      rw [hXb, abs_mul, abs_of_pos hjZ]
      nlinarith [hs, habs]

open Classical in
/-- **The `k = 0` branch count**: `2·⌊(q−1)/(14·max(i',j'))⌋` at coprime
moduli — for EVERY reduced pair, locked or sparse. -/
theorem branch_zero_count (g : ℤ) (i' j' : ℕ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hcop : Nat.Coprime i' j')
    (hgcd : Int.gcd g (q : ℤ) = 1) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
        ∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
          14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
          (j' : ℤ) * wa - (i' : ℤ) * wb = 0).card
      = 2 * ((q - 1) / (14 * max i' j')) := by
  have hM : 1 ≤ max i' j' := le_trans hi (le_max_left i' j')
  have hfilter : ((Finset.Ioo 0 q).filter fun p : ℕ =>
      ∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
        14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb = 0)
      = ((Finset.Ioo 0 q).filter fun p : ℕ =>
          14 * (max i' j') * ((g * (p : ℤ) % (q : ℤ)).toNat) < q ∨
          14 * (max i' j') * (q - ((g * (p : ℤ) % (q : ℤ)).toNat)) < q) := by
    apply Finset.filter_congr
    intro p _
    rw [branch_zero_iff g i' j' q p hi hj hcop,
      witness_mod_bridge (g * (p : ℤ)) (max i' j') q hM hq]
  rw [hfilter]
  have htrans := card_mod_filter_eq g q hq hgcd
    (fun r => 14 * (max i' j') * r < q ∨ 14 * (max i' j') * (q - r) < q)
  rw [htrans]
  exact narrowFailN_count (max i' j') q hM hq

/-- **The `k = 1` branch in Bézout normal form**: with `j'·X₀ − i'·Y₀ = 1`,
the branch is exactly the two-band condition on `Z = g·p − t·q`. -/
theorem branch_one_iff (g : ℤ) (i' j' : ℕ) (X₀ Y₀ : ℤ) (q p : ℕ)
    (hi : 1 ≤ i') (hj : 1 ≤ j') (hcop : Nat.Coprime i' j')
    (hbez : (j' : ℤ) * X₀ - (i' : ℤ) * Y₀ = 1) :
    (∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
        14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
        (j' : ℤ) * wa - (i' : ℤ) * wb = 1)
      ↔ (∃ t : ℤ, 14 * |(i' : ℤ) * (g * (p : ℤ) - t * q) - X₀ * q| < q ∧
          14 * |(j' : ℤ) * (g * (p : ℤ) - t * q) - Y₀ * q| < q) := by
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (0 : ℤ) < (j' : ℤ) := by exact_mod_cast hj
  have hcopZ : IsCoprime (i' : ℤ) (j' : ℤ) := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    simpa using hcop
  constructor
  · rintro ⟨wa, wb, ha, hb, hk⟩
    -- (wa − X₀, wb − Y₀) solves the homogeneous equation ⟹ = t·(i', j')
    have hhom : (j' : ℤ) * (wa - X₀) = (i' : ℤ) * (wb - Y₀) := by
      have h1 : (j' : ℤ) * wa - (i' : ℤ) * wb = (j' : ℤ) * X₀ - (i' : ℤ) * Y₀ := by
        rw [hbez]; exact hk
      linarith [h1]
    have hdvd : (i' : ℤ) ∣ (j' : ℤ) * (wa - X₀) := ⟨wb - Y₀, hhom⟩
    obtain ⟨t, ht⟩ := hcopZ.dvd_of_dvd_mul_left hdvd
    have hwb : wb - Y₀ = (j' : ℤ) * t := by
      have h1 : (i' : ℤ) * (wb - Y₀) = (i' : ℤ) * ((j' : ℤ) * t) := by
        calc (i' : ℤ) * (wb - Y₀) = (j' : ℤ) * (wa - X₀) := hhom.symm
          _ = (j' : ℤ) * ((i' : ℤ) * t) := by rw [ht]
          _ = (i' : ℤ) * ((j' : ℤ) * t) := by ring
      exact mul_left_cancel₀ (by linarith) h1
    refine ⟨t, ?_, ?_⟩
    · have hXa : (i' : ℤ) * (g * (p : ℤ) - t * q) - X₀ * q
          = (g * i') * (p : ℤ) - wa * q := by
        have hwa : wa = X₀ + (i' : ℤ) * t := by linarith [ht]
        rw [hwa]; ring
      rw [hXa]
      exact ha
    · have hXb : (j' : ℤ) * (g * (p : ℤ) - t * q) - Y₀ * q
          = (g * j') * (p : ℤ) - wb * q := by
        have hwb' : wb = Y₀ + (j' : ℤ) * t := by linarith [hwb]
        rw [hwb']; ring
      rw [hXb]
      exact hb
  · rintro ⟨t, ha, hb⟩
    refine ⟨X₀ + (i' : ℤ) * t, Y₀ + (j' : ℤ) * t, ?_, ?_, by linarith [hbez]⟩
    · have hXa : (g * i') * (p : ℤ) - (X₀ + (i' : ℤ) * t) * q
          = (i' : ℤ) * (g * (p : ℤ) - t * q) - X₀ * q := by ring
      rw [hXa]
      exact ha
    · have hXb : (g * j') * (p : ℤ) - (Y₀ + (j' : ℤ) * t) * q
          = (j' : ℤ) * (g * (p : ℤ) - t * q) - Y₀ * q := by ring
      rw [hXb]
      exact hb

/-- **The branch interval avoids multiples of `q`**: `Z = m·q` inside both
Bézout bands would force `j'·X₀ − i'·Y₀ = 0 ≠ 1`.  Hence the `p ↦ Z`
correspondence on the `k = 1` branch needs no `p ≡ 0` boundary correction. -/
theorem branch_no_qmultiple (i' j' : ℕ) (X₀ Y₀ m : ℤ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hj : 1 ≤ j')
    (hbez : (j' : ℤ) * X₀ - (i' : ℤ) * Y₀ = 1)
    (ha : 14 * |(i' : ℤ) * (m * q) - X₀ * q| < q)
    (hb : 14 * |(j' : ℤ) * (m * q) - Y₀ * q| < q) :
    False := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- both bands at Z = m·q force i'·m = X₀ and j'·m = Y₀
  have hfa : (i' : ℤ) * (m * q) - X₀ * q = ((i' : ℤ) * m - X₀) * q := by ring
  have hfb : (j' : ℤ) * (m * q) - Y₀ * q = ((j' : ℤ) * m - Y₀) * q := by ring
  rw [hfa, abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))] at ha
  rw [hfb, abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ))] at hb
  have hXa : (i' : ℤ) * m - X₀ = 0 := by
    by_contra hne
    have h1 : 1 ≤ |(i' : ℤ) * m - X₀| := Int.one_le_abs hne
    nlinarith [ha]
  have hXb : (j' : ℤ) * m - Y₀ = 0 := by
    by_contra hne
    have h1 : 1 ≤ |(j' : ℤ) * m - Y₀| := Int.one_le_abs hne
    nlinarith [hb]
  -- substitute into the Bézout identity: 1 = j'·i'·m − i'·j'·m = 0
  have hX : X₀ = (i' : ℤ) * m := by linarith
  have hY : Y₀ = (j' : ℤ) * m := by linarith
  rw [hX, hY] at hbez
  have : (0 : ℤ) = 1 := by linarith [hbez]
  omega

open Classical in
/-- **Mirror branch symmetry**: the involution `p ↦ q − p` negates the Bézout
residue (`w ↦ g·v' − w` on each witness), so the `k = 1` and `k = −1` branches
have EQUAL cardinality — the total sparse count is `N⁰ + 2·N⁺`. -/
theorem branch_mirror_card (g : ℤ) (i' j' : ℕ) (q : ℕ) :
    ((Finset.Ioo 0 q).filter fun p : ℕ =>
        ∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
          14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
          (j' : ℤ) * wa - (i' : ℤ) * wb = 1).card
      = ((Finset.Ioo 0 q).filter fun p : ℕ =>
        ∃ wa wb : ℤ, 14 * |(g * i') * (p : ℤ) - wa * q| < q ∧
          14 * |(g * j') * (p : ℤ) - wb * q| < q ∧
          (j' : ℤ) * wa - (i' : ℤ) * wb = -1).card := by
  apply Finset.card_bij (fun p _ => q - p)
  · -- maps into the k = −1 filter
    intro p hp
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp
    obtain ⟨⟨hp0, hpq⟩, wa, wb, ha, hb, hk⟩ := hp
    rw [Finset.mem_filter, Finset.mem_Ioo]
    have hcast : ((q - p : ℕ) : ℤ) = (q : ℤ) - p := by
      push_cast [le_of_lt hpq]
      ring
    refine ⟨⟨by omega, by omega⟩, g * i' - wa, g * j' - wb, ?_, ?_, by linarith [hk]⟩
    · have hrw : (g * i') * ((q - p : ℕ) : ℤ) - (g * i' - wa) * q
          = -((g * i') * (p : ℤ) - wa * q) := by
        rw [hcast]; ring
      rw [hrw, abs_neg]
      exact ha
    · have hrw : (g * j') * ((q - p : ℕ) : ℤ) - (g * j' - wb) * q
          = -((g * j') * (p : ℤ) - wb * q) := by
        rw [hcast]; ring
      rw [hrw, abs_neg]
      exact hb
  · -- injective
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    omega
  · -- surjective
    intro r hr
    rw [Finset.mem_filter, Finset.mem_Ioo] at hr
    obtain ⟨⟨hr0, hrq⟩, wa, wb, ha, hb, hk⟩ := hr
    refine ⟨q - r, ?_, by omega⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    have hcast : ((q - r : ℕ) : ℤ) = (q : ℤ) - r := by
      push_cast [le_of_lt hrq]
      ring
    refine ⟨⟨by omega, by omega⟩, g * i' - wa, g * j' - wb, ?_, ?_, by linarith [hk]⟩
    · have hrw : (g * i') * ((q - r : ℕ) : ℤ) - (g * i' - wa) * q
          = -((g * i') * (r : ℤ) - wa * q) := by
        rw [hcast]; ring
      rw [hrw, abs_neg]
      exact ha
    · have hrw : (g * j') * ((q - r : ℕ) : ℤ) - (g * j' - wb) * q
          = -((g * j') * (r : ℤ) - wb * q) := by
        rw [hcast]; ring
      rw [hrw, abs_neg]
      exact hb

/-! ## Axiom audit -/
#print axioms witness_unique
#print axioms branch_zero_iff
#print axioms branch_zero_count
#print axioms branch_one_iff
#print axioms branch_no_qmultiple
#print axioms branch_mirror_card

end LRC14Concrete
end LonelyRunner
