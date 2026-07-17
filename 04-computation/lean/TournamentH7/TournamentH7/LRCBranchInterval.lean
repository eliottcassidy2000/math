/-
  TournamentH7.LRCBranchInterval — THE GENERAL WITNESS CROSS BOUND AND THE
  BRANCH INTERVAL COUNT (death-star-2026-07-17-S50, HYP-7255).

  Two deliverables, one from free-thinking over an old thread, one the last
  pair-layer Lean piece:

  * `witness_cross_bound` — [unification, from re-mining THM-949's witness
    ladder] for ANY two positive speeds — NO ratio structure whatsoever —
    joint failure bounds the witness cross: `14·|w_a·v_b − w_b·v_a| < v_a + v_b`.
    The exact identity `v_b·X_a − v_a·X_b = (v_a·w_b − v_b·w_a)·q` needs no
    hypothesis at all; the ladder (n_j ≥ 3·n_i on ratio [3,13]) was the
    inequality shadow of this, and THM-966/967's locks are its structured
    sharpenings.
  * `cross_lock_of_sum_le13` — corollary: `v_a + v_b ≤ 13` locks
    `v_a·w_b = v_b·w_a` OUTRIGHT — an unstructured lock for small speed sums
    (for the canonical family: every pair with `a + b ≤ 13`).
  * `branch_Z_mem_Icc` — the k = 1 branch interval in closed form: with
    `j'·X₀ − i'·Y₀ = 1` and `i', j' ≤ 13`, the two Bézout bands on `Z` hold
    IFF `A·q/C + 1 ≤ Z ≤ (B·q − 1)/D` where `A = 14X₀−1, B = 14Y₀+1,
    C = 14i', D = 14j'` (ℤ ediv = floor).  The two non-binding constraints are
    DOMINATED via the Bézout identity alone:
    `(14X₀−1)·j' − (14Y₀−1)·i' = 14 + i' − j' ≥ 0` and
    `(14X₀+1)·j' − (14Y₀+1)·i' = 14 + j' − i' ≥ 0`.
  * `branch_interval_card` — hence the branch count is
    `((B·q−1)/D − A·q/C).toNat` — matching the S49 recon formula verified
    232/232 across all 29 sparse pairs of {1,…,13}.  The pair layer of the
    canonical family is now Lean-closed at the interval level.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCSparseBranch

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **THE GENERAL WITNESS CROSS BOUND**: any two positive speeds, no ratio
structure — joint failure bounds the witness cross by the speed sum. -/
theorem witness_cross_bound (va vb wa wb : ℤ) (q p : ℕ) (hq : 0 < q)
    (hva : 0 < va) (hvb : 0 < vb)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q) :
    14 * |wa * vb - wb * va| < va + vb := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- the exact identity: v_b·X_a − v_a·X_b = (v_a·w_b − v_b·w_a)·q
  have hkey : (va * wb - vb * wa) * (q : ℤ)
      = vb * (va * (p : ℤ) - wa * q) - va * (vb * (p : ℤ) - wb * q) := by ring
  have habs : |va * wb - vb * wa| * (q : ℤ)
      ≤ vb * |va * (p : ℤ) - wa * q| + va * |vb * (p : ℤ) - wb * q| := by
    have h1 : |(va * wb - vb * wa) * (q : ℤ)|
        ≤ |vb * (va * (p : ℤ) - wa * q)| + |va * (vb * (p : ℤ) - wb * q)| := by
      rw [hkey]
      exact abs_sub _ _
    rw [abs_mul, abs_of_nonneg (by positivity : (0 : ℤ) ≤ (q : ℤ)),
      abs_mul, abs_mul, abs_of_pos hvb, abs_of_pos hva] at h1
    exact h1
  have hscale : 14 * (|va * wb - vb * wa| * (q : ℤ)) < (va + vb) * q := by
    have h1 : vb * (14 * |va * (p : ℤ) - wa * q|) < vb * q :=
      mul_lt_mul_of_pos_left ha hvb
    have h2 : va * (14 * |vb * (p : ℤ) - wb * q|) < va * q :=
      mul_lt_mul_of_pos_left hb hva
    nlinarith [habs]
  have hcomm : |wa * vb - wb * va| = |va * wb - vb * wa| := by
    rw [← abs_neg]
    congr 1
    ring
  rw [hcomm]
  by_contra hcon
  push Not at hcon
  nlinarith [hscale, hqZ, hcon]

/-- **The unstructured small-sum lock**: `v_a + v_b ≤ 13` forces exact witness
proportionality — no ratio hypothesis at all. -/
theorem cross_lock_of_sum_le13 (va vb wa wb : ℤ) (q p : ℕ) (hq : 0 < q)
    (hva : 0 < va) (hvb : 0 < vb) (hsum : va + vb ≤ 13)
    (ha : 14 * |va * (p : ℤ) - wa * q| < q)
    (hb : 14 * |vb * (p : ℤ) - wb * q| < q) :
    va * wb = vb * wa := by
  have hcross := witness_cross_bound va vb wa wb q p hq hva hvb ha hb
  have hk0 : |wa * vb - wb * va| = 0 := by
    by_contra hne
    have h1 : 1 ≤ |wa * vb - wb * va| := by
      rcases Nat.eq_zero_or_pos (|wa * vb - wb * va|).toNat with h0 | hpos
      · exfalso
        apply hne
        have := Int.toNat_of_nonneg (abs_nonneg (wa * vb - wb * va))
        omega
      · have := Int.toNat_of_nonneg (abs_nonneg (wa * vb - wb * va))
        omega
    linarith [hcross, hsum]
  have : wa * vb - wb * va = 0 := abs_eq_zero.mp hk0
  linarith [this]

/-- **The branch interval, closed Icc form**: with the Bézout identity and
`i', j' ≤ 13`, the two branch bands on `Z` hold iff
`A·q/C + 1 ≤ Z ≤ (B·q − 1)/D` (ℤ ediv; `A = 14X₀−1, B = 14Y₀+1, C = 14i',
D = 14j'`).  The non-binding constraints are dominated via Bézout alone. -/
theorem branch_Z_mem_Icc (i' j' : ℕ) (X₀ Y₀ Z : ℤ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hi13 : i' ≤ 13) (hj : 1 ≤ j') (hj13 : j' ≤ 13)
    (hbez : (j' : ℤ) * X₀ - (i' : ℤ) * Y₀ = 1) :
    (14 * |(i' : ℤ) * Z - X₀ * q| < q ∧ 14 * |(j' : ℤ) * Z - Y₀ * q| < q)
      ↔ ((14 * X₀ - 1) * q / (14 * i') + 1 ≤ Z ∧
         Z ≤ ((14 * Y₀ + 1) * q - 1) / (14 * j')) := by
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hiZ : (0 : ℤ) < (i' : ℤ) := by exact_mod_cast hi
  have hjZ : (0 : ℤ) < (j' : ℤ) := by exact_mod_cast hj
  have hi13Z : (i' : ℤ) ≤ 13 := by exact_mod_cast hi13
  have hj13Z : (j' : ℤ) ≤ 13 := by exact_mod_cast hj13
  have hCpos : (0 : ℤ) < 14 * (i' : ℤ) := by linarith
  have hDpos : (0 : ℤ) < 14 * (j' : ℤ) := by linarith
  -- unfold the four strict inequalities of the two abs bands
  have habs_iff : ∀ (c w : ℤ), 0 < c →
      (14 * |c * Z - w * q| < q ↔
        (14 * w - 1) * q < 14 * c * Z ∧ 14 * c * Z < (14 * w + 1) * q) := by
    intro c w hc
    constructor
    · intro h
      rcases abs_cases (c * Z - w * q) with ⟨he, _⟩ | ⟨he, _⟩ <;>
        rw [he] at h <;> constructor <;> nlinarith [h]
    · rintro ⟨h1, h2⟩
      rcases abs_cases (c * Z - w * q) with ⟨he, _⟩ | ⟨he, _⟩ <;>
        rw [he] <;> nlinarith [h1, h2]
  rw [habs_iff (i' : ℤ) X₀ hiZ, habs_iff (j' : ℤ) Y₀ hjZ]
  constructor
  · rintro ⟨⟨ha1, _⟩, ⟨_, hb2⟩⟩
    constructor
    · -- (14X₀−1)q < 14i'Z ⟹ Aq/C + 1 ≤ Z
      have hlt : (14 * X₀ - 1) * q / (14 * i') < Z := by
        rw [Int.ediv_lt_iff_lt_mul hCpos]
        calc (14 * X₀ - 1) * (q : ℤ) < 14 * i' * Z := ha1
          _ = Z * (14 * i') := by ring
      omega
    · -- 14j'Z < (14Y₀+1)q ⟹ Z ≤ ((14Y₀+1)q − 1)/D
      rw [Int.le_ediv_iff_mul_le hDpos]
      calc Z * (14 * (j' : ℤ)) = 14 * j' * Z := by ring
        _ ≤ (14 * Y₀ + 1) * q - 1 := by omega
  · rintro ⟨hlo, hhi⟩
    -- recover the two binding inequalities
    have ha1 : (14 * X₀ - 1) * (q : ℤ) < 14 * i' * Z := by
      have h1 : (14 * X₀ - 1) * q / (14 * i') < Z := by omega
      rw [Int.ediv_lt_iff_lt_mul hCpos] at h1
      calc (14 * X₀ - 1) * (q : ℤ) < Z * (14 * i') := h1
        _ = 14 * i' * Z := by ring
    have hb2 : 14 * (j' : ℤ) * Z < (14 * Y₀ + 1) * q := by
      have h1 := (Int.le_ediv_iff_mul_le hDpos).mp hhi
      calc 14 * (j' : ℤ) * Z = Z * (14 * j') := by ring
        _ ≤ (14 * Y₀ + 1) * q - 1 := h1
        _ < (14 * Y₀ + 1) * q := by omega
    -- dominate the two non-binding constraints via the Bézout identity
    refine ⟨⟨ha1, ?_⟩, ⟨?_, hb2⟩⟩
    · -- 14i'Z < (14X₀+1)q  ⟸  j'-side upper + (14X₀+1)j' ≥ (14Y₀+1)i'
      have hdom : ((14 * X₀ + 1) * (j' : ℤ)) = (14 * Y₀ + 1) * i' + (14 + j' - i') := by
        have h := hbez
        nlinarith [h]
      have h1 : (j' : ℤ) * (14 * i' * Z) < (j' : ℤ) * ((14 * X₀ + 1) * q) := by
        have h2 : (i' : ℤ) * (14 * j' * Z) < (i' : ℤ) * ((14 * Y₀ + 1) * q) :=
          mul_lt_mul_of_pos_left hb2 hiZ
        nlinarith [h2, hqZ]
      exact lt_of_mul_lt_mul_left h1 (by linarith)
    · -- (14Y₀−1)q < 14j'Z  ⟸  i'-side lower + (14X₀−1)j' ≥ (14Y₀−1)i'
      have h1 : (i' : ℤ) * (14 * j' * Z) > (i' : ℤ) * ((14 * Y₀ - 1) * q) := by
        have h2 : (j' : ℤ) * (14 * i' * Z) > (j' : ℤ) * ((14 * X₀ - 1) * q) :=
          mul_lt_mul_of_pos_left ha1 hjZ
        nlinarith [h2, hqZ, hbez]
      exact lt_of_mul_lt_mul_left h1 (by linarith)

/-- **The branch interval count**: the Bézout-band `Z`-set is the integer
interval `[A·q/C + 1, (B·q−1)/D]`, of size `((B·q−1)/D − A·q/C).toNat` —
the S49 recon formula (232/232), now in-kernel. -/
theorem branch_interval_card (i' j' : ℕ) (X₀ Y₀ : ℤ) (q : ℕ) (hq : 0 < q)
    (hi : 1 ≤ i') (hi13 : i' ≤ 13) (hj : 1 ≤ j') (hj13 : j' ≤ 13)
    (hbez : (j' : ℤ) * X₀ - (i' : ℤ) * Y₀ = 1) :
    ((Finset.Icc ((14 * X₀ - 1) * q / (14 * i') + 1)
        (((14 * Y₀ + 1) * q - 1) / (14 * j'))).filter fun Z : ℤ =>
      14 * |(i' : ℤ) * Z - X₀ * q| < q ∧ 14 * |(j' : ℤ) * Z - Y₀ * q| < q).card
      = ((((14 * Y₀ + 1) * q - 1) / (14 * j'))
          - ((14 * X₀ - 1) * q / (14 * i'))).toNat := by
  have hfull : ((Finset.Icc ((14 * X₀ - 1) * q / (14 * i') + 1)
      (((14 * Y₀ + 1) * q - 1) / (14 * j'))).filter fun Z : ℤ =>
      14 * |(i' : ℤ) * Z - X₀ * q| < q ∧ 14 * |(j' : ℤ) * Z - Y₀ * q| < q)
      = Finset.Icc ((14 * X₀ - 1) * q / (14 * i') + 1)
          (((14 * Y₀ + 1) * q - 1) / (14 * j')) := by
    apply Finset.filter_true_of_mem
    intro Z hZ
    rw [Finset.mem_Icc] at hZ
    exact (branch_Z_mem_Icc i' j' X₀ Y₀ Z q hq hi hi13 hj hj13 hbez).mpr hZ
  rw [hfull, Int.card_Icc]
  congr 1
  ring

/-! ## Axiom audit -/
#print axioms witness_cross_bound
#print axioms cross_lock_of_sum_le13
#print axioms branch_Z_mem_Icc
#print axioms branch_interval_card

end LRC14Concrete
end LonelyRunner
