/-
  TournamentH7.LRCAggregatedDispersal -- THM-673(A) in Lean (klein-2026-07-09-S220).

  THE AGGREGATED DISPERSAL LEMMA: summed over all moduli q ∈ (V, 2V], the
  per-modulus resonance mass (total weight of relation values divisible by q)
  is at most K times the full weight, where every relation value is ≤ K·V:

      Σ_{q ∈ (V,2V]} resonanceMass J w q  ≤  K · Σ_{n ∈ J} w n.

  V-INDEPENDENT: each relation value N ≤ K·V has at most K divisors in the
  window (death-star's dvd_Ioc_card_le cofactor injection), so its weight is
  counted at most K times across all moduli.  This is the proved leg (A) of
  THM-673 (klein-S211: "each (j, m≠0) pins exactly one q"), the m ≠ 0 half of
  the aggregated modular supply route.

  Corollaries: THE THIN-MODULUS PIGEONHOLE (some window modulus carries at
  most the average mass -- with the total W₁ = K·Σw V-free, thin moduli abound
  as V grows), and the conditional composition socket into the certificate
  pipe (thin modulus + the off-line/off-peak transfer hypothesis ⟹ Mreach ≥
  1/14 via mreach_ge_of_bonf_pos) -- the Lean shape of THM-671 part 6(i),
  with the transfer as the named remaining mathematical hypothesis.

  Kernel-pure: no native_decide, no sorry.
-/

import Mathlib
import TournamentH7.LRCDiscreteBonferroni

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The per-modulus resonance mass: total weight of the relation values
divisible by `q`. -/
def resonanceMass (J : Finset ℕ) (w : ℕ → ℕ) (q : ℕ) : ℕ :=
  ∑ n ∈ J.filter (fun n => q ∣ n), w n

/-- **THE AGGREGATED DISPERSAL (THM-673(A))**: over the whole window the
resonance mass totals at most `K` times the full weight -- each relation value
is counted once per window divisor, and `dvd_Ioc_card_le` caps those at `K`. -/
theorem aggregated_dispersal (J : Finset ℕ) (w : ℕ → ℕ) (K V : ℕ)
    (hJ : ∀ n ∈ J, 0 < n ∧ n ≤ K * V) :
    ∑ q ∈ Finset.Ioc V (2 * V), resonanceMass J w q ≤ K * ∑ n ∈ J, w n := by
  classical
  have hswap : ∑ q ∈ Finset.Ioc V (2 * V), resonanceMass J w q
      = ∑ n ∈ J, ((Finset.Ioc V (2 * V)).filter fun q => q ∣ n).card * w n := by
    unfold resonanceMass
    calc ∑ q ∈ Finset.Ioc V (2 * V), ∑ n ∈ J.filter (fun n => q ∣ n), w n
        = ∑ q ∈ Finset.Ioc V (2 * V), ∑ n ∈ J, if q ∣ n then w n else 0 := by
          refine Finset.sum_congr rfl fun q _ => ?_
          rw [Finset.sum_filter]
      _ = ∑ n ∈ J, ∑ q ∈ Finset.Ioc V (2 * V), if q ∣ n then w n else 0 :=
          Finset.sum_comm
      _ = ∑ n ∈ J, ((Finset.Ioc V (2 * V)).filter fun q => q ∣ n).card * w n := by
          refine Finset.sum_congr rfl fun n _ => ?_
          rw [← Finset.sum_filter, Finset.sum_const, smul_eq_mul]
  rw [hswap]
  calc ∑ n ∈ J, ((Finset.Ioc V (2 * V)).filter fun q => q ∣ n).card * w n
      ≤ ∑ n ∈ J, K * w n := by
        refine Finset.sum_le_sum fun n hn => ?_
        exact Nat.mul_le_mul_right (w n)
          (dvd_Ioc_card_le n K V (hJ n hn).1 (hJ n hn).2)
    _ = K * ∑ n ∈ J, w n := by rw [Finset.mul_sum]

/-- **THE THIN-MODULUS PIGEONHOLE**: some modulus in the window carries at
most the average resonance mass: `V · mass(q) ≤ K · (total weight)`.  The
right side is V-free, so the per-modulus mass of the thin modulus vanishes
as the window grows -- the a-priori existence half of THM-671 part 6(i). -/
theorem exists_thin_modulus (J : Finset ℕ) (w : ℕ → ℕ) (K V : ℕ) (hV : 0 < V)
    (hJ : ∀ n ∈ J, 0 < n ∧ n ≤ K * V) :
    ∃ q ∈ Finset.Ioc V (2 * V),
      V * resonanceMass J w q ≤ K * ∑ n ∈ J, w n := by
  classical
  by_contra hcon
  push_neg at hcon
  have hne : (Finset.Ioc V (2 * V)).Nonempty := by
    refine ⟨V + 1, ?_⟩
    simp only [Finset.mem_Ioc]
    omega
  have hcard : (Finset.Ioc V (2 * V)).card = V := by
    rw [Nat.card_Ioc]
    omega
  have hstrict : ∑ q ∈ Finset.Ioc V (2 * V), (K * ∑ n ∈ J, w n)
      < ∑ q ∈ Finset.Ioc V (2 * V), V * resonanceMass J w q :=
    Finset.sum_lt_sum_of_nonempty hne fun q hq => hcon q hq
  rw [Finset.sum_const, hcard, smul_eq_mul, ← Finset.mul_sum] at hstrict
  have hagg := aggregated_dispersal J w K V hJ
  have := Nat.mul_le_mul_left V hagg
  omega

/-- **The composition socket (THM-671 part 6(i), conditional shape)**: if the
off-line/off-peak TRANSFER holds -- thin resonance mass at a window modulus
forces a positive odd-depth Bonferroni functional there -- then the window
certifies loneliness outright.  The transfer hypothesis is the named remaining
mathematical content (THM-677 Addenda 3-4, THM-680/681); everything else in
the chain is discharged. -/
theorem mreach_ge_of_thin_transfer (v : Fin 13 → ℤ) (D : ℕ) (hD : D % 2 = 1)
    (J : Finset ℕ) (w : ℕ → ℕ) (K V : ℕ) (hV : 0 < V)
    (hJ : ∀ n ∈ J, 0 < n ∧ n ≤ K * V)
    (htransfer : ∀ q ∈ Finset.Ioc V (2 * V),
      V * resonanceMass J w q ≤ K * ∑ n ∈ J, w n → 0 < bonf D v q) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  obtain ⟨q, hq, hthin⟩ := exists_thin_modulus J w K V hV hJ
  have hqpos : 0 < q := by
    have := (Finset.mem_Ioc.mp hq).1
    omega
  exact mreach_ge_of_bonf_pos v D hD q hqpos (htransfer q hq hthin)

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14Concrete.aggregated_dispersal
#print axioms LonelyRunner.LRC14Concrete.exists_thin_modulus
#print axioms LonelyRunner.LRC14Concrete.mreach_ge_of_thin_transfer
