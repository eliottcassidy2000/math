/-
  TournamentH7.LRCTightAPFreeRider — THE CRT TIGHT-AP FREE-RIDER (klein-2026-07-05-S132,
  HYP-4095; mac-mini-S46/S47 handoff).  The extremal half of hcomp's compressed peel:
  when the peeled base is the dilated AP `c·{1,…,12}` (the ENTIRE n=13 tight locus,
  mac-mini-S47), the killer rides free at some base optimum `k/(13c)`.

  THE MECHANISM: at `t = k/(13c)` with `13 ∤ k`, every base runner `c·j` has residue
  `c·(jk mod 13) ∈ [c, 12c]` — margin `1/13`, the base optimum.  The killer needs its
  residue `v*·k mod 13c` in the same band; primitivity (`gcd(c, v*) = 1`, supplied by
  `hcomp_of_primitive`) makes the band REACHABLE:
   * `13 ∤ v*`: `v*` is invertible mod `13c` (Bezout), aim at `m* = c` (or `c+1` when
     `13 ∣ c`, keeping `13 ∤ m*` and hence `13 ∤ k`);
   * `13 ∣ v*` (`v* = 13u`): the Bezout certificate of `gcd(c, 13u) = 1` splits into
     `IsCoprime c u` and `IsCoprime 13 c`; steer `u·k ≡ s (mod c)` with
     `13s ∈ [c, c+12] ⊆ [c, 12c]` and pin `k ≡ 1 (mod 13)` by explicit CRT.
  Kernel-pure: `lonely14_of_ratio` + per-runner `residue_key`.
-/
import TournamentH7.LRCOneSwapLadders
import TournamentH7.LRCSpread13

namespace LonelyRunner
namespace TightAPFreeRider

open LRC14

/-- The canonical tight-base-plus-killer family: `c·{1,…,12}` and the killer `v*`. -/
def apKiller (c vstar : ℤ) : Fin 13 → ℤ :=
  fun i => if (i : ℕ) < 12 then c * ((i : ℕ) + 1) else vstar

/-- **Killer targeting**: primitivity steers the killer's residue mod `13c` into the
base band `[c, 12c]` at a base optimum (`13 ∤ k`). -/
theorem killer_target (c vstar : ℤ) (hc : 2 ≤ c) (hgcd : Int.gcd c vstar = 1) :
    ∃ k r qq : ℤ, ¬ ((13 : ℤ) ∣ k) ∧ vstar * k = qq * (13 * c) + r ∧
      c ≤ r ∧ r ≤ 12 * c := by
  by_cases h13 : (13 : ℤ) ∣ vstar
  · -- v* = 13u
    obtain ⟨u, rfl⟩ := h13
    obtain ⟨x, y, hxy⟩ := Int.isCoprime_iff_gcd_eq_one.mpr hgcd
    -- hxy : x * c + y * (13 * u) = 1
    have hcu : IsCoprime c u := ⟨x, 13 * y, by linear_combination hxy⟩
    have hc13 : IsCoprime (13 : ℤ) c := ⟨y * u, x, by linear_combination hxy⟩
    set s : ℤ := (c + 12) / 13 with hs
    have hmod := Int.emod_nonneg (c + 12) (by norm_num : (13:ℤ) ≠ 0)
    have hmod' := Int.emod_lt_of_pos (c + 12) (by norm_num : (0:ℤ) < 13)
    have hde := Int.ediv_add_emod (c + 12) 13
    have hs13 : c ≤ 13 * s ∧ 13 * s ≤ c + 12 := by omega
    obtain ⟨a, b, hab⟩ := hcu
    -- hab : a * c + b * u = 1
    obtain ⟨x₁, y₁, hxy₁⟩ := hc13
    -- hxy₁ : x₁ * 13 + y₁ * c = 1
    refine ⟨y₁ * c + 13 * x₁ * (b * s), 13 * s,
      u * y₁ - 13 * x₁ * a * s - y₁ * s, ?_, ?_, by omega, by omega⟩
    · rintro ⟨d, hd⟩
      have h131 : (13 : ℤ) ∣ 1 := ⟨x₁ + d - x₁ * (b * s), by linear_combination hd - hxy₁⟩
      norm_num at h131
    · linear_combination (169 * x₁ * s) * hab + (13 * s) * hxy₁
  · -- 13 ∤ v*: v* invertible mod 13c
    have hvne : vstar ≠ 0 := by
      intro h0
      rw [h0, Int.gcd_zero_right] at hgcd
      omega
    have hg13 : Int.gcd 13 vstar = 1 := by
      have hp : Nat.Prime 13 := by norm_num
      have hdvd' : (Int.gcd 13 vstar) ∣ 13 := by
        have h := Int.gcd_dvd_left (a := (13 : ℤ)) (b := vstar)
        exact_mod_cast h
      rcases (Nat.Prime.eq_one_or_self_of_dvd hp _ hdvd') with h1 | h13g
      · exact h1
      · exfalso
        apply h13
        have := Int.gcd_dvd_right (a := (13:ℤ)) (b := vstar)
        rw [h13g] at this
        exact_mod_cast this
    have hco : IsCoprime vstar (13 * c) := by
      apply IsCoprime.mul_right
      · exact (Int.isCoprime_iff_gcd_eq_one.mpr hg13).symm
      · exact (Int.isCoprime_iff_gcd_eq_one.mpr hgcd).symm
    obtain ⟨x, y, hxy⟩ := hco
    -- hxy : x * vstar + y * (13 * c) = 1
    set mstar : ℤ := if (13 : ℤ) ∣ c then c + 1 else c with hmstar
    have hm13 : ¬ ((13 : ℤ) ∣ mstar) := by
      rw [hmstar]
      split_ifs with hdc
      · rintro ⟨e, he⟩
        obtain ⟨f, hf⟩ := hdc
        have : (13 : ℤ) ∣ 1 := ⟨e - f, by omega⟩
        norm_num at this
      · exact hdc
    have hmrange : c ≤ mstar ∧ mstar ≤ 12 * c := by
      rw [hmstar]
      split_ifs <;> omega
    refine ⟨x * mstar, mstar, -(y * mstar), ?_, ?_, hmrange.1, hmrange.2⟩
    · rintro ⟨d, hd⟩
      apply hm13
      exact ⟨vstar * d + y * mstar * c, by linear_combination vstar * hd - mstar * hxy⟩
    · linear_combination mstar * hxy

/-- **THE TIGHT-AP FREE-RIDER**: the dilated-AP base `c·{1,…,12}` plus any killer
attached primitively (`gcd(c, v*) = 1`) is lonely at a base optimum `k/(13c)` — the
base runners sit at margin `1/13`, the killer is steered into the same band by CRT.
The extremal (tight-locus) half of hcomp's compressed peel, pure arithmetic. -/
theorem tight_ap_free_rider (c vstar : ℤ) (hc : 2 ≤ c) (hgcd : Int.gcd c vstar = 1) :
    ∃ t : ℝ, Lonely 14 (apKiller c vstar) t := by
  obtain ⟨k, r, qq, hk13, hEqK, hrlo, hrhi⟩ := killer_target c vstar hc hgcd
  refine ⟨_, lonely14_of_ratio (apKiller c vstar) k (13 * c) (by omega) ?_⟩
  intro i m
  by_cases hi : (i : ℕ) < 12
  · have hval : apKiller c vstar i = c * (((i : ℕ) : ℤ) + 1) := by
      simp [apKiller, hi]
    rw [hval]
    set j : ℤ := ((i : ℕ) : ℤ) + 1 with hj
    have hj1 : 1 ≤ j := by omega
    have hj12 : j ≤ 12 := by omega
    have hp13 : Prime (13 : ℤ) := by
      rw [Int.prime_iff_natAbs_prime]
      norm_num
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
  · have hval : apKiller c vstar i = vstar := by
      simp [apKiller, hi]
    rw [hval]
    exact residue_key k (13 * c) c m vstar qq r (by omega) (by omega) hEqK hrlo (by omega)

#print axioms killer_target
#print axioms tight_ap_free_rider

end TightAPFreeRider
end LonelyRunner
