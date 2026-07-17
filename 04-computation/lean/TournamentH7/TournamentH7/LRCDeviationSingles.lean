/-
  TournamentH7.LRCDeviationSingles — THE SINGLES RUNG of the deviation ledger
  (death-star-2026-07-17-S36, HYP-7173; the first unconditional bound in THM-940's
  discrete identification).

  THM-940 wrote `B5 = (q−1)·2052/16807 + Σ (−1)^k Σ_{|T|=k} D_T` exactly.  This module
  computes the k = 1 rung in closed form via the unit bijection `p ↦ (v i·p) mod q`:

  * `bandSize_eq` — the safe band is the integer interval `[(q+13)/14, 13q/14]`
    (ℕ-division; `omega` carries everything downstream).
  * `jointFail_singleton_eq` — for `gcd(v i, q) = 1` the singleton failure count is
    EXACT: `N_{i} = (q − 1) − bandSize q` (multiplication by a unit permutes the
    nonzero residues; `0` always fails the band).
  * `deviation_singleton_bounds` — for the gcd case at any `q ≥ 14`:
        −13/7 ≤ D_{i} ≤ 0.
    The singles NEVER threaten the equilibrium budget: their total ledger
    contribution lies in `[−169/7, 0]` — o(q).  **The deviation debt is carried
    entirely by |T| ≥ 2.**
  * `deviation_singleton_of_dvd` — at `14 ∣ q` every singleton deviation is the
    CONSTANT `−13/7` exactly (7·bandSize = 6q + 7 on the nose).

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCB5SubsetExpansion

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The number of safe-band residues mod `q`. -/
def bandSize (q : ℕ) : ℕ :=
  ((Finset.range q).filter fun r => q ≤ 14 * r ∧ 14 * r ≤ 13 * q).card

/-- The band is the integer interval `[(q+13)/14, 13q/14]` (ℕ-division). -/
theorem bandSize_eq (q : ℕ) (hq : 14 ≤ q) :
    bandSize q = 13 * q / 14 - (q + 13) / 14 + 1 := by
  unfold bandSize
  have hset : (Finset.range q).filter (fun r => q ≤ 14 * r ∧ 14 * r ≤ 13 * q)
      = Finset.Icc ((q + 13) / 14) (13 * q / 14) := by
    ext r
    simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Icc]
    omega
  rw [hset, Nat.card_Icc]
  omega

/-- **The unit bijection count**: for `gcd(v i, q) = 1`, the singleton joint-failure
count is exactly `(q − 1) − bandSize q`. -/
theorem jointFail_singleton_eq (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13)
    (hq : 0 < q) (hgcd : Int.gcd (v i) (q : ℤ) = 1) :
    jointFail v q {i} = (q - 1) - bandSize q := by
  set f : ℕ → ℕ := fun p => ((v i * (p : ℤ)) % (q : ℤ)).toNat with hf
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hmod_lt : ∀ p : ℕ, ((v i * (p : ℤ)) % (q : ℤ)) < (q : ℤ) :=
    fun p => Int.emod_lt_of_pos _ hqZ
  have hmod_nonneg : ∀ p : ℕ, (0 : ℤ) ≤ ((v i * (p : ℤ)) % (q : ℤ)) :=
    fun p => Int.emod_nonneg _ (by omega)
  have hfval : ∀ p : ℕ, ((f p : ℕ) : ℤ) = (v i * (p : ℤ)) % (q : ℤ) := by
    intro p
    rw [hf]
    exact Int.toNat_of_nonneg (hmod_nonneg p)
  have hband_iff : ∀ p : ℕ, inBand v q p i ↔ (q ≤ 14 * f p ∧ 14 * f p ≤ 13 * q) := by
    intro p
    unfold inBand
    constructor
    · intro ⟨h1, h2⟩
      constructor
      · have h1' : (q : ℤ) ≤ 14 * ((f p : ℕ) : ℤ) := by rw [hfval]; exact h1
        exact_mod_cast h1'
      · have h2' : 14 * ((f p : ℕ) : ℤ) ≤ 13 * (q : ℤ) := by rw [hfval]; exact h2
        exact_mod_cast h2'
    · intro ⟨h1, h2⟩
      constructor
      · have h1' : (q : ℤ) ≤ 14 * ((f p : ℕ) : ℤ) := by exact_mod_cast h1
        rwa [hfval] at h1'
      · have h2' : 14 * ((f p : ℕ) : ℤ) ≤ 13 * (q : ℤ) := by exact_mod_cast h2
        rwa [hfval] at h2'
  have hjf : jointFail v q {i}
      = ((Finset.Ioo 0 q).filter fun p => ¬ inBand v q p i).card := by
    unfold jointFail
    congr 1
    apply Finset.filter_congr
    intro p _
    simp
  rw [hjf]
  have htarget : ((Finset.Ioo 0 q).filter fun r =>
      ¬ (q ≤ 14 * r ∧ 14 * r ≤ 13 * q)).card = (q - 1) - bandSize q := by
    have hbandIoo : bandSize q
        = ((Finset.Ioo 0 q).filter fun r => q ≤ 14 * r ∧ 14 * r ≤ 13 * q).card := by
      unfold bandSize
      congr 1
      ext r
      simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Ioo]
      omega
    have hpart := Finset.filter_card_add_filter_neg_card_eq_card
      (s := Finset.Ioo 0 q) (p := fun r => q ≤ 14 * r ∧ 14 * r ≤ 13 * q)
    have hcard : (Finset.Ioo 0 q).card = q - 1 := by
      simp [Nat.card_Ioo]
    omega
  rw [← htarget]
  obtain ⟨a, b, hab⟩ := Int.isCoprime_iff_gcd_eq_one.mpr hgcd
  have hva : v i * a ≡ 1 [ZMOD (q : ℤ)] :=
    Int.modEq_iff_dvd.mpr ⟨b, by linarith⟩
  apply Finset.card_bij (fun p _ => f p)
  · -- into the target
    intro p hp
    rw [Finset.mem_filter] at hp
    obtain ⟨hpIoo, hpfail⟩ := hp
    rw [Finset.mem_Ioo] at hpIoo
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, fun hcon => hpfail ((hband_iff p).mpr hcon)⟩
    · rcases Nat.eq_zero_or_pos (f p) with h0 | hpos
      · exfalso
        have hzero : (v i * (p : ℤ)) % (q : ℤ) = 0 := by
          have hv := hfval p
          rw [h0] at hv
          simpa using hv.symm
        have hdvd : (q : ℤ) ∣ v i * (p : ℤ) := Int.dvd_of_emod_eq_zero hzero
        have hp_eq : (p : ℤ) = a * (v i * (p : ℤ)) + (b * (p : ℤ)) * (q : ℤ) := by
          linear_combination (-(p : ℤ)) * hab
        have hqp : (q : ℤ) ∣ (p : ℤ) := by
          rw [hp_eq]
          exact dvd_add (Dvd.dvd.mul_left hdvd a) (Dvd.intro_left _ rfl)
        have hppos : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hpIoo.1
        have hple := Int.le_of_dvd hppos hqp
        have hplt : (p : ℤ) < (q : ℤ) := by exact_mod_cast hpIoo.2
        omega
      · exact hpos
    · have h1 : ((f p : ℕ) : ℤ) < (q : ℤ) := by rw [hfval]; exact hmod_lt p
      exact_mod_cast h1
  · -- injective
    intro p₁ hp₁ p₂ hp₂ hEq
    rw [Finset.mem_filter, Finset.mem_Ioo] at hp₁ hp₂
    have hmodeq : v i * (p₁ : ℤ) ≡ v i * (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      show (v i * (p₁ : ℤ)) % (q : ℤ) = (v i * (p₂ : ℤ)) % (q : ℤ)
      have hc : ((f p₁ : ℕ) : ℤ) = ((f p₂ : ℕ) : ℤ) := by exact_mod_cast hEq
      rw [hfval, hfval] at hc
      exact hc
    have hpq : (p₁ : ℤ) ≡ (p₂ : ℤ) [ZMOD (q : ℤ)] := by
      calc (p₁ : ℤ) = 1 * (p₁ : ℤ) := (one_mul _).symm
        _ ≡ (v i * a) * (p₁ : ℤ) [ZMOD (q : ℤ)] := (hva.symm).mul_right _
        _ = a * (v i * (p₁ : ℤ)) := by ring
        _ ≡ a * (v i * (p₂ : ℤ)) [ZMOD (q : ℤ)] := hmodeq.mul_left _
        _ = (v i * a) * (p₂ : ℤ) := by ring
        _ ≡ 1 * (p₂ : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (p₂ : ℤ) := one_mul _
    have h1 : (p₁ : ℤ) % (q : ℤ) = (p₁ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₁.1.2)
    have h2 : (p₂ : ℤ) % (q : ℤ) = (p₂ : ℤ) :=
      Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hp₂.1.2)
    have hfin : (p₁ : ℤ) = (p₂ : ℤ) := by
      have hpq' : (p₁ : ℤ) % (q : ℤ) = (p₂ : ℤ) % (q : ℤ) := hpq
      rw [h1, h2] at hpq'
      exact hpq'
    exact_mod_cast hfin
  · -- surjective onto the nonzero failing residues
    intro r hr
    rw [Finset.mem_filter, Finset.mem_Ioo] at hr
    obtain ⟨⟨hr0, hrq⟩, hrfail⟩ := hr
    set p : ℕ := ((a * (r : ℤ)) % (q : ℤ)).toNat with hpdef
    have hpmod_nonneg : (0 : ℤ) ≤ (a * (r : ℤ)) % (q : ℤ) := Int.emod_nonneg _ (by omega)
    have hpmod_lt : (a * (r : ℤ)) % (q : ℤ) < (q : ℤ) := Int.emod_lt_of_pos _ hqZ
    have hpval : ((p : ℕ) : ℤ) = (a * (r : ℤ)) % (q : ℤ) := by
      rw [hpdef]
      exact Int.toNat_of_nonneg hpmod_nonneg
    have hp_ar : (p : ℤ) ≡ a * (r : ℤ) [ZMOD (q : ℤ)] := by
      show (p : ℤ) % (q : ℤ) = (a * (r : ℤ)) % (q : ℤ)
      rw [hpval]
      exact Int.emod_emod_of_dvd _ dvd_rfl
    have hvp : v i * (p : ℤ) ≡ (r : ℤ) [ZMOD (q : ℤ)] := by
      calc v i * (p : ℤ)
          ≡ v i * (a * (r : ℤ)) [ZMOD (q : ℤ)] := hp_ar.mul_left _
        _ = (v i * a) * (r : ℤ) := by ring
        _ ≡ 1 * (r : ℤ) [ZMOD (q : ℤ)] := hva.mul_right _
        _ = (r : ℤ) := one_mul _
    have hfp : f p = r := by
      have hfpZ : ((f p : ℕ) : ℤ) = (r : ℤ) := by
        rw [hfval]
        have hrsmall : (r : ℤ) % (q : ℤ) = (r : ℤ) :=
          Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hrq)
        have hvp' : (v i * (p : ℤ)) % (q : ℤ) = (r : ℤ) % (q : ℤ) := hvp
        rw [hrsmall] at hvp'
        exact hvp'
      exact_mod_cast hfpZ
    refine ⟨p, ?_, hfp⟩
    rw [Finset.mem_filter, Finset.mem_Ioo]
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · rcases Nat.eq_zero_or_pos p with h0 | hpos
      · exfalso
        have hr' : f p = r := hfp
        rw [h0] at hr'
        have hf0 : f 0 = 0 := by
          rw [hf]
          simp
        omega
      · exact hpos
    · have hlt : ((p : ℕ) : ℤ) < (q : ℤ) := by rw [hpval]; exact hpmod_lt
      exact_mod_cast hlt
    · intro hcon
      exact hrfail (by
        have hb := (hband_iff p).mp hcon
        rwa [hfp] at hb)

/-- **The singleton deviation bounds**: for the gcd case at any `q ≥ 14`,
`−13/7 ≤ D_{i} ≤ 0`.  The singles never threaten the equilibrium budget. -/
theorem deviation_singleton_bounds (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13)
    (hq : 14 ≤ q) (hgcd : Int.gcd (v i) (q : ℤ) = 1) :
    -(13 : ℚ) / 7 ≤ deviation v q {i} ∧ deviation v q {i} ≤ 0 := by
  have hq0 : 0 < q := by omega
  unfold deviation
  rw [jointFail_singleton_eq v q i hq0 hgcd, Finset.card_singleton, bandSize_eq q hq]
  set s : ℕ := 13 * q / 14 - (q + 13) / 14 + 1 with hs
  have hkey : 6 * q - 6 ≤ 7 * s ∧ 7 * s ≤ 6 * q + 7 := by
    constructor <;> omega
  have hsle : s ≤ q - 1 := by omega
  have h1 : ((q - 1 - s : ℕ) : ℚ) = (q : ℚ) - 1 - (s : ℚ) := by
    have hq1 : 1 ≤ q := by omega
    push_cast [Nat.cast_sub hsle, Nat.cast_sub hq1]
    ring
  rw [h1, pow_one]
  have hk1 : (6 : ℚ) * (q : ℚ) - 6 ≤ 7 * (s : ℚ) := by
    have h := hkey.1
    have hq6 : 6 ≤ 6 * q := by omega
    have hcast : ((6 * q - 6 : ℕ) : ℚ) ≤ ((7 * s : ℕ) : ℚ) := by exact_mod_cast h
    rw [Nat.cast_sub hq6] at hcast
    push_cast at hcast
    linarith
  have hk2 : (7 : ℚ) * (s : ℚ) ≤ 6 * (q : ℚ) + 7 := by
    have := hkey.2
    exact_mod_cast this
  constructor
  · linarith
  · linarith

/-- **The exact constant at `14 ∣ q`**: every singleton deviation is `−13/7`. -/
theorem deviation_singleton_of_dvd (v : Fin 13 → ℤ) (q : ℕ) (i : Fin 13)
    (hq : 14 ≤ q) (hdvd : 14 ∣ q) (hgcd : Int.gcd (v i) (q : ℤ) = 1) :
    deviation v q {i} = -(13 : ℚ) / 7 := by
  have hq0 : 0 < q := by omega
  unfold deviation
  rw [jointFail_singleton_eq v q i hq0 hgcd, Finset.card_singleton, bandSize_eq q hq]
  set s : ℕ := 13 * q / 14 - (q + 13) / 14 + 1 with hs
  have hexact : 7 * s = 6 * q + 7 := by omega
  have hsle : s ≤ q - 1 := by omega
  have h1 : ((q - 1 - s : ℕ) : ℚ) = (q : ℚ) - 1 - (s : ℚ) := by
    have hq1 : 1 ≤ q := by omega
    push_cast [Nat.cast_sub hsle, Nat.cast_sub hq1]
    ring
  rw [h1, pow_one]
  have hexactQ : (7 : ℚ) * (s : ℚ) = 6 * (q : ℚ) + 7 := by exact_mod_cast hexact
  linear_combination (-1 / 7 : ℚ) * hexactQ

/-! ## Axiom audit -/
#print axioms jointFail_singleton_eq
#print axioms deviation_singleton_bounds
#print axioms deviation_singleton_of_dvd

end LRC14Concrete
end LonelyRunner
