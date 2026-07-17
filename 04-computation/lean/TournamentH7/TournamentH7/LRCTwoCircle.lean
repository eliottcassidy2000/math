/-
  TournamentH7.LRCTwoCircle — THE TWO-CIRCLE DEEP CERTIFICATE, PART I
  (death-star-2026-07-17-S52, HYP-7265; the wagner-circle hint decoded).

  THE CLAIM (recon-exact over 1185 moduli, every multiplier, zero mismatches):
  on the canonical family v = (1,…,13), the deep set (bandCount ≥ 6) is
  EXACTLY two resonance circles:
    (I)  the integer circle  `84·p < q ∨ 84·(q−p) < q`
    (II) the half circle     `84·|2p − q| < q`.

  This module delivers Part I — the certificate machinery and the forward
  half plus the first exclusion cases:

  * `divisor_descent` — if speed `k` fails with witness `w`, then speed `k/g`
    fails with witness `w/g` for `g = gcd(k, w)`: failures descend to the
    reduced ray.  (The workhorse for constraining witnesses by which SMALLER
    speeds are safe.)
  * `circleI_deep` / `circleII_deep` — the ⟸ direction, in full: on circle I
    the six smallest speeds fail (nested); on circle II the six even speeds
    fail (coherently, witnesses `k·(2p−q)`-shifted).  Both give
    `bandCount ≥ 6` for the canonical family.
  * `bandCount_ge_of_card` / `not_inBand_of_witness` — the assembly tools:
    any six witnessed failures defeat the cap.

  Part II (the ⟹ direction: deep ⟹ one of the two circles, by case analysis
  on the smallest failing speed k₀ — k₀ = 1 nests via the hub lock; k₀ = 2
  forces the six evens via parity locks + the 13-branch ray estimate;
  k₀ ∈ {3,…,6, 8} die by parity/divisor locks; k₀ = 7 dies by the residue
  gem `l·w₇ ≡ ±1 (mod 7)` forcing ≤ 2 residue classes; k₀ ≥ 9 is trivial)
  is the named next step, fully structured in the S52 session log.

  Kernel-pure: no `sorry`, no `native_decide`, no new axioms.
-/
import Mathlib
import TournamentH7.LRCRelationLock

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- The canonical family: speeds 1,…,13. -/
def canonical : Fin 13 → ℤ := fun i => (i : ℤ) + 1

/-- **Divisor descent**: failures descend to the reduced ray — if speed `k`
fails with witness `w` and `d` divides both, speed `k/d` fails with witness
`w/d`. -/
theorem divisor_descent (k w d : ℤ) (q p : ℕ) (hd : 0 < d)
    (hdk : d ∣ k) (hdw : d ∣ w)
    (h : 14 * |k * (p : ℤ) - w * q| < q) :
    14 * |(k / d) * (p : ℤ) - (w / d) * q| < q := by
  obtain ⟨k', hk⟩ := hdk
  obtain ⟨w', hw⟩ := hdw
  have hkd : k / d = k' := by rw [hk]; exact Int.mul_ediv_cancel_left k' (by omega)
  have hwd : w / d = w' := by rw [hw]; exact Int.mul_ediv_cancel_left w' (by omega)
  rw [hkd, hwd]
  have hexp : k * (p : ℤ) - w * q = d * (k' * p - w' * q) := by
    rw [hk, hw]; ring
  rw [hexp, abs_mul, abs_of_pos hd] at h
  have habs : (0 : ℤ) ≤ |k' * (p : ℤ) - w' * q| := abs_nonneg _
  nlinarith [h, habs, hd]

/-- Circle I puts the six smallest canonical speeds in the band:
`84·p < q` (low side) makes speeds 1..6 fail with witness 0. -/
theorem circleI_low_fails (q p : ℕ) (k : ℤ) (hk1 : 1 ≤ k) (hk6 : k ≤ 6)
    (hlow : 84 * (p : ℤ) < q) (hp : 0 < p) :
    14 * |k * (p : ℤ) - 0 * q| < q := by
  have hpZ : (0 : ℤ) < (p : ℤ) := by exact_mod_cast hp
  have habs : |k * (p : ℤ) - 0 * (q : ℤ)| = k * p := by
    rw [zero_mul, sub_zero]
    exact abs_of_pos (by positivity)
  rw [habs]
  nlinarith [hlow, hpZ]

/-- Circle I high side: `84·(q−p) < q` makes speeds 1..6 fail with witness k. -/
theorem circleI_high_fails (q p : ℕ) (k : ℤ) (hk1 : 1 ≤ k) (hk6 : k ≤ 6)
    (hhigh : 84 * ((q : ℤ) - p) < q) (hpq : p < q) :
    14 * |k * (p : ℤ) - k * q| < q := by
  have hdZ : (0 : ℤ) < (q : ℤ) - p := by
    have : (p : ℤ) < q := by exact_mod_cast hpq
    linarith
  have habs : |k * (p : ℤ) - k * (q : ℤ)| = k * ((q : ℤ) - p) := by
    have h1 : k * (p : ℤ) - k * q = -(k * ((q : ℤ) - p)) := by ring
    rw [h1, abs_neg]
    exact abs_of_pos (by positivity)
  rw [habs]
  nlinarith [hhigh, hdZ]

/-- Circle II puts the six even canonical speeds in the band: the even speed
`2k` fails with witness `k·p·2 − k·(2p−q)`-shifted — concretely, witness `k`
against the residue `k·(2p − q)`. -/
theorem circleII_fails (q p : ℕ) (k : ℤ) (hk1 : 1 ≤ k) (hk6 : k ≤ 6)
    (hhalf : 84 * |2 * (p : ℤ) - q| < q) :
    14 * |(2 * k) * (p : ℤ) - k * q| < q := by
  have habs : |(2 * k) * (p : ℤ) - k * (q : ℤ)| = k * |2 * (p : ℤ) - q| := by
    have h1 : (2 * k) * (p : ℤ) - k * q = k * (2 * p - q) := by ring
    rw [h1, abs_mul, abs_of_pos (by linarith : (0 : ℤ) < k)]
  rw [habs]
  nlinarith [hhalf, abs_nonneg (2 * (p : ℤ) - q), hk1, hk6]

/-- Six failing indices (any decidable selector) give the bandCount bound. -/
theorem bandCount_ge_of_card (v : Fin 13 → ℤ) (q p : ℕ)
    (P : Fin 13 → Prop) [DecidablePred P]
    (hcard : 6 ≤ (Finset.univ.filter P).card)
    (h : ∀ i : Fin 13, P i → ¬ inBand v q p i) :
    6 ≤ bandCount v q p := by
  unfold bandCount
  have hsub : (Finset.univ.filter P)
      ⊆ Finset.univ.filter fun i : Fin 13 => ¬ inBand v q p i := by
    intro i hi
    rw [Finset.mem_filter] at hi ⊢
    exact ⟨Finset.mem_univ _, h i hi.2⟩
  calc 6 ≤ (Finset.univ.filter P).card := hcard
    _ ≤ _ := Finset.card_le_card hsub

/-- A band failure in witness form defeats `inBand`. -/
theorem not_inBand_of_witness (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13)
    (hq : 0 < q) (w : ℤ) (h : 14 * |v i * (p : ℤ) - w * q| < q) :
    ¬ inBand v q p i := by
  intro hband
  obtain ⟨h1, h2⟩ := hband
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  -- the residue r = v i * p mod q sits in [q/14, 13q/14]; the witness puts
  -- v i * p within q/14 of w*q — the two are incompatible
  have hr := Int.emod_emod_of_dvd (v i * (p : ℤ)) (dvd_refl (q : ℤ))
  set r : ℤ := (v i * (p : ℤ)) % (q : ℤ) with hrdef
  have hdm := Int.ediv_add_emod (v i * (p : ℤ)) (q : ℤ)
  set m : ℤ := (v i * (p : ℤ)) / (q : ℤ) with hmdef
  -- v i * p − w q = (m − w) q + r
  have hexp : v i * (p : ℤ) - w * q = ((m - w) * q) + r := by
    have : (q : ℤ) * m + r = v i * (p : ℤ) := hdm
    linarith
  rw [hexp] at h
  rcases le_or_gt (m - w) (-1) with hle | hge
  · -- v i p − w q ≤ −q + r ≤ −q + 13q/14 < 0, abs ≥ q − r ≥ q/14
    have hub : ((m - w) * q) + r ≤ -(q : ℤ) + r := by nlinarith [hqZ]
    have habs2 : |((m - w) * q) + r| ≥ (q : ℤ) - r := by
      have h3 : ((m - w) * q) + r < 0 := by nlinarith [h2, hqZ]
      rw [abs_of_neg h3]
      nlinarith [hub]
    nlinarith [h, h1, h2, habs2]
  · rcases le_or_gt 1 (m - w) with hge1 | heq
    · have hlb : ((m - w) * q) + r ≥ (q : ℤ) + r := by nlinarith [hqZ]
      have habs : |((m - w) * q) + r| ≥ (q : ℤ) + r := by
        rw [abs_of_pos (by nlinarith [h1, hqZ] : (0:ℤ) < ((m - w) * q) + r)]
        exact hlb
      nlinarith [h, h1]
    · have heq0 : m - w = 0 := by omega
      rw [heq0, zero_mul, zero_add] at h
      have hrnn : (0 : ℤ) ≤ r := Int.emod_nonneg _ (by omega)
      rw [abs_of_nonneg hrnn] at h
      linarith [h1]

/-- **Circle I ⟹ deep**: on the integer circle the canonical family has
`bandCount ≥ 6` — the six smallest speeds fail. -/
theorem circleI_deep (q p : ℕ) (hq : 0 < q) (hp : 0 < p) (hpq : p < q)
    (hcirc : 84 * (p : ℤ) < q ∨ 84 * ((q : ℤ) - p) < q) :
    6 ≤ bandCount canonical q p := by
  apply bandCount_ge_of_card canonical q p (fun i => (i : ℕ) < 6) (by decide)
  intro i hi
  have hk1 : (1 : ℤ) ≤ (i : ℤ) + 1 := by
    have := Int.natCast_nonneg (i : ℕ)
    omega
  have hk6 : (i : ℤ) + 1 ≤ 6 := by
    have : ((i : ℕ) : ℤ) < 6 := by exact_mod_cast hi
    omega
  rcases hcirc with hlow | hhigh
  · exact not_inBand_of_witness canonical q p i hq 0
      (by simpa [canonical] using circleI_low_fails q p ((i : ℤ) + 1) hk1 hk6 hlow hp)
  · exact not_inBand_of_witness canonical q p i hq ((i : ℤ) + 1)
      (by simpa [canonical] using circleI_high_fails q p ((i : ℤ) + 1) hk1 hk6 hhigh hpq)

/-- **Circle II ⟹ deep**: on the half circle the canonical family has
`bandCount ≥ 6` — the six even speeds fail. -/
theorem circleII_deep (q p : ℕ) (hq : 0 < q)
    (hhalf : 84 * |2 * (p : ℤ) - q| < q) :
    6 ≤ bandCount canonical q p := by
  apply bandCount_ge_of_card canonical q p
    (fun i => (i : ℕ) % 2 = 1 ∧ (i : ℕ) < 12) (by decide)
  intro i ⟨hodd, hlt⟩
  -- index i has speed i+1 = 2k with k = (i+1)/2 ∈ [1,6]
  set k : ℤ := ((i : ℕ) + 1 : ℕ) / 2 with hk_def
  have hk1 : (1 : ℤ) ≤ k := by
    rw [hk_def]
    have : 1 ≤ ((i : ℕ) + 1) / 2 := by omega
    exact_mod_cast this
  have hk6 : k ≤ 6 := by
    rw [hk_def]
    have : ((i : ℕ) + 1) / 2 ≤ 6 := by omega
    exact_mod_cast this
  have hspeed : canonical i = 2 * k := by
    rw [hk_def]
    unfold canonical
    have h2 : ((i : ℕ) + 1) % 2 = 0 := by omega
    have : ((((i : ℕ) + 1) / 2 : ℕ) : ℤ) * 2 = ((i : ℕ) + 1 : ℕ) := by
      have := Nat.div_mul_cancel (Nat.dvd_of_mod_eq_zero h2)
      exact_mod_cast this
    push_cast at this ⊢
    linarith
  exact not_inBand_of_witness canonical q p i hq k
    (by rw [hspeed]; exact circleII_fails q p k hk1 hk6 hhalf)

/-! ## Axiom audit -/
#print axioms divisor_descent
#print axioms bandCount_ge_of_card
#print axioms not_inBand_of_witness
#print axioms circleI_low_fails
#print axioms circleI_high_fails
#print axioms circleII_fails
#print axioms circleI_deep
#print axioms circleII_deep

end LRC14Concrete
end LonelyRunner
