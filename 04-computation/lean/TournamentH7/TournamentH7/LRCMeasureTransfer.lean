/-
  TournamentH7.LRCMeasureTransfer -- THM-685's transfer pipeline in Lean
  (klein-2026-07-10-S242, HYP-5935): rational safe-interval CERTIFICATES yield
  strictly-live rulers at EVERY sufficiently large modulus.  Fully integer --
  no measure theory, no reals, no Fourier.

  The measure program (THM-685..692) ends in explicit safe intervals
  [x/D, y/D] on which every speed stays strictly inside the band (the
  first-window slivers, the 7-torsion slivers, the engine-emitted component
  lists).  This file is the bridge that makes any such certificate
  operational:

    SafeIvStrict v D x y   :=  exists j, D < 14(vx - jD)  and  14(vy - jD) < 13D
      -- one integer floor-witness certifies speed v strictly in-band on the
         WHOLE interval (linearity: the fractional part v*alpha - j is
         monotone on the interval, pinned at both endpoints)
    strict_band_of_cert    :  certificate + any grid point c/q inside the
                              interval  =>  q < 14*((v*c) % q) < 13*q
                              (the scaled strict rounding identity)
    exists_grid_point      :  D < (y-x)*q  =>  the interval contains a grid
                              point c/q with 1 <= c < q  (c = xq/D + 1)
    strictlyLive_of_cert   :  13 certificates  =>  exists p, StrictlyLive v q p
                              -- at EVERY modulus q > D/(y-x), prime or not
    strictWitness_of_cert  :  ... => StrictWitness v   (kps's LRCStrictRuler
                              chain: strict witness => measure floor => lonely)

  DEMO: the deep well {1..12, 182} (the covering-min extremizer) with the
  explicit certificate [93/1274, 96/1274]: StrictWitness -- hence a strict
  safe time with clearance > 1/14 -- obtained from ONE certificate, valid
  through a strictly-live ruler at EVERY q >= 425.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCStrictRuler

namespace LonelyRunner
namespace MeasureTransfer

open LRC14Grand

/-- Speed `v` stays STRICTLY inside the safe band on the rational interval
`[x/D, y/D]`: a single integer floor-witness `j` pins both endpoints
(`v·(x/D) - j` and `v·(y/D) - j` both lie strictly in `(1/14, 13/14)`, scaled
by `14D`). -/
def SafeIvStrict (v D x y : ℤ) : Prop :=
  ∃ j : ℤ, D < 14 * (v * x - j * D) ∧ 14 * (v * y - j * D) < 13 * D

/-- **The scaled strict rounding identity**: a certified speed is strictly
in-band at every grid point `c/q` of the certificate interval. -/
theorem strict_band_of_cert {v D x y q c : ℤ} (hv : 0 < v) (hD : 0 < D)
    (hq : 0 < q) (hcert : SafeIvStrict v D x y)
    (hc1 : x * q ≤ c * D) (hc2 : c * D ≤ y * q) :
    q < 14 * ((v * c) % q) ∧ 14 * ((v * c) % q) < 13 * q := by
  obtain ⟨j, hj1, hj2⟩ := hcert
  set r : ℤ := v * c - j * q with hr
  have key1 : v * (x * q) ≤ v * (c * D) :=
    mul_le_mul_of_nonneg_left hc1 hv.le
  have key2 : v * (c * D) ≤ v * (y * q) :=
    mul_le_mul_of_nonneg_left hc2 hv.le
  have hlo : D * q < D * (14 * r) := by
    have step1 : D * q < 14 * (v * x - j * D) * q := by nlinarith [hj1, hq]
    have step2 : 14 * (v * x - j * D) * q ≤ D * (14 * r) := by
      nlinarith [key1]
    linarith
  have hhi : D * (14 * r) < D * (13 * q) := by
    have step1 : D * (14 * r) ≤ 14 * (v * y - j * D) * q := by
      nlinarith [key2]
    have step2 : 14 * (v * y - j * D) * q < 13 * D * q := by nlinarith [hj2, hq]
    nlinarith [step1, step2]
  have hqr : q < 14 * r := by
    by_contra h
    push_neg at h
    nlinarith [hlo, hD, h]
  have hr13 : 14 * r < 13 * q := by
    by_contra h
    push_neg at h
    nlinarith [hhi, hD, h]
  have hr0 : 0 ≤ r := by omega
  have hrq : r < q := by omega
  have hmod : (v * c) % q = r := by
    have hvc : v * c = r + q * j := by rw [hr]; ring
    rw [hvc, Int.add_mul_emod_self_left, Int.emod_eq_of_lt hr0 hrq]
  rw [hmod]
  exact ⟨hqr, hr13⟩

/-- **An interval longer than the grid spacing contains a grid point**:
if `D < (y-x)·q` then some `c/q` with `1 ≤ c < q` lies in `[x/D, y/D]`.
Constructive: `c = ⌊xq/D⌋ + 1`. -/
theorem exists_grid_point {x y D q : ℤ} (hD : 0 < D) (hq : 0 < q)
    (hx : 0 ≤ x) (hy : y < D) (hlen : D < (y - x) * q) :
    ∃ c : ℤ, 1 ≤ c ∧ c < q ∧ x * q ≤ c * D ∧ c * D ≤ y * q := by
  refine ⟨x * q / D + 1, ?_, ?_, ?_, ?_⟩
  all_goals
    have hdm := Int.ediv_add_emod (x * q) D
    have hm0 := Int.emod_nonneg (x * q) hD.ne'
    have hmD := Int.emod_lt_of_pos (x * q) hD
  · have : 0 ≤ x * q / D := Int.ediv_nonneg (by positivity) hD.le
    omega
  · -- c < q:  c·D ≤ xq + D ≤ yq < Dq
    have hcd : (x * q / D + 1) * D ≤ x * q + D := by nlinarith [hdm, hm0]
    have hyq : y * q < D * q := by nlinarith [hy, hq]
    nlinarith [hcd, hlen, hyq]
  · -- xq ≤ c·D  (strictly, in fact)
    nlinarith [hdm, hmD]
  · -- c·D ≤ yq
    have hcd : (x * q / D + 1) * D ≤ x * q + D := by nlinarith [hdm, hm0]
    nlinarith [hcd, hlen]

/-- **THE TRANSFER (THM-685, strict form)**: thirteen safe-interval
certificates yield a strictly-live ruler at EVERY modulus `q > D/(y-x)` --
prime or not.  The continuum measure program hands over exactly here. -/
theorem strictlyLive_of_cert (v : Fin 13 → ℤ) (D x y q : ℤ)
    (hv : ∀ i, 0 < v i) (hD : 0 < D) (hq : 0 < q) (hx : 0 ≤ x) (hy : y < D)
    (hlen : D < (y - x) * q) (hcert : ∀ i, SafeIvStrict (v i) D x y) :
    ∃ p, StrictlyLive v q p := by
  obtain ⟨c, hc1, hcq, hg1, hg2⟩ := exists_grid_point hD hq hx hy hlen
  exact ⟨c, by omega, hcq,
    fun i => strict_band_of_cert (hv i) hD hq (hcert i) hg1 hg2⟩

/-- The certificate-to-witness composition: any rational safe interval of
positive length gives a STRICT WITNESS (kps's chain consumes it from here:
strict witness ⟹ positive safe measure ⟹ lonely). -/
theorem strictWitness_of_cert (v : Fin 13 → ℤ) (D x y q : ℤ)
    (hv : ∀ i, 0 < v i) (hD : 0 < D) (hq : 0 < q) (hx : 0 ≤ x) (hy : y < D)
    (hlen : D < (y - x) * q) (hcert : ∀ i, SafeIvStrict (v i) D x y) :
    StrictWitness v := by
  obtain ⟨p, hp⟩ := strictlyLive_of_cert v D x y q hv hD hq hx hy hlen hcert
  exact strictWitness_of_strictlyLive hp

/-! ## Demo: the deep well, certified once, live at every large modulus -/

/-- The deep well `{1,…,12,182}` — the covering-min extremizer `M = 14/183`. -/
def deepWell : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 182]

/-- The explicit certificate: every deep-well speed is strictly in-band on
`[93/1274, 96/1274]` (floor-witness `j = 0` for speeds 1–12, `j = 13` for
182).  Thirteen pairs of integer inequalities. -/
theorem deepWell_cert : ∀ i, SafeIvStrict (deepWell i) 1274 93 96 := by
  intro i
  fin_cases i <;>
    first
      | exact ⟨0, by decide⟩
      | exact ⟨13, by decide⟩

/-- **The standing transfer, demonstrated**: the deep well has a strictly-live
ruler at EVERY modulus `q ≥ 425` — one certificate, infinitely many rulers. -/
theorem deepWell_strictlyLive (q : ℤ) (hq : 425 ≤ q) :
    ∃ p, StrictlyLive deepWell q p := by
  refine strictlyLive_of_cert deepWell 1274 93 96 q ?_ (by norm_num)
    (by omega) (by norm_num) (by norm_num) (by norm_num; omega) deepWell_cert
  intro i
  fin_cases i <;> decide

/-- The deep well has a strict witness — the full pipeline, kernel-pure. -/
theorem deepWell_strictWitness : StrictWitness deepWell :=
  strictWitness_of_cert deepWell 1274 93 96 425
    (by intro i; fin_cases i <;> decide)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    deepWell_cert

end MeasureTransfer
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.MeasureTransfer.strict_band_of_cert
#print axioms LonelyRunner.MeasureTransfer.exists_grid_point
#print axioms LonelyRunner.MeasureTransfer.strictlyLive_of_cert
#print axioms LonelyRunner.MeasureTransfer.strictWitness_of_cert
#print axioms LonelyRunner.MeasureTransfer.deepWell_strictWitness
