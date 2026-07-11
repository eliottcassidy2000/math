/-
  TournamentH7.LRCTestPointCore -- the test-point theorems' arithmetic cores
  in Lean (klein-2026-07-10-S243, HYP-5940): THM-690 (13-torsion), THM-691(A)
  (q*), THM-692 (7-torsion middle-chain rigidity), plus THE FATTENING LEMMA
  that turns a q*-test-point into a SafeIvStrict certificate consumable by
  LRCMeasureTransfer (S242).

    net_value_strictly_in_band : q <= 13, 1 <= r < q  =>  q < 14r < 13q
      -- every nonzero q*-net value is STRICTLY inside the band: the P-side
    qstar_p_nonzero  :  gcd(a,q) = 1, 1 <= p <= 13, p != q  =>  q does not
      divide p*a  (coprimality moves the divisibility to p; p < 2q kills it)
    residue_in_range :  the two combined -- (p*a) % q lands in [1, q-1]
    qstar_cert (THE FATTENING) :  the test point a/q fattens to the EXPLICIT
      certificate interval [ (4732a - q)/4732q , (4732a + q)/4732q ], giving
      SafeIvStrict p (4732q) (4732a - q) (4732a + q) for every p in [1,13]
      with q not dividing p*a -- the uniform slack W(14r - q) >= 4732 >
      14*13*13 = 2366 covers all q <= 13, p <= 13 at once.  Composes DIRECTLY
      into MeasureTransfer.strictlyLive_of_cert: a q*-test-point is a
      certificate factory for the P-part.
    pigeonhole_missed_class :  |E| < q  =>  some residue class mod q carries
      no co-offset (the E-side of THM-690/691A)
    middle_chain_rigidity :  THM-692's core -- if the middle classes' min/max
      co-offsets satisfy BOTH sign chains, they coincide pointwise, which
      distinct residues forbid: coverage cannot hold on both sides of a
      7-torsion point.

  Kernel-pure: no native_decide, no sorry.
-/
import Mathlib
import TournamentH7.LRCMeasureTransfer

namespace LonelyRunner
namespace TestPointCore

open MeasureTransfer

/-- **The P-side band bound**: every nonzero value on a `q`-net with `q ≤ 13`
is STRICTLY inside the safe band. -/
theorem net_value_strictly_in_band {q r : ℤ} (hq : q ≤ 13) (hr1 : 1 ≤ r)
    (hrq : r < q) : q < 14 * r ∧ 14 * r < 13 * q := by omega

/-- **The P-side nonvanishing**: for `a` coprime to `q` and a speed
`p ∈ [1,13]` with `p ≠ q`, the residue `(p·a) % q` is nonzero — `q ∤ pa`
would force `q ∣ p`, impossible below `2q`. -/
theorem qstar_p_nonzero {q p a : ℤ} (hq8 : 8 ≤ q) (hcop : IsCoprime a q)
    (hp1 : 1 ≤ p) (hp13 : p ≤ 13) (hpq : p ≠ q) : ¬ (q ∣ p * a) := by
  intro hdvd
  have hqp : q ∣ p := by
    have := (IsCoprime.symm hcop).dvd_of_dvd_mul_right hdvd
    exact this
  obtain ⟨k, hk⟩ := hqp
  have hq0 : 0 < q := by omega
  rcases le_or_gt k 0 with h | h
  · nlinarith
  · have : 1 ≤ k := h
    rcases eq_or_lt_of_le this with rfl | h2
    · omega
    · nlinarith

/-- The residue of a certified speed lands in `[1, q−1]`. -/
theorem residue_in_range {q p a : ℤ} (hq8 : 8 ≤ q) (hcop : IsCoprime a q)
    (hp1 : 1 ≤ p) (hp13 : p ≤ 13) (hpq : p ≠ q) :
    1 ≤ (p * a) % q ∧ (p * a) % q < q := by
  have hq0 : (0 : ℤ) < q := by omega
  have h0 : 0 ≤ (p * a) % q := Int.emod_nonneg _ (by omega)
  have hlt : (p * a) % q < q := Int.emod_lt_of_pos _ hq0
  have hne : (p * a) % q ≠ 0 := by
    intro h
    exact qstar_p_nonzero hq8 hcop hp1 hp13 hpq (Int.dvd_of_emod_eq_zero h)
  exact ⟨by omega, hlt⟩

/-- **THE FATTENING LEMMA**: the `q*`-test-point `a/q` fattens into the
explicit SafeIvStrict certificate interval
`[(4732a − q)/(4732q), (4732a + q)/(4732q)]` — for EVERY speed `p ∈ [1,13]`
whose residue is nonzero.  The uniform slack `4732 > 14·13·13` absorbs all
`q ≤ 13, p ≤ 13` at once.  A test point is a certificate factory. -/
theorem qstar_cert {q p a : ℤ} (hq8 : 8 ≤ q) (hq13 : q ≤ 13)
    (hp1 : 1 ≤ p) (hp13 : p ≤ 13)
    (hr : 1 ≤ (p * a) % q) :
    SafeIvStrict p (4732 * q) (4732 * a - q) (4732 * a + q) := by
  refine ⟨p * a / q, ?_, ?_⟩
  all_goals
    have hq0 : (0 : ℤ) < q := by omega
    have hdm := Int.ediv_add_emod (p * a) q
    have hlt : (p * a) % q < q := Int.emod_lt_of_pos _ hq0
  · -- 4732q < 14 * (p(4732a − q) − j·4732q) = 14(4732 r − p q)
    nlinarith [hdm, hr, hlt, hp1, hp13]
  · -- 14 * (p(4732a + q) − j·4732q) = 14(4732 r + p q) < 13·4732q
    nlinarith [hdm, hr, hlt, hp1, hp13]

/-- **The E-side pigeonhole (THM-690/691A)**: fewer than `q` co-offsets miss
some residue class mod `q`. -/
theorem pigeonhole_missed_class (E : Finset ℤ) (q a : ℤ) (hq : 0 < q)
    (hcard : (E.card : ℤ) < q) :
    ∃ c, 0 ≤ c ∧ c < q ∧ ∀ e ∈ E, (e * a) % q ≠ c := by
  classical
  have himg : (E.image fun e => (e * a) % q).card ≤ E.card :=
    Finset.card_image_le
  have hqnat : (Finset.Ico (0 : ℤ) q).card = q.toNat := by
    rw [Int.card_Ico]
    simp
  have hproper : (E.image fun e => (e * a) % q).card < (Finset.Ico (0 : ℤ) q).card := by
    rw [hqnat]
    omega
  obtain ⟨c, hc1, hc2⟩ := Finset.exists_mem_notMem_of_card_lt_card hproper
  refine ⟨c, (Finset.mem_Ico.mp hc1).1, (Finset.mem_Ico.mp hc1).2, ?_⟩
  intro e he heq
  exact hc2 (Finset.mem_image.mpr ⟨e, he, heq⟩)

/-- **The middle-chain rigidity (THM-692's core)**: if the middle classes'
extreme co-offsets satisfy BOTH sign chains (coverage on both sides of the
7-torsion point), they coincide pointwise — impossible when adjacent classes
carry distinct residues mod 7.  Coverage holds on at most one side. -/
theorem middle_chain_rigidity (L R : Fin 6 → ℤ)
    (hres : ∀ j : Fin 5, L j.succ % 7 ≠ R j.castSucc % 7)
    (hdown : ∀ j : Fin 5, L j.succ ≤ R j.castSucc)
    (hup : ∀ j : Fin 5, R j.castSucc ≤ L j.succ) : False := by
  have heq : L (0 : Fin 5).succ = R (0 : Fin 5).castSucc :=
    le_antisymm (hdown 0) (hup 0)
  exact hres 0 (by rw [heq])

/-! ## Demo: the q*-certificate factory, instantiated

`P = {1, 13}` has `q* = 12`; at `a = 1` the residues are `1` and `1`
(`13 ≡ 1 mod 12`), both in `[1, 11]`, so both speeds get SafeIvStrict
certificates on `[(4732 − 12)/56784, (4732 + 12)/56784]` — decidably. -/

theorem demo_cert_speed1 :
    SafeIvStrict 1 (4732 * 12) (4732 * 1 - 12) (4732 * 1 + 12) :=
  qstar_cert (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by decide)

theorem demo_cert_speed13 :
    SafeIvStrict 13 (4732 * 12) (4732 * 1 - 12) (4732 * 1 + 12) :=
  qstar_cert (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by decide)

end TestPointCore
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.TestPointCore.net_value_strictly_in_band
#print axioms LonelyRunner.TestPointCore.qstar_p_nonzero
#print axioms LonelyRunner.TestPointCore.residue_in_range
#print axioms LonelyRunner.TestPointCore.qstar_cert
#print axioms LonelyRunner.TestPointCore.pigeonhole_missed_class
#print axioms LonelyRunner.TestPointCore.middle_chain_rigidity
