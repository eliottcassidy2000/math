/-
  TournamentH7.LRCUniformCell — THE UNIFORM CELL LEMMA CORE
  (mac-mini-2026-07-06-S3, HYP-4252; generalizes LRCKStratification's
  (3, 38) instance to every Farey cell (c, q) of THM-622's reduction).

  For ANY attained value M(W) = c/q (reduced) with maximizer m/(q·k):
  the scaling identity descends mod-(qk) distances of k-multiples to mod-q
  distances of quotients (dInt_scale_cell); binders are k-multiples
  (binder_dvd_cell); the attaining grid is a q-multiple (grid_div_cell);
  up/down binder quotients sum to 0 mod q (pair_sum_dvd_cell); and — the
  two new pieces — every quotient binder is a UNIT mod q (binder_unit_cell:
  gcd(v', q) = 1, since a common prime would divide c against gcd(c,q) = 1),
  which pins the witness dilation class per pair shape
  (witness_determined_cell, by coprime cancellation — no prime-by-prime
  case split needed).

  Pure integer arithmetic; no analysis.  Draft with proofs + verification:
  03-artifacts/drafts/uniform-cell-lemma-macmini-S3.md;
  04-computation/lrc_uniform_cell_lemma_macmini_S3.py (0 violations,
  436 families / 268 distinct attained cells).
-/
import Mathlib.RingTheory.Int.Basic
import Mathlib.Tactic

namespace LonelyRunner
namespace UniformCell

/-- Integer circle distance at modulus q (canonical residue form). -/
def dInt (q x : ℤ) : ℤ := min (x % q) (q - x % q)

/-- THE SCALING IDENTITY, cell-parametric: distances at modulus q·k of
    k-multiples are k times the mod-q distances of the quotients. -/
theorem dInt_scale_cell (q k x : ℤ) (hk : 0 < k) :
    dInt (q * k) (k * x) = k * dInt q x := by
  unfold dInt
  have h1 : k * x % (q * k) = k * (x % q) := by
    rw [show (q : ℤ) * k = k * q by ring]
    exact Int.mul_emod_mul_of_pos x q hk
  rw [h1, show (q : ℤ) * k - k * (x % q) = k * (q - x % q) by ring]
  rcases le_total (x % q) (q - x % q) with h | h
  · rw [min_eq_left (by nlinarith), min_eq_left h]
  · rw [min_eq_right (by nlinarith), min_eq_right h]

/-- BINDER DIVISIBILITY, cell-parametric: a runner binding at ±c·k on the
    q·k-grid under a dilation coprime to k is divisible by k. -/
theorem binder_dvd_cell (v m k c q : ℤ) (hco : IsCoprime k m)
    (h : (q * k : ℤ) ∣ (v * m - c * k) ∨ (q * k : ℤ) ∣ (v * m + c * k)) :
    k ∣ v := by
  have hkdvd : (k : ℤ) ∣ q * k := ⟨q, by ring⟩
  have hkvm : k ∣ v * m := by
    rcases h with h | h
    · have hk1 : k ∣ v * m - c * k := hkdvd.trans h
      have h2 : v * m = (v * m - c * k) + c * k := by ring
      rw [h2]
      exact dvd_add hk1 ⟨c, by ring⟩
    · have hk1 : k ∣ v * m + c * k := hkdvd.trans h
      have h2 : v * m = (v * m + c * k) - c * k := by ring
      rw [h2]
      exact dvd_sub hk1 ⟨c, by ring⟩
  exact hco.dvd_of_dvd_mul_right hkvm

/-- GRID DIVISIBILITY, cell-parametric: gcd(c, q) = 1 and q ∣ c·q* force
    q ∣ q* — the attaining grid of a c/q-value is a q-multiple grid. -/
theorem grid_div_cell (c q qs : ℤ) (hco : IsCoprime q c)
    (h : (q : ℤ) ∣ c * qs) : (q : ℤ) ∣ qs :=
  hco.dvd_of_dvd_mul_left h

/-- THE BINDING-PAIR DIVISIBILITY, cell-parametric: an up-binder (≡ +c) and a
    down-binder (≡ −c) under the same dilation coprime to q sum to 0 mod q. -/
theorem pair_sum_dvd_cell (vi vj m c q : ℤ)
    (hi : (q : ℤ) ∣ (vi * m - c)) (hj : (q : ℤ) ∣ (vj * m + c))
    (hm : IsCoprime (q : ℤ) m) :
    (q : ℤ) ∣ (vi + vj) := by
  have hsum : (q : ℤ) ∣ ((vi + vj) * m) := by
    have hadd := dvd_add hi hj
    have heq : (vi * m - c) + (vj * m + c) = (vi + vj) * m := by ring
    rwa [heq] at hadd
  exact hm.dvd_of_dvd_mul_right hsum

/-- BINDER UNITS (the new piece): a runner binding at ±c on the q-grid is a
    UNIT mod q — any common divisor of v and q divides c, so gcd(c, q) = 1
    forces gcd(v, q) = 1.  (At q = 38 this is binder parity PLUS 19 ∤ v;
    in general it collapses the k = 1 pair shapes to the φ(q)/2 unit pairs.) -/
theorem binder_unit_cell (v m c q : ℤ)
    (h : (q : ℤ) ∣ (v * m - c) ∨ (q : ℤ) ∣ (v * m + c))
    (hcq : IsCoprime c q) :
    IsCoprime v q := by
  rw [Int.isCoprime_iff_gcd_eq_one]
  have hg1 : (↑(Int.gcd v q) : ℤ) ∣ v := Int.gcd_dvd_left v q
  have hg2 : (↑(Int.gcd v q) : ℤ) ∣ q := Int.gcd_dvd_right v q
  have hgvm : (↑(Int.gcd v q) : ℤ) ∣ v * m := hg1.mul_right m
  have hgc : (↑(Int.gcd v q) : ℤ) ∣ c := by
    rcases h with h | h
    · have h1 : (↑(Int.gcd v q) : ℤ) ∣ v * m - c := hg2.trans h
      have h2 : (c : ℤ) = v * m - (v * m - c) := by ring
      rw [h2]
      exact dvd_sub hgvm h1
    · have h1 : (↑(Int.gcd v q) : ℤ) ∣ v * m + c := hg2.trans h
      have h2 : (c : ℤ) = (v * m + c) - v * m := by ring
      rw [h2]
      exact dvd_sub h1 hgvm
  have hgcd1 : Int.gcd v q ∣ Int.gcd c q := Int.dvd_gcd hgc hg2
  rw [Int.isCoprime_iff_gcd_eq_one] at hcq
  rw [hcq] at hgcd1
  exact Nat.dvd_one.mp hgcd1

/-- WITNESS DETERMINISM, cell-parametric: two dilations binding the same
    runner at +c on the q-grid agree mod q — by coprime cancellation through
    binder_unit_cell (no prime-by-prime split). -/
theorem witness_determined_cell (v m m' c q : ℤ)
    (hm : (q : ℤ) ∣ (v * m - c)) (hm' : (q : ℤ) ∣ (v * m' - c))
    (hcq : IsCoprime c q) :
    (q : ℤ) ∣ (m - m') := by
  have hunit : IsCoprime v q := binder_unit_cell v m c q (Or.inl hm) hcq
  have hdiff : (q : ℤ) ∣ v * (m - m') := by
    have h := dvd_sub hm hm'
    have heq : (v * m - c) - (v * m' - c) = v * (m - m') := by ring
    rwa [heq] at h
  exact hunit.symm.dvd_of_dvd_mul_left hdiff

end UniformCell
end LonelyRunner
