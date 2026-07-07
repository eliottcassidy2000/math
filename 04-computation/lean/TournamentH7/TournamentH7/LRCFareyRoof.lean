/-
  TournamentH7.LRCFareyRoof — THE FAREY ROOF, pointwise part (THM-637, opus-S134/S135).

  For x in the open Farey-k cell (p/q, p'/q'), the config {frac(i·x) : 1 ≤ i ≤ k} of the
  AP {1..k} satisfies:

    (A)  every point is ≥ q·x − p          [lemmaA: first return from above]
    (B)  every point is ≤ 1 + q'·x − p'    [lemmaB: first return from below]
    (C)  above every config point there is another config point within
         roof := (q·x − p) + (p' − q'·x)   [lemmaC ⟹ no circular gap exceeds the roof]
    (0)  the open interval (q'·x − p', q·x − p) ∋ 0, of length roof, contains no config
         point (zero_gap_empty = A + B in lifted form)

  Hence maxgap = gap∋0 = roof: the AP max-gap function is the "Farey roof", linear on each
  Farey cell with node values 1/q.  This is the pointwise engine behind the diameter floor
  (kps-S59 HYP-4797 subset lemma + monad-S2 HYP-4817 exact crossing:
  μ_{1/7}(E) ≥ μ_{1/7}(AP_{D+1}) ≥ m_P for primitive diameter D ≤ 75 — the k=13 leg of
  hlarge on bounded diameter) and the exact canon constants (93/440, 477/1078, …).

  Proofs are the S135 self-contained arguments — pure determinant arithmetic plus one strict
  inequality per case; no continued fractions, no three-distance citation, no measure theory:
    · no_middle_fraction:  p·i < a·q ∧ a·q' < p'·i  ⟹  i = q(p'i−aq') + q'(aq−pi) ≥ q+q' > k.
    · lemmaA hard case: fract(i·x) < q·x−p forces m := pi−aq ≥ 1 with m·q'+i < q, while
      q·(ip'−aq') = m·q'+i forces q ∣ m·q'+i — contradiction (Int.le_of_dvd).
    · lemmaB: mirrored with b := ⌊i·x⌋+1, m' := bq'−p'i, and q'·(bq−ip) = m'·q+i.
    · lemmaC: the three exhibited indices i+q, i−q', i+q−q'.

  All Farey-neighbor facts are DERIVED from the determinant q·p'−p·q' = 1 and q+q' > k.
  Hypotheses are stated in cleared (division-free) form:  p < q·x  and  q'·x < p'.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner
namespace FareyRoof

variable {p q p' q' k i : ℤ} {x : ℝ}

/-- No fraction with denominator ≤ k lies strictly inside a Farey-k cell:
from `p·i < a·q` and `a·q' < p'·i` (both integer-strict) and the determinant,
`i = q·(p'·i − a·q') + q'·(a·q − p·i) ≥ q + q' > k ≥ i`. -/
theorem no_middle_fraction
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1)
    (hsum : k < q + q') (_hi1 : 1 ≤ i) (hik : i ≤ k)
    (a : ℤ) (h1 : p * i < a * q) (h2 : a * q' < p' * i) : False := by
  have key : q * (p' * i - a * q') + q' * (a * q - p * i) = (q * p' - p * q') * i := by
    ring
  rw [hdet, one_mul] at key
  have h1' : 1 ≤ a * q - p * i := by omega
  have h2' : 1 ≤ p' * i - a * q' := by omega
  have e1 : q * 1 ≤ q * (p' * i - a * q') :=
    mul_le_mul_of_nonneg_left h2' (le_of_lt hq)
  have e2 : q' * 1 ≤ q' * (a * q - p * i) :=
    mul_le_mul_of_nonneg_left h1' (le_of_lt hq')
  omega

/-- Coprimality of `q` and `p` from the determinant. -/
theorem coprime_q_p (hdet : q * p' - p * q' = 1) : IsCoprime q p :=
  ⟨p', -q', by linarith⟩

/-- Coprimality of `q'` and `p'` from the determinant. -/
theorem coprime_q'_p' (hdet : q * p' - p * q' = 1) : IsCoprime q' p' :=
  ⟨-p, q, by linarith⟩

/-! ### Lemma A — first return from above -/

/-- **Lemma A (THM-637).** For `x` strictly inside the cell and `1 ≤ i ≤ k`:
`q·x − p ≤ Int.fract (i·x)`. -/
theorem lemmaA
    (hq : 0 < q) (hq' : 0 < q')
    (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hi1 : 1 ≤ i) (hik : i ≤ k)
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p') :
    (q : ℝ) * x - p ≤ Int.fract ((i : ℝ) * x) := by
  set a : ℤ := ⌊(i : ℝ) * x⌋ with ha
  have hfract : Int.fract ((i : ℝ) * x) = (i : ℝ) * x - a := rfl
  have hfl : (a : ℝ) ≤ (i : ℝ) * x := Int.floor_le _
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hq'R : (0 : ℝ) < (q' : ℝ) := by exact_mod_cast hq'
  have hiR : (1 : ℝ) ≤ (i : ℝ) := by exact_mod_cast hi1
  have hqx : (0 : ℝ) < (q : ℝ) * x - p := by linarith
  -- Step 1: a·q ≤ p·i, else (a, i) is a middle fraction.
  have haq : a * q ≤ p * i := by
    by_contra hcon
    push_neg at hcon
    refine no_middle_fraction hq hq' hdet hsum hi1 hik a (by omega) ?_
    have hr : (a : ℝ) * q' < (p' : ℝ) * i := by
      have h1 : (a : ℝ) * q' ≤ ((i : ℝ) * x) * q' :=
        mul_le_mul_of_nonneg_right hfl (le_of_lt hq'R)
      have h2 : ((i : ℝ) * x) * q' < (p' : ℝ) * i := by
        have h3 : (i : ℝ) * ((q' : ℝ) * x) < (i : ℝ) * p' :=
          mul_lt_mul_of_pos_left hx' (by linarith)
        nlinarith [h3]
      linarith
    exact_mod_cast hr
  rcases eq_or_lt_of_le haq with heq | hlt
  · -- Equality: q ∣ i, i = q·t, a = p·t, fract = t·(q·x − p) ≥ q·x − p.
    have hdvd : q ∣ i := by
      have hqi : q ∣ p * i := ⟨a, by linarith⟩
      exact (coprime_q_p hdet).dvd_of_dvd_mul_left hqi
    obtain ⟨t, ht⟩ := hdvd
    have ht1 : 1 ≤ t := by
      by_contra hcon
      push_neg at hcon
      have ht0 : t ≤ 0 := by omega
      have : q * t ≤ q * 0 := mul_le_mul_of_nonneg_left ht0 (le_of_lt hq)
      omega
    have hat : a = p * t := by
      have h0 : a * q = p * (q * t) := by rw [← ht]; exact heq
      have h1 : q * (a - p * t) = 0 := by linear_combination h0
      have h2 := mul_eq_zero.mp h1
      rcases h2 with h | h
      · omega
      · omega
    rw [hfract, ht, hat]
    have htR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast ht1
    have hmul : 1 * ((q : ℝ) * x - p) ≤ (t : ℝ) * ((q : ℝ) * x - p) :=
      mul_le_mul_of_nonneg_right htR (le_of_lt hqx)
    push_cast
    nlinarith [hmul]
  · -- Strict: m := p·i − a·q ≥ 1; the divisibility contradiction.
    by_contra hcon
    push_neg at hcon
    rw [hfract] at hcon
    set m : ℤ := p * i - a * q with hm
    have hm1 : 1 ≤ m := by omega
    have hmR : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm1
    -- i·(q·x − p) + m < q·(q·x − p)   [multiply hcon by q, substitute a·q = p·i − m]
    have hstep0 : (i : ℝ) * ((q : ℝ) * x - p) + (m : ℝ) < (q : ℝ) * ((q : ℝ) * x - p) := by
      have h1 : (q : ℝ) * ((i : ℝ) * x - a) < (q : ℝ) * ((q : ℝ) * x - p) :=
        mul_lt_mul_of_pos_left hcon hqR
      have h2 : (q : ℝ) * ((i : ℝ) * x - a) = (i : ℝ) * ((q : ℝ) * x - p) + (m : ℝ) := by
        push_cast [hm]
        ring
      linarith
    -- i < q
    have hiq : i < q := by
      by_contra hge
      push_neg at hge
      have hgeR : (q : ℝ) ≤ (i : ℝ) := by exact_mod_cast hge
      have : (q : ℝ) * ((q : ℝ) * x - p) ≤ (i : ℝ) * ((q : ℝ) * x - p) :=
        mul_le_mul_of_nonneg_right hgeR (le_of_lt hqx)
      linarith
    -- q'·(q·x − p) < 1
    have hsmall : (q' : ℝ) * ((q : ℝ) * x - p) < 1 := by
      have h1 : (q : ℝ) * ((q' : ℝ) * x) < (q : ℝ) * p' :=
        mul_lt_mul_of_pos_left hx' hqR
      have hdetR : (q : ℝ) * p' - (p : ℝ) * q' = 1 := by exact_mod_cast hdet
      nlinarith [h1]
    -- m·q' + i < q  over ℤ
    have hmq' : m * q' + i < q := by
      have hd : (m : ℝ) < ((q : ℝ) - i) * ((q : ℝ) * x - p) := by nlinarith [hstep0]
      have hA1 : (m : ℝ) * q' < ((q : ℝ) - i) * ((q' : ℝ) * ((q : ℝ) * x - p)) := by
        have := mul_lt_mul_of_pos_left hd hq'R
        nlinarith [this]
      have hiqR : (i : ℝ) < q := by exact_mod_cast hiq
      have hA2 : ((q : ℝ) - i) * ((q' : ℝ) * ((q : ℝ) * x - p)) < ((q : ℝ) - i) * 1 :=
        mul_lt_mul_of_pos_left hsmall (by linarith)
      have h5 : (m : ℝ) * q' < (q : ℝ) - i := by linarith
      have h6 : (m * q' : ℤ) < q - i := by exact_mod_cast h5
      omega
    have hmq'pos : 0 < m * q' + i := by
      have : 0 < m * q' := mul_pos (by omega) hq'
      omega
    -- the identity q·(i·p' − a·q') = m·q' + i
    have hkey : q * (i * p' - a * q') = m * q' + i := by
      have hdet' : q * p' = p * q' + 1 := by linarith
      calc q * (i * p' - a * q') = i * (q * p') - a * q * q' := by ring
        _ = i * (p * q' + 1) - a * q * q' := by rw [hdet']
        _ = q' * (p * i - a * q) + i := by ring
        _ = m * q' + i := by rw [hm]; ring
    have hdvd : q ∣ m * q' + i := ⟨i * p' - a * q', hkey.symm⟩
    have := Int.le_of_dvd hmq'pos hdvd
    omega

/-! ### Lemma B — first return from below -/

/-- **Lemma B (THM-637).** For `x` strictly inside the cell and `1 ≤ i ≤ k`:
`Int.fract (i·x) ≤ 1 + q'·x − p'`. -/
theorem lemmaB
    (hq : 0 < q) (hq' : 0 < q')
    (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hi1 : 1 ≤ i) (hik : i ≤ k)
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p') :
    Int.fract ((i : ℝ) * x) ≤ 1 + (q' : ℝ) * x - p' := by
  set a : ℤ := ⌊(i : ℝ) * x⌋ with ha
  have hfract : Int.fract ((i : ℝ) * x) = (i : ℝ) * x - a := rfl
  have hqu : (i : ℝ) * x < (a : ℝ) + 1 := Int.lt_floor_add_one _
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hq'R : (0 : ℝ) < (q' : ℝ) := by exact_mod_cast hq'
  have hiR : (1 : ℝ) ≤ (i : ℝ) := by exact_mod_cast hi1
  have hq'x : (0 : ℝ) < (p' : ℝ) - q' * x := by linarith
  -- Step 1: p'·i ≤ (a+1)·q', else (a+1, i) is a middle fraction.
  have hbq' : p' * i ≤ (a + 1) * q' := by
    by_contra hcon
    push_neg at hcon
    refine no_middle_fraction hq hq' hdet hsum hi1 hik (a + 1) ?_ hcon
    have hr : (p : ℝ) * i < ((a : ℝ) + 1) * q := by
      have h1 : (i : ℝ) * p < (i : ℝ) * ((q : ℝ) * x) :=
        mul_lt_mul_of_pos_left hx (by linarith)
      have h2 : ((i : ℝ) * x) * q < ((a : ℝ) + 1) * q :=
        mul_lt_mul_of_pos_right hqu hqR
      nlinarith [h1, h2]
    have hr' : (p * i : ℤ) < (a + 1) * q := by exact_mod_cast hr
    omega
  rcases eq_or_lt_of_le hbq' with heq | hlt
  · -- Equality: q' ∣ i, i = q'·t, a+1 = p'·t, fract = 1 − t·(p' − q'·x).
    have hdvd : q' ∣ i := by
      have hqi : q' ∣ p' * i := ⟨a + 1, by linear_combination heq⟩
      exact (coprime_q'_p' hdet).dvd_of_dvd_mul_left hqi
    obtain ⟨t, ht⟩ := hdvd
    have ht1 : 1 ≤ t := by
      by_contra hcon
      push_neg at hcon
      have ht0 : t ≤ 0 := by omega
      have : q' * t ≤ q' * 0 := mul_le_mul_of_nonneg_left ht0 (le_of_lt hq')
      omega
    have hat : a + 1 = p' * t := by
      have h0 : p' * (q' * t) = (a + 1) * q' := by rw [← ht]; exact heq
      have h1 : q' * (p' * t - (a + 1)) = 0 := by linear_combination h0
      have h2 := mul_eq_zero.mp h1
      rcases h2 with h | h
      · omega
      · omega
    rw [hfract]
    have htR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast ht1
    have hiRt : (i : ℝ) = (q' : ℝ) * t := by exact_mod_cast ht
    have haRt : (a : ℝ) + 1 = (p' : ℝ) * t := by exact_mod_cast hat
    have hval : (i : ℝ) * x - a = 1 - (t : ℝ) * ((p' : ℝ) - q' * x) := by
      rw [hiRt]
      linear_combination -haRt
    rw [hval]
    have hmul : 1 * ((p' : ℝ) - q' * x) ≤ (t : ℝ) * ((p' : ℝ) - q' * x) :=
      mul_le_mul_of_nonneg_right htR (le_of_lt hq'x)
    linarith
  · -- Strict: m' := (a+1)·q' − p'·i ≥ 1; mirrored divisibility contradiction.
    by_contra hcon
    push_neg at hcon
    rw [hfract] at hcon
    set b : ℤ := a + 1 with hb
    set m' : ℤ := b * q' - p' * i with hm'
    have hm'1 : 1 ≤ m' := by omega
    have hm'R : (1 : ℝ) ≤ (m' : ℝ) := by exact_mod_cast hm'1
    -- i·(p' − q'·x) + m' < q'·(p' − q'·x)
    have hstep0 : (i : ℝ) * ((p' : ℝ) - q' * x) + (m' : ℝ) < (q' : ℝ) * ((p' : ℝ) - q' * x) := by
      have h0 : (b : ℝ) - (i : ℝ) * x < (p' : ℝ) - (q' : ℝ) * x := by
        push_cast [hb]
        linarith
      have h1 : (q' : ℝ) * ((b : ℝ) - (i : ℝ) * x) < (q' : ℝ) * ((p' : ℝ) - q' * x) :=
        mul_lt_mul_of_pos_left h0 hq'R
      have h2 : (q' : ℝ) * ((b : ℝ) - (i : ℝ) * x)
          = (i : ℝ) * ((p' : ℝ) - q' * x) + (m' : ℝ) := by
        push_cast [hm', hb]
        ring
      linarith
    have hiq' : i < q' := by
      by_contra hge
      push_neg at hge
      have hgeR : (q' : ℝ) ≤ (i : ℝ) := by exact_mod_cast hge
      have : (q' : ℝ) * ((p' : ℝ) - q' * x) ≤ (i : ℝ) * ((p' : ℝ) - q' * x) :=
        mul_le_mul_of_nonneg_right hgeR (le_of_lt hq'x)
      linarith
    -- q·(p' − q'·x) < 1
    have hsmall : (q : ℝ) * ((p' : ℝ) - q' * x) < 1 := by
      have h1 : (q' : ℝ) * p < (q' : ℝ) * ((q : ℝ) * x) :=
        mul_lt_mul_of_pos_left hx hq'R
      have hdetR : (q : ℝ) * p' - (p : ℝ) * q' = 1 := by exact_mod_cast hdet
      nlinarith [h1]
    -- m'·q + i < q'  over ℤ
    have hm'q : m' * q + i < q' := by
      have hd : (m' : ℝ) < ((q' : ℝ) - i) * ((p' : ℝ) - q' * x) := by nlinarith [hstep0]
      have hA1 : (m' : ℝ) * q < ((q' : ℝ) - i) * ((q : ℝ) * ((p' : ℝ) - q' * x)) := by
        have := mul_lt_mul_of_pos_left hd hqR
        nlinarith [this]
      have hiq'R : (i : ℝ) < q' := by exact_mod_cast hiq'
      have hA2 : ((q' : ℝ) - i) * ((q : ℝ) * ((p' : ℝ) - q' * x)) < ((q' : ℝ) - i) * 1 :=
        mul_lt_mul_of_pos_left hsmall (by linarith)
      have h5 : (m' : ℝ) * q < (q' : ℝ) - i := by linarith
      have h6 : (m' * q : ℤ) < q' - i := by exact_mod_cast h5
      omega
    have hm'qpos : 0 < m' * q + i := by
      have : 0 < m' * q := mul_pos (by omega) hq
      omega
    -- the mirrored identity q'·(b·q − i·p) = m'·q + i
    have hkey : q' * (b * q - i * p) = m' * q + i := by
      have hdet' : p' * q = p * q' + 1 := by linarith
      calc q' * (b * q - i * p) = b * q * q' - i * (p * q') := by ring
        _ = b * q * q' - i * (p' * q - 1) := by rw [hdet']; ring_nf
        _ = q * (b * q' - p' * i) + i := by ring
        _ = m' * q + i := by rw [hm']; ring
    have hdvd : q' ∣ m' * q + i := ⟨b * q - i * p, hkey.symm⟩
    have := Int.le_of_dvd hm'qpos hdvd
    omega

/-! ### Lemma C — the gap bound via three exhibited indices -/

/-- **Lemma C (THM-637).** Above every config point `frac(i·x)` there is a config point
within `roof := (q·x − p) + (p' − q'·x)`: in lifted form, there exist `b ∈ [1,k]`, `n : ℤ`
with `0 < (b·x + n) − i·x ≤ roof`.  Hence no circular gap exceeds the roof. -/
theorem lemmaC
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1)
    (hqk : q ≤ k) (hq'k : q' ≤ k)
    (hi1 : 1 ≤ i) (hik : i ≤ k)
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p') :
    ∃ b n : ℤ, 1 ≤ b ∧ b ≤ k ∧
      0 < ((b : ℝ) * x + n) - (i : ℝ) * x ∧
      ((b : ℝ) * x + n) - (i : ℝ) * x ≤ ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x) := by
  have hroof1 : (0 : ℝ) < (q : ℝ) * x - p := by linarith
  have hroof2 : (0 : ℝ) < (p' : ℝ) - q' * x := by linarith
  rcases le_or_gt (i + q) k with hcase1 | hcase1
  · -- b = i + q, n = −p: distance q·x − p
    refine ⟨i + q, -p, by omega, hcase1, ?_, ?_⟩
    · push_cast
      linarith [hroof1]
    · push_cast
      linarith [hroof2]
  rcases le_or_gt (q' + 1) i with hcase2 | hcase2
  · -- b = i − q', n = p': distance p' − q'·x
    refine ⟨i - q', p', by omega, by omega, ?_, ?_⟩
    · push_cast
      linarith [hroof2]
    · push_cast
      linarith [hroof1]
  · -- window case k − q < i ≤ q': b = i + q − q', n = p' − p: distance = roof exactly
    refine ⟨i + q - q', p' - p, by omega, by omega, ?_, ?_⟩
    · push_cast
      linarith [hroof1, hroof2]
    · push_cast
      linarith

/-! ### The roof theorem, pointwise assembly -/

/-- The fract value at `i = q` is exactly the left roof leg: `⌊q·x⌋ = p`. -/
theorem fract_q_mul
    (hq : 0 < q) (hq' : 0 < q')
    (hdet : q * p' - p * q' = 1)
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p') :
    Int.fract ((q : ℝ) * x) = (q : ℝ) * x - p := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hq'R : (1 : ℝ) ≤ (q' : ℝ) := by exact_mod_cast hq'
  have hsmall : (q' : ℝ) * ((q : ℝ) * x - p) < 1 := by
    have h1 : (q : ℝ) * ((q' : ℝ) * x) < (q : ℝ) * p' := mul_lt_mul_of_pos_left hx' hqR
    have hdetR : (q : ℝ) * p' - (p : ℝ) * q' = 1 := by exact_mod_cast hdet
    nlinarith [h1]
  have hub : (q : ℝ) * x - p < 1 := by nlinarith [hsmall]
  have hfloor : ⌊(q : ℝ) * x⌋ = p := by
    apply Int.floor_eq_iff.mpr
    constructor
    · linarith
    · push_cast; linarith
  unfold Int.fract
  rw [hfloor]

/-- **The 0-gap (THM-637 roof theorem, lower half):** the open interval
`(q'·x − p', q·x − p)`, which contains `0` and has length `roof`, contains no config point:
for all `1 ≤ i ≤ k` and every integer lift `n`, `i·x + n` is not strictly inside. Combined
with `lemmaC` (no gap exceeds the roof) this gives `maxgap = gap∋0 = roof` pointwise. -/
theorem zero_gap_empty
    (hq : 0 < q) (hq' : 0 < q')
    (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hi1 : 1 ≤ i) (hik : i ≤ k)
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p')
    (n : ℤ) :
    ¬((q' : ℝ) * x - p' < (i : ℝ) * x + n ∧ (i : ℝ) * x + n < (q : ℝ) * x - p) := by
  rintro ⟨hlo, hhi⟩
  have hA := lemmaA hq hq' hdet hsum hi1 hik hx hx'
  have hB := lemmaB hq hq' hdet hsum hi1 hik hx hx'
  have hfr0 : (0 : ℝ) ≤ Int.fract ((i : ℝ) * x) := Int.fract_nonneg _
  have hfr1 : Int.fract ((i : ℝ) * x) < 1 := Int.fract_lt_one _
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hq'R : (0 : ℝ) < (q' : ℝ) := by exact_mod_cast hq'
  have hq1R : (1 : ℝ) ≤ (q : ℝ) := by exact_mod_cast hq
  have hq'1R : (1 : ℝ) ≤ (q' : ℝ) := by exact_mod_cast hq'
  -- interval endpoint bounds: q·x − p < 1 and p' − q'·x < 1
  have hsmall1 : (q' : ℝ) * ((q : ℝ) * x - p) < 1 := by
    have h1 : (q : ℝ) * ((q' : ℝ) * x) < (q : ℝ) * p' := mul_lt_mul_of_pos_left hx' hqR
    have hdetR : (q : ℝ) * p' - (p : ℝ) * q' = 1 := by exact_mod_cast hdet
    nlinarith [h1]
  have hsmall2 : (q : ℝ) * ((p' : ℝ) - q' * x) < 1 := by
    have h1 : (q' : ℝ) * p < (q' : ℝ) * ((q : ℝ) * x) := mul_lt_mul_of_pos_left hx hq'R
    have hdetR : (q : ℝ) * p' - (p : ℝ) * q' = 1 := by exact_mod_cast hdet
    nlinarith [h1]
  have hub : (q : ℝ) * x - p < 1 := by nlinarith [hsmall1]
  have hlb : (p' : ℝ) - q' * x < 1 := by nlinarith [hsmall2]
  -- i·x + n = fract(i·x) + z with z := n + ⌊i·x⌋
  set z : ℤ := n + ⌊(i : ℝ) * x⌋ with hz
  have hzval : (i : ℝ) * x + n = Int.fract ((i : ℝ) * x) + (z : ℝ) := by
    have hf : Int.fract ((i : ℝ) * x) = (i : ℝ) * x - ⌊(i : ℝ) * x⌋ := rfl
    rw [hf]
    push_cast [hz]
    ring
  rw [hzval] at hlo hhi
  -- z ∈ {0, −1}
  have hz0 : z ≤ 0 := by
    have : (z : ℝ) < 1 := by linarith
    have : z < 1 := by exact_mod_cast this
    omega
  have hzm1 : -1 ≤ z := by
    have : (-2 : ℝ) < (z : ℝ) := by linarith
    have : (-2 : ℤ) < z := by exact_mod_cast this
    omega
  interval_cases z
  · -- z = −1: fract > 1 + q'·x − p' contradicts lemmaB
    push_cast at hlo
    linarith
  · -- z = 0: fract < q·x − p contradicts lemmaA
    push_cast at hhi
    linarith

end FareyRoof
end LonelyRunner
