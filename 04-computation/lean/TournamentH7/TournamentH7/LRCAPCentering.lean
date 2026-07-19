/-
  TournamentH7.LRCAPCentering — the centering witness for arithmetic-progression loneliness
  (boxeph-2026-07-18-S118).

  Closes the spread case (`d ≥ 17`) that S117 left open in the AP-rigidity face of HYP-4382.

  IDEA (rescale by the inverse of the common difference).  For a primitive 12-term AP
  `C = {a + d·k : k = 0..11}` with `d` ODD, set `q = 2a + 11d` and `p ≡ d⁻¹ (mod q)`.  Then at
  `t = p/q` each speed reduces to `(a+d·k)·p ≡ a·p + k (mod q)` — a run of 12 CONSECUTIVE residues —
  and `2·a·p ≡ -11 (mod q)` forces that run to be `{(q-11)/2, …, (q+11)/2}`, 12 consecutive integers
  SYMMETRIC about `q/2`.  The nearest multiples of `q` are `0` and `q`, both at distance `(q-11)/2`
  from the run, so every runner sits `≥ (q-11)/(2q)` from the integers:
      `min_k ‖(a+d·k)·(p/q)‖ = (q-11)/(2q) = 1/2 − 11/(2q) > 1/13   ⟺   q > 13`.
  Since `q = 2a+11d = 13 ⟺ a=d=1`, every primitive AP other than `{1,…,12}` has `M > 1/13`; hence the
  only tight 12-term APs are the dilates `c·{1,…,12}`.  (For `d` even, `gcd(a,d)=1` makes all terms odd,
  so `t = 1/2` already gives `M = 1/2`.)

  This file formalizes the reusable HEART:
   * `centered_block_far` — the integer core: a residue in the centered band `[(q-11)/2, (q+11)/2]`
     is at distance `≥ (q-11)/2` from every multiple of `q`.
   * `centered_block_witness` — the real witness: if every speed's numerator `N i` reduces into the
     centered band mod `q`, then `t = 1/q` (applied to the reduced numerators) is a common
     `(q-11)/(2q)`-lonely time — the `≥` half that powers the AP rigidity.

  Companion to `LRCMod13Blocking` (the `p = 13` middle-band sieve).  Same shape: an omega-checked
  integer inequality lifted to a real circle-distance bound by `abs_div`/`gcongr`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- **Integer core (PROVED).**  For `q ≥ 13`, a residue `r` in the centered band
`(q-11)/2 ≤ r ≤ (q+11)/2` is at distance `≥ (q-11)/2` from every multiple `q·k`:
    `q - 11 ≤ 2·|r − q·k|`.
The band sits inside `(0,q)`, so the two nearest multiples are `0` (`k=0`, distance `r ≥ (q-11)/2`)
and `q` (`k=1`, distance `q − r ≥ (q-11)/2`); every other multiple is farther. -/
theorem centered_block_far (q r k : ℤ) (hq : 13 ≤ q)
    (hlo : q - 11 ≤ 2 * r) (hhi : 2 * r ≤ q + 11) :
    q - 11 ≤ 2 * |r - q * k| := by
  rcases le_or_gt k 0 with hk | hk
  · -- k ≤ 0 ⟹ q·k ≤ 0 ⟹ r − q·k ≥ r ≥ (q-11)/2 ≥ 0
    have hqk : q * k ≤ 0 := by nlinarith [hq, hk]
    rw [abs_of_nonneg (by omega)]
    omega
  · -- k ≥ 1 ⟹ q·k ≥ q ⟹ r − q·k ≤ r − q < 0, so |·| = q·k − r ≥ q − r ≥ (q-11)/2
    have hqk : q ≤ q * k :=
      by nlinarith [mul_nonneg (by omega : (0:ℤ) ≤ q) (by omega : (0:ℤ) ≤ k - 1)]
    rw [abs_of_nonpos (by omega)]
    omega

/-- **The centering witness (PROVED).**  Let `q ≥ 13` and let each speed contribute a reduced
numerator `N i` (in the AP application `N i = (a + d·k)·d⁻¹`, so `N i ≡ a·p + k (mod q)`).  If every
`N i` lands in the centered band mod `q`, then at `t = 1/q` every runner sits at circle-distance
`≥ (q-11)/(2q)` from the integers:
    `∀ i m, (q-11)/(2q) ≤ ‖N i / q − m‖`.
Hence `M(C) ≥ (q-11)/(2q) = 1/2 − 11/(2q)`, which exceeds `1/13` for every `q > 13`. -/
theorem centered_block_witness {ι : Type*} (q : ℤ) (hq : 13 ≤ q)
    (N : ι → ℤ)
    (hcen : ∀ i, q - 11 ≤ 2 * (N i % q) ∧ 2 * (N i % q) ≤ q + 11) :
    ∀ i, ∀ m : ℤ, ((q : ℝ) - 11) / (2 * q) ≤ |(N i : ℝ) / q - m| := by
  intro i m
  obtain ⟨hlo, hhi⟩ := hcen i
  have hqpos : (0 : ℝ) < q := by exact_mod_cast (by omega : (0:ℤ) < q)
  have hqne : (q : ℝ) ≠ 0 := ne_of_gt hqpos
  -- the real value equals the integer numerator over q
  have hreal : (N i : ℝ) / q - m = ((N i - q * m : ℤ) : ℝ) / q := by
    rw [eq_div_iff hqne]; push_cast; field_simp
  -- decompose N i − q·m into `residue − q·(shift)`, then apply the integer core
  -- (q is a variable, so feed omega the div/mod identity and the product expansion)
  have hdecomp : N i - q * m = (N i % q) - q * (m - N i / q) := by
    have hdm : q * (N i / q) + N i % q = N i := Int.ediv_add_emod (N i) q
    have hexp : q * (m - N i / q) = q * m - q * (N i / q) := by ring
    omega
  have hint : q - 11 ≤ 2 * |N i - q * m| := by
    rw [hdecomp]; exact centered_block_far q (N i % q) (m - N i / q) hq hlo hhi
  -- lift to reals: |N i − q·m| ≥ (q-11)/2
  have hintR : ((q : ℝ) - 11) / 2 ≤ |((N i - q * m : ℤ) : ℝ)| := by
    have h2 : ((q : ℝ) - 11) ≤ 2 * ((|N i - q * m| : ℤ) : ℝ) := by exact_mod_cast hint
    rw [Int.cast_abs] at h2
    linarith
  -- divide by q > 0
  rw [hreal, abs_div, show |(q : ℝ)| = q from abs_of_pos hqpos, le_div_iff₀ hqpos]
  have e : ((q : ℝ) - 11) / (2 * q) * q = ((q : ℝ) - 11) / 2 := by
    rw [div_mul_eq_mul_div, mul_div_mul_right _ _ hqne]
  linarith [hintR, e]

end LonelyRunner

#print axioms LonelyRunner.centered_block_far
#print axioms LonelyRunner.centered_block_witness
