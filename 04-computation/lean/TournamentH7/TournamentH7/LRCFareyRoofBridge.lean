/-
  TournamentH7.LRCFareyRoofBridge — bridging the Farey roof (opus-S135, THM-637) to
  the tail-diameter good set (mac-mini-S42 / monad-S2 HYP-4817).  boxeph-2026-07-07.

  `LRCFareyRoof.zero_gap_empty` proves that on the open Farey-k cell `(p/q, p'/q')`
  no config point of the AP `{1,…,k}` lies strictly inside the open interval
  `(q'x − p', qx − p)`, which contains `0` and has length
  `roof := (qx − p) + (p' − q'x)`.

  `LRCTailDiameter.Good θ E = {x | ∃ a, ∀ e ∈ E, θ < fract(e·x − a)}` is the good set
  (some closed θ-arc is empty of the orbit); its measure `muGood` is the density
  `μ_θ`, and the whole diameter route is GREEN conditional on the AP₇₆ certificate
  `muGood (1/7) {0..75} ≥ 2314528732/40290957525`.

  THIS FILE supplies the missing pointwise link:

    roof(x) > θ   ⟹   x ∈ Good θ (Finset.Icc 1 k).

  Consequence: `{x in a Farey-k cell : roof > θ} ⊆ Good θ (AP_k)`, so
  `muGood θ (AP_k) ≥ meas{ roof > θ }` — reducing the certificate to a PURE
  real-superlevel (Farey-sum) measure computation, with no orbit reasoning left.

  Proof idea: put the θ-arc's left end at `a := (q'x − p') + (roof − θ)/2`, the
  midpoint offset that makes the closed θ-arc `[a, a+θ]` sit STRICTLY inside the
  empty roof-interval.  Any orbit point in `(a, a+θ]` would then land strictly
  inside `(q'x − p', qx − p)`, contradicting `zero_gap_empty`.

  Kernel-pure: no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LRCFareyRoof
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace FareyRoofBridge

open TailDiameter

/-- **The roof→good bridge.**  On the open Farey-`k` cell `(p/q, p'/q')` (encoded in
cleared form `p < q·x`, `q'·x < p'`, determinant `q·p' − p·q' = 1`, `k < q + q'`),
if the roof `(q·x − p) + (p' − q'·x)` exceeds `θ ≥ 0`, then `x` is in the AP-`k`
good set: some closed θ-arc is empty of the orbit `{frac(i·x) : 1 ≤ i ≤ k}`. -/
theorem good_of_roof_gt {p q p' q' k : ℤ} {x θ : ℝ}
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1) (hsum : k < q + q')
    (hx : (p : ℝ) < q * x) (hx' : (q' : ℝ) * x < p')
    (hroof : θ < ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x)) :
    x ∈ Good θ (Finset.Icc (1 : ℤ) k) := by
  -- roof interval endpoints  a0 = q'x − p' < 0 < b0 = qx − p
  set a0 : ℝ := (q' : ℝ) * x - p' with ha0
  set b0 : ℝ := (q : ℝ) * x - p with hb0
  have hroof' : θ < b0 - a0 := by rw [ha0, hb0]; linarith [hroof]
  -- witness: left end of the θ-arc, pushed in by half the slack
  refine ⟨a0 + (b0 - a0 - θ) / 2, ?_⟩
  intro e he
  rcases Finset.mem_Icc.mp he with ⟨he1, hek⟩
  set a : ℝ := a0 + (b0 - a0 - θ) / 2 with ha
  by_contra hcon
  have hcon' : Int.fract ((e : ℝ) * x - a) ≤ θ := not_lt.mp hcon
  have hf0 : (0 : ℝ) ≤ Int.fract ((e : ℝ) * x - a) := Int.fract_nonneg _
  -- n := −⌊e·x − a⌋ realizes  e·x + n = a + fract(e·x − a)
  set n : ℤ := -⌊(e : ℝ) * x - a⌋ with hn
  have hfl : ((⌊(e : ℝ) * x - a⌋ : ℤ) : ℝ) + Int.fract ((e : ℝ) * x - a)
      = (e : ℝ) * x - a := Int.floor_add_fract _
  have hval : (e : ℝ) * x + (n : ℝ) = a + Int.fract ((e : ℝ) * x - a) := by
    rw [hn]; push_cast; linarith [hfl]
  -- a + fract ∈ (a0, b0): left end from slack>0 and fract≥0; right end from fract≤θ
  have hlo : a0 < (e : ℝ) * x + (n : ℝ) := by rw [hval, ha]; nlinarith [hf0, hroof']
  have hhi : (e : ℝ) * x + (n : ℝ) < b0 := by rw [hval, ha]; nlinarith [hcon', hroof']
  exact FareyRoof.zero_gap_empty hq hq' hdet hsum he1 hek hx hx' n ⟨hlo, hhi⟩

/-- Set form: within a fixed Farey-`k` cell, the roof-superlevel set is contained
in the AP-`k` good set.  (The cell hypotheses are carried as a predicate on `x`.) -/
theorem roof_superlevel_subset_good
    {p q p' q' k : ℤ} {θ : ℝ}
    (hq : 0 < q) (hq' : 0 < q') (hdet : q * p' - p * q' = 1) (hsum : k < q + q') :
    {x : ℝ | (p : ℝ) < q * x ∧ (q' : ℝ) * x < p' ∧
        θ < ((q : ℝ) * x - p) + ((p' : ℝ) - q' * x)} ⊆ Good θ (Finset.Icc (1 : ℤ) k) := by
  rintro x ⟨hx, hx', hroof⟩
  exact good_of_roof_gt hq hq' hdet hsum hx hx' hroof

end FareyRoofBridge
end LonelyRunner
