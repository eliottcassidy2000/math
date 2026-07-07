/-
  TournamentH7.LRCSmallModFloor — THE SMALL-MODULUS 1/q FLOOR CERTIFICATE
  (kind-pasteur-2026-07-06-S44, HYP-4597).

  Companion to `LRCMod25Floor`.  The case-2 residual of (G) — the pair-blocking
  rigidity — reduces to a covering system: every non-AP pair-blocker is loose at
  some small modulus `q`.  For `q ≤ 12` the clearance band `2q/25 < 1`, so clearing
  at `q` needs only that every rotated speed avoids residue `0` — i.e. a rotation `c`
  with `q ∤ (v i · c)` for all `i`.  At `c = 1` this is just "no speed is a multiple
  of `q`", and it forces the far stronger floor `M ≥ 1/q` (not merely `2/25`).

  For `q ∈ {7,9,10,11,12}` we have `1/q > 2/25`, so:

      no speed divisible by `q`  ⟹  `M ≥ 1/q > 2/25`  (the family is LOOSE).

  This closes, by one margin certificate each, every 12-speed family that misses a
  multiple of some `q ∈ {7,…,12}` — the bulk of the covering system.  The remaining
  families carry a multiple of every such `q` (highly divisible) and are handled by
  the `q ≥ 13` layer of the covering.

  A direct instance of `rational_point_margin` (`s = q`, `μ = 1`); kernel-pure.
-/
import Mathlib
import TournamentH7.LRCHarmonicGate

namespace LonelyRunner
namespace SmallModFloor

open HarmonicGate

/-- **The `1/q` floor from a zero-avoiding rotation.**  If some `c` keeps every
`v i · c` off residue `0` mod `q` (`1 ≤ (v i·c) % q ≤ q−1`), then at `t = c/q` every
runner is `≥ 1/q` from every integer.  Instance of `rational_point_margin` (`μ=1`). -/
theorem zero_avoid_floor {ι : Type*} (v : ι → ℤ) (q c : ℤ) (hq : 0 < q)
    (hcov : ∀ i, 1 ≤ (v i * c) % q ∧ (v i * c) % q ≤ q - 1) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / q ≤ |(v i : ℝ) * ((c : ℝ) / q) - m| := by
  intro i m
  have h := rational_point_margin v c q 1 hq hcov i m
  simpa using h

/-- **No multiple of `q` ⟹ `1/q` floor** (at `c = 1`).  If no speed is divisible by
`q` (each `v i % q ∈ [1, q−1]`), the family has margin `≥ 1/q` at `t = 1/q`. -/
theorem no_multiple_floor {ι : Type*} (v : ι → ℤ) (q : ℤ) (hq : 0 < q)
    (hnd : ∀ i, 1 ≤ v i % q ∧ v i % q ≤ q - 1) :
    ∀ i, ∀ m : ℤ, (1 : ℝ) / q ≤ |(v i : ℝ) * ((1 : ℝ) / q) - m| := by
  have h := zero_avoid_floor v q 1 hq (by simpa using hnd)
  simpa using h

/-- **The `q = 12` loose certificate** (`1/12 > 2/25`).  No speed divisible by 12 ⟹
every runner is `≥ 1/12 > 2/25` from every integer at `t = 1/12` — LOOSE. -/
theorem loose_of_no_multiple_12 {ι : Type*} (v : ι → ℤ)
    (hnd : ∀ i, 1 ≤ v i % 12 ∧ v i % 12 ≤ 11) :
    ∃ t : ℝ, ∀ i, ∀ m : ℤ, (2 : ℝ) / 25 ≤ |(v i : ℝ) * t - m| := by
  refine ⟨(1 : ℝ) / 12, fun i m => ?_⟩
  have h := no_multiple_floor v 12 (by norm_num) (by simpa using hnd) i m
  have h2 : (2 : ℝ) / 25 ≤ (1 : ℝ) / 12 := by norm_num
  exact le_trans h2 (by simpa using h)

#print axioms zero_avoid_floor
#print axioms no_multiple_floor
#print axioms loose_of_no_multiple_12

end SmallModFloor
end LonelyRunner
