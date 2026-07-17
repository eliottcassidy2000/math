/- Thm878ClockTable.lean — klein-2026-07-16-S317.
   THM-878 rung three, the kernel-pure table certificate.

   For the tight AP {1..13}, the pair-overlap functional at x = a/q is
     S₂(a/q) = Σ_{d=1}^{12} (13−d) · max(0, 1/7 − ‖d·a/q‖),
   and THM-878's deficit is D(q) = (1/φ(q)) Σ_{a primitive} S₂(a/q) − 6/7 ≥ 0,
   vanishing exactly at the clock moduli q ∈ {7, 13, 14}.

   Clearing denominators: 7q·S₂(a/q) = Σ_d (13−d)·(q − 7·min(r, q−r)) with
   r = d·a mod q (ℕ-subtraction implements the tent clamp exactly), so
     D(q) = 0  ⟺  clockSum q = 6·q·φ(q),
     D(q) ≥ 0  ⟺  clockSum q ≥ 6·q·φ(q),
   both PURE ℕ facts, decided by the kernel below for 2 ≤ q ≤ 60 (the analytic
   q ≥ 15 tail is THM-878's one-line a = 1 bound on paper; the table covers the
   entire finite window with margin).  No sorries, no native_decide. -/
import Mathlib.Data.Nat.Totient
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

namespace LonelyRunner
namespace LRC14
namespace ClockTable

/-- `7q · S₂(a/q)` for the tight AP `{1..13}`, as a pure ℕ expression:
    the pair at difference `d` contributes `(13−d)·(q − 7·min(r, q−r))`,
    `r = d·a mod q` (ℕ-subtraction clamps the tent at zero). -/
def s2Scaled (q a : ℕ) : ℕ :=
  ∑ d ∈ Finset.Icc 1 12, (13 - d) * (q - 7 * min ((d * a) % q) (q - (d * a) % q))

/-- The primitive-class sum `7q·φ(q)·A(q)`. -/
def clockSum (q : ℕ) : ℕ :=
  ∑ a ∈ (Finset.range q).filter (fun a => Nat.gcd a q = 1), s2Scaled q a

/-- **THM-878, the clock table (kernel-decided)**: on `2 ≤ q ≤ 60` the
    Ramanujan primitive-mean deficit vanishes exactly at the clock moduli
    `q ∈ {7, 13, 14}`. -/
theorem clock_table : ∀ q ∈ Finset.Icc 2 60,
    (clockSum q = 6 * q * Nat.totient q ↔ q = 7 ∨ q = 13 ∨ q = 14) := by
  decide

/-- **The FT-floor face of the table**: the deficit is nonnegative on the whole
    window (LEM-020's pointwise floor, averaged). -/
theorem deficit_nonneg_table : ∀ q ∈ Finset.Icc 2 60,
    6 * q * Nat.totient q ≤ clockSum q := by
  decide

/-- The three clock values, exactly: at the moduli the mean sits ON the floor. -/
theorem clock_values :
    clockSum 7 = 6 * 7 * Nat.totient 7 ∧
    clockSum 13 = 6 * 13 * Nat.totient 13 ∧
    clockSum 14 = 6 * 14 * Nat.totient 14 := by
  decide

end ClockTable
end LRC14
end LonelyRunner
