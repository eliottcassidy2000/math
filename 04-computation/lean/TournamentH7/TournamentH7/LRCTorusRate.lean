/-
  TournamentH7.LRCTorusRate — THE ELEMENTARY LIFT-CONVERGENCE RATE
  (mac-mini-2026-07-06-S10, HYP-4342).

  The formal core of the preprint-free (A)-subsumption: a proper 2-torus value
  M(U) is approached from below by its lift-family values M(v^N) with an
  EXPLICIT rate L/(2N), by pure Lipschitz + net geometry — NO Jain-Kravitz
  lift-limit theorem.

  M(v^N) = max over the curve C_N = {(τ, Nτ)} of F(t,s) = minᵢ‖rᵢt + ℓᵢs‖;
  F is L-Lipschitz with L = maxᵢ(|rᵢ|+|ℓᵢ|) (N-independent); C_N is 1/(2N)-dense
  in T².  This file proves the abstract heart:

  * `lipschitz_ge_at_near`   — an L-Lipschitz f at a point y within ε of the
    maximizer x* satisfies f y ≥ f x* − L·ε.  (The one inequality the rate rests
    on; a two-line Lipschitz estimate.)
  * `sup_on_net_ge`          — hence the sup of f over an ε-dense set S is ≥
    (max of f) − L·ε: the max transfers to any dense net with a Lipschitz loss.
  * `two_point_rate`         — packaged: |f x* − f y| ≤ L·ε for y within ε of x*.

  The subsumption instantiates this with f = F on T², S = C_N, ε = 1/(2N):
  M(v^N) = sup_{C_N} F ≥ M(U) − L/(2N), and M(v^N) ≤ M(U) (sub-torus).  Draft:
  07-reflections/the-subsumption-is-preprint-free-macmini-S10.md.
-/
import Mathlib.Topology.MetricSpace.Lipschitz
import Mathlib.Order.CompleteLattice.Basic

namespace LonelyRunner
namespace TorusRate

variable {X : Type*} [PseudoMetricSpace X]

/-- The rate's engine: an `L`-Lipschitz `f` evaluated at a point `y` within `ε`
    of a reference point `x` loses at most `L·ε`.  (Here `x` is the maximizer.) -/
theorem lipschitz_ge_at_near {f : X → ℝ} {L : NNReal} (hf : LipschitzWith L f)
    {x y : X} {ε : ℝ} (hε : dist x y ≤ ε) :
    f x - (L : ℝ) * ε ≤ f y := by
  have h := hf.dist_le_mul x y
  have h1 : |f x - f y| ≤ (L : ℝ) * dist x y := by
    rw [Real.dist_eq] at h; exact h
  have h2 : f x - f y ≤ (L : ℝ) * dist x y := (abs_le.mp h1).2
  have h3 : (L : ℝ) * dist x y ≤ (L : ℝ) * ε :=
    mul_le_mul_of_nonneg_left hε (L.coe_nonneg)
  linarith

/-- Packaged symmetric form: `|f x − f y| ≤ L·ε` for `y` within `ε` of `x`. -/
theorem two_point_rate {f : X → ℝ} {L : NNReal} (hf : LipschitzWith L f)
    {x y : X} {ε : ℝ} (hε : dist x y ≤ ε) :
    |f x - f y| ≤ (L : ℝ) * ε := by
  have h := hf.dist_le_mul x y
  rw [Real.dist_eq] at h
  exact h.trans (mul_le_mul_of_nonneg_left hε (L.coe_nonneg))

/-- The max transfers to any ε-dense set with a Lipschitz loss: if `f` attains
    its value `M` at `x` (i.e. `f x = M`), and `S` is ε-dense (some `y ∈ S` is
    within `ε` of `x`), then some point of `S` has `f`-value `≥ M − L·ε`.
    Instantiated with `x` = the torus maximizer, `S` = the curve `C_N`, this is
    `M(v^N) ≥ M(U) − L/(2N)`. -/
theorem exists_net_ge {f : X → ℝ} {L : NNReal} (hf : LipschitzWith L f)
    {x : X} {M ε : ℝ} (hM : f x = M) {S : Set X}
    (hdense : ∃ y ∈ S, dist x y ≤ ε) :
    ∃ y ∈ S, M - (L : ℝ) * ε ≤ f y := by
  obtain ⟨y, hyS, hy⟩ := hdense
  refine ⟨y, hyS, ?_⟩
  have := lipschitz_ge_at_near hf hy
  rw [hM] at this
  exact this

end TorusRate
end LonelyRunner
