/-
  TournamentH7.LRCSmoothBridge — the continuum reformulation bridge on the smooth surrogate
  (kind-pasteur-2026-07-09-S112).

  The density floor gives `E_x[W] = ∫₀¹ W > 0` for the smooth uncovered-measure surrogate
  `W = Σ_i (gap_i − 1/7)_+` (opus-S170; THM-530/661; kps-S108 verified `≈ 0.13–0.14 > 0`).  The
  reformulation bridge turns that into a lonely time.  This file formalizes the **genuinely new core**
  of that bridge — the step that is *why the smooth surrogate works* — the **measure-to-point**
  passage:

      a nonnegative `W` with **positive integral** is **strictly positive at a point**.

  This is the desingularization of kps-S107: the sharp indicator's good set is bounded by
  grid-invisible rational pinches (measure-zero nodes, `LRCPinch`), so *sampling* it on any ruler grid
  is unreliable.  The smooth surrogate sidesteps grids entirely — a positive **measure** `∫W > 0`
  hands you a positive **point** `W(x) > 0` directly (`exists_pos_of_integral_pos`), and a positive
  point of `W = Σ(gap−1/7)_+` forces a genuine wide gap `gap_i > 1/7` (`exists_gap_gt_of_smoothW_pos`)
  — i.e. a good period at the continuum, no drift, no window.

  Composed with the reformulation (`good x ⟹ ∃ lonely instant`, the shared Part-A hypothesis, now
  stated at the continuum) and kps-S99b's `Mreach_ge_of_lonely_instant`, this yields `Mreach ≥ 1/14`
  (`mreach_ge_of_smooth_surrogate`).  Self-contained apart from `LRCReachWitness`.
-/
import Mathlib
import TournamentH7.LRCReachWitness

namespace LonelyRunner
namespace LRC14Concrete

open scoped BigOperators

/-- **Measure → point (the smooth-surrogate core).**  A nonnegative function whose interval integral
over `[0,1]` is strictly positive is strictly positive at some point.  Contrapositive: if `W ≤ 0`
everywhere then (being `≥ 0`) `W ≡ 0`, so `∫ W = 0`.  This is the desingularization at the heart of
the bridge — a positive *measure* yields a positive *point* with no grid sampling, so the
grid-invisible pinches of the sharp good set (kps-S107, `LRCPinch`) are bypassed. -/
theorem exists_pos_of_integral_pos {W : ℝ → ℝ} (hW : ∀ x, 0 ≤ W x)
    (hint : 0 < ∫ x in (0:ℝ)..1, W x) : ∃ x, 0 < W x := by
  by_contra h
  push_neg at h
  have hz : W = fun _ => (0 : ℝ) := funext fun x => le_antisymm (h x) (hW x)
  rw [hz] at hint
  simp at hint

/-- **A positive value of the surrogate forces a wide gap.**  With `W x = Σ_{i} (gap i x − 1/7)_+`
(sum of positive parts of the tooth gaps minus the threshold `1/7`), if `W x > 0` then some gap
exceeds `1/7`: `∃ i, 1/7 < gap i x`.  Otherwise every positive part is `0` and the sum vanishes.
This is the continuum "good period at `x`": a genuine clearing gap wider than `1/7` (half of it,
`> 1/14`, is where a lonely phase sits). -/
theorem exists_gap_gt_of_smoothW_pos {ι : Type*} (s : Finset ι) (gap : ι → ℝ → ℝ) (x : ℝ)
    (hpos : 0 < ∑ i ∈ s, max (gap i x - 1 / 7) 0) :
    ∃ i ∈ s, (1 : ℝ) / 7 < gap i x := by
  by_contra h
  push_neg at h
  have hz : ∀ i ∈ s, max (gap i x - 1 / 7) 0 = 0 := by
    intro i hi
    exact max_eq_right (by linarith [h i hi])
  rw [Finset.sum_congr rfl hz] at hpos
  simp at hpos

/-- **`nearInt y ≥ c` on the central band `[c, 1−c]`.**  A helper: on `[c, 1−c] ⊆ (0,1)` the fractional
part is the identity, so `nearInt y = min(y, 1−y) ≥ c`. -/
theorem nearInt_ge_of_mem {c y : ℝ} (hc : 0 < c) (h1 : c ≤ y) (h2 : y ≤ 1 - c) :
    c ≤ nearInt y := by
  have hy0 : (0 : ℝ) ≤ y := le_trans (le_of_lt hc) h1
  have hy1 : y < 1 := by linarith
  have hfr : Int.fract y = y := Int.fract_eq_self.mpr ⟨hy0, hy1⟩
  unfold nearInt
  rw [hfr]
  exact le_min h1 (by linarith)

/-- **The drift-free observer (the geometric heart of the reconstruction).**  If the cluster positions
`q i` all lie in an arc `[a, a+w]` of length `w ≤ 1 − 2c` — equivalently the complementary circular
**gap is `≥ 2c`** — then the gap midpoint `φ = a + (1+w)/2` clears every runner: `nearInt(φ − q i) ≥ c`
for all `i`.  With `c = 1/14` this is `w ≤ 6/7` (gap `≥ 1/7`) ⟹ observer at `≥ 1/14`.

This is the `good period ⟹ lonely phase` step **exactly, with no drift** — the `Vmax → ∞` continuum
limit of the slow-fast embedding (klein-S205 / kps-S106 carry the finite-`Vmax` drift `e·φ/Vmax`; here
it has vanished).  Mac-mini's local refutation of hembed was precisely the finite-`Vmax` drift; at the
continuum the observer is clean.  What remains of the reformulation is only the Kronecker realization —
that the cluster `{e_i·τ}` can be simultaneously confined to such an arc and that `φ` lifts to a genuine
`τ`; that residue is the named hypothesis in the bridge below. -/
theorem observer_of_confined {ι : Type*} (q : ι → ℝ) (a w c : ℝ) (hc : 0 < c)
    (hw : w ≤ 1 - 2 * c) (hconf : ∀ i, q i ∈ Set.Icc a (a + w)) :
    ∀ i, c ≤ nearInt ((a + (1 + w) / 2) - q i) := by
  intro i
  obtain ⟨hlo, hhi⟩ := hconf i
  exact nearInt_ge_of_mem hc (by linarith) (by linarith)

/-- **Observer at the LRC(14) threshold (drift-free).**  Specialising `observer_of_confined` to
`c = 1/14`: a cluster confined to an arc of length `≤ 6/7` (circular gap `≥ 1/7`) admits a phase `φ`
clearing every runner by `≥ 1/14` — the exact LRC(14) loneliness margin, with no drift.  This is the
positive-measure good period turned into an explicit lonely phase in the `Vmax → ∞` continuum. -/
theorem observer_at_threshold {ι : Type*} (q : ι → ℝ) (a w : ℝ) (hw : w ≤ 6 / 7)
    (hconf : ∀ i, q i ∈ Set.Icc a (a + w)) :
    ∃ φ : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt (φ - q i) :=
  ⟨a + (1 + w) / 2, observer_of_confined q a w (1 / 14) (by norm_num) (by linarith) hconf⟩

/-- **The continuum reformulation bridge on the smooth surrogate.**  Given
  * the smooth surrogate `W ≥ 0` with positive integral `∫₀¹ W > 0` (the density floor
    `E_x[W] ≥ m_P > 0`, THM-530/661; kps-S108), and
  * the reformulation `(∃ x, 0 < W x) → ∃ lonely instant` — a positive point of the surrogate is a
    good period, and a good period reconstructs a lonely time (the shared Part-A hypothesis, stated
    here at the continuum, where there is no drift and no window),

the runner set is lonely: `Mreach v ≥ 1/14`.  The measure-to-point core `exists_pos_of_integral_pos`
is the bridge's new content; the reformulation is the sole remaining hypothesis. -/
theorem mreach_ge_of_smooth_surrogate (v : Fin 13 → ℤ) (W : ℝ → ℝ)
    (hW : ∀ x, 0 ≤ W x) (hint : 0 < ∫ x in (0:ℝ)..1, W x)
    (hrefl : (∃ x, 0 < W x) → ∃ τ : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt ((v i : ℝ) * τ)) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  Mreach_ge_of_lonely_instant v (hrefl (exists_pos_of_integral_pos hW hint))

/-- **The bridge, packaged from the density floor.**  Same conclusion, taking the positivity of the
integral through the density-floor constants: `0 < m_P ≤ E_x[W] = ∫₀¹ W`.  This is the exact shape the
Bonferroni floor (THM-530: `m_P = 14249/252252`) and its transport to the smooth surrogate (THM-661)
supply, so the only user-facing hypothesis left is the reformulation. -/
theorem mreach_ge_of_density_floor (v : Fin 13 → ℤ) (W : ℝ → ℝ) (mP : ℝ)
    (hW : ∀ x, 0 ≤ W x) (hmP : 0 < mP) (hfloor : mP ≤ ∫ x in (0:ℝ)..1, W x)
    (hrefl : (∃ x, 0 < W x) → ∃ τ : ℝ, ∀ i, (1 : ℝ) / 14 ≤ nearInt ((v i : ℝ) * τ)) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  mreach_ge_of_smooth_surrogate v W hW (lt_of_lt_of_le hmP hfloor) hrefl

end LRC14Concrete
end LonelyRunner
