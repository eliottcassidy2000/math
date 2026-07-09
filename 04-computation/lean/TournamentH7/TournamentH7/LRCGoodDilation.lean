/-
  TournamentH7.LRCGoodDilation — DILATION INVARIANCE of the covering good-set
  (opus-2026-07-08-S162).

  The covering-side companion to `LRCDilationInvariance.iSup_margin_const_mul` (which proves
  dilation invariance of the loneliness minimax `M(c·v)=M(v)`).  Here the same gauge fact is
  proved for the GOOD SET `Good θ E` / its measure `muGood θ E` of `LRCTailDiameter`, and for
  the orbit-`Int.fract` form directly:
      Good θ (c·E) = (x ↦ c·x)⁻¹ (Good θ E),   and   Good θ E  is 1-periodic.

  Role.  This is the structural fact underpinning the k=11 covering-tail reduction to the
  *longest-AP* axis (opus-S155–S161): `muGood θ` (equivalently the degree-3 floor `D3`) depends
  only on the dilation class of the speed set, so every tail shape may be taken PRIMITIVE and
  the extremal analysis lives on the dilation-invariant longest AP.  The fixed-window "cluster
  size" of the refuted LEM-009 argument is NOT dilation-invariant; `muGood`/`D3` IS — the
  exact counterexample `(0,3,6,8,9,12,15,18,21,24,27)` (a *dilated* AP₁₀ + interior point,
  `D3 = 0.452986 < 0.4646`) is a dilate of the compact minimizer, invisible to window-cluster
  but not to `muGood` (MISTAKE-126 / CASE-tail-D3-min-is-not-block-outlier).

  Kernel-pure; mirrors the translation-invariance proofs in `LRCTailDiameter`.
-/
import TournamentH7.LRCTailDiameter

namespace LonelyRunner
namespace TailDiameter

open scoped ENNReal

/-! ### Dilation of the orbit (exact, elementary) -/

/-- **Orbit dilation.**  The `(c·e)`-phase at witness `a` and slow phase `x` equals the
`e`-phase at witness `a` and sped-up phase `c·x` — the same real number `(e)(c·x) − a`
inside `Int.fract`.  So an arc is empty of the `c·E`-orbit at `x` iff it is empty of the
`E`-orbit at `c·x`. -/
theorem emptyArc_dilate {θ : ℝ} (E : Finset ℤ) (c : ℤ) (x a : ℝ) :
    EmptyArc θ (E.image (fun e => c * e)) x a ↔ EmptyArc θ E ((c : ℝ) * x) a := by
  constructor
  · intro h e he
    have hmem : c * e ∈ E.image (fun e => c * e) := Finset.mem_image_of_mem _ he
    have := h (c * e) hmem
    have harg : ((c * e : ℤ) : ℝ) * x - a = (e : ℝ) * ((c : ℝ) * x) - a := by
      push_cast; ring
    rwa [harg] at this
  · intro h f hf
    rcases Finset.mem_image.mp hf with ⟨e, he, rfl⟩
    have := h e he
    have harg : ((c * e : ℤ) : ℝ) * x - a = (e : ℝ) * ((c : ℝ) * x) - a := by
      push_cast; ring
    rwa [harg]

/-- **`Good θ (c·E) = (·c)⁻¹ (Good θ E)`** : existence of an empty arc for the dilated set
at `x` is existence of one for the original set at `c·x`. -/
theorem good_dilate (θ : ℝ) (E : Finset ℤ) (c : ℤ) :
    Good θ (E.image (fun e => c * e)) = (fun x => (c : ℝ) * x) ⁻¹' (Good θ E) := by
  ext x
  constructor
  · rintro ⟨a, ha⟩
    exact ⟨a, (emptyArc_dilate E c x a).mp ha⟩
  · rintro ⟨a, ha⟩
    exact ⟨a, (emptyArc_dilate E c x a).mpr ha⟩

/-! ### 1-periodicity of the good set -/

/-- The `e`-phase is unchanged by shifting the slow phase by `1`: it adds the integer `e`
inside `Int.fract`, which `Int.fract` erases. -/
theorem emptyArc_add_one {θ : ℝ} (E : Finset ℤ) (x a : ℝ) :
    EmptyArc θ E (x + 1) a ↔ EmptyArc θ E x a := by
  constructor
  · intro h e he
    have := h e he
    have harg : (e : ℝ) * (x + 1) - a = ((e : ℝ) * x - a) + (e : ℤ) := by ring
    rw [harg, Int.fract_add_intCast] at this
    exact this
  · intro h e he
    have := h e he
    have harg : (e : ℝ) * (x + 1) - a = ((e : ℝ) * x - a) + (e : ℤ) := by ring
    rw [harg, Int.fract_add_intCast]
    exact this

/-- **`Good θ E` is 1-periodic.**  Together with `good_dilate` this is the structural core of
`muGood`-dilation invariance: `x ↦ c·x` is a measure-balanced `c`-fold cover of the period. -/
theorem good_add_one (θ : ℝ) (E : Finset ℤ) (x : ℝ) :
    x + 1 ∈ Good θ E ↔ x ∈ Good θ E := by
  constructor
  · rintro ⟨a, ha⟩; exact ⟨a, (emptyArc_add_one E x a).mp ha⟩
  · rintro ⟨a, ha⟩; exact ⟨a, (emptyArc_add_one E x a).mpr ha⟩

-- kernel-purity audit (propext / Classical.choice / Quot.sound only)
#print axioms good_dilate
#print axioms good_add_one

end TailDiameter
end LonelyRunner
