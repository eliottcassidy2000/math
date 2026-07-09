/-
  TournamentH7.LRCComposedRealization — the Lean core of LEM-014, the P-SEPARATED COMPOSED
  REALIZATION (death-star-2026-07-09-S2, HYP-5750; boxeph-2026-07-09-S1's handoff H1).

  CONTEXT.  klein-S205's `LRCDriftEmbed.minReach_ge_of_driftGap` realizes a good period as a real
  lonely time, but it binds EVERY runner to the ruler (`v_i = Vmax − e_i`, `|e_i| ≤ spread`) — the
  "all-13 fold".  mac-mini-S64 showed that fold is vacuous on covering sets (they contain small
  speeds, so the full-set spread ≈ Vmax and the drift budget dies).  boxeph's LEM-014 fixed the
  MATH with the P-SEPARATED composition — cluster runners cleared by drift absorption around the
  observer anchor, SLOW runners cleared by carrying an eroded safety margin from a nearby good
  slow time `x*` (the robust density-floor point) — and verified it exactly, including `P ≠ ∅`
  (k = 10 margin +0.038, k = 8 margin +0.069).  This file is its Lean transcription:

    · `nearInt_ge_of_close`   — the SLOW LEG: `nearInt` moves by at most `|y − x|`
      (the ∀m-characterization makes this a two-line triangle inequality; the ε-eroded
      `G_P^ε` hypothesis of LEM-014 is exactly "safety at `x` with margin `|v|·Δ`").
    · `nearInt_clear_of_driftGap_single` — the CLUSTER LEG: klein-S205's per-runner
      drift-absorbed clearance, factored out of the all-13 fold so a SUBSET of runners
      can consume it (same algebra, same anchor, credit klein-S205).
    · `minReach_ge_of_composed_realization` / `Mreach_ge_of_composed_realization` — the
      COMPOSITION: a per-runner disjunction [cluster-bound at grid j] ∨ [ε-eroded slow
      safety at x], concluding loneliness at the explicit LEM-014 time
      `τ = (j + a + g/2)/Vmax`.

  The measure-side inputs — that a robust good `x` with eroded `G_P` safety EXISTS — are the
  analytic hypotheses; they are supplied by LEM-014's robust floor (shell bound
  `E[W] ≤ 6δ + (6/7)·μ_{1/7+δ}` + `G_P` erosion), kps-S113's `LRCDensityFloorCert`, and the
  θ-transfer thread (boxeph-S2 HYP-5722, monad THM-669/670).  Instantiating `j = round(Vmax·x)`
  and `Δ = 3/(2·Vmax)` (so `|τ − x| ≤ (1/2 + φ)/Vmax ≤ Δ`) recovers LEM-014's exact constants
  (`δ = 3s/Vmax` drift budget, `ε = 20/Vmax` for slow speeds ≤ 13).

  Scope note (the delineation, death-star-S1): with `P = ∅` this file is NOT needed —
  `LRCPureClusterCorner.mreach_ge_of_pure_cluster` (via kps-S28) closes ratio ≤ 13 outright.
  The composition is the ratio > 13 instrument.

  Kernel-pure target: no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LRCCriterionC
import TournamentH7.LRCGapReach
import TournamentH7.LRCDriftEmbed
import TournamentH7.LRCAPTight
import TournamentH7.LRCHembedScaleSep
import TournamentH7.LRCReachWitness

namespace LonelyRunner
namespace LRC14Concrete

open Set

/-- **The slow-leg transport: `nearInt` is 1-Lipschitz, in the one-sided form the composition
consumes.**  If the safety at `x` carries a margin covering the move to `y` — `c + |y − x| ≤
nearInt x` — then the safety survives at `y`: `c ≤ nearInt y`.  Two lines from the
∀m-characterization: for every integer `m`, `|y − m| ≥ |x − m| − |y − x| ≥ nearInt x − |y − x|`.
This is LEM-014's `G_P^ε` erosion made pointwise. -/
theorem nearInt_ge_of_close {c x y : ℝ} (h : c + |y - x| ≤ nearInt x) : c ≤ nearInt y := by
  apply le_nearInt_of_forall_int
  intro m
  have h1 : nearInt x ≤ |x - (m : ℝ)| := nearInt_le_abs_sub x m
  have h2 : |x - (m : ℝ)| ≤ |x - y| + |y - (m : ℝ)| := abs_sub_le x y m
  have h3 : |x - y| = |y - x| := abs_sub_comm x y
  linarith

/-- **The cluster-leg clearance, per runner** (klein-S205's drift-absorbed argument, factored out
of the all-13 fold of `minReach_ge_of_driftGap` so that a SUBSET of runners can consume it).
A ruler-bound runner `vi = Vmax − e` with `|e| ≤ s`, whose tooth `e·j/Vmax` avoids the open gap
`(a, a+g)` in every integer translate, is cleared strictly beyond `1/14` by the gap-midpoint fast
phase `φ = a + g/2`, provided the drift margin `1/7 + 2·(s·φ/Vmax) < g` holds. -/
theorem nearInt_clear_of_driftGap_single
    (Vmax j e : ℤ) (vi a g s : ℝ)
    (hV : 0 < Vmax)
    (hbind : vi = (Vmax : ℝ) - (e : ℝ))
    (hsp : |(e : ℝ)| ≤ s)
    (ha : 0 ≤ a) (hgpos : 0 < g)
    (hfree : ∀ n : ℤ, ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)) ∉ Ioo a (a + g))
    (hgap : 1 / 7 + 2 * (s * (a + g / 2) / (Vmax : ℝ)) < g) :
    (1 : ℝ) / 14 < nearInt (vi * (((j : ℝ) + (a + g / 2)) / (Vmax : ℝ))) := by
  have hVR : (0 : ℝ) < (Vmax : ℝ) := by exact_mod_cast hV
  set φ : ℝ := a + g / 2 with hφ
  have hφ0 : 0 ≤ φ := by rw [hφ]; linarith
  -- the midpoint clears every tooth-translate by g/2
  have hmargin : ∀ n : ℤ, g / 2 ≤ |φ - ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))| := by
    intro n
    by_contra hcon
    push_neg at hcon
    rw [abs_lt] at hcon
    obtain ⟨h1, h2⟩ := hcon
    exact hfree n ⟨by rw [hφ] at h2; linarith, by rw [hφ] at h1; linarith⟩
  -- the drift is φ-proportional and spread-bounded
  have hdrift : |(e : ℝ) * φ / (Vmax : ℝ)| ≤ s * φ / (Vmax : ℝ) := by
    rw [abs_div, abs_mul, abs_of_pos hVR, abs_of_nonneg hφ0]
    gcongr
  -- algebra: vi·τ = (φ − tooth − drift) + j
  have halg : vi * (((j : ℝ) + φ) / (Vmax : ℝ))
      = (φ - (e : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e : ℝ) * φ / (Vmax : ℝ)) + (j : ℝ) := by
    rw [hbind]; field_simp; ring
  rw [halg, nearInt_add_int]
  apply GapReach.nearInt_gt_of_forall_int
  intro n
  have h1 := hmargin n
  have h2 := hdrift
  have h3 : |φ - ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))|
        - |(e : ℝ) * φ / (Vmax : ℝ)|
      ≤ |φ - (e : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e : ℝ) * φ / (Vmax : ℝ) - (n : ℝ)| := by
    have hh := abs_sub_abs_le_abs_sub
      (φ - ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ))) ((e : ℝ) * φ / (Vmax : ℝ))
    rwa [show (φ - ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)))
            - (e : ℝ) * φ / (Vmax : ℝ)
          = φ - (e : ℝ) * (j : ℝ) / (Vmax : ℝ) - (e : ℝ) * φ / (Vmax : ℝ) - (n : ℝ)
        by ring] at hh
  linarith

/-- **The P-separated composed realization (LEM-014's Lean core).**  At the explicit LEM-014 time
`τ = (j + a + g/2)/Vmax`, every runner clears `1/14`, provided EACH runner satisfies one of:

  * **CLUSTER leg** — it is ruler-bound (`v i = Vmax − e`, `|e| ≤ s`) and its tooth avoids the
    open gap `(a, a+g)` at the grid point `j` in every integer translate (drift margin `hgap`
    shared across the cluster); or
  * **SLOW leg** — it carries the ε-eroded safety at a nearby slow time `x`:
    `1/14 + |v i|·Δ ≤ nearInt(v i · x)` with `|τ − x| ≤ Δ`  (LEM-014's `G_P^ε` event).

The measure-side existence of such `(x, j, a, g)` is LEM-014's robust floor; this theorem is the
pointwise composition it feeds. -/
theorem minReach_ge_of_composed_realization
    (v : Fin 13 → ℤ) (Vmax j : ℤ) (a g s Δ x : ℝ)
    (hV : 0 < Vmax)
    (ha : 0 ≤ a) (hag : a + g ≤ 1) (hgpos : 0 < g)
    (hgap : 1 / 7 + 2 * (s * (a + g / 2) / (Vmax : ℝ)) < g)
    (hclose : |((j : ℝ) + (a + g / 2)) / (Vmax : ℝ) - x| ≤ Δ)
    (hsplit : ∀ i : Fin 13,
      (∃ e : ℤ, (v i : ℝ) = (Vmax : ℝ) - (e : ℝ) ∧ |(e : ℝ)| ≤ s ∧
        ∀ n : ℤ, ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)) ∉ Ioo a (a + g))
      ∨ (1 / 14 + |(v i : ℝ)| * Δ ≤ nearInt ((v i : ℝ) * x))) :
    (1 : ℝ) / 14 ≤ minReach v (((j : ℝ) + (a + g / 2)) / (Vmax : ℝ)) := by
  rw [le_minReach_iff]
  intro i
  rcases hsplit i with ⟨e, hbind, hsp, hfree⟩ | hslow
  · -- cluster leg: strict clearance from the factored drift argument
    exact le_of_lt
      (nearInt_clear_of_driftGap_single Vmax j e ((v i : ℝ)) a g s hV hbind hsp ha
        hgpos hfree hgap)
  · -- slow leg: transport the eroded safety from x to τ
    have hmove : |(v i : ℝ) * (((j : ℝ) + (a + g / 2)) / (Vmax : ℝ)) - (v i : ℝ) * x|
        ≤ |(v i : ℝ)| * Δ := by
      rw [← mul_sub, abs_mul]
      exact mul_le_mul_of_nonneg_left hclose (abs_nonneg _)
    exact nearInt_ge_of_close (x := (v i : ℝ) * x) (by linarith)

/-- **`Mreach ≥ 1/14` from the composed realization** — the `hpartA`-shaped conclusion for the
genuinely multi-scale (speeds-ratio > 13) shapes, the residual class of the delineation
(death-star-S1).  Composes with `Mreach_ge_of_minReach` (klein-S205). -/
theorem Mreach_ge_of_composed_realization
    (v : Fin 13 → ℤ) (Vmax j : ℤ) (a g s Δ x : ℝ)
    (hV : 0 < Vmax)
    (ha : 0 ≤ a) (hag : a + g ≤ 1) (hgpos : 0 < g)
    (hgap : 1 / 7 + 2 * (s * (a + g / 2) / (Vmax : ℝ)) < g)
    (hclose : |((j : ℝ) + (a + g / 2)) / (Vmax : ℝ) - x| ≤ Δ)
    (hsplit : ∀ i : Fin 13,
      (∃ e : ℤ, (v i : ℝ) = (Vmax : ℝ) - (e : ℝ) ∧ |(e : ℝ)| ≤ s ∧
        ∀ n : ℤ, ((e : ℝ) * (j : ℝ) / (Vmax : ℝ) + (n : ℝ)) ∉ Ioo a (a + g))
      ∨ (1 / 14 + |(v i : ℝ)| * Δ ≤ nearInt ((v i : ℝ) * x))) :
    (1 : ℝ) / 14 ≤ Mreach v :=
  Mreach_ge_of_minReach v _
    (minReach_ge_of_composed_realization v Vmax j a g s Δ x hV ha hag hgpos hgap hclose hsplit)

/-! ## Axiom audit -/
#print axioms nearInt_ge_of_close
#print axioms nearInt_clear_of_driftGap_single
#print axioms minReach_ge_of_composed_realization
#print axioms Mreach_ge_of_composed_realization

end LRC14Concrete
end LonelyRunner
