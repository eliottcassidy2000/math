/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: opus (LRC multi-agent project, 2026-07-02-S40)
-/
import TournamentH7.LRCWitnessCert

/-!
# Module 6, per-level budgets: `cert_ladder'` (THM-606 in its sharp form)

kps's `cert_ladder` uses one uniform drift budget `μ`, which forces consecutive reference
ratios `> 1 + 1/μ` and rejects tuples that THM-606 certifies (e.g. the S36 depth-3 tuple
`(50, 2000, 90000)`, ratios 40 and 45).  This file generalizes to **per-level budgets**
`L.μ` with the last level at `μ = 0` (THM-606: inflation does not accumulate; each level
pays only its own drift).  The induction is kps's `ladder_walk` verbatim with two changes:
the window-width hypothesis relaxes to `0 ≤ δ` (the nil case tolerates a degenerate window,
so `μ_d = 0` is allowed), and the head level absorbs drift `≤ L.μ` into its own band
`h + L.μ`.  The cons case re-derives `0 < δ` from its own separation.

Regression instance: `depth3_perLevel_row` certifies the S36 tuple with the ORIGINAL S36
certificates — the configuration the uniform-μ `SepChain` rejects.
-/

namespace LonelyRunner
namespace WitnessCert

/-- One ladder level with its own drift budget. -/
structure Level' where
  offs : List ℤ
  c    : ℚ
  V    : ℤ
  μ    : ℚ

/-- Chained separations with per-level budgets: the window passed down by level `L` has
width `L.μ / L.V`. -/
def SepChain' : ℚ → List Level' → Prop
  | _, [] => True
  | δ, L :: rest => 0 < L.V ∧ 0 ≤ L.μ ∧ 1 < (L.V : ℚ) * (δ - L.μ / L.V) ∧
      SepChain' (L.μ / L.V) rest

instance instDecidableSepChain' : ∀ (δ : ℚ) (Ls : List Level'), Decidable (SepChain' δ Ls)
  | _, [] => by unfold SepChain'; infer_instance
  | δ, L :: rest => by
      have := instDecidableSepChain' (L.μ / L.V) rest
      unfold SepChain'
      infer_instance

/-- The per-level ladder walk: from any window `[a, a+δ]` (possibly degenerate) inside
`[lo, hi]`, a common time making every level `h`-safe. -/
theorem ladder_walk' {h : ℚ} (h0 : 0 ≤ h) {lo hi : ℚ} :
    ∀ (Ls : List Level') (a δ : ℚ), 0 ≤ δ →
    lo ≤ a → a + δ ≤ hi →
    SepChain' δ Ls →
    (∀ L ∈ Ls, ∀ o ∈ L.offs, 0 ≤ o ∧ arcSafe (h + L.μ) o L.c lo hi) →
    ∃ τ : ℝ, ((a : ℚ) : ℝ) ≤ τ ∧ τ ≤ ((a : ℚ) : ℝ) + ((δ : ℚ) : ℝ) ∧
      ∀ L ∈ Ls, ∀ o ∈ L.offs,
        (h : ℝ) ≤ ‖((((L.V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  intro Ls
  induction Ls with
  | nil =>
      intro a δ hδ _ _ _ _
      refine ⟨((a : ℚ) : ℝ), le_refl _, ?_, ?_⟩
      · have : (0:ℝ) ≤ ((δ : ℚ) : ℝ) := by exact_mod_cast hδ
        linarith
      · intro L hL; cases hL
  | cons L rest ih =>
      intro a δ hδ halo hahi hsep hsafe
      obtain ⟨hVpos, hμnonneg, hVlen, hsep'⟩ := hsep
      have hVR : (0 : ℝ) < (L.V : ℝ) := by exact_mod_cast hVpos
      have hVQpos : (0 : ℚ) < (L.V : ℚ) := by exact_mod_cast hVpos
      set δ' : ℚ := L.μ / L.V with hδ'def
      have hδ'nonneg : 0 ≤ δ' := div_nonneg hμnonneg hVQpos.le
      obtain ⟨j, hj1, hj2⟩ := exists_int_in_long_interval
        (a := (L.V : ℝ) * ((a : ℚ) : ℝ) - ((L.c : ℚ) : ℝ))
        (b := (L.V : ℝ) * (((a : ℚ) : ℝ) + ((δ : ℚ) : ℝ) - ((δ' : ℚ) : ℝ)) - ((L.c : ℚ) : ℝ))
        (by
          have h1 : (1 : ℚ) < (L.V : ℚ) * (δ - δ') := by rw [hδ'def]; exact hVlen
          have h1R : (1 : ℝ) < (L.V : ℝ) * (((δ : ℚ) : ℝ) - ((δ' : ℚ) : ℝ)) := by
            exact_mod_cast h1
          nlinarith)
      set t0Q : ℚ := ((j : ℚ) + L.c) / (L.V : ℚ) with ht0Q
      have ht0cast : ((t0Q : ℚ) : ℝ) = ((j : ℝ) + ((L.c : ℚ) : ℝ)) / (L.V : ℝ) := by
        rw [ht0Q]; push_cast; ring
      have hVt0 : (L.V : ℝ) * ((t0Q : ℚ) : ℝ) = (j : ℝ) + ((L.c : ℚ) : ℝ) := by
        rw [ht0cast]; field_simp
      have ht0a : ((a : ℚ) : ℝ) ≤ ((t0Q : ℚ) : ℝ) := by
        rw [ht0cast, le_div_iff₀ hVR]; linarith
      have ht0b : ((t0Q : ℚ) : ℝ) ≤ ((a : ℚ) : ℝ) + ((δ : ℚ) : ℝ) - ((δ' : ℚ) : ℝ) := by
        rw [ht0cast, div_le_iff₀ hVR]; linarith
      have ht0aQ : lo ≤ t0Q := by
        have h1 : ((lo : ℚ) : ℝ) ≤ ((a : ℚ) : ℝ) := by exact_mod_cast halo
        have h2 : ((lo : ℚ) : ℝ) ≤ ((t0Q : ℚ) : ℝ) := le_trans h1 ht0a
        exact_mod_cast h2
      have ht0bQ : t0Q + δ' ≤ hi := by
        have h1 : ((a : ℚ) : ℝ) + ((δ : ℚ) : ℝ) ≤ ((hi : ℚ) : ℝ) := by exact_mod_cast hahi
        have h2 : ((t0Q + δ' : ℚ) : ℝ) ≤ ((hi : ℚ) : ℝ) := by
          push_cast
          push_cast at ht0b
          linarith
        exact_mod_cast h2
      obtain ⟨τ, hτ1, hτ2, hτrest⟩ := ih t0Q δ' hδ'nonneg ht0aQ ht0bQ hsep'
        (fun L' hL' => hsafe L' (List.mem_cons_of_mem _ hL'))
      refine ⟨τ, ?_, ?_, ?_⟩
      · exact le_trans ht0a hτ1
      · have h3 : ((δ' : ℚ) : ℝ) + ((t0Q : ℚ) : ℝ) ≤ ((a : ℚ) : ℝ) + ((δ : ℚ) : ℝ) := by
          push_cast at ht0b ⊢
          linarith
        linarith [hτ2]
      · intro L' hL' o ho
        rcases List.mem_cons.mp hL' with rfl | hL'
        · obtain ⟨ho0, hosafe⟩ := hsafe L' (List.mem_cons_self ..) o ho
          have hτlo' : ((lo : ℚ) : ℝ) ≤ τ := by
            have h1 : ((lo : ℚ) : ℝ) ≤ ((t0Q : ℚ) : ℝ) := by exact_mod_cast ht0aQ
            linarith [hτ1]
          have hτhi' : τ ≤ ((hi : ℚ) : ℝ) := by
            have h1 : ((t0Q + δ' : ℚ) : ℝ) ≤ ((hi : ℚ) : ℝ) := by exact_mod_cast ht0bQ
            push_cast at h1
            linarith [hτ2]
          have hsafeτ := norm_ge_of_arcSafe ho0
            (by linarith : (0:ℚ) ≤ h + L'.μ) hosafe hτlo' hτhi'
          have hcast : (((h + L'.μ : ℚ)) : ℝ) = (h : ℝ) + ((L'.μ : ℚ) : ℝ) := by
            push_cast; ring
          rw [hcast] at hsafeτ
          have hkey : ((L'.V - o : ℤ) : ℝ) * τ
              = ((((L'.c : ℚ) : ℝ) - (o : ℝ) * τ) + (L'.V : ℝ) * (τ - ((t0Q : ℚ) : ℝ))) + (j : ℝ) := by
            push_cast
            push_cast at hVt0
            nlinarith [hVt0]
          rw [hkey, coe_add_int]
          have hdrift : |((((L'.c : ℚ) : ℝ) - (o : ℝ) * τ) + (L'.V : ℝ) * (τ - ((t0Q : ℚ) : ℝ)))
              - (((L'.c : ℚ) : ℝ) - (o : ℝ) * τ)| ≤ ((L'.μ : ℚ) : ℝ) := by
            have h1 : ((((L'.c : ℚ) : ℝ) - (o : ℝ) * τ) + (L'.V : ℝ) * (τ - ((t0Q : ℚ) : ℝ)))
                - (((L'.c : ℚ) : ℝ) - (o : ℝ) * τ) = (L'.V : ℝ) * (τ - ((t0Q : ℚ) : ℝ)) := by ring
            rw [h1]
            have h2 : (0:ℝ) ≤ τ - ((t0Q : ℚ) : ℝ) := by linarith [hτ1]
            have h3 : τ - ((t0Q : ℚ) : ℝ) ≤ ((δ' : ℚ) : ℝ) := by linarith [hτ2]
            rw [abs_of_nonneg (by nlinarith)]
            have hVδ' : (L'.V : ℝ) * ((δ' : ℚ) : ℝ) = ((L'.μ : ℚ) : ℝ) := by
              rw [hδ'def]; push_cast; field_simp
            nlinarith
          have hflip : ‖(((o : ℝ) * τ - ((L'.c : ℚ) : ℝ) : ℝ) : UnitAddCircle)‖
              = ‖((((L'.c : ℚ) : ℝ) - (o : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
            have hne : (((L'.c : ℚ) : ℝ) - (o : ℝ) * τ : ℝ)
                = -(((o : ℝ) * τ - ((L'.c : ℚ) : ℝ))) := by ring
            rw [hne]
            have hcoe : ((-(((o : ℝ) * τ - ((L'.c : ℚ) : ℝ))) : ℝ) : UnitAddCircle)
                = -((((o : ℝ) * τ - ((L'.c : ℚ) : ℝ) : ℝ) : UnitAddCircle)) := rfl
            rw [hcoe, norm_neg]
          have hlip := norm_ge_norm_sub_abs (((L'.c : ℚ) : ℝ) - (o : ℝ) * τ)
            (((((L'.c : ℚ) : ℝ) - (o : ℝ) * τ) + (L'.V : ℝ) * (τ - ((t0Q : ℚ) : ℝ))))
          rw [hflip] at hsafeτ
          linarith
        · exact hτrest L' hL' o ho

/-- **Module 6, sharp form**: per-level drift budgets (THM-606).  Accepts every THM-607
sharp-region tuple; the uniform-`μ` `cert_ladder` is the constant-budget special case. -/
theorem cert_ladder' {h : ℚ} (h0 : 0 ≤ h) {lo hi : ℚ}
    (P : List ℤ) (Ls : List Level')
    (hPpos : ∀ s ∈ P, 0 ≤ s) (hPsafe : ∀ s ∈ P, arcSafe h s 0 lo hi)
    (hlohi : lo ≤ hi)
    (hsep : SepChain' (hi - lo) Ls)
    (hLsafe : ∀ L ∈ Ls, ∀ o ∈ L.offs, 0 ≤ o ∧ arcSafe (h + L.μ) o L.c lo hi) :
    ∃ τ : ℝ,
      (∀ s ∈ P, (h : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      ∀ L ∈ Ls, ∀ o ∈ L.offs,
        (h : ℝ) ≤ ‖((((L.V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ := by
  obtain ⟨τ, hτ1, hτ2, hτL⟩ := ladder_walk' h0 Ls lo (hi - lo)
    (by linarith) (le_refl _) (by linarith) hsep hLsafe
  refine ⟨τ, ?_, hτL⟩
  intro s hs
  have hτlo : ((lo : ℚ) : ℝ) ≤ τ := hτ1
  have hτhi : τ ≤ ((hi : ℚ) : ℝ) := by
    have : ((lo : ℚ) : ℝ) + ((hi - lo : ℚ) : ℝ) = ((hi : ℚ) : ℝ) := by push_cast; ring
    linarith [hτ2, this.symm.le, this.le]
  have := norm_ge_of_arcSafe (hPpos s hs) h0 (hPsafe s hs) hτlo hτhi
  simpa using this

/-- **The S36 regression row**: the tuple `(50, 2000, 90000)` with the ORIGINAL S36
certificates and per-level budgets `(253/9000, 11/450, 0)` — certified here, rejected by
the uniform-`μ` `SepChain` (ratios 40 and 45 < 1 + 1/μ ≈ 41). -/
def packLevels' : List Level' :=
  [⟨[0, 1, 2], 399/4000, 50, 253/9000⟩,
   ⟨[0, 1, 3], 221/1000, 2000, 11/450⟩,
   ⟨[0, 1, 2, 4, 7], 143/2000, 90000, 0⟩]

theorem depth3_perLevel_row :
    ∃ τ : ℝ,
      (∀ s ∈ ([1, 2] : List ℤ), ((1/14 : ℚ) : ℝ) ≤ ‖(((s : ℝ) * τ : ℝ) : UnitAddCircle)‖) ∧
      ∀ L ∈ packLevels', ∀ o ∈ L.offs,
        ((1/14 : ℚ) : ℝ) ≤ ‖((((L.V - o : ℤ) : ℝ) * τ : ℝ) : UnitAddCircle)‖ :=
  cert_ladder' (h := 1/14) (by norm_num) (lo := 7/20) (hi := 3/8) [1, 2] packLevels'
    (by decide) (by native_decide) (by norm_num) (by native_decide) (by native_decide)

end WitnessCert
end LonelyRunner
