/-
  TournamentH7.LRCReachTransport  (kind-pasteur-2026-07-11-S127 cont.50)

  DILATION PRESERVES LOOSENESS. For a 13-speed family `v` and an integer `c ≥ 1`,

        reach(c·v) ≥ reach(v),      reach X = sSup (margin X '' [0,1]),

  by the scaled witness: if `t₀` attains `reach v`, then `t₀/c ∈ [0,1]` attains the SAME margin for
  `c·v` (because `(c·vᵢ)·(t₀/c) = vᵢ·t₀`). This is the *easy* direction of the dilation-invariance of
  `M` (THM-531) — the one that matters for the endgame: it TRANSPORTS looseness up dilation, so a
  bounded structural *core* that is loose (`reach > 1/14`) makes ALL its dilates loose at every scale.

  This is the formal underpinning of MISTAKE-140 (boxeph-S19): "min M grows with diameter" is a
  sampling artifact because dilation carries any small-M / large-M structure to every scale — the DC
  class stratifies by STRUCTURE, not diameter. `loose_dilate` is the machine-checked transport that
  reduces the unbounded near-dilate slice to the bounded cores (the finite base).
-/
import TournamentH7.LRCDecorrelation13

open TournamentH7.LRCWitness TournamentH7.LRCDecorr13

namespace TournamentH7.LRCReachTransport

/-- **Dilation does not decrease reach.**  For `c ≥ 1`, `reach(v) ≤ reach(c·v)`. -/
theorem reach_dilate_ge (v : Fin 13 → ℤ) (c : ℤ) (hc : 1 ≤ c) :
    sSup (margin v '' Set.Icc (0 : ℝ) 1) ≤ sSup (margin (fun i => c * v i) '' Set.Icc (0 : ℝ) 1) := by
  have hc0 : (0 : ℝ) < (c : ℝ) := by exact_mod_cast (by omega : 0 < c)
  obtain ⟨t0, ht0mem, ht0max⟩ := exists_max_margin v
  have hbddV : BddAbove (margin v '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨s, _, rfl⟩; exact margin_le_half13 v s⟩
  have hsupV : sSup (margin v '' Set.Icc (0 : ℝ) 1) = margin v t0 := by
    apply le_antisymm
    · refine csSup_le ((show (Set.Icc (0 : ℝ) 1).Nonempty from ⟨t0, ht0mem⟩).image (margin v)) ?_
      rintro _ ⟨s, hs, rfl⟩; exact ht0max s hs
    · exact le_csSup hbddV ⟨t0, ht0mem, rfl⟩
  set s : ℝ := t0 / (c : ℝ) with hs
  have ht0_01 : t0 ∈ Set.Icc (0 : ℝ) 1 := ht0mem
  have hs_mem : s ∈ Set.Icc (0 : ℝ) 1 := by
    refine ⟨div_nonneg (by linarith [ht0_01.1]) (le_of_lt hc0), ?_⟩
    rw [hs, div_le_one hc0]
    nlinarith [ht0_01.2, (by exact_mod_cast hc : (1 : ℝ) ≤ (c : ℝ))]
  have hcs : (c : ℝ) * s = t0 := by rw [hs]; field_simp
  -- `margin v t0 ≤ margin (c·v) s` because each dilated residue equals the original at `t0`
  have hkey : margin v t0 ≤ margin (fun i => c * v i) s := by
    rw [le_margin_iff]
    intro i m
    have h1 : margin v t0 ≤ distZ ((v i : ℝ) * t0) := Finset.inf'_le _ (Finset.mem_univ i)
    have h3 : ((c * v i : ℤ) : ℝ) * s - (m : ℝ) = (v i : ℝ) * t0 - (m : ℝ) := by
      push_cast; rw [← hcs]; ring
    show margin v t0 ≤ |((c * v i : ℤ) : ℝ) * s - (m : ℝ)|
    rw [h3]
    exact h1.trans ((le_distZ_iff (distZ ((v i : ℝ) * t0)) _).1 le_rfl m)
  have hbddCV : BddAbove (margin (fun i => c * v i) '' Set.Icc (0 : ℝ) 1) :=
    ⟨1 / 2, by rintro _ ⟨u, _, rfl⟩; exact margin_le_half13 _ u⟩
  have hmem : margin (fun i => c * v i) s ∈ margin (fun i => c * v i) '' Set.Icc (0 : ℝ) 1 :=
    ⟨s, hs_mem, rfl⟩
  calc sSup (margin v '' Set.Icc (0 : ℝ) 1) = margin v t0 := hsupV
    _ ≤ margin (fun i => c * v i) s := hkey
    _ ≤ sSup (margin (fun i => c * v i) '' Set.Icc (0 : ℝ) 1) := le_csSup hbddCV hmem

/-- **Looseness transports up dilation.**  If a core `v` is loose (`reach v > 1/14`) then every
integer dilate `c·v` (`c ≥ 1`) is loose — at every scale, without recomputation.  Formal underpinning
of MISTAKE-140: the near-dilate slice reduces to the bounded loose cores. -/
theorem loose_dilate (v : Fin 13 → ℤ) (c : ℤ) (hc : 1 ≤ c)
    (hv : (1 : ℝ) / 14 < sSup (margin v '' Set.Icc (0 : ℝ) 1)) :
    (1 : ℝ) / 14 < sSup (margin (fun i => c * v i) '' Set.Icc (0 : ℝ) 1) :=
  lt_of_lt_of_le hv (reach_dilate_ge v c hc)

#print axioms reach_dilate_ge
#print axioms loose_dilate

end TournamentH7.LRCReachTransport
