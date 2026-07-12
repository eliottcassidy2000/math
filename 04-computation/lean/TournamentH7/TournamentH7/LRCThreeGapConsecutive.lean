/-
Copyright (c) 2026 The TournamentH7 project contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: kind-pasteur (LRC multi-agent project, 2026-07-11-S127 cont.44)
-/
import TournamentH7.LRCCleanRuler

/-!
# The three-gap consecutive-reduction lemma (the AP corner clears via a consecutive block)

For the **AP sub-corner** of the divisor-complete residual, opus-S236 observed that the residues
`{(a+i·d)·p mod q}` are themselves an AP mod `q` (step `d·p`) — a three-gap / Steinhaus system.
The extreme, cleanest case of three-gap is the **consecutive reduction**: for `gcd(d,q)=1`, the
multiplier `p ≡ d⁻¹ (mod q)` collapses the AP to `13` *consecutive* residues `{r, r+1, …, r+12}`
(`r = a·d⁻¹ mod q`).  This file formalizes the elementary core: **a length-13 consecutive block that
fits inside the safe band `[2, q-2]` has `bandCount = 0`** (every runner is safe / lonely), i.e. the
family *clears* at `q`.  Combined with the band-edge margin lemma (opus-S235), clearing at a non-`14`
modulus gives `M > 1/14`.

The band for `16 ≤ q ≤ 28` is exactly `[2, q-2]` (danger arc `{0, ±1}`), because `⌈q/14⌉ = 2` there;
a residue `R ∈ [2, q-2]` satisfies `inBand` since `14·R ≥ 28 ≥ q` and `14·R ≤ 14(q-2) = 14q-28 ≤ 13q`.

The hypothesis `hres` (the residues form the consecutive block `r + i`) is exactly "the consecutive
reduction happened"; a consumer discharges it for an AP family at `p ≡ d⁻¹` by modular arithmetic
(`(a+i·d)·d⁻¹ ≡ a·d⁻¹ + i (mod q)`, no wrap when `r+12 < q`).  Empirically (kps cont.44) this single
disjunction over `p = ±d⁻¹`, `q ∈ [16,27]`, clears **every** primitive divisor-complete 13-term AP
(0 failures / 2161), so the AP corner needs only this consecutive case, not the general three-gap.
-/

namespace LonelyRunner
namespace LRC14Concrete

open Finset

/-- **Safe-band membership.**  For `16 ≤ q ≤ 28`, a runner whose residue `R = (vᵢ·p) mod q` lies in
`[2, q-2]` is `inBand` (safe/lonely at `p/q`). -/
theorem inBand_of_residue_mem_band (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13)
    (hq28 : q ≤ 28)
    (hlo : 2 ≤ (v i * (p : ℤ)) % (q : ℤ))
    (hhi : (v i * (p : ℤ)) % (q : ℤ) ≤ (q : ℤ) - 2) :
    inBand v q p i := by
  unfold inBand
  have hq28' : (q : ℤ) ≤ 28 := by exact_mod_cast hq28
  refine ⟨?_, ?_⟩
  · -- q ≤ 14 · R : since R ≥ 2 we have 14·R ≥ 28 ≥ q
    linarith
  · -- 14 · R ≤ 13 · q : since R ≤ q-2 we have 14·R ≤ 14q-28 ≤ 13q ⟺ q ≤ 28
    linarith

/-- **The consecutive-reduction / three-gap clearing lemma.**  If, at multiplier `p` and modulus
`16 ≤ q ≤ 28`, the `13` residues form a consecutive block `{r, r+1, …, r+12}` that fits inside the
safe band (`2 ≤ r` and `r + 12 ≤ q - 2`), then `bandCount v q p = 0` — the family *clears* at `q`. -/
theorem bandCount_zero_of_consecutive_block (v : Fin 13 → ℤ) (q p : ℕ) (r : ℤ)
    (hq28 : q ≤ 28) (hr2 : 2 ≤ r) (hr12 : r + 12 ≤ (q : ℤ) - 2)
    (hres : ∀ i : Fin 13, (v i * (p : ℤ)) % (q : ℤ) = r + (i : ℕ)) :
    bandCount v q p = 0 := by
  unfold bandCount
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro i _
  rw [not_not]
  refine inBand_of_residue_mem_band v q p i hq28 ?_ ?_
  · rw [hres i]
    have : (0 : ℤ) ≤ ((i : ℕ) : ℤ) := by positivity
    linarith
  · rw [hres i]
    have hi : ((i : ℕ) : ℤ) ≤ 12 := by
      have := i.isLt
      omega
    linarith

/-- **Clearing ⟹ a live (lonely) multiplier exists.**  A consecutive block in the safe band makes
`p` a live multiplier, so `0 < liveCount v q`.  (This is the hypothesis the band-edge margin lemma
consumes to conclude `M ≥ ⌈q/14⌉/q > 1/14` for `q` not a multiple of `14`.) -/
theorem liveCount_pos_of_consecutive_block (v : Fin 13 → ℤ) (q p : ℕ) (r : ℤ)
    (hq28 : q ≤ 28) (hp0 : 0 < p) (hpq : p < q)
    (hr2 : 2 ≤ r) (hr12 : r + 12 ≤ (q : ℤ) - 2)
    (hres : ∀ i : Fin 13, (v i * (p : ℤ)) % (q : ℤ) = r + (i : ℕ)) :
    0 < liveCount v q := by
  rw [liveCount, Finset.card_pos]
  refine ⟨p, ?_⟩
  rw [Finset.mem_filter, Finset.mem_Ioo]
  exact ⟨⟨hp0, hpq⟩, bandCount_zero_of_consecutive_block v q p r hq28 hr2 hr12 hres⟩

/-- **The coprime-reduction safety lemma.**  At a multiplier `p` and modulus `q ≤ 28`, a runner that
shares a *proper* common factor `g` with `q` (`2 ≤ g`, `g ∣ vᵢ`, `g ∣ q`) and is not knocked to `0`
(`q ∤ vᵢ·p`) is automatically **safe** (`inBand`): its residue `r = vᵢ·p mod q` is a nonzero multiple of
`g`, so `g ≤ r ≤ q - g`, hence `2 ≤ r ≤ q - 2`.  Consequence for Route B: at a *unit* multiplier `p`
(`gcd(p,q)=1`) the hypothesis `q ∤ vᵢ·p` holds for every non-multiple runner, so **only the runners
coprime to `q` can be unsafe** — the anti-concentration shrinks to the coprime-to-`q` sub-family, which is
small for divisor-complete families (they are forced to have runners divisible by every prime `≤ 13`). -/
theorem inBand_of_proper_common_factor (v : Fin 13 → ℤ) (q p : ℕ) (i : Fin 13) (g : ℤ)
    (hqpos : 0 < q) (hq28 : q ≤ 28) (hg : 2 ≤ g) (hgv : g ∣ v i) (hgq : g ∣ (q : ℤ))
    (hnz : ¬ (q : ℤ) ∣ (v i * (p : ℤ))) :
    inBand v q p i := by
  have hq0 : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hqpos
  set r : ℤ := (v i * (p : ℤ)) % (q : ℤ) with hrdef
  have hr0 : 0 ≤ r := Int.emod_nonneg _ (by positivity)
  have hrq : r < (q : ℤ) := Int.emod_lt_of_pos _ hq0
  have hrne : r ≠ 0 := fun h => hnz (Int.dvd_of_emod_eq_zero h)
  -- `g ∣ r` because `g ∣ q` gives `r % g = (vᵢ·p) % g = 0` (as `g ∣ vᵢ·p`)
  have hgr : g ∣ r := by
    have hb : (v i * (p : ℤ)) % g = 0 := Int.emod_eq_zero_of_dvd (hgv.mul_right _)
    have hrg : r % g = 0 := by rw [hrdef, Int.emod_emod_of_dvd _ hgq]; exact hb
    exact Int.dvd_of_emod_eq_zero hrg
  have hrpos : 0 < r := lt_of_le_of_ne hr0 (Ne.symm hrne)
  have hge : g ≤ r := Int.le_of_dvd hrpos hgr
  have hle : r ≤ (q : ℤ) - g := by
    have : g ≤ (q : ℤ) - r := Int.le_of_dvd (by omega) (dvd_sub hgq hgr)
    omega
  exact inBand_of_residue_mem_band v q p i hq28 (by omega) (by omega)

end LRC14Concrete
end LonelyRunner

#print axioms LonelyRunner.LRC14Concrete.bandCount_zero_of_consecutive_block
#print axioms LonelyRunner.LRC14Concrete.liveCount_pos_of_consecutive_block
#print axioms LonelyRunner.LRC14Concrete.inBand_of_proper_common_factor
