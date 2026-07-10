/-
  TournamentH7.LRCDepth4Dispatch — the depth-4 dyadic dispatch (mac-mini LEM-021, S65 cont.13).

  The first live layer of the 2-adic parity tower (LEM-020): at `q = 16` every EVEN residue
  clears the band `[2, 14]` automatically (the "even free pass": `2u·m, 4u·m, 8u·m mod 16` land
  in `{2,6,10,14} ∪ {4,12} ∪ {8}` for odd `u, m` — and `14·2 = 28 > 16`; depth 4 is the UNIQUE
  such layer since `32 > 28`).  Only odd-unit residues can fail, and they kill odd multipliers
  in ±inverse pairs.  Hence (LEM-021): if no speed is `≡ 0 (mod 16)` and the odd speeds MISS one
  of the four unit ±classes `{±1, ±3, ±5, ±7} mod 16`, the missed class's inverse multiplier
  `m ∈ {1, 11, 13, 7}` puts every residue `(vᵢ·m) mod 16` in `[2, 14]`, so `t = m/16` is a
  lonely instant with clearance `2/16 = 1/8 > 1/14`.

  PRODUCER for kps-S114's consumer `mreach_ge_of_pairsum_band` (q = 16, p = m): this file emits
  the (q, p) pair from pure residue hypotheses — the depth-4 member of the dispatch family
  (detuned 1/13, common-residue 8/17, all-odd 1/2, depth-4 1/8).  Decides 18.8% of covering
  sets unconditionally (census in lrc14_depth4_dispatch_macmini_S65cont12).
-/
import Mathlib
import TournamentH7.LRCPairSumDispatch

namespace LonelyRunner
namespace LRC14Concrete

/-- **The depth-4 dyadic dispatch (LEM-021, sufficient direction).**  If no speed is
`≡ 0 (mod 16)` and every speed avoids the residues `c` and `16 − c (mod 16)` for one of the
four unit classes `c ∈ {1,3,5,7}` (with its explicit inverse multiplier `m = c⁻¹ mod 16`),
then `Mreach v ≥ 1/14`, witnessed at the dyadic instant `t = m/16`. -/
theorem mreach_ge_of_depth4 (v : Fin 13 → ℤ) (c m : ℤ)
    (hcm : (c = 1 ∧ m = 1) ∨ (c = 3 ∧ m = 11) ∨ (c = 5 ∧ m = 13) ∨ (c = 7 ∧ m = 7))
    (h16 : ∀ i, v i % 16 ≠ 0)
    (hmiss : ∀ i, v i % 16 ≠ c ∧ v i % 16 ≠ 16 - c) :
    (1 : ℝ) / 14 ≤ Mreach v := by
  apply mreach_ge_of_pairsum_band v m 16 (by norm_num)
  intro i
  have h1 := h16 i
  have h2 := (hmiss i).1
  have h3 := (hmiss i).2
  have hr0 : 0 ≤ v i % 16 := Int.emod_nonneg _ (by norm_num)
  have hr16 : v i % 16 < 16 := Int.emod_lt_of_pos _ (by norm_num)
  rw [Int.mul_emod]
  rcases hcm with ⟨hc, hm⟩ | ⟨hc, hm⟩ | ⟨hc, hm⟩ | ⟨hc, hm⟩ <;> subst hc <;> subst hm <;>
    (generalize hgen : v i % 16 = r at h1 h2 h3 hr0 hr16 ⊢
     interval_cases r <;> omega)

end LRC14Concrete
end LonelyRunner
