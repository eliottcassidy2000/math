/-
  TournamentH7.LRCDilation — WLOG gcd = 1 by dilation invariance
  (mac-mini-2026-07-03-S24, HYP-4043).

  Loneliness is invariant under scaling the speed vector by a nonzero integer `g`:
  `g • v` is lonely at `t` iff `v` is lonely at `g • t` (both sides are, verbatim,
  `∀ i m, 1/n ≤ |g·vᵢ·t − m|`).  Hence the EXISTENCE of a lonely time is gcd-invariant,
  and every LRC(n) statement may assume `gcd = 1`.

  This closes the completeness gap HYP-4043: `CoveringFarLonely` as stated admits gcd>1
  dilations of the tight AP (e.g. `{2,4,…,26} = 2·{1,…,13}`) which are lonely only at an
  isolated point (measure-zero good region), unreachable by the positive-length far-peel.
  After the reduction below, such a family becomes `{1,…,13}` — a window family handled by
  the census — so the far-peel only ever sees gcd=1 families, all of positive good-region.

  Kernel-pure (`propext, Classical.choice, Quot.sound`); no `sorry`, no `native_decide`.
-/
import Mathlib
import TournamentH7.LonelyRunner

namespace LonelyRunner

/-- **Dilation identity**: `g • v` is lonely at `t` iff `v` is lonely at `g • t`.
Both unfold to `∀ i m, 1/n ≤ |g·vᵢ·t − m|`. -/
theorem lonely_smul {ι : Type*} (n : ℕ) (v : ι → ℤ) (g : ℤ) (t : ℝ) :
    Lonely n (fun i => g * v i) t ↔ Lonely n v ((g : ℝ) * t) := by
  unfold Lonely
  refine forall_congr' fun i => forall_congr' fun m => ?_
  have h : ((g * v i : ℤ) : ℝ) * t = ((v i : ℤ) : ℝ) * ((g : ℝ) * t) := by
    push_cast; ring
  rw [h]

/-- **WLOG gcd = 1**: existence of a lonely time is invariant under dilation by a nonzero
integer `g`.  So `∃ t, Lonely n (g • v) t ↔ ∃ t, Lonely n v t`. -/
theorem exists_lonely_smul {ι : Type*} (n : ℕ) (v : ι → ℤ) (g : ℤ) (hg : g ≠ 0) :
    (∃ t : ℝ, Lonely n (fun i => g * v i) t) ↔ (∃ t : ℝ, Lonely n v t) := by
  have hgR : (g : ℝ) ≠ 0 := by exact_mod_cast hg
  constructor
  · rintro ⟨t, ht⟩
    exact ⟨(g : ℝ) * t, (lonely_smul n v g t).mp ht⟩
  · rintro ⟨s, hs⟩
    refine ⟨s / g, ?_⟩
    rw [lonely_smul]
    have hgs : (g : ℝ) * (s / g) = s := by field_simp
    rwa [hgs]

/-- The reduced-speed vector `v / g` when `g` divides every entry. -/
theorem exists_lonely_of_dvd {ι : Type*} [Fintype ι] (n : ℕ) (v : ι → ℤ) (g : ℤ)
    (hg : g ≠ 0) (hdvd : ∀ i, g ∣ v i) :
    (∃ t : ℝ, Lonely n v t) ↔ (∃ t : ℝ, Lonely n (fun i => v i / g) t) := by
  have hrw : (fun i => g * (v i / g)) = v := by
    funext i; exact Int.mul_ediv_cancel' (hdvd i)
  calc (∃ t : ℝ, Lonely n v t)
      ↔ (∃ t : ℝ, Lonely n (fun i => g * (v i / g)) t) := by rw [hrw]
    _ ↔ (∃ t : ℝ, Lonely n (fun i => v i / g) t) := exists_lonely_smul n _ g hg

-- Note: `LonelyRunner.lonely_scale` (LonelyRunner.lean:286) is the sibling `t/c` form of
-- `lonely_smul`.  The NEW content here is the EXISTENCE / gcd-reduction (`exists_lonely_smul`,
-- `exists_lonely_of_dvd`) — the WLOG-gcd=1 step that closes the HYP-4043 far-peel gap.
#print axioms lonely_smul
#print axioms exists_lonely_smul
#print axioms exists_lonely_of_dvd

end LonelyRunner
