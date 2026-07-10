/-
  TournamentH7.LRCParityBank -- the parity law wired into the 966 bank
  (klein-2026-07-09-S229; the S227/S228 law made a standing bank invariant).

  Three upgrades to kps-S115's LRCCovering966, ALL FREE (zero new
  native_decide -- the law does the work):

   * bandOK_mirror: the band certificate is mirror-symmetric in pure ℤ
     (p ↦ q − p), by the sum-to-zero-mod-q route (cast-free);
   * coveringWitnesses_twin_band: every one of the 966 stored witnesses
     yields its TWIN certificate (l, q − p, q) by the law;
   * coveringWitnesses_lonely_twin: 966 second loneliness proofs through the
     dispatch -- every bank entry now carries TWO independent witnesses;
   * coveringWitnesses_parity: the LM-EVEN validation invariant holds at
     every bank modulus (from coveringPrim's q = 2 covering clause) -- the
     standing bug-detector: any future bank evaluation finding an odd live
     count contradicts a kernel-pure theorem.

  Kernel-pure: no native_decide, no sorry (the bank's own native_decides are
  upstream, untouched).
-/

import Mathlib
import TournamentH7.LRCCovering966
import TournamentH7.LRCParityPairing

namespace LonelyRunner
namespace LRC14Concrete

/-- **The ℤ-band mirror**: a band certificate at `p` yields one at `q − p`.
Sum-to-zero route: the two residues add to `0 (mod q)`, both lie in `[0, q)`,
and the band forces the first positive. -/
theorem bandOK_mirror (l : List ℤ) (p q : ℤ) (h : bandOK (l, p, q)) :
    bandOK (l, q - p, q) := by
  obtain ⟨hq0, hband0⟩ := h
  have hq : (0 : ℤ) < q := hq0
  have hband : ∀ i : Fin 13,
      q ≤ 14 * ((vOf l i * p) % q) ∧ 14 * ((vOf l i * p) % q) ≤ 13 * q := hband0
  refine ⟨hq, fun i => ?_⟩
  show q ≤ 14 * ((vOf l i * (q - p)) % q) ∧ 14 * ((vOf l i * (q - p)) % q) ≤ 13 * q
  obtain ⟨h1, h2⟩ := hband i
  set v : ℤ := vOf l i with hv
  set x : ℤ := (v * p) % q with hx
  set y : ℤ := (v * (q - p)) % q with hy
  have hx0 : 0 ≤ x := Int.emod_nonneg _ (ne_of_gt hq)
  have hxq : x < q := Int.emod_lt_of_pos _ hq
  have hy0 : 0 ≤ y := Int.emod_nonneg _ (ne_of_gt hq)
  have hyq : y < q := Int.emod_lt_of_pos _ hq
  have hxpos : 0 < x := by
    rcases lt_or_eq_of_le hx0 with h | h
    · exact h
    · exfalso; rw [← h] at h1; omega
  have hsum : (x + y) % q = 0 := by
    have htot : v * p + v * (q - p) = v * q := by ring
    calc (x + y) % q
        = (v * p + v * (q - p)) % q := by
          rw [hx, hy, Int.add_emod, Int.emod_emod_of_dvd _ dvd_rfl,
              Int.emod_emod_of_dvd _ dvd_rfl, ← Int.add_emod]
      _ = (v * q) % q := by rw [htot]
      _ = 0 := Int.mul_emod_left _ _
  obtain ⟨k, hk⟩ : q ∣ (x + y) := Int.dvd_of_emod_eq_zero hsum
  have hyeq : y = q - x := by
    have hk01 : k = 0 ∨ k = 1 := by
      rcases lt_trichotomy k 0 with h | h | h
      · exfalso; nlinarith
      · exact Or.inl h
      · have : k = 1 := by nlinarith
        exact Or.inr this
    rcases hk01 with rfl | rfl
    · exfalso; simp at hk; omega
    · omega
  clear_value v x y
  constructor <;> omega

/-- **The twin bank, free**: every stored witness yields its mirror witness
by the law -- 966 second certificates with zero new computation. -/
theorem coveringWitnesses_twin_band :
    ∀ e ∈ coveringWitnesses, bandOK (e.1, e.2.2 - e.2.1, e.2.2) :=
  fun e he => bandOK_mirror e.1 e.2.1 e.2.2 (coveringWitnesses_band e he)

/-- **966 second loneliness proofs**: the twin witnesses drive through the
dispatch just like the originals. -/
theorem coveringWitnesses_lonely_twin :
    ∀ e ∈ coveringWitnesses, (1 : ℝ) / 14 ≤ Mreach (vOf e.1) := by
  intro e he
  obtain ⟨hq, hband⟩ := coveringWitnesses_twin_band e he
  exact mreach_ge_of_pairsum_band (vOf e.1) (e.2.2 - e.2.1) e.2.2 hq hband

/-- A covering list has an even entry at some `Fin 13` index (from the `q = 2`
covering clause and the length). -/
theorem even_speed_of_coveringPrim (l : List ℤ) (h : coveringPrim l) :
    ∃ i : Fin 13, (2 : ℤ) ∣ vOf l i := by
  obtain ⟨hlen, _, _, hcov, _⟩ := h
  obtain ⟨x, hxmem, hx2⟩ := hcov 2 (by simp)
  obtain ⟨n, hn, hnx⟩ := List.mem_iff_getElem.mp hxmem
  refine ⟨⟨n, by omega⟩, ?_⟩
  have : vOf l ⟨n, by omega⟩ = x := by
    unfold vOf
    simp only
    rw [List.getD_eq_getElem l 0 (by omega)]
    exact hnx
  rw [this]
  exact Int.dvd_of_emod_eq_zero hx2

/-- **The standing bank invariant**: the live count is EVEN at every bank
modulus -- any future evaluation finding it odd contradicts this kernel-pure
theorem (the bug-detector, wired). -/
theorem coveringWitnesses_parity :
    ∀ e ∈ coveringWitnesses, ∀ q : ℕ, liveCount (vOf e.1) q % 2 = 0 :=
  fun e he q =>
    liveCount_even_of_even_speed (vOf e.1) q
      (even_speed_of_coveringPrim e.1 (coveringWitnesses_valid e he))

end LRC14Concrete
end LonelyRunner

-- kernel-purity audit (fleet convention)
#print axioms LonelyRunner.LRC14Concrete.bandOK_mirror
#print axioms LonelyRunner.LRC14Concrete.coveringWitnesses_twin_band
#print axioms LonelyRunner.LRC14Concrete.coveringWitnesses_lonely_twin
#print axioms LonelyRunner.LRC14Concrete.coveringWitnesses_parity
