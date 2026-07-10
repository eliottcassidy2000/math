/-
  TournamentH7.LRCDangerCount — the danger-set and gcd-count ingredients of the pair-sum ledger
  (kind-pasteur-2026-07-09-S119).

  The pair-sum liveness ledger (mac-mini's C1, feeding `HasLiveRuler`/`mreach_ge_of_blocked_lt`) bounds
  the number of BLOCKED multipliers by counting the **danger set** `D` — the residues mod `q` that lie
  OUTSIDE the safe band `[q/14, 13q/14]`, i.e. `14r < q` or `13q < 14r` — and lifting that count to each
  runner through the gcd.  This file formalises the two number-theoretic ingredients:

  * `dangerCard_eq` — the danger set has exactly `2⌊(q−1)/14⌋ + 1` residues (the `2m+1` of the ledger,
    `m = ⌈q/14⌉−1 = ⌊(q−1)/14⌋`).
  * `blocked_card_coprime` — when a speed `v` is coprime to `q`, multiplication by `v` permutes the
    residues, so it blocks *exactly* `|D|` multipliers.  This is the `g = 1` case of the ledger's
    `|B_l| = g(2⌊m/g⌋+1)` (the generic / prime-modulus case, mac-mini's C3).

  These bound the per-runner blocked count; the union bound with ± class-merging (mac-mini's THM-672/674)
  assembles them into `blocked < q−1`, which the consumer (kps-S117) turns into loneliness.
  Self-contained (imports Mathlib only).
-/
import Mathlib

namespace LonelyRunner
namespace LRC14Ledger

open Finset

/-- The number of **danger residues** mod `q`: those outside the safe band `[q/14, 13q/14]`. -/
def dangerCard (q : ℕ) : ℕ :=
  ((Finset.range q).filter (fun r => 14 * r < q ∨ 13 * q < 14 * r)).card

/-- **The danger-set size is `2⌊(q−1)/14⌋ + 1`** (the ledger's `2m+1`).  The near-`0` arc
`{r : 14r < q}` has `⌊(q−1)/14⌋ + 1` residues; its reflection `{r : 13q < 14r}` near `q` has
`⌊(q−1)/14⌋`; they are disjoint. -/
theorem dangerCard_eq (q : ℕ) (hq : 0 < q) : dangerCard q = 2 * ((q - 1) / 14) + 1 := by
  unfold dangerCard
  rw [Finset.filter_or, Finset.card_union_of_disjoint]
  · -- the two arcs, as explicit ranges
    have hA : (Finset.range q).filter (fun r => 14 * r < q) = Finset.range ((q - 1) / 14 + 1) := by
      ext r; simp only [Finset.mem_filter, Finset.mem_range]; omega
    have hB : ((Finset.range q).filter (fun r => 13 * q < 14 * r)).card = (q - 1) / 14 := by
      have key : (Finset.range q).filter (fun r => 13 * q < 14 * r)
          = (Finset.range q) \ Finset.range (13 * q / 14 + 1) := by
        ext r; simp only [Finset.mem_filter, Finset.mem_sdiff, Finset.mem_range]; omega
      have hsub : Finset.range (13 * q / 14 + 1) ⊆ Finset.range q := by
        intro r hr; simp only [Finset.mem_range] at hr ⊢; omega
      rw [key, Finset.card_sdiff, Finset.card_range, Finset.inter_eq_left.mpr hsub,
        Finset.card_range]
      omega
    rw [hA, Finset.card_range, hB]; omega
  · -- the two arcs are disjoint
    rw [Finset.disjoint_left]
    intro r hrA hrB
    simp only [Finset.mem_filter] at hrA hrB
    omega

/-- The **blocked multipliers** of speed `v` at modulus `q`: those `p ∈ {0,…,q−1}` whose residue
`(v·p) mod q` is a danger residue (leaves the band). -/
def blockedCard (v q : ℕ) : ℕ :=
  ((Finset.range q).filter (fun p => 14 * (v * p % q) < q ∨ 13 * q < 14 * (v * p % q))).card

/-- **Coprime case of the gcd-count (`g = 1`).**  If `gcd(v, q) = 1`, multiplication by `v` permutes
the residues mod `q`, so `v` blocks exactly `|D|` multipliers: `blockedCard v q = dangerCard q =
2⌊(q−1)/14⌋+1`.  This is the ledger's `|B_l| = g(2⌊m/g⌋+1)` at `g = 1` — the generic / prime case. -/
theorem blocked_card_coprime (v q : ℕ) (hq : 0 < q) (hcop : Nat.Coprime v q) :
    blockedCard v q = dangerCard q := by
  unfold blockedCard dangerCard
  have hinj : Set.InjOn (fun p => v * p % q) (Finset.range q) := by
    intro p hp p' hp' h
    rw [Finset.mem_coe, Finset.mem_range] at hp hp'
    have hme : v * p ≡ v * p' [MOD q] := h
    have hpp : p ≡ p' [MOD q] := Nat.ModEq.cancel_left_of_coprime hcop.symm hme
    rwa [Nat.ModEq, Nat.mod_eq_of_lt hp, Nat.mod_eq_of_lt hp'] at hpp
  have himg : (Finset.range q).image (fun p => v * p % q) = Finset.range q := by
    apply Finset.eq_of_subset_of_card_le
    · intro x hx
      rw [Finset.mem_image] at hx
      obtain ⟨p, _, rfl⟩ := hx
      exact Finset.mem_range.mpr (Nat.mod_lt _ hq)
    · rw [Finset.card_image_of_injOn hinj]
  refine Finset.card_bij (fun p _ => v * p % q) ?_ ?_ ?_
  · -- blocked ↦ danger
    intro p hp
    rw [Finset.mem_filter, Finset.mem_range] at hp
    exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (Nat.mod_lt _ hq), hp.2⟩
  · -- injective on the blocked set
    intro p hp p' hp' h
    rw [Finset.mem_filter, Finset.mem_range] at hp hp'
    exact hinj (Finset.mem_coe.mpr (Finset.mem_range.mpr hp.1))
      (Finset.mem_coe.mpr (Finset.mem_range.mpr hp'.1)) h
  · -- surjective onto the danger set
    intro r hr
    rw [Finset.mem_filter, Finset.mem_range] at hr
    have hrimg : r ∈ (Finset.range q).image (fun p => v * p % q) := by
      rw [himg]; exact Finset.mem_range.mpr hr.1
    rw [Finset.mem_image] at hrimg
    obtain ⟨p, hp, hpr⟩ := hrimg
    refine ⟨p, ?_, hpr⟩
    rw [Finset.mem_filter, Finset.mem_range]
    rw [Finset.mem_range] at hp
    exact ⟨hp, by rw [hpr]; exact hr.2⟩

end LRC14Ledger
end LonelyRunner
