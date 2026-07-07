/-
  TournamentH7.LRCTenSwapLadder — the k/(10k+7) one-swap ladder, closed by formula.
  (kind-pasteur-2026-07-07-S54, HYP-4717.)

  The single-scale census (kps-S54) found the FIRST EXCITED rung above the AP `{1..13}`:
  the family `{1,…,9, 11,12,13, 10k} = AP[10 → 10k]` has `M = k/(10k+7)` (first rung `2/27`
  at `k = 2`; `→ 1/10`), an infinite family of LRC(14) lonely certificates, sibling to the
  residue-liar `{1,…,11,13,12k}` (`k/(12k+5)`, `LRCResidueLiar`).  The two families flank the
  Goddyn-Wong tight family `{1,…,11,13,24}` (the `k=2` boundary of the residue-liar, `M=1/14`).

  Lonely at `t = (3k+2)/(10k+7)`: at that time each runner `v` sits at `r_v/(10k+7)` with
  `k ≤ r_v ≤ (10k+7) − k` (residue table linear in `k`; binding LOWER at `v=7` (`r=k`) and
  UPPER at `v=10k` (`r=9k+7`)), so the lattice distance is `≥ k`, and `14k ≥ 10k+7 ⟺ k ≥ 2`.
  Kernel-pure (propext, Classical.choice, Quot.sound); no sorry, no native_decide.
-/
import Mathlib
import TournamentH7.LonelyRunner
import TournamentH7.LRCSpread13
import TournamentH7.LRCResidueLiar

namespace LonelyRunner
namespace LRC14

/-- The one-swap ladder family `{1,…,9, 11,12,13, 10k}` (13 speeds, AP with `10 → 10k`). -/
def tenSwap (k : ℤ) : Fin 13 → ℤ :=
  ![1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 10 * k]

/-- **THE k/(10k+7) LADDER IS LONELY** at `t = (3k+2)/(10k+7)` for every `k ≥ 2`
(margin `k/(10k+7) ≥ 2/27 > 1/14`): an infinite family of LRC(14) lonely certificates in
closed form, kernel-pure, no `native_decide`. -/
theorem tenSwap_lonely (k : ℤ) (hk : 2 ≤ k) :
    Lonely 14 (tenSwap k) (((3 * k + 2 : ℤ) : ℝ) / ((10 * k + 7 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (tenSwap k) (3 * k + 2) (10 * k + 7) (by omega)
  intro i m
  have key : ∀ val qq r : ℤ, val * (3 * k + 2) = qq * (10 * k + 7) + r →
      k ≤ r → r ≤ (10 * k + 7) - k →
      (10 * k + 7) ≤ 14 * |val * (3 * k + 2) - m * (10 * k + 7)| := by
    intro val qq r hEq hlo hhi
    have hm : k ≤ |val * (3 * k + 2) - m * (10 * k + 7)| :=
      lattice_dist_ge qq r (by omega) hEq hlo hhi m
    nlinarith [hm, hk]
  fin_cases i <;> simp only [tenSwap]
  · exact key 1 0 (3 * k + 2) (by ring) (by omega) (by omega)
  · exact key 2 0 (6 * k + 4) (by ring) (by omega) (by omega)
  · exact key 3 0 (9 * k + 6) (by ring) (by omega) (by omega)
  · exact key 4 1 (2 * k + 1) (by ring) (by omega) (by omega)
  · exact key 5 1 (5 * k + 3) (by ring) (by omega) (by omega)
  · exact key 6 1 (8 * k + 5) (by ring) (by omega) (by omega)
  · exact key 7 2 k (by ring) (by omega) (by omega)
  · exact key 8 2 (4 * k + 2) (by ring) (by omega) (by omega)
  · exact key 9 2 (7 * k + 4) (by ring) (by omega) (by omega)
  · exact key 11 3 (3 * k + 1) (by ring) (by omega) (by omega)
  · exact key 12 3 (6 * k + 3) (by ring) (by omega) (by omega)
  · exact key 13 3 (9 * k + 5) (by ring) (by omega) (by omega)
  · exact key (10 * k) (3 * k - 1) (9 * k + 7) (by ring) (by omega) (by omega)

/-- The first excited rung: `{1,…,9,11,12,13,20}` (`k=2`) is lonely at `8/27`, margin
`2/27 > 1/14` — the closest non-tight single-scale family to the `1/14` floor. -/
theorem tenSwap20_lonely : Lonely 14 (tenSwap 2) ((8 : ℝ) / 27) := by
  have h := tenSwap_lonely 2 (by norm_num)
  norm_num at h
  exact h

/-- The **Goddyn-Wong tight family** `GW = {1,…,11, 13, 24} = AP[12 → 24]` (mac-mini THM-612,
verified kps-S54) — the SECOND `M = 1/14` family besides the AP `{1..13}`, correcting the
kps-S53 "AP is the unique tight family" claim (MISTAKE-100 recurrence).  It is the `k=2`
boundary of the residue-liar `{1..11,13,12k}`: for `k ≥ 3` that family has `M = k/(12k+5)`,
but at `k = 2` it drops to the TIGHT value `1/14`. -/
def gw : Fin 13 → ℤ := ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24]

/-- **GW is lonely at `t = 1/14`** (margin exactly `1/14`, with equality — a tight family).
At `t = 1/14` every runner's residue mod 14 lies in `[1,13]` (`24 ≡ 10`), so the lattice
distance is `≥ 1` and `14 ≤ 14·1`. -/
theorem gw_lonely : Lonely 14 gw (((1 : ℤ) : ℝ) / ((14 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio gw 1 14 (by norm_num)
  intro i m
  have key : ∀ val qq r : ℤ, val * 1 = qq * 14 + r → 1 ≤ r → r ≤ 14 - 1 →
      (14 : ℤ) ≤ 14 * |val * 1 - m * 14| := by
    intro val qq r hEq hlo hhi
    have hm : (1 : ℤ) ≤ |val * 1 - m * 14| :=
      lattice_dist_ge qq r (by norm_num) hEq hlo hhi m
    nlinarith [hm]
  fin_cases i <;> simp only [gw]
  · exact key 1 0 1 (by ring) (by omega) (by omega)
  · exact key 2 0 2 (by ring) (by omega) (by omega)
  · exact key 3 0 3 (by ring) (by omega) (by omega)
  · exact key 4 0 4 (by ring) (by omega) (by omega)
  · exact key 5 0 5 (by ring) (by omega) (by omega)
  · exact key 6 0 6 (by ring) (by omega) (by omega)
  · exact key 7 0 7 (by ring) (by omega) (by omega)
  · exact key 8 0 8 (by ring) (by omega) (by omega)
  · exact key 9 0 9 (by ring) (by omega) (by omega)
  · exact key 10 0 10 (by ring) (by omega) (by omega)
  · exact key 11 0 11 (by ring) (by omega) (by omega)
  · exact key 13 0 13 (by ring) (by omega) (by omega)
  · exact key 24 1 10 (by ring) (by omega) (by omega)

/-- The **j=13 ladder** `{1,…,12, 13k} = AP[13 → 13k]` — the cleanest one-swap ladder: it
**contains the AP itself as `k = 1`** (`{1..12,13} = {1..13}`), and has `M = k/(13k+1)` for
**every `k ≥ 1`** (`→ 1/13`).  Witness `t = k/(13k+1)`; residues `r_v = v·k` for `v ≤ 12`
(binding lower at `v=1`, `r=k`) and `r = 12k+1` for `v=13k` (binding upper), so lattice
distance `≥ k`, and `14k ≥ 13k+1 ⟺ k ≥ 1`. -/
def thirteenLadder (k : ℤ) : Fin 13 → ℤ :=
  ![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 * k]

/-- **THE j=13 LADDER IS LONELY** at `t = k/(13k+1)` for every `k ≥ 1` — an infinite family
that includes the arithmetic progression `{1..13}` (at `k=1`, `M = 1/14`) and every dilated
"top-swap" above it (`{1..12,26}` at `k=2`, `M = 2/27`, …).  Kernel-pure. -/
theorem thirteenLadder_lonely (k : ℤ) (hk : 1 ≤ k) :
    Lonely 14 (thirteenLadder k) (((k : ℤ) : ℝ) / ((13 * k + 1 : ℤ) : ℝ)) := by
  apply lonely14_of_ratio (thirteenLadder k) k (13 * k + 1) (by omega)
  intro i m
  have key : ∀ val qq r : ℤ, val * k = qq * (13 * k + 1) + r →
      k ≤ r → r ≤ (13 * k + 1) - k →
      (13 * k + 1) ≤ 14 * |val * k - m * (13 * k + 1)| := by
    intro val qq r hEq hlo hhi
    have hm : k ≤ |val * k - m * (13 * k + 1)| :=
      lattice_dist_ge qq r (by omega) hEq hlo hhi m
    nlinarith [hm, hk]
  fin_cases i <;> simp only [thirteenLadder]
  · exact key 1 0 k (by ring) (by omega) (by omega)
  · exact key 2 0 (2 * k) (by ring) (by omega) (by omega)
  · exact key 3 0 (3 * k) (by ring) (by omega) (by omega)
  · exact key 4 0 (4 * k) (by ring) (by omega) (by omega)
  · exact key 5 0 (5 * k) (by ring) (by omega) (by omega)
  · exact key 6 0 (6 * k) (by ring) (by omega) (by omega)
  · exact key 7 0 (7 * k) (by ring) (by omega) (by omega)
  · exact key 8 0 (8 * k) (by ring) (by omega) (by omega)
  · exact key 9 0 (9 * k) (by ring) (by omega) (by omega)
  · exact key 10 0 (10 * k) (by ring) (by omega) (by omega)
  · exact key 11 0 (11 * k) (by ring) (by omega) (by omega)
  · exact key 12 0 (12 * k) (by ring) (by omega) (by omega)
  · exact key (13 * k) (k - 1) (12 * k + 1) (by ring) (by omega) (by omega)

/-- The **arithmetic progression `{1..13}` is lonely at `t = 1/14`** (margin exactly `1/14`),
recovered as the `k=1` member of the j=13 ladder — the canonical LRC(14) tight family. -/
theorem ap_lonely : Lonely 14 (thirteenLadder 1) ((1 : ℝ) / 14) := by
  have h := thirteenLadder_lonely 1 (by norm_num)
  simpa using h

#print axioms tenSwap_lonely
#print axioms tenSwap20_lonely
#print axioms gw_lonely
#print axioms thirteenLadder_lonely
#print axioms ap_lonely

end LRC14
end LonelyRunner
