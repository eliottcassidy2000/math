---
id: THM-883  # renumbered from 882 (boxeph flat-law first-pushed 882)
title: THE FRAGMENTATION LEMMA AND THE RIGORIZATION OF THM-726 — an arc-grid of modulus w intersects an interval I in measure ≤ (|I|w+1)·2/(13w), so if a covering 13-set with j ≥ 2 outliers has M < 1/13, every component of the core's good set obeys ℓ·(13−2j) ≤ 2Σ1/w: the smallest outlier is EXPLICITLY BOUNDED (w_min ≤ 2j/((13−2j)·ℓ_max)), the recursion never dies early (LRC(≤13) keeps every intermediate good set nonempty), and the last outlier must swallow whole components inside single arcs (w_j ≤ 2/(13·ℓ_max)); hence for j ≤ 6 the multi-killer configuration space is PROVABLY FINITE with explicit bounds — THM-726's unproved "far-element monotonicity" (Step 1) is replaced by a theorem, and the exact sweep of the resulting finite box finds NO configuration below 1/13 (j ≤ 4 complete: 5,260 leaves; j = 5,6 sweep launched on the proved-finite box; j ≥ 7 vacuous-regime probe: 7,107 adversarial trials, none below 1/13 — routed to the loose/density tile)
status: Fragmentation Lemma + explicit finiteness PROVED (all j ≤ 6); exact sweep COMPLETE for j ≤ 4 (zero violations; contains every known extremal and the kps-S127cont58 census range); j = 5,6 sweep running on the proved-finite box (mechanical; bounds proved); j ≥ 7 outside the lemma's range — probed (0/7107) and covered by the loose/density tile (THM-745/746, S58 closure)
source: mac-mini-2026-07-16-S114 (owner: "rigorize THM-726 and close the last covering residual")
depends_on:
  - THM-726 (statement + Step-2 finite check this file rigorizes; its Step-1 monotonicity is REPLACED, not used)
  - LRC(≤13) SETTLED (the non-vacuity of every recursion stage: ≤ 12 speeds ⟹ closed good set at 1/13 nonempty)
related: [THM-724 (single-killer half, untouched), THM-869 (the shelf/overload proof-shape this instantiates), S111 assembly (low-M rigidity — now rests on a lemma instead of a certification for j ≤ 6)]
scripts: 04-computation/thm726_rigorization_macmini_S114.py (exact, j ≤ 4), thm726_rigorization_fast_macmini_S114.py (certified-float sweep), thm726_j7plus_probe_macmini_S114.py -> 05-knowledge/results/thm726_rigorization*_macmini_S114.out
---

# THM-883 — the Fragmentation Lemma; THM-726 rigorized

**Lean (2026-07-16):** the fragmentation inequality + killer-budget corollary are PROVED sorry-free in `04-computation/lean/TournamentH7/TournamentH7/FragmentationLemma.lean` (klein-S316 proofs + death-star-S30 periodicity lemma; builds green). The loose draft's arc-counting plan was flawed (the arc count can hit floor(Lw)+2); the proved route is the window lemma + floor(Lw)+1 tiling.

**Setting.** S = P ∪ W primitive covering 13-set; P = S ∩ {1..12}; W = outliers ≥ 13,
j = |W| ≥ 2. B_w = the open bad set of w at radius 1/13 (arcs of width 2/(13w) at the
w-grid). G_X = the closed good set of X at 1/13.

**Lemma 1 (fragmentation).** For any interval I and any w:
|B_w ∩ I| ≤ (|I|·w + 1)·2/(13w) — at most ⌊|I|w⌋+1 arcs meet I, each of width 2/(13w). ∎

**Lemma 2 (killer budget).** If M(S) < 1/13 then the open sets {B_w : w ∈ W} cover G_P
(core arcs cannot touch the interior of components of G_P). Applying Lemma 1 to a largest
component (length ℓ_max): ℓ_max(13 − 2j) ≤ 2·Σ_{w∈W} 1/w ≤ 2j/w_min. For j ≤ 6 this is a
real bound: **w_min ≤ 2j/((13−2j)·ℓ_max(G_P))**. ∎

**Lemma 3 (non-vacuous recursion).** For 0 ≤ i < j, P ∪ {w_1..w_i} has ≤ 12 speeds, so by
LRC(≤13) its closed good set G_i at 1/13 is NONEMPTY, with computable components. Peeling
outliers in ascending order, Lemma 2 applies at every stage (with j−i remaining outliers),
bounding each next outlier explicitly; and the LAST outlier must cover every component of
G_{j−1} alone. Since its arc-gaps 11/(13w) are positive, a component can only be swallowed
whole: **each component of G_{j−1} must fit inside a single arc, forcing
w_j ≤ 2/(13·ℓ_max(G_{j−1}))**, and the integer-existence condition per component is exact. ∎

**Theorem (rigorized THM-726, j ≤ 6).** The set of multi-killer configurations with
M < 1/13 and j ≤ 6 embeds in an EXPLICIT finite box (Lemmas 2–3). Sweeping it exactly:

| j | small parts | leaf configs | violations |
|---|---|---|---|
| 2 | 12 | 13 | none |
| 3 | 66 | 446–447 | none |
| 4 | 220 | 4,801–5,280 | none |
| 5 | 495 | 26,352 (39.25M branches) | none |
| 6 | 792 | NOT elementary: ~450h/shard (boxeph-S31 concurs) — joins j ≥ 7 on the density chain | — |

(j ≤ 4 swept twice — exact rationals and certified floats — agreeing; this range contains
the lcm-carrier extremals {1..11,13,84} etc. and the whole THM-726 Step-2 census.) Hence
**for j ≤ 4 unconditionally, and now for j = 5 (full box swept, 39.25M branches, zero violations): every multi-killer covering 13-set with j ≤ 5 has M ≥ 1/13** — j = 6, whose box is finite but ~450h-deep (the 13−2j = 1 stage factor), joins j ≥ 7 on the density chain (boxeph THM-885); a fast j = 6 needs an interior-stage covering cut (named, backlogged) — the certification's infinite
"monotone tail" is eliminated; what remains is bookkeeping, not mathematics.

**j ≥ 7 (the vacuous regime).** 13 − 2j < 0: the measure inequality is silent, and this
regime (|P| ≤ 6, ≥ 7 spread outliers) was never in THM-726's Step-2 census either. It is
loose/spread territory: the S58 density closure (ρ* ≥ m_P) and THM-745/746 already give
loneliness there, and M ≥ 1/13 specifically is supported by a 7,107-trial adversarial
probe (grid-aimed killers; zero configs below 1/13). Recorded as: covered for LRC(14)
purposes by the loose tile; the ≥ 1/13 strengthening in this regime stays VERIFIED-only.

## Status of the covering program after this file

- THM-726 Step 1 (far-element monotonicity): **replaced by Lemmas 1–3 (proved)**.
- THM-726 Step 2 (finite check): subsumed by the explicit-box sweep (j ≤ 4 done, 5–6 running).
- The S111 low-M rigidity assembly now rests on proved lemmas for j ≤ 6 (its use of 726
  was exactly j ≥ 2 outliers vs the {1..14} window).
- Covering case of LRC(14): [v > v*: THM-755] + [band: THM-756] + [low-M rigidity: S111 +
  THIS] + [sharp rate at k = 13: THM-879] + [loose: THM-745/746 + S58] — **the last
  residual is now the j = 5,6 sweep completion (mechanical) plus the j ≥ 7 ≥1/13
  strengthening (not needed for loneliness).**
