---
id: THM-547
title: Boundary-collar closure for LRC(14) — every cluster with second-largest co-offset ≤ 14 and a single far element w > 14 satisfies p0(E) ≤ cap_k, by THM-546 (sharpened) for w > w* plus a FINITE check for 14 < w ≤ w*, with w* = 54/90/103 for k = 8/9/10
status: REDUCTION PROVED (rigorous: Plat(E')≤Qb(k−1) and V(E')≤V_max are finite maxima over E'⊆[0,14]; THM-546 sharpened gives |Δ_w|≤(6/49)V(E')/w; for w>w* the sum is <cap_k). The residual FINITE check 14<w≤w* is feasible (spot-verified, no violations) but not yet exhaustively run. Closes codex HYP-2675's boundary-collar branch modulo that finite computation. LRC(14) NOT proved (the true-wide branch, second-largest>14, remains — the Ruzsa/Plünnecke program).
source: mac-mini-2026-06-20-S2
depends_on:
  - THM-546   # the far-element decorrelation bound, sharpened to |Δ_w|≤(6/49)V(E')/w
  - HYP-2642  # the far-element recursion p0(E)=Plat(E')+Δ_w
  - HYP-2675  # codex's split of span(E)>14 into boundary collar (2nd-largest≤14) + true-wide
related:
  - HYP-2657  # QR reality enabling the signed (6/49) constant
  - HYP-2637  # Freiman-dimension penalty (the true-wide complement)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case = 13 speeds).
---

# THM-547 — Boundary-collar closure

## Statement

Let `E` be a primitive cluster (co-offsets, `0∈E`, `|E|=k`) with **second-largest element ≤ 14**
and `w := max(E) > 14` (codex HYP-2675's *boundary collar*). Then `E' := E∖{w} ⊆ {0,…,14}`, and
> `p0(E) ≤ cap_k`  for all such `E`,
via the explicit dichotomy below, with `cap_8 = 2243/5880`, `cap_9 = 1979/4004`, `cap_10 ≥ 4/7`.

## Proof (reduction)

The far-element recursion (HYP-2642) gives `p0(E) = Plat(E') + Δ_w`, `Plat(E') := p0(E')+(1/7)p1(E')`.
Because `E' ⊆ {0,…,14}` is **bounded**, two finite maxima control it (both computed exactly over the
`C(14,k−1)` admissible `E'`):
- `Plat(E') ≤ Qb(k−1) := max_{E'⊆[0,14], |E'|=k−1} Plat(E')`, and
- `V(E') ≤ V_max(k−1) := max_{E'⊆[0,14]} V(E')`  (`V` = arc-complexity, THM-546).

THM-546 sharpened (the signed Abel bound) gives `|Δ_w| ≤ (6/49)·V(E')/w ≤ (6/49)·V_max/w`. Hence
> `p0(E) ≤ Qb(k−1) + (6/49)·V_max/w`.

Define the **collar cutoff** `w*(k) := (6/49)·V_max(k−1) / (cap_k − Qb(k−1))`. The computed values:

| k | #E' | Qb(k−1) | cap_k | margin | V_max | **w*** |
|---|-----|---------|-------|--------|-------|--------|
| 8 | 3003 | 0.19660 | 0.38146 | 0.18486 | 81 | **54** |
| 9 | 3432 | 0.36210 | 0.49426 | 0.13216 | 97 | **90** |
| 10 | 3003 | 0.44789 | 0.57143 | 0.12354 | 104 | **103** |

(`margin = cap_k − Qb(k−1) > 0` in every row — already a finite-check fact: the bounded `(k−1)`-extremal
`Qb` sits strictly below `cap_k`.) Two cases:

1. **`w > w*`:** `(6/49)V_max/w < margin`, so `p0(E) < Qb(k−1) + margin = cap_k`. **Closed, rigorously**,
   by THM-546 + the two finite maxima.
2. **`14 < w ≤ w*`:** a **finite check** — `E' ⊆ {0,…,14}` of size `k−1`, `w ∈ {15,…,w*}`. Feasible
   (`≈ C(14,k−1)·(w*−14)` clusters, ~`2.6·10^5` at `k=9`); spot-verified with no violation. This is the
   natural extension of kps-S19's proved finite half (`max(E) ≤ 14`) to `max(E) ≤ w*` while keeping the
   second-largest `≤ 14`. ∎ (modulo running case 2 exhaustively)

## Remarks

- **The sharpening matters.** With the *full* signed sum (5–76× tighter than `(6/49)V`, HYP-2657/S2),
  `w*` drops by that factor to `~10–20`, shrinking case 2 to little more than kps-S19's existing check.
- **cap_10.** The table uses the conservative floor `cap_10 ≥ 4/7` (THM-535); the true cap is larger, so
  the real margin is wider and `w*` smaller — the row closes *a fortiori*.
- **What remains:** the **true-wide** branch (second-largest `> 14`, ≥ 2 far elements). There `E'` is itself
  wide, `V(E')` unbounded, and the closed-form bound is loose; that branch needs the Ruzsa/Plünnecke +
  Freiman-dimension program (small doubling ⇒ low-dim GAP ⇒ scale-invariance/dimension penalty; large
  doubling ⇒ dissociated ⇒ `p0 ≈ M7`). See HYP-2678.

## Net

With kps-S19 (finite half, `max ≤ 14`) + THM-547 (boundary collar, `2nd-largest ≤ 14`), **two of the
three regions of the LRC(14) sector crux are closed** (the collar modulo a feasible finite check). Only the
true-wide region (`2nd-largest > 14`) is open.
