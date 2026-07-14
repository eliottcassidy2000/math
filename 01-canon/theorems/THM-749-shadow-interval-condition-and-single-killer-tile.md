---
id: THM-749
title: The shadow-interval condition and the single-killer shadow theorem — a covering-case TILE. (A) An EXACT sufficient loneliness condition at a bounded-height rational t=a/k+δ (k≤13), rigorous, in closed rational form; (B) every single-killer covering 13-set {1,…,12,182m} is 1/14-lonely by an elementary k=13 shadow witness (a THIRD proof of the covering-min class); (C) honest scope — this is ONE tile (single-killer + tight/packed), NOT uniform: it provably misses near-AP-with-far families (verified counterexample {1,2,3,4,5,7,8,9,10,11,12,13,182}, M=2/23, no k≤13 shadow), which are covered by kps THM-734.
status: PROVED (A: elementary interval arithmetic; B: elementary, all m; C: verified counterexample + honest scope). Part of the covering-case triangulation (klein-S302). Corrects the "uniform shadow closure" framing: the shadow is a tile, the covering case is the UNION of tiles (shadow + kps near-AP + opus density).
source: mac-mini-2026-07-14-S97/S99 (canonizing the shadow route with the correct scope)
depends_on:
  - THM-744   # klein's k=2 shadow (tight clusters, ratio<13) — the tight sub-tile
  - THM-734   # kps near-AP closure — the tile that covers what the shadow misses
  - THM-736   # mac-mini Farey-disc single-killer (a different proof of the same class as B)
  - THM-724   # balance single-killer (a third proof of the same class)
related:
  - THM-731   # klein disc_v certificate (the density route; the shadow is the elementary alternative for single-killer)
  - HYP-6625  # the exact condition + census (S97)
  - HYP-6660  # the correction: low-M escapees exist; complete tiling map (S99)
external: LRC(≤13) SETTLED.
---

# THM-749 — The shadow-interval condition and the single-killer shadow tile

**One line.** For a covering 13-set, an explicit **middle** lonely time can often be exhibited at a
bounded-height rational `t = a/k + δ` (`k ≤ 13`); part (A) gives the **exact closed-form δ-interval**
(a rigorous sufficient loneliness condition), part (B) uses it to prove **every single-killer covering
family is 1/14-lonely** in ~6 lines, and part (C) delimits the tile honestly: it does **not** reach
near-AP-with-far families, which are a separate tile (kps THM-734). The covering case is the **union**
of tiles, not a uniform shadow closure.

## Setup

`c := 1/14`. For a speed `w`, `‖wt‖ ≥ 1/14` is *safe*. Fix a reduced fraction `a/k`
(`gcd(a,k)=1`, `1 ≤ a ≤ k−1`) with `1/14 ≤ a/k ≤ 13/14` (so the observer speed `1` is middle-safe at
`a/k`). Consider `t = a/k + δ`, `δ > 0` small. For a speed `c`, write `r = (ca) \bmod k` and the
**signed residue** `s = r` if `r ≤ k/2`, else `s = r − k` (so `|s| ≤ k/2`).

## (A) The exact shadow-interval condition (rigorous, sufficient)

`ct \bmod 1 = r/k + cδ` (since `ca/k` is an integer plus `r/k`). Hence for `δ > 0` small:

- **`r = 0` (`k ∣ c`, the *shadow*):** `‖ct‖ = ‖cδ‖ = cδ`, safe `⟺ δ ∈ [1/(14c),\, 13/(14c)]`.
- **`r ≠ 0`:** `‖ct‖` starts at `|s|/k ≥ 1/k ≥ 1/13 > 1/14` (safe at `δ=0`) and drifts linearly. It
  leaves `[1/14, 13/14]` first at
  `U_c = (13/14 − |s|/k)/c` if `s > 0` (value drifts up past `13/14`), or
  `U_c = (|s|/k − 1/14)/c` if `s < 0` (value at `1 − |s|/k` drifts up toward `1`, i.e. `‖ct‖`
  drifts **down** past `1/14`).

> **Condition (A).** `t = a/k + δ` is lonely for all `δ` in
> `I(a,k) := \big[\ \max_{c:\,k\mid c} \tfrac{1}{14c},\ \ \min_c U_c\ \big]`
> (with `U_c = 13/(14c)` for `k∣c`). **If `I(a,k) ≠ ∅`, then `{1}∪S` is 1/14-lonely**
> (`M ≥ 1/14`), witnessed at any `t = a/k + δ`, `δ ∈ I(a,k)`.

This is a **rigorous sufficient** condition, in exact rational arithmetic, decidable per `(a,k)`.
(Verified against true loneliness on every family tested: the interval midpoint is exactly lonely,
mac-mini-S97.) The **binding upper bound** is typically a *drift-down killer*: a speed `c ≡ −a^{-1}
\pmod k` (signed residue `s=−1`) with `c` large, forcing `U_c = (14−k)/(14kc)` small.

## (B) Single-killer, proved via the `k=13` shadow (all `m`)

A **single-killer** covering 13-set has exactly one outlier `≥ 13`; the 12 speeds `≤ 12` are then
forced to be `{1,…,12}`, and covering the missing moduli `13, 14` forces the outlier
`= 182m` (`182 = \mathrm{lcm}(13,14)`), `m ≥ 1`. So the class is exactly `\{1,…,12,182m\}`.

**Claim.** For every `m ≥ 1`, `t = 1/13 + δ` with `δ ∈ [\,1/(2548m),\ 1/2184\,]` is lonely; hence
`M(\{1,…,12,182m\}) ≥ 1/14`.

*Proof.* Take `a=1, k=13` (`a/k = 1/13 > 1/14`, middle-safe). Then:
- `182m` is `13`-divisible (`182 = 14·13`): the shadow, `δ ≥ 1/(14·182m) = 1/(2548m)`.
- Runner `i` (`1 ≤ i ≤ 12`) has `r = i`; for `7 ≤ i ≤ 12` the signed residue `s = i−13 < 0`, giving
  `U_i = (|s|/13 − 1/14)/i = (169 − 14i)/(182 i)`, minimized at `i = 12`: `U_{12} = 1/2184`. For
  `i ≤ 6` (`s>0`) the bound is loose, and runner `1` at `1/13` is middle-safe.
So `I(1,13) = [\,1/(2548m),\ 1/2184\,]`, nonempty because `1/(2548m) ≤ 1/2548 < 1/2184` for all
`m ≥ 1`. By (A), lonely. ∎

This is a **third** elementary proof of the covering-min class (with THM-724 balance and THM-736
Farey-disc), and — unlike the disc_v machinery it parallels — a single short interval computation.

## (C) Honest scope — a TILE, not a uniform closure

The shadow route closes **single-killer** (B) and **tight/packed** clusters (klein THM-744:
`\max(C) < 6·\min\text{-even}`, refined `6`-odd/`13`-even). It does **not** close all covering
families. **Verified counterexample:**

> `S = \{1,2,3,4,5,7,8,9,10,11,12,13,182\}` (= drop-6 near-AP `\{1..13\}\setminus\{6\}` `∪ \{182\}`)
> is covering, `M(S) = 2/23 ≈ 0.087` (low — near the covering-min region, **not** loose), and its
> lonely times sit entirely at `k ∈ \{17,23,25\}`. It has **no `k ≤ 13` shadow witness**:
> `13` lies in the core, so there are **two** `13`-carriers `\{13,182\}` with ratio `182/13 = 14 > 13`,
> collapsing the `k=13` shadow window; and the near-AP small speeds occupy the residues that block
> every smaller `k`.

Such **near-AP-with-far** families are covered by **kps THM-734** (all 13-sets with `≥ 11` speeds in
`\{1,…,14\}`). Hence the covering case is the **union of tiles**, not a uniform shadow closure:

| stratum | tile | route |
|---|---|---|
| single-killer `\{1..12,182m\}` | this THM (B), THM-724, THM-736 | shadow `k=13` / balance / Farey-disc |
| tight/packed (ratio `<13`) | THM-744 | shadow `k=2` (parity split) |
| near-AP (`≥10` in `\{1..14\}`) | THM-733/734/738 | kps Bonferroni / disc_v |
| spread, high-diameter (loose, `M ≥ 0.12`) | THM-745/746 | opus density floor + bounded-`W` |

The **low-M / binding** region is fully covered by *shadow ∪ kps near-AP*; the residual outside both
is empirically all **loose** (`M ≥ 0.12`, diameter `≥ 24`; mac-mini-S99: `0` low-M families among
`2500` sampled outside those two tiles), the domain of opus's density floor. There is **no low-M
gap** outside the tiles.

## Honest status

- **Proved:** (A) the exact sufficient condition; (B) single-killer, all `m`.
- **Corrects:** the "uniform shadow closure" framing (mac-mini-S97 overclaim; klein-S300's
  `120/120` equivalence holds on packed clusters, fails on near-AP-with-far). The shadow is one tile.
- **Open (elsewhere):** the near-AP tile beyond kps's `j≤3`, and the rigorous closure of the spread
  loose stratum (opus's `W>W0` floor + bounded-`W` finite check) — neither is a shadow question.

*Artifacts:* `04-computation/lrc14_shadow_residue_condition_macmini_S97.py`,
`lrc14_lowM_escapee_search_macmini_S99.py`, `lrc14_residual_map_macmini_S99.py` (+outs).
Credits: klein THM-744/S299/S300 (the shadow route + grid localization), kps THM-734 (the near-AP
tile), opus THM-745/746 (the loose tile), THM-724/736 (the other single-killer proofs).
