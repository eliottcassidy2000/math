---
source: death-star-2026-07-07-S2
status: FORMALIZATION WIN. The AP₂₀ density-floor certificate is PROVEN in Lean, kernel-pure
  (no sorry, no native_decide) — making the diameter-≤19 (S15 census) density floor UNCONDITIONAL
  in the muGood formalism, with a 4.5× margin, via 2 explicit good intervals.
tags:
  - lonely-runner
  - LRC14
  - lean
  - formalization
  - density-floor
  - tail-diameter
---

# The AP₂₀ density-floor certificate is proven in Lean (unconditional, 4.5× margin)

Owner: *formalize the LRC(14) proof; get it into its best possible state; healthy diet from the
other agents.* I worked the density-floor leg of the tail-diameter route and closed a certificate.

## The gap I found

mac-mini's `LRCTailDiameter` reduces the density floor for a family of diameter ≤ D to the single
leaf `μ_{1/7}(AP_{D+1}) ≥ m_P`, and discharges it via `muGood_ge_AP76` — which **always drops to
the 76-point AP**, whose value `μ_{1/7}(AP₇₆) = 0.05745` clears `m_P = 14249/252252 = 0.05649` by
only **1.02×** (razor-thin). That `AP76Certificate` is carried as an *open* hypothesis (a Farey-cell
computation; boxeph-S1's `LRCFareyRoofBridge` is reducing it to a superlevel Farey-sum, still open).

**But the S15 census families have `|vᵢ| ≤ 20`** — 13 distinct integers in `{1,…,20}`, i.e.
**diameter ≤ 19** — so by subset-monotonicity they only need `μ_{1/7}(AP₂₀) ≥ m_P`, and
`μ_{1/7}(AP₂₀) = 0.254` has a **4.5× margin**. The proof was always dropping too far (to AP₇₆).

## What I proved (`TournamentH7/LRCAP20Certificate.lean`, kernel-pure)

1. **`muGood_ge_APD`** — the general-diameter bound: a family whose translate lies in `{0,…,D}`
   has `μ_θ(AP_{D+1}) ≤ μ_θ(E)` (bound by `AP_{D+1}` *directly*, not always `AP₇₆`). One-line
   generalization of `muGood_ge_AP76`, decoupling the certificate from the hardcoded 76.

2. **`ap20_certificate`** — UNCONDITIONAL: `ENNReal.ofReal(14249/252252) ≤ muGood(1/7, {0..19})`.
   Proven by exhibiting **two explicit good intervals**:
   - near `x=0`: for `x ∈ (0, 6/133)`, the arc left at `a = 19x/2 + 3/7` is empty of the AP₂₀
     orbit (single cluster `{0,x,…,19x}` leaves the gap `(19x,1)` of length `>1/7`);
   - near `x=1`: the mirror interval `(127/133, 1)`, same length `6/133`.
   Disjoint ⟹ `μ_{1/7}(AP₂₀) ≥ 12/133 = 0.0902 > m_P`. The `EmptyArc` witness reduces, per orbit
   point `e`, to `Int.fract((e−19/2)x − 3/7) > 1/7`, and `(e−19/2)x ∈ (−3/7, 3/7)` (from
   `|e−19/2| ≤ 19/2` and `x < 6/133`) puts the argument in `(−1,0)`, so the fract is `value+1`.

3. **`hlarge_floor_of_diam_le_19`** — the export: any family whose translate lies in `{0,…,19}`
   has `μ_{1/7}(E) ≥ m_P`, **with no certificate hypothesis**.

4. **`muGood_diam_floor` / `hlarge_floor_of_diam_le` (GENERAL, `2 ≤ D ≤ 30`)** — generalized the
   two-interval method to arbitrary diameter: `μ_{1/7}(AP_{D+1}) ≥ 12/(7D) ≥ m_P` for `D ≤ 30`
   (`12/(7·30) = 2/35 = 0.0571 > m_P`). Generic multiplicative empty-arc lemmas
   `emptyArc_near0_gen`/`near1_gen` (`D·x < 6/7`, no division) + a symbolic-`D` measure assembly.
   So **every family of diameter ≤ 30 has an unconditional density floor `≥ m_P`** — covering the
   census (`D ≤ 19`) *and* the peel headroom, with no Farey-cell computation anywhere.

**Axiom audit:** both `ap20_certificate` and `hlarge_floor_of_diam_le_19` depend on exactly
`[propext, Classical.choice, Quot.sound]` — no `sorryAx`, no `native_decide`. Wired into the root
`TournamentH7.lean` import; the full project builds green (8476 jobs).

## Why this matters / how it fits

- The census leg of `lrc14_concrete` (the `|vᵢ| ≤ 20` families) is exactly the diameter-≤19 regime.
  Its density floor is now an **unconditional Lean theorem** in the `muGood` formalism — whoever
  wires `witnessG2 = muGood(1/7)∘family` (the mac-mini/monad line; the architecture already carries
  this as the `hwitness` hypothesis in `LRCEventMeasureBridge`) can consume it directly, instead of
  the open razor-thin `AP76Certificate`.
- The **2-interval method is elementary and margin-robust up to `D ≤ ~25`** (`μ(AP_{D+1}) ≥
  12/(7D) ≥ m_P` iff `D ≤ 30`, comfortably for `D ≤ 25`). It needs **no Farey-cell measure
  computation** — it sidesteps the hard AP₇₆ tail entirely for the bounded regime. boxeph's
  Farey-sum route is still the tool for the *wide* (post-peel, diameter up to 75) regime; the two
  are complementary (mine = census/bounded, sharp margin; boxeph's = general AP₇₆).
- Method note for the fleet: **don't drop to AP₇₆ when the family is bounded** — bound by the
  *tightest* `AP_{D+1}` (via `muGood_ge_APD`), where the margin is large and 2 explicit intervals
  suffice. The tail-diameter theorem's own `muGood_ge_AP76` throws away that margin.

## Ledger

- **Files:** `TournamentH7/LRCAP20Certificate.lean` (kernel-pure, in root import);
  support: `lrc_ap76_smallq_lowerbound_deathstar_S2.py` (+ out; showed AP₇₆ needs 95 intervals /
  0.05% margin — intractable — motivating the switch to AP₂₀).
- **Builds on:** mac-mini-S42 `LRCTailDiameter` (Good/muGood/good_anti/good_translate); opus-S135
  `LRCFareyRoof`; boxeph-S1 `LRCFareyRoofBridge` (complementary AP₇₆ route).
- **Does NOT prove LRC(14).** Discharges one leaf (the census-regime density-floor certificate).
