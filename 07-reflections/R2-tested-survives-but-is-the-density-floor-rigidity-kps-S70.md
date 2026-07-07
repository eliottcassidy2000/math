---
source: kind-pasteur-2026-07-07-S70
status: HONEST NEGATIVE ON A PROOF + evidence-strong survival + a rigorous local-minimality
  crumb + a reduction. Owner: prove R2 (spread APs are the global PA_2-minimizer). I could
  not prove it; I tested it hard (it survives), pinned down exactly why it is hard, and
  reduced (A') to it. Reporting honestly rather than fabricating a proof (MISTAKES culture).
tags:
  - lonely-runner
  - LRC14
  - anchor-floor
  - rigidity
  - honest-negative
  - density-floor
---

# R2 survives a hard test — but it is the density-floor rigidity, and I could not prove it

**kind-pasteur-2026-07-07-S70 (HYP-5077).** The task was to prove **R2**: among all `k`-element
integer families, `PA_2(E) = P(max(gap@0, gap@1/2) > 1/7)` is minimized by a spread AP. I did
not prove it. What I did: tested it hard (it survives), established the practically important
fact it was standing in for, secured one rigorous sub-result, and diagnosed precisely why a
proof is hard. This is the honest outcome, and it is worth more than a fabricated proof.

## The reframe (a clean rigorous identity)

`gap@0(E,x) > 1/7 ⟺ max_i frac(e_i x) − min_i frac(e_i x) < 6/7` — the origin gap exceeds
`1/7` iff the phases fit in a `6/7`-arc (whose complementary `1/7`-arc holds `0`). So

> `PA_2(E) = P(` the phases avoid a `>1/7` arc at `0` **or** at `1/2` `)`,

a two-point avoidance functional. Minimizing it wants families whose phases sit near **both**
`0` and `1/2` for most `x`. The spread AP (Steinhaus/three-gap config) is the natural
candidate — maximally spread, hits both anchors often — which is *why* R2 is plausible.

## What the strong adversary found

The first cut was weak (the classic MISTAKE-102 trap: local descent from random seeds got
stuck at `0.86` for `k=8`, never near the AP min `0.72`). Strengthened with a proper
local-minimum test around the exact winner, plus a structured large-diameter adversary
(spread/dilated/defected APs, GW, interleaves, two-block, geometric, Sidon, dilated-large):

- **Q1 — the fact `(A')` actually needs — `inf_E PA_2(E) ≥ T_k` HOLDS robustly.** Zero families
  dipped below `T_k` (`0.6185` at `k=8`, `0.0565` at `k=13`) over the entire adversary. So the
  2-anchor route is intact independent of R2.
- **R2 survives.** The spread AP is a **rigorous exact 1-move local minimum** at `k=8`: all 16
  single-element neighbors of `{5,7,…,19}` have exact rational `PA_2 ≥ 1164298/1616615` (nearest
  neighbor `0.7656` vs the winner `0.7202`). Local minimum at `k=13` too (`d=2, a=15`,
  `PA_2 = 0.3188`). No structured or large-diameter family beat it.

## Why I could not prove it (the honest diagnosis)

R2 is a **global extremal rigidity**, and three independent facts say it is the hard kind:

1. **It is the same statement, one anchor-family down, as the density-floor `μ` AP-minimality.**
   `μ_{1/7}(E) = P(maxgap > 1/7) ≥ μ_{1/7}(AP)` is the fleet's load-bearing `(A')`, verified
   exhaustively only to `k ≤ 10` (opus-S131) and **open in general**. R2 is its 2-anchor
   analog. Proving R2 in the regime that matters is not obviously easier than the thing it
   feeds.
2. **The landscape has no majorization structure.** My S63 refuted Schur-monotonicity of exactly
   this class of functional under step transfers (77% of equalizing moves go *uphill* — the
   resonance-rugged landscape). So there is **no clean rearrangement/majorization proof** that
   the most-even family (the AP) minimizes; the natural proof template is dead.
3. **It is the σ-even measure core.** In the S67 grading, R2 is a `σ`-even (measure) extremal
   statement — exactly the part that resists the covering/parity/moment tools, and closes only
   by finitization. My S68 (bounded diameter) and S69 (exact spread-AP inf) *finitized* their
   slices; R2 is the un-finitized global-minimality that remains.

So the honest verdict: **R2 is true on all evidence but is a genuine open rigidity of
density-floor hardness; I did not prove it.**

## Where this leaves `(A')` (the useful reduction)

```
(A')  inf_E PA_2(E) ≥ T_k   ⟺
    bounded primitive diameter D ≤ D_0(k):  PA_2(E) ≥ PA_2(AP_{D+1}) ≥ T_k   [S68, ALL families, PROVED]
    unbounded diameter:                     spread-AP inf ≥ T_k              [S69, EXACT]
                                            + spread AP is the minimizer      [R2, the single open rigidity]
```

Two honest silver linings:
- **Q1 is robust:** `inf PA_2 ≥ T_k` survived a strong adversary, so the practical bound `(A')`
  needs is well-supported even before R2 is proved.
- **R2 might be strictly easier than the density floor:** it is a bound on two fixed anchor gaps,
  not the max over all `k` gaps (boxeph's original point). A concrete, cleaner future target than
  the full `μ` AP-minimality — or bypass it entirely by closing the unbounded regime with an
  explicit Erdős–Turán decorrelation bound plus a stability neighborhood of the spread AP (both
  routine analysis, both currently un-executed).

## Ledger

- Rigorous: the range reframe (identity); the exact 1-move local minimality of the spread AP at
  `k=8` (`1164298/1616615`, all 16 neighbors ≥ it).
- Evidence: R2 survives a strong structured + local + descent adversary; `inf PA_2 ≥ T_k` (Q1)
  holds on everything tested (k=8, 13).
- NOT proved: R2 (global minimality) — diagnosed as the density-floor rigidity (σ-even, no
  majorization, same class as μ AP-minimality).
- Files: `lrc_R2_adversary_kps_S70.py`, `lrc_R2_strong_kps_S70.py` (+outs).
- Builds on: boxeph-S1 (2-anchor; R2 = its open claim), my S68/S69 (finitizations), S63
  (majorization refuted), S67 (σ-grading), opus-S131 (μ AP-min = same class), MISTAKE-102.
- Does NOT prove LRC(14) or R2. It is an honest test-and-diagnose: R2 stands as the single
  remaining rigidity, evidence-strong, density-floor-hard.
