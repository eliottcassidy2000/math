---
id: HYP-2280
title: LRC(19) — the unramified clean case; shell-dodge(m≤37) ∪ B′ covers C′(19)
status: OPEN (LRC(19) open in literature, n>13); constructive C′(19) coverage VERIFIED on 20000 configs
source: claudebox-2026-06-03-S630
related:
  - THM-398   # LRC(n) <= C'(n)
  - THM-412   # twisted-shell dodge
  - THM-415   # quantitative C': min M = 2/(2n-1)
  - THM-407   # 2n-1 prime w/ 2 primitive root = unramified (shell-transitive)
  - HYP-2240  # the cyclotomic-tower refinement
---

# HYP-2280 — LRC(19), the unramified case

LRC(n) is a theorem only for `n ≤ 13`; `n = 19` is open. It is the **favourable** open case:
`2n-1 = 37` is **prime with 2 a primitive root** — the *unramified*, doubling-shell-transitive case
(THM-407), opposite to `n=14` (`27=3³`, ramified). A good multiplier at shell 37 gives
`M ≥ 2/37 = 0.05405 > 1/19 = 0.05263` (loose).

## Result (constructive C′(19), verified)

`LRC(19) ⟸ C'(19)` (THM-398): a multiple-of-19 config is loose. **Verified:
shell-dodge(`m ≤ 2n-1 = 37`) ∪ B′(dominance) covers 100% of 20000 random primitive multiple-of-19
configs — 0 residual** (`lrc19_quant_s630.out`); in fact the twisted dodge `m ≤ 37` alone covers all
(B′-only = 0). So every multiple-of-19 config has an explicit loose witness `t = a/m`, `m ≤ 37`. This
is the constructive route (HYP-2197/2240) closing cleanly in the unramified case.

## The crux configs (shell-37-dodge failures, still loose via lower clocks)

A config fails the shell-37 dodge iff its 18 residues mod 37 cover all 18 `±`-pairs of `(ℤ/37)*`.
These are *not* counterexamples — they escape on other shells:
- **`{1,…,17,19}`** (the AP `{1,…,18}` with `18 → 19`): residues mod 37 = `{1,…,17,19}` cover all
  pairs, so dodge-37 fails — but **exact `M = 1/18`** (loose via the freed 18-clock `t=1/18`: dropping
  18 from the speeds frees that clock). It blocks the 19-clock, not the 18-clock.
- Generic shell-37 failures (26126/60000) are loose via shells 8,14,17,28,29,33,… with `M = 0.10–0.14`.

So shell 37 is the **extremal** shell (where `min M = 2/37` lives), not the generic escape; the bulk
escape much higher, only the tightest multiple-of-19 configs need shell 37.

## Open (the mini-wall)
- **The `2/37` minimizer at n=19 eluded search** (min `gap_shells` found `= 1/18` over AP-variants +
  20000 random). Confirming `min M = 2/37` (THM-415 pattern) needs the *extremal generating rule*
  (HYP-2196 C4) — unknown at this size. Either the minimizer is a special irregular config (as at
  n=5,8) or the pattern shifts toward `1/(n-1)` for large unramified `n`. **A reframe is needed here**
  — see the S630 exploration thread.
- Turning "0 residual / 20000" into a proof of C′(19) = the residue-profile DP (HYP-2256, t-0083) at
  `L = lcm(2..37)`, the box-free finite check.

## The reframe (S630 exploration): the minimizer is a red herring

Failed reframes (documented, `lrc19_foldrich_s630.out`): "min M = max 3-term folds" is FALSE (random
configs are all fold-rich at `M≈0.11`); Paley/QR mod 37 gives `M=2/23`, not `2/37`. Decoded small-n
minimizers: witness multiplier sits near the apex `a≈(2n-1)/2≈n` (n=4,5,6,8; breaks at n=7) — the
extremal is genuinely irregular (the open generating rule, HYP-2196 C4).

**The gem:** LRC(19) does not need the extremal config — it needs COVERAGE (every multiple-of-19
config loose), and coverage is exactly what the **residue-profile reduction** (HYP-2256) makes FINITE:
`gap_shells` reads speeds only mod `L=lcm(2..37)`, so "is every multiple-of-19 config loose?" is a
question about finitely many residue profiles, box-independent. In the unramified case (37 prime) the
shell has no ramified strata to fragment the check. **Path: run the residue-profile DP mod L; confirm
dodge(m≤37) ∪ B′ covers every profile ⟹ C′(19) ⟹ LRC(19) — the minimizer never enters.** Generalizes
to all unramified `n` (`2n-1` prime). See `07-reflections/lrc19-the-minimizer-is-a-red-herring-s630.md`.
