# The concordance decision, the E_n dual, and how the minimal test saved us

*opus-2026-07-14-S308. Owner directive: compute the E_n dual Smith diagram,
prove the concordance law, and validate above n=7 minimally — "make sure we are
not misguided." The last clause was the important one.*

## 1. The decision: the concordance law is FALSE beyond small n

S307's headline said the harmonic potential of G_n orders classes exactly like
the score-variance axis (exact through n=6; 132 float-suspect pairs at n=7).
Before attempting a proof, the suspect pairs were DECIDED (float64 LU + exact
integer-residual iterative refinement; zero pairs left undecided):

- **n = 7: the 132 discordances are REAL.** Exact concordance fails. But every
  one sits between ADJACENT axis levels (gap exactly 8) — the corrected
  candidate was "coarse concordance" (strict order at level-distance ≥ 16).
- **n = 8 (the first above-7 test, 6,880-class certified network): coarse
  concordance FAILS TOO.** 65,921 discordant pairs (0.31%), now reaching
  level-distances 16, 24, and 32: {8: 61,954; 16: 3,811; 24: 150; 32: 6}.

The trajectory — exact (n ≤ 6) → adjacent-only (n = 7) → spreading (n = 8) —
says the strong forms are SMALL-n ARTIFACTS. Had we proved-by-effort what the
small cases suggested, we would have canonized a false law. The one-session
sequence [decide the suspect pairs exactly → jump one n with the cheapest
certified instrument] is the minimal misguidedness test, and it worked.

## 2. What survives (and is worth proving)

- **Level-MEAN and level-MEDIAN monotonicity hold at every computed n**
  (5, 6, 7 verified; n=8 check running at write time): the potential's level
  averages descend strictly along the axis. The degradation is in the TAILS:
  the adjacent-level overlap/mean-gap ratio runs −0.611 (n=5, separated),
  −0.245 (n=6, separated), **+0.726 (n=7, interpenetrating)** — the level
  intervals began to merge at exactly the n where exact concordance broke.
- The max-principle endpoints are trivially exact: the transitive (source) has
  the strict maximum potential; the distributed rail the minimum.
- The honest statement of the phenomenon: **the axis is a mean-field harmonic
  coordinate — exact for level aggregates, increasingly wrong for individual
  classes as the level populations diversify.** The right provable target is
  the mean/median law, not the pairwise law.

## 3. The E_n dual Smith diagram (exact, n = 4..8*)

The cycle-space bijection makes every tiling an even graph; the dual network:
nodes = E_n classes (A002854-certified), conductances = wiggly multiplicities,
source = the EMPTY even graph's class (fiber 1 — the dual of the transitive),
sink = the maximal even graphs, axis = edge count.

| n | E_n classes | R^E_n | dual concordance |
|---|---|---|---|
| 4 | 3 | **2/5** | perfect |
| 5 | 7 | **5183/14960** | perfect |
| 6 | 16 | exact ≈ 0.109086 | ONE adjacent overlap |
| 7 | 54 | exact ≈ 0.134011 | perfect |
| 8 | 243 | (rerun pending — classifier tuning) | — |

- The dual side shows the SAME phenomenon (an adjacent-level overlap appears
  at n=6) — the mean-field reading applies to both quotients of Q_m.
- **R^E is non-monotone in n** (0.400, 0.346, 0.109, 0.134): a parity
  structure — at odd n the sink is K_n itself (all degrees even), at even n
  the maximal even graphs are further from complete; the dual geometry
  alternates.
- **The reciprocity probe is an honest negative**: R^G·R^E = 0.171, 0.089,
  0.012, 0.009 — no BSST-style R·R′ = 1 law. G_n and E_n are ALGEBRAIC
  cut/cycle duals (Ihara/Bass), not electrical-planar duals; the products
  simply track both resistances shrinking.
- Curiosity logged: R^E_6 ≈ 0.109086 vs R^G_6 ≈ 0.109041 — equal to three
  digits, not equal exactly. No mechanism proposed.

## 4. The method note (what "do as little as possible" bought)

Total compute for the decision: one exact-refined 456-node solve (minutes), one
reuse of the already-certified n=8 classifier (minutes of new code), and two
small exact dual networks. No new enumeration beyond what existed. The
discriminating power came from choosing the RIGHT next instances (the decided
suspect pairs; one step up in n on the cheapest certified instrument), not from
scale. The standing lesson joins the MISTAKE-140/147 genus from the positive
side: **before proving an empirical law, spend one session trying to break it
at the cheapest scale you have not yet certified.**
