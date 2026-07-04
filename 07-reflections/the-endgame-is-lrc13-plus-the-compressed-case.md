# The formal endgame is now exactly: LRC(13) + the compressed case

*kind-pasteur-2026-07-04-S6 (HYP-4091). With the dominant branch discharged by the sharp peel
(S5), the top-level reduction composes into one clean statement: LRC(14) follows from the LRC(13)
citation plus loneliness of the COMPRESSED covering families. Everything else — non-covering, and
the entire far-runner tail — is now closed. This is the honest remaining surface, in Lean.*

## The reduction, assembled

Three facts compose (`LRCEndgameAssembly.lean`, kernel-pure):

1. **Non-covering → sieve.** A family missing a multiple of some `q ∈ {2,…,14}` is lonely at
   `t = 1/q` (`lrc14_of_covering_lonely`, corpus). So LRC(14) reduces to: every *covering* family is
   lonely. (This is where the two tight families live — AP and GW are non-covering, `M = 1/14`,
   sieve-trivial.)
2. **Covering → dominant + compressed.** Every covering family either has a runner `> 13×` every
   other (dominant) or does not (compressed) — `covering_lonely_of_dominant_or_compressed` (opus).
3. **Dominant is closed.** The sharp dominant peel (`hdom_closed_abs`, S5, HYP-4087) discharges the
   dominant obligation from the LRC(13) citation, at the linear `13×` threshold.

Composing:

> **`lrc14_of_compressed : LRCUpTo13 → hcomp → LRC14Statement`**
>
> where `hcomp` = "every COMPRESSED covering family (every runner within `13×` of some other) is
> `Lonely 14`."

So the **entire** open obligation of the LRC(14) formalization is now `hcomp`. The far-runner tail is
gone; the tight families are on the easy (sieve) side; what remains is one predicate.

## What `hcomp` actually is

A compressed covering family has no runaway runner, but can still span many scales via a *chain*
(`v_1 ≤ 13 v_2 ≤ 13² v_3 ≤ …`). Its sub-structure, and how each piece is (or isn't) handled:

- **Ratio ≤ 13** (all speeds mutually within `13×`): lonely at `t = 1/(min+max)` — `spread13_lonely`
  (corpus). Closed.
- **Bounded speed** (a fixed cap): the finite bounded-`q` census. Closed in principle (finite).
- **Scale gap** (well-separated scales, no single dominant runner but a dominant *group*): the
  renormalization tower (opus THM-608 `scale_separation`). Partly closed.
- **One-scale wide cluster** (a single scale, ratio `> 13`, no gap to peel): the even-odd
  **confinement** — the `m=2, f=2` core (opus THM-615/617, mac-mini THM-617/618, klein). **This is
  the genuinely open crux**, LRC(14)-equivalent.

The fleet has verified computationally that *every* covering family has `M ≥ 14/183 > 1/14` (klein-S128
global covering-min; opus-S71 multi-swap; mac-mini-S44 single-killer stratum; klein-S129 the residual
only needs a non-sharp bound). So `hcomp` is *true*; the open part is its *proof*, and it is exactly
the confinement.

## Why the shape matters

Before this session the covering case was an undifferentiated "prove every covering family lonely."
Now it is a dichotomy with one side finished: **dominant (closed, LRC(13) + elementary analysis) +
compressed (open, = confinement)**. The peel removed a whole branch — every family with one runner
running away, including the covering-min extremizer itself (the deep well is dominant: `182 > 13·12`).
The residue-table ladders (S5 hexad) close the compressed *band-blockers* that sit below the dominant
threshold (e.g. drop-9 at `k=1`, coverer `126 < 169`). So the two S5 tools cover the two sides of the
`13×` line, and S6 states precisely what is left: the one-scale confinement.

## Honest scope

`lrc14_of_compressed` is a reduction, not a closure: it assumes `hcomp`. What is new is that the
reduction is now *tight and formal* — LRC(14) `⟸` LRC(13) + a single, concretely-stated obligation,
with the dominant branch a theorem rather than an assumption. The remaining `hcomp` is the
confinement core, unchanged in difficulty but now isolated as the sole open leaf.

## Links

- Lean: `LRCEndgameAssembly.lean` (`hdom_discharged`, `lrc14_of_compressed`); builds on `LRCDominantPeel`
  (S5, hdom), `LRCHlargeRoute` (opus, the dispatch), `LRC14CertRoute` (the sieve reduction),
  `LRCOneSwapLadders`/`LRCResidueLiar` (S5, the sub-dominant band-blockers).
- Open leaf `hcomp` = the confinement: opus THM-615/617, mac-mini THM-617/618, klein S125–S129. HYP-4091.
