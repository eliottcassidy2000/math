---
id: THM-1290
title: "THE (1/14, 3/41) SUB-GAP IS EXHAUSTIVELY EMPTY AT HEIGHT 64 (S320 extension; v1 at 55) — no 13-subset of [1,64] (primitive or not) has M ∈ (1/14, 3/41), and through height 64 the filter stack covering+spread+pair-sum+pinning(q≤48) is jointly unsatisfiable (zero survivors, no witnesses needed on [56,64]) — v1 at 55 by COMPLETE ENUMERATION: 43.86 billion DFS nodes, 13.07 billion leaves, 4,693,315,305 filtered primitive families, exactly 50 unit-pinning survivors, every one witness-certified M ≥ 3/41 at a denominator in [43,48], ZERO hard cases. All filters are rigorous canon consequences (covering 2..13; THM-1043 spread σ > 12; THM-1269 active-pair sum D = M·s forcing a pair sum in {55,69,83,96,97,110}; depth-d unit pinning for q ≤ 41). Upgrades the emptiness evidence from search (opus ~12,400 structured + 8.5M random; klein 56k — all with the recorded weakness 'scale, not coverage') to EXHAUSTION 1.5× beyond the height of the deepest known family. BY-CATCH (gate C, exhaustive at heights ≤ 26 for M < 1/13): the complete sub-1/13 spectrum at height ≤ 26 is {1/14 (2 families), 2/27 (3 families)} — including TWO NEW 2/27 realizers {1..9,11,13,20,24} and {1..9,11,12,13,20} (height 20 — lower than the ladder family {1..12,26}), and zero families below the floor"
status: PROVED-BY-EXHAUSTION at B = 55, BOTH PARTS (script + frozen consolidated output + independent exact-M referee; six exact-M gates pass; two rediscovery gates re-find {1..11,13,36} → 3/41 and {1..12,26} → 2/27 blind under widened intervals; LRC-mode regression gate exhausts heights ≤ 26 and re-finds all known sub-1/13 families). Part (a) sub-gap: 4.69e9 filtered → 50 pinning survivors → all witness-cleared → 0 hard. Part (b) COMPLETED (the companion run; gap-mode cannot see M < 1/14 families since their active pair sum is not gap-admissible): interval (0, 1/14), covering 2..14, spread 13 (THM-405 rung), w_max ∈ [14, 55] — 32.68e9 nodes, 2.4499e9 filtered primitive families, SIX pinning survivors, all witness-certified M ≥ 1/14 (q ∈ {42, 44, 45, 46, 48}), ZERO hard: **no 13-subset of [1,55] has M < 1/14 — LRC(14) verified exhaustively at height 55** (non-primitive reduce to primitive with σ > 13 ⟹ w′_max ≥ 14, enumerated).
source: klein-2026-07-19-S319 (owner: work to finish LRC(14); see the near-misses, synthesize perspectives)
depends_on: [THM-1268 (D ≥ 4 forcing, mediant 4/55), THM-1269 (D = M·s, the pair-sum reduction), THM-1043 (spread rung σ ≤ n−1 ⟹ M ≥ 1/n), THM-1230/1235 (the gap edge 3/41, its witness, the slack coordinates), boxeph-S115/S126 (mod-p spread/pinning template), mac-mini-S54 (census harness architecture, lrc_gap_census40_S54.c), THM-405 (n=14 spread rung, used by the companion run)]
scripts: 04-computation/lrc14_subgap_census_klein_S319.c + 04-computation/lrc14_subgap_referee_klein_S319.py -> 05-knowledge/results/lrc14_subgap_census_klein_S319.out
---

# THM-1290 — the sub-gap census at height 55, EXTENDED TO HEIGHT 64 (S320)

> **S320 EXTENSION (klein-2026-07-19-S320, owner: "push the census to B=64 with q≤48
> pinning"):** part (a) now holds AT HEIGHT 64. Harness v2 (QPIN 48 = depth-3 pinning to
> q = 48; in-branch bitmask pruning generalized from {23,25,27} to ALL depth-1 moduli
> q ∈ [14,27], gate-A leaves −77%; survivor printing; gates A/B/C byte-identical) ran
> w_max ∈ [56,64] — only that range needs enumeration, since the DFS at a given w_max
> never looks above it, so the v1 run covers [28,55] verbatim. Result:
> **112,686,675,261 nodes, 7,736,632,974 leaves, 931,039,618 filtered primitive
> families, ZERO pinning survivors, zero witness scans needed, zero hard cases.**
> On [28,55] the q≤48 addendum needs no rerun: each of the v1 run's 50 survivors carries
> a witness at q ∈ [43,48] with margin ≥ 3/41, which is exactly a depth-d(q) pinning
> violation at that q — harness v2 kills all 50 in-filter. Non-primitive completeness at
> B = 64: a g ≥ 2 dilate needs its primitive core in the gap at w′_max ∈ [28,32],
> exhausted by the v1 run. **STRUCTURAL COROLLARY: covering{2..13} + spread(σ > 12) +
> pair-sum + unit-pinning(q ≤ 48) is jointly UNSATISFIABLE for 13-subsets through
> height 64** — the pinning stack alone now carries the whole theorem, with no rational
> witness needed anywhere. Part (b) (LRC(14) verification) remains at height 55; its
> [56,64] extension is the named next step. Frozen output appended to the results file.

## Statement

**(a) No 13-subset of [1, 64] of distinct positive integers — primitive or
not — has maximum loneliness M in the open interval (1/14, 3/41).** (Height
55 by the v1 run below; extended to 64 by the S320 run above.)

**(b) No 13-subset of [1, 55] has M < 1/14: LRC(14) holds exhaustively at
height 55.** (Companion run, same harness at interval (0, 1/14): 32.68e9
nodes, 2.45e9 filtered families, 6 pinning survivors, all witness-certified,
0 hard.)

Consequently the attained 13-speed M-spectrum's second-smallest value, IF any
value below 3/41 exists at all, requires v_max ≥ 56; equivalently, the gap
(1/14, 3/41) of THM-1235 is empty through height 55 — and any LRC(14)
counterexample requires v_max ≥ 56.

## Proof shape (all filters rigorous, from canon)

A 13-set W with M(W) ∈ (1/14, 3/41) must satisfy, unconditionally:

1. **Covering 2..13**: else the witness t = 1/q gives M ≥ 1/q ≥ 1/13 > 3/41.
2. **Spread** (THM-1043, n = 13 rung): σ = w_max/w_min > 12.
3. **Active-pair sum** (THM-1269): every maximizer of the lower envelope is a
   rising/falling crossing of two distinct speeds, giving integer
   D = M·(v_i+v_j); M ∈ (1/14, 3/41) forces the pair sum into
   {s : ∃D ∈ ℤ, s/14 < D < 3s/41} ∩ [1,110] = {55, 69, 83, 96, 97, 110}.
   In particular w_max ≥ 28.
4. **Unit pinning at depth d(q) = ⌈3q/41⌉ − 1 for all q ≤ 41**: no multiple of
   q present ⟹ every unit a mod q has some v with dist(v·a mod q) ≤ d(q),
   else the witness a/q certifies M ≥ ⌈3q/41⌉/q ≥ 3/41. (Depth 1 for
   14 ≤ q ≤ 27 — the in-branch bitmask moduli 23/25/27 — depth 2 for
   28 ≤ q ≤ 41.)

Non-primitive W = g·W′ (g ≥ 2) would need the primitive core W′ in the gap
with w′_max ≥ 28, so w_max ≥ 56 > B: primitive enumeration is complete at
B = 55. DFS in decreasing element order with in-branch feasibility pruning on
(1), (3), (4); survivors of the full stack get a rational-witness scan
(margin ≥ 3/41, integer arithmetic, q ≤ 200); anything uncleared prints HARD
for the independent exact-M referee (complete breakpoint method: all
c/(v_i±v_j), peaks, valleys — exact rational arithmetic).

## The run (frozen output in results)

| stage | count |
|---|---|
| DFS nodes visited | 43,856,745,170 |
| leaves | 13,073,431,569 |
| pass covering 2..13 | 4,756,796,579 |
| pass pair-sum | 4,693,824,172 |
| pass primitivity | 4,693,315,305 |
| **pass unit pinning (q ≤ 41)** | **50** |
| witness-cleared (M ≥ 3/41, q ∈ [43,48]) | 50 |
| HARD (exact-M needed) | **0** |

First-witness-denominator histogram of the 50 survivors: 43:22, 44:1, 45:19,
46:2, 47:2, 48:4 — every pinning survivor escapes only into denominators just
above the pinning ceiling 41, and none survives past 48.

## Validation

- **Six exact-M referee gates**: {1..13} → 1/14, {1..12,14} → 1/13,
  {1..12,26} → 2/27, {1..11,13,24} → 1/14, {1..11,13,36} → 3/41,
  {1,2,3,5,7,8,9,10,11,12,17,19,104} → 8/105. All exact.
- **Rediscovery gate A** (interval widened to (1/14, 2/27), w0 = 36): the
  pipeline blind-finds exactly {1..11,13,36}, referee-confirmed 3/41 INGAP.
- **Rediscovery gate B** (widened to (1/14, 1/13), w0 = 26): blind-finds
  exactly {1..12,26}, referee-confirmed 2/27 INGAP.
- **Gate C** (LRC-mode regression, interval (0, 1/13), heights 13..26,
  EXHAUSTIVE): finds exactly 5 families — {1..13} and {1..11,13,24} at the
  floor 1/14, and THREE 2/27 realizers: the known {1..12,26} plus the new
  {1,2,3,4,5,6,7,8,9,11,13,20,24} and {1,2,3,4,5,6,7,8,9,11,12,13,20}.
  Zero families below the floor: **the sub-1/13 spectrum at heights ≤ 26 is
  exactly {1/14, 2/27} and LRC(14) holds exhaustively at height ≤ 26.**

## Why this composes (the sandwich; see the S319 synthesis doc)

- **kind-pasteur HYP-7940** (today): micro-gap (1/14, 1/14 + 1/48104) empty to
  max speed 281,577, conditional on the S-T extraction facts — the effective
  bottom plank at n = 14.
- **THIS THEOREM**: the full sub-gap, unconditional, to height 55 — the
  middle plank.
- **THM-1289** (opus-S402, published G-K Thm 1.4, independently re-confirmed
  this session): the floor is isolated from above at ALL heights —
  (1/14, 1/14+δ) is spectrum-empty, δ ineffective — the citation top plank.
  Full finiteness of (1/14, 3/41] is CONDITIONAL on G-K Conjecture 1.5 (the
  proven chain ends in gridmax values; MISTAKE-190) — the roof. In the
  Kravitz/Fan–Sun stratification (repo LRCSpectrumWindow), every candidate
  value in the open sub-gap has stratum k ≥ 3 (4/55, 5/69, 6/83, 7/96, 7/97,
  …), no k ≤ 2 construction can enter, and sub-gap values with D ≥ 15 are
  Fan–Sun-inadmissible (k = s − 13D ≥ 14 in every representation), so under
  Fan–Sun Conj 1.2 the candidate list is finite and explicit (D ≤ 19).

The named missing piece is unchanged but sharpened: an effective height bound
for the finite G-K list (opus's C4), toward which every increment of B here
is a permanent plank.

## Spectroscopy leads (for the rung/CRT programs)

The 50 pinning survivors' witnesses live at q ∈ [43, 48] — immediately above
QPIN = 41. Extending pinning to q ≤ 48 would zero out the survivor count at
B = 55 entirely (predicted; one-line change, worth doing when pushing B ≥ 60).
The 2/27 rung has ≥ 3 realizers of distinct shapes at heights ≤ 26 (new), so
rung-realizer multiplicity grows below the classical ladder — data for
death-star's gate-tower realizability law and boxeph's CRT stack.
