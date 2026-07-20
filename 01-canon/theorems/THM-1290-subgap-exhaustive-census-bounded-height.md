---
id: THM-1290
title: "[CLAIM STUB — census running] THE (1/14, 3/41) SUB-GAP IS EXHAUSTIVELY EMPTY AT BOUNDED HEIGHT — every 13-subset of [1,B] (primitive or not) has M ∉ (1/14, 3/41), by complete enumeration with rigorous necessary-condition pruning: covering 2..13 (else M ≥ 1/q ≥ 1/13), spread w_max > 12·w_min (THM-1043), an active pair summing into the sparse admissible set {s : ∃D ∈ ℤ, s/14 < D < 3s/41} = {55, 69, 83, 96, 97, 110, 111, ...} (THM-1269's D = M·s), unit-direction pinning at depth 1 for all q ≤ 27 and depth 2 for 28 ≤ q ≤ 41 (spread-lemma template), survivors certified out by rational witness (margin ≥ 3/41) or exact pinch-set M. Target B = 36 first, push toward 55."
status: CLAIMED (klein-2026-07-19-S319), census in flight. What is KNOWN now — the necessary-condition filters are all rigorous consequences of canon (THM-1043, THM-1268/1269, the spread/pinning lemma template of boxeph-S115/S126 and mac-mini S54); non-primitive families reduce to primitive ones at height ≤ B/2 < 28 ≤ (forced v_max), so primitive enumeration suffices for B ≤ 55; prior emptiness evidence is SEARCH only (opus ~12,400 structured + 8.5M random, klein 56k). What NEEDS EVIDENCE — the enumeration itself (script + frozen out + referee re-verification of all HARD survivors and a sample of witness-certified ones + rediscovery gates on {1..12,26} → 2/27 and {1..11,13,36} → 3/41). Nothing below is citable until the run and referee complete.
source: klein-2026-07-19-S319 (owner: work to finish LRC(14); see the near-misses, synthesize perspectives)
depends_on: [THM-1268 (D ≥ 4 forcing, mediant 4/55), THM-1269 (D = M·s, the pair-sum reduction), THM-1043 (spread rung), THM-1230/1235 (the gap edge 3/41 and its witness), boxeph-S115/S126 (mod-p spread/pinning template), mac-mini-S54 (the census harness architecture this adapts, lrc_gap_census40_S54.c)]
scripts: 04-computation/lrc14_subgap_census_klein_S319.c + 04-computation/lrc14_subgap_referee_klein_S319.py -> 05-knowledge/results/lrc14_subgap_census_klein_S319.out
---

# THM-1290 — the sub-gap census at bounded height (claim stub)

## Why this census is new

THM-1268/1269 record that all existing emptiness evidence for (1/14, 3/41) is
search: ~12,400 structured families (opus-S397), 8.5M random families
(opus-S398, "scale, not coverage" — it missed both known near-floor families),
56k mixed (klein). No exhaustive bounded-height statement exists for THIS
interval. The n=12 shadow gap (1/13, 2/25) has one (mac-mini S54, height 48);
this file transfers that architecture to the n=14 sub-gap.

## Why bounded height composes forward (the sandwich)

opus-S401's HYP-7930 (G-K accumulation revival, citation-pending) retypes
(1/14, 3/41] as a FINITE list of attained values at ALL heights, ineffectively.
The named missing piece is an effective height bound near the floor. An
exhaustive census at height B is exactly the statement that composes with any
future effective bound H₀ ≤ B: cage-style effectivization from below +
accumulation finiteness from above + this census in the middle. It is also the
n=14 twin of HYP-7920's height-258,276 rigidity at n=13.

## The a-priori structure (rigorous, from canon)

A 13-set W ⊆ [1,B] with M(W) ∈ (1/14, 3/41) must satisfy:

1. **Covering 2..13**: no multiple of q ⟹ witness t = 1/q ⟹ M ≥ 1/q ≥ 1/13.
2. **Spread** (THM-1043): σ ≤ 12 ⟹ M ≥ 1/13; so w_max > 12·w_min.
3. **Active pair sum** (THM-1269): at any maximizer, the rising/falling pair
   gives D = M·s with s = v_i + v_j, D ∈ ℤ₊; M ∈ (1/14, 3/41) forces
   s ∈ {55, 69, 83, 96, 97, 110, 111, 124, 125, ...} (integers whose interval
   (s/14, 3s/41) contains an integer). In particular v_max ≥ 28, and at
   B = 36 the family must contain a pair summing to 55 or 69 exactly.
4. **Pinning** (spread-lemma template): for every q with 3q/41 < d+1, i.e.
   depth d = 1 for 14 ≤ q ≤ 27 and d = 2 for 28 ≤ q ≤ 41: either q | some v,
   or for every unit a mod q some v has dist(v·a mod q) ≤ d.

Non-primitive W = g·W′ have M(W) = M(W′) and w′_max ≥ 28 (by 3), so
g·28 ≤ B < 56 forces g = 1: primitive enumeration suffices for B ≤ 55.

## Deliverables when the run completes

- Exhaustive emptiness of (1/14, 3/41) at height B (this theorem).
- By the same run (the filters are implied by M < 3/41 ⊃ M ≤ 1/14): every
  enumerated 13-subset of [1,B] has exact M ∉ (1/14, 3/41) — any family with
  M < 1/14 would surface at the exact-M stage and be a counterexample to
  LRC(14) itself; none is expected.
- The HARD-survivor spectroscopy (which q certifies, residue anatomy) as data
  for the rung-realizability program (death-star tower) and the mod-p CRT
  stack (boxeph).
