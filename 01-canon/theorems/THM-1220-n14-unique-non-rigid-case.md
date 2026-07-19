---
id: THM-1220
title: n = 14 IS THE UNIQUE NON-RIGID CASE — AP UNIQUENESS HOLDS AT n = 10,11,12,13,15,16,17,18 AND FAILS ONLY AT 14 — searching single substitutions of {1,…,n−1} for families attaining the floor 1/n exactly, the AP is the unique extremal at every n tested EXCEPT n = 14, where 12 → 24 ties it ({1,…,11,13,24}, THM-1120). At n = 13 the same substitution gives M = 2/25 > 1/13 — it MISSES the floor and is instead klein's stability minimiser (THM-1004). So the rigidity that klein proved at n = 13 (M = 1/13 or M ≥ 2/25, equality only at the AP) has no analogue at n = 14: there is no "gap above a unique extremal" because the extremal set is not a single point. Separately, 26 adversarial runs found NO family with M strictly in (1/14, 1/13), smallest M observed 6/61 ≈ 0.0984
status: the uniqueness sweep is exhaustive over SINGLE substitutions with replacement r ≤ 3n at n = 10…18 (r ≤ 60 for n ≤ 15) — not over multi-element perturbations, so "unique" is at Hamming radius 1. The exceptional tie at n = 14 is exact. The absence of families in (1/14, 1/13) is adversarial-search evidence, not proof
source: opus-2026-07-19-S393 (owner: branch 2 sharper search; n=12 AP uniqueness; mine 3/4/(1/12) threads)
depends_on: [THM-1004/1005 (klein-S313d/e — Hamming-1 and -2 rigidity at n=13, prior art), THM-1002 (klein — pair-sum denominator bound), THM-1120 (the second tight family at n=14), THM-1215 (the branch-2 target), THM-1125 (the essential-region criterion explaining the swap)]
scripts: 04-computation/stability_gap_opus_S393.py -> 05-knowledge/results/stability_gap_opus_S393.out
---

# THM-1220 — why 14 is the odd one out

## Prior art mined

klein-S313 established the n = 13 picture, which is the exact analogue of
what this programme has been doing at n = 14:

- **THM-1002** — the pair-sum denominator bound: M(A) = val/q with q | vᵢ+vⱼ.
  (Same territory as boxeph-S120's located maximizer and my THM-1205.)
- **THM-1004** — **Hamming-1 rigidity**: replacing any single element of
  {1,…,12} forces M ≥ 2/25, equality **exactly at {1,…,11,24}**.
- **THM-1005** — **Hamming-2 rigidity**: same bound at radius 2. Radius ≥ 3 open.

So at n = 13 there is a genuine **stability gap**: M = 1/13 (the AP) or
M ≥ 2/25, nothing between — and the AP is the unique extremal.

## The parallel that breaks

The substitution 12 → 24 appears at both n:

| n | AP | M | 12 → 24 family | M | verdict |
|---|---|---|---|---|---|
| 13 | {1,…,12} | 1/13 = floor | {1,…,11,24} | **2/25 > 1/13** | misses the floor |
| 14 | {1,…,13} | 1/14 = floor | {1,…,11,13,24} | **1/14 = floor** | **ties** |

At n = 13 the swap is the *stability minimiser*; at n = 14 it is a *second
extremal*. Same substitution, different status.

## The uniqueness sweep

Searching all single substitutions of {1,…,n−1} for families attaining 1/n
exactly (r ≤ 3n):

| n | 10 | 11 | 12 | 13 | **14** | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|
| ties | 0 | 0 | 0 | 0 | **1 (12→24)** | 0 | 0 | 0 | 0 |
| swappable speeds | — | — | none | none | **{12}** | none | none | none | none |

**n = 14 is the only non-rigid case in the range.** It is not a parity
effect — 10, 12, 16, 18 are even and rigid.

## Why this matters

The rigidity/stability-gap route is exactly how n = 13 was handled: prove
the AP is the unique extremal, then prove a definite gap above it
(klein's 2/25 at radius ≤ 2). **That route cannot transfer to n = 14**,
because the hypothesis it rests on is false there — {1,…,11,13,24} sits on
the floor alongside the AP, so there is no gap above a unique extremal to
establish.

That is a concrete structural explanation for why 14 is the open case,
rather than an accident of difficulty: the standard rigidity argument has
its hypothesis fail at exactly this n.

## Branch 2, sharpened

Adversarial search minimising distance into the open interval (1/14, 1/13)
found **0 families inside** across 26 runs, smallest M observed 6/61 ≈
0.0984. So the near-floor ladder remains 1/14 (twice), then 1/13, with
nothing between — consistent with an n = 14 stability gap of the klein type
*above the two-element extremal set* rather than above a single point.
Evidence, not proof.
