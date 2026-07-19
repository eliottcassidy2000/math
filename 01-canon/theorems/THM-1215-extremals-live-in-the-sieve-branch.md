---
id: THM-1215
title: THE EXTREMAL BEHAVIOUR OF LRC(14) LIVES ENTIRELY IN THE D=1 SIEVE BRANCH — the D≥2 branch is not tight — every extremal and near-extremal family has D = 1 with gap exactly 1/q₀: {1,…,13} and {1,…,11,13,24} give 1/14 at pair (1,13), {1,…,12,14} and {1,…,12,15} give 1/13 at (1,12), {1,…,11,13,25} gives 1/12 at (1,11). So the smallest gaps come from pushing q₀ up to 14 within the sieve branch, NOT from the hard stratum. Conversely, adversarial minimisation restricted to the hard stratum (q₀ > 14, where the sieve is too weak) bottoms at g = 2/19 ≈ 0.1053 — a margin factor of **1.474× above 1/14** — with the rescue explicit: at V = {2,8,9,10,11,12,14,18,21,32,44,66,78} the sieve at q₀ = 15 would give only 1/15 < 1/14, while the D = 2 pair (8,11) delivers 2/19
status: the extremal/near-extremal (D,s) values are exact; the hard-stratum minimum 2/19 is from adversarial hill-climbing over 30 restricted starts, which is MEASURED not proved — and hill-climbing has misled me repeatedly (MISTAKE-152/154/156/157), so the 1.474× margin is evidence, not a bound. The hard-stratum constructor forces 8…14 among the speeds, which is one family shape and not the whole stratum
source: opus-2026-07-19-S392 (owner: work the D≥2 branch)
depends_on: [THM-1210 (the D=1 ⟺ sieve identification and the two-mechanism split), THM-1205 (the active-pair ratio), THM-1105 (the hard stratum), boxeph-S120 / HYP-7745]
scripts: 04-computation/dge2_branch_opus_S392.py, dge2_hardmin_opus_S392.py -> 05-knowledge/results/
---

# THM-1215 — where the tightness actually lives

## The near-extremal ladder is entirely D = 1

Reading the known small-gap families in the (D, s) coordinates of THM-1205:

| family | g | D | s | active pair |
|---|---|---|---|---|
| **{1,…,13}** | **1/14** | **1** | **14** | **(1,13)** |
| **{1,…,11,13,24}** | **1/14** | **1** | **14** | **(1,13)** |
| {1,…,12,14} | 1/13 | 1 | 13 | (1,12) |
| {1,…,12,15} | 1/13 | 1 | 13 | (1,12) |
| {1,…,11,13,25} | 1/12 | 1 | 12 | (1,11) |
| {2,…,14} | 1/8 | 2 | 16 | (2,14) |

Every extremal and near-extremal is **D = 1**, with gap exactly **1/q₀** — the
classical sieve (THM-1210). The ladder toward tightness is *raising q₀ toward
14 inside the sieve branch*, and 14 is where it stops, because q₀ = 15 would
give 1/15 < 1/14 and the sieve would no longer suffice.

That corrects an expectation I had going in. I assumed D ≥ 2 — the hard
stratum — was the dangerous branch, since that is where the sieve fails. It
is the opposite: D ≥ 2 is where the **rescue** happens.

## The D ≥ 2 branch has a margin

Adversarial minimisation of g restricted to q₀ > 14 (so no small-modulus
sieve exists) bottoms at

> **g = 2/19 ≈ 0.10526**, a factor **1.474×** above 1/14

at V = {2, 8, 9, 10, 11, 12, 14, 18, 21, 32, 44, 66, 78}, with q₀ = 15,
active pair **(8,11)**, D = 2, s = 19, ratio s/D = 9.5 ≤ 14.

The rescue is explicit there: the sieve at q₀ = 15 would deliver only
1/15 ≈ 0.0667, below the threshold — and the D = 2 pair supplies 2/19
instead. That is the two-mechanism picture of THM-1210 caught in the act.

## What this reduces LRC(14) to

Combining with THM-1210:

1. **q₀ ≤ 14** — the sieve (D = 1) gives g = 1/q₀ ≥ 1/14. **Closed**, and it
   contains every extremal family, at q₀ = 14 exactly.
2. **q₀ > 14** — the sieve is too weak; a D ≥ 2 pair must supply ≥ 1/14.
   Measured margin 1.474× at the adversarial minimum.

So the whole conjecture rests on branch 2, and the target is sharp: **show
that when q₀ > 14, some pair achieves D/(vᵢ+vⱼ) ≥ 1/14.** The measurements
say this holds with room to spare, which is a materially better position
than the knife-edge that branch 1 sits on — the difficulty is no longer at
the extremal families, which are now fully accounted for.

## Status

Structural claims exact; the margin measured. I want to be explicit that
1.474× comes from hill-climbing over 30 restricted starts, and that method
has misled me four times in this programme. It is evidence that branch 2 is
not tight, not a proof that it cannot be.
