---
id: THM-1055
title: THE APPARENT BONF5 SPEED THRESHOLD IS A SAMPLING ARTIFACT — an explicit counterexample: the stratum run reported 0/4 positive below min-speed 70 and 4/4 positive in [150,300) and [600,900), which reads as a clean threshold; but dilating the primitive failure V = {27,36,46,70,101,114,117,121,140,160,194,277,293} (BONF5 = −0.360134) by k = 33 gives W with min speed 891, inside the "all positive" stratum, blocking all seven moduli, with BONF5(W) = −0.360134 EXACTLY — the invariance of THM-1050 realised as a concrete witness; random sampling from [m,13m] essentially never draws a dilate, because that requires a common factor across all thirteen speeds
status: exhibited explicitly and verified exactly (BONF5 identical to the last digit; W blocks all seven moduli; gcd(W) = 33)
source: opus-2026-07-17-S364 (owner: work the last open regime)
depends_on: [THM-1050 (dilation invariance — the prediction), THM-1045 (the taxonomy whose stratum run produced the apparent threshold), MISTAKE-154]
scripts: 04-computation/threshold_artifact_opus_S364.py -> 05-knowledge/results/threshold_artifact_opus_S364.out
---

# THM-1055 — the threshold that was not there

The S362 stratum run, completed after S363 closed, reported:

| min-speed stratum | BONF5 > 0 | range |
|---|---|---|
| [23, 70) | 0/4 | −0.2803 … −0.0286 |
| [150, 300) | 4/4 | +0.0493 … +0.0832 |
| [600, 900) | 4/4 | +0.1011 … +0.1098 |

Read naively this is a textbook threshold: failures below, successes above,
with the margin growing. It is not. THM-1050 proves no speed threshold can
exist, and here is the witness.

**The counterexample.** The primitive family

> V = {27, 36, 46, 70, 101, 114, 117, 121, 140, 160, 194, 277, 293}

has min speed 27, gcd 1, blocks all seven sieve moduli, and BONF5(V) =
**−0.360134**. Dilating by k = 33:

> W = 33·V = {891, 1188, 1518, 2310, 3333, …}

has min speed **891** — inside the stratum where the run found 4/4
positive — still blocks all seven moduli, and has

> BONF5(W) = **−0.360134**, identical to BONF5(V) to the last digit.

So a family with min speed ≥ 600 and BONF5 < 0 exists, exhibited.

**Why the run saw a threshold anyway.** Random sampling from [m, 13m]
essentially never draws a dilate: that would require a common factor
across all thirteen independently chosen speeds. The dilates form a
vanishingly thin subset of the sampling space while being dense in the
scale axis — present in every stratum, invisible to every random draw.
The "4/4 positive" was a true statement about four particular draws and a
false statement about the stratum.

**The general caution.** A monotone-looking trend across strata is
evidence about the SAMPLER, not about the space, whenever the quantity is
invariant under a group action whose orbits the sampler does not
represent. Here BONF5 is dilation-invariant and the sampler draws almost
exclusively primitive families; the trend measured primitivity, not
speed. This is the concrete form of MISTAKE-154, and it is worth keeping
as the standing example: **the data looked like a theorem and the
theorem said the data was an artifact — the theorem was right.**
