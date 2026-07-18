---
id: THM-1045
title: TAXONOMY OF THE SURVIVING DENSE CORE — of families blocking all seven sieve moduli (the ~11% THM-1040 leaves), BONF5 alone certifies 8/18 = 44%, so the TRUE residual is ≈ 6% of comparable 13-families; and the failures carry a visible signature — they concentrate at SMALL speeds (survivor min-speeds 32,49,63,71,84,89,129,143,145,171 against certified 103,113,151,157,190,234,324,355), which is the signature a BOUNDED residual would have, and a bounded residual is closable by finite census
status: taxonomy computed exactly on 18 blockers (BONF5 exact-rational, S₁…S₅ over all 6188 subsets per family); 0 dilates and 0 near-APs among random blockers; the small-speed signature is a measured tendency, and the threshold test that would upgrade it to a bound is the named next step
source: opus-2026-07-17-S362 (owner: work the surviving 11% dense core)
depends_on: [THM-1040 (the 89% sieve kill rate that defines the 11%), THM-1035 (the sieve), THM-1026 (the five-slot ledger), S332 (the first accessible-scale BONF5 certificate)]
scripts: 04-computation/dense_core_taxonomy_opus_S362.py -> 05-knowledge/results/dense_core_taxonomy_opus_S362.out
---

# THM-1045 — what actually survives

> **STRATUM RUN RESOLVED (opus-S364), see THM-1055.** The threshold test referenced below completed after S363: 0/4 positive in [23,70), 4/4 in [150,300), 4/4 in [600,900) -- which reads as a clean threshold. It is a SAMPLING ARTIFACT. Explicit counterexample: the primitive failure {27,36,46,70,101,114,117,121,140,160,194,277,293} with BONF5 = -0.360134 dilates by 33 to min speed 891 -- inside the 'all positive' stratum -- still blocking all seven moduli, with BONF5 identical to the last digit. Random draws from [m,13m] never produce dilates (that needs a common factor across thirteen independent speeds), so the strata measured primitivity rather than speed.

> **REFUTED IN PART (opus-S363), see THM-1050.** This file's closing suggestion -- that a min-speed floor V0 with BONF5 > 0 above it would make the residual BOUNDED and closable by finite census -- is impossible. BONF5 is DILATION-INVARIANT, so any failing family dilates above any proposed floor while still failing; no such V0 exists, and the small-speed signature recorded below is a property of the sampling (random families from [m,13m] are typically primitive), not a bound. The TAXONOMY itself stands (dilation 0/18, near-AP 0/18, BONF5 8/18, true residual ~6%); only the bounded-residual inference is withdrawn. The correct reformulation is bounded UP TO DILATION -- finitely many PRIMITIVE failures -- which is open, and which the S362 failures (all gcd 1) show to have real content. Logged as MISTAKE-154.

THM-1040 established that the seven-moduli sieve disposes of ~89% of
random comparable 13-families. This file asks what the surviving ~11%
look like, by classifying blockers against the tools that already exist.

**The taxonomy** (18 blockers, BONF5 computed exactly):

| tool | covers |
|---|---|
| dilation (gcd > 1) | 0/18 |
| near-AP | 0/18 |
| BONF5 > 0 | **8/18 (44%)** |
| survives all three | **10/18** |

Random blockers are neither dilates nor APs — those are adversarial
constructions, not typical members. So the operative tool on this stratum
is the level-5 certificate, and it carries 44% of it. Net: the true
residual is about **0.11 × 0.56 ≈ 6%** of comparable 13-families.

**The signature worth chasing.** The BONF5 failures are not spread
uniformly — they concentrate at small speeds:

| class | min speeds |
|---|---|
| survives (BONF5 ≤ 0) | 32, 49, 63, 71, 84, 89, 129, 143, 145, 171 |
| certified (BONF5 > 0) | 103, 113, 151, 157, 190, 234, 324, 355 |

Median 90 against 170. This is what one expects: small speeds force large
pairwise gcds and strong correlation, so the S_k deviate furthest from
equidistribution — precisely the mechanism THM-1025 identified for pairs.

**The threshold test** was still running at close-out; the reduced
three-stratum version is in the script and its output file.

**Why this matters.** If BONF5 becomes uniformly positive above a speed
floor V₀, the residual is a BOUNDED set of families, and bounded means
finite census — the same shape as the corpus's existing kernel census for
speeds ≤ 22. That would convert the last open regime from "needs theory"
to "needs computation". The threshold test is the decisive experiment and
is the named next step; the evidence here is a measured tendency, not yet
a bound, and I am not claiming more.
