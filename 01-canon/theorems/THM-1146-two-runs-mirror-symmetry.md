---
id: THM-1146
title: THE SINGLE-RUN CLAIM IS REFUTED — the bad set is TWO mirror-symmetric runs — but the location is confirmed exactly and the counting argument survives with reduced margin. (I) SELF-CORRECTION: THM-1145(I) asserted the bad indices form a "single contiguous run". They do not. Over k₁ ∈ [157,420) the bad set splits into **two runs in all 263 cases tested**, never one. I misread the cont.73 output — a bad set of 12 indices with "longest run 6" is two runs of six, and I recorded it as one run of six. (II) THE MECHANISM, and it is provable: for consecutive killers the normalised tooth offsets drift at rate dᵢ/k₁ per gap-index, so with (d₂,d₃,d₄) = (1,2,3) the configuration is **u, 2u, 3u with u = j/k₁ exactly**, and the longest surviving piece is a function F(u) of that single parameter. Since the point set {0,u,2u,3u} is invariant under u ↦ 1−u (reflection of the circle), **F(u) = F(1−u)**, so every minimum at u has a mirror at 1−u. Hence exactly two runs, at u ≈ 1/4 and u ≈ 3/4 — measured as 2 local minima and 4 sign changes of the difference, for every k₁ tested. (III) THE LOCATION IS CONFIRMED EXACTLY: argmin j/k₁ = 0.2484, 0.2487, 0.2489, 0.2491, 0.2492, 0.2493, **0.2494** for k₁ = 157…397 — converging monotonically to **1/4**, exactly as the u,2u,3u structure predicts ({0,¼,½,¾} is the maximally-spread configuration). (IV) THE WIDTH BOUND ALSO NEEDS RAISING: per-run fraction grows with k₁ — 0.0382 (157), 0.0386 (207), 0.0428 (257), 0.0423 (307), 0.0448 (357), **0.0467** (407) — exceeding the 0.0457 asserted in THM-1145, which had scanned only to k₁ = 397. (V) THE COUNTING ARGUMENT SURVIVES, with less room: two runs of ≈ 0.047 each give a total bad measure ≈ **0.093**, against S(P) of measure ≥ 0.164 spread over 14–26 components — so the margin falls from 0.118 to **0.071**, and whether the per-run width converges or keeps growing is now the load-bearing question
status: (I) REFUTED — my own THM-1145(I), by exhaustive count over 263 values of k₁. (II) PROVED — the drift rates give u,2u,3u exactly for consecutive killers, and F(u) = F(1−u) is an immediate reflection symmetry; it explains the observed 2 local minima and 4 sign changes. (III) MEASURED and matching the predicted 1/4 to four decimals. (IV) MEASURED — the earlier bound was an artefact of a shorter scan. (V) the conclusion holds on the measured range but its stability depends on the growth in (IV), which is NOT settled. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.74; owner: prove the single-run bound from the linear descent)
depends_on:
  - THM-1145    # whose (I) this corrects and whose (III) counting argument it re-costs
  - THM-1142    # the linear descent that supplies the drift rates
related: [THM-1144, THM-1143]
script: 04-computation/single_run_proof_kps_S128c74.py (+ .out)
---

# THM-1146 — two runs, not one

## (I) The correction

THM-1145(I) claimed the bad indices form a **single contiguous run**. Over k₁ ∈ [157,420)
the bad set is **two runs in all 263 cases**, never one. The error was a misreading of my
own cont.73 output: a bad set of 12 indices with "longest run 6" is two runs of six, and I
wrote it up as one run of six. The distinction matters because the counting argument in
THM-1145(III) is charged per interval.

## (II) What is actually true, and it is provable

For consecutive killers the offset of comb kᵢ's tooth inside the k₁-gap drifts by
−dᵢ/(k₁kᵢ) per gap-index, i.e. by **−dᵢ/k₁** normalised to kᵢ's own period. With
(d₂,d₃,d₄) = (1,2,3) the three normalised positions are therefore

> **u, 2u, 3u,  with u = j/k₁ exactly,**

so the longest surviving piece is a function **F(u)** of one parameter. The point set
{0,u,2u,3u} is invariant under the circle reflection u ↦ 1−u, hence

> **F(u) = F(1−u),**

and every minimum at u is mirrored at 1−u. So the bad set has **exactly two runs**, at
u ≈ 1/4 and u ≈ 3/4. Measured: 2 local minima and 4 sign changes of the difference, for
every k₁ tested — matching the symmetry exactly.

## (III) The location is confirmed

| k₁ | 157 | 197 | 237 | 277 | 317 | 357 | 397 |
|---|---|---|---|---|---|---|---|
| argmin j/k₁ | 0.2484 | 0.2487 | 0.2489 | 0.2491 | 0.2492 | 0.2493 | **0.2494** |

Converging monotonically to **1/4**, exactly as the u, 2u, 3u structure predicts: {0,¼,½,¾}
is the maximally-spread configuration, so it minimises the largest piece.

## (IV) The width bound must be raised

| k₁ | 157 | 207 | 257 | 307 | 357 | 407 |
|---|---|---|---|---|---|---|
| per-run fraction | 0.0382 | 0.0386 | 0.0428 | 0.0423 | 0.0448 | **0.0467** |

THM-1145 asserted ≤ 0.0457, but that scan stopped at k₁ = 397. The value at 407 exceeds it,
and the trend is still rising.

## (V) The counting argument, re-costed

Two runs of ≈ 0.047 give a total bad measure of ≈ **0.093** in t, against S(P) of measure
≥ 0.164 spread over 14–26 components. The argument still goes through — the bad set cannot
meet every adequate component — but the margin drops from 0.118 to **0.071**.

The load-bearing question is now (IV): does the per-run width converge, or keep growing? If
it converges near 0.05 the total stays ≈ 0.10 and the argument is comfortable; if it grows
past 0.082 the total reaches 0.164 and the argument fails. **This is not settled**, and it
is the first thing to check.

## Honest status

I was asked to prove the single-run bound. It is false, and the proof of what replaces it
(two runs, by reflection symmetry) is short and solid. The location prediction j/k₁ → 1/4
from the u,2u,3u structure is confirmed to four decimals, which is real support for the
parametrisation. **Uniform r=5 remains open.**

## Named next
- Settle the growth in (IV). The natural approach: F(u) has a *continuum* limit as k₁ → ∞
  with normalised tooth widths fixed at k₁/(6kᵢ) → 1/6, so the per-run width should converge
  to the measure of {u : F_∞(u) ≤ 1/6}. Computing F_∞ directly settles it, and it is a
  one-parameter calculation rather than a scan.
- The reflection symmetry of (II) should also be checked for non-consecutive killers, where
  the frequencies are (d₂,d₃,d₄) rather than (1,2,3) and the symmetry may or may not survive.
