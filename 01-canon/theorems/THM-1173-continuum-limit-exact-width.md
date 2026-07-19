---
id: THM-1173
title: THE CONTINUUM LIMIT, COMPUTED EXACTLY — the bad-run width CONVERGES to 1/21 and the counting argument holds with margin 0.0688. (I) THE LIMIT MODEL, derived: normalising the k₁-gap [j/k₁+1/(14k₁), (j+1)/k₁−1/(14k₁)] to [0,1], a tooth of kᵢ = k₁+dᵢ has half-width k₁/(12kᵢ) → **1/12** and centre (7/6)·frac(−dᵢu) − 1/12 with u = j/k₁, while the threshold 1/(7k₄) normalises to k₁/(6k₄) → **1/6**. So the limit is three teeth of width 1/6 in [0,1] at centres determined by one parameter, and BAD means the longest surviving piece is ≤ 1/6. (II) THE PROFILE IS EXACTLY LINEAR: for u ∈ [0,1/4] the d=3 tooth has centre 13/12 − (7/2)u, hence left edge 1 − (7/2)u, and the free piece [0, 1−(7/2)u] is the longest — so **F_∞(u) = 1 − (7/2)u**, verified as an exact rational identity at all 26 tested points. (III) THE ENDPOINTS ARE EXACT: entry where 1 − (7/2)u = 1/6, i.e. **u = 5/21**, and exit at **u = 2/7 = 6/21**; both give F_∞ = **1/6 exactly**, with the values just inside and just outside straddling the threshold. So the per-run width is **exactly 1/21**. (IV) IT CONVERGES — settling THM-1172(V)'s load-bearing question in the favourable direction. The finite per-run fractions 0.0382, 0.0386, 0.0428, 0.0423, 0.0448, 0.0467 for k₁ = 157…407 rise monotonically toward 1/21 = 0.047619 from BELOW; they do not grow past it. (V) THE COUNTING ARGUMENT HOLDS: two mirror runs give total bad measure **2/21 = 0.0952381**, against S(P) of measure ≥ 0.164 spread over 14–26 components — margin **0.0688**. A fine exact grid confirms 0.0952834 against 2/21, difference 4.5·10⁻⁵ (grid resolution)
status: (I) PROVED — the normalisations are exact limits. (II) PROVED and verified as a rational identity: the d=3 tooth's left edge is 1−(7/2)u and the piece to its left is longest on [0,1/4]. (III) PROVED — both endpoints evaluate to exactly 1/6 in rational arithmetic. (IV) the convergence is established for the continuum value, and the finite values are measured approaching it from below; the monotonicity of that approach is measured, not proved. (V) holds on these constants, but S(P) ≥ 0.164 is the atlas figure for eight-speed cores and the whole argument is for CONSECUTIVE killers (d = 1,2,3). **Uniform r=5 remains OPEN** — non-consecutive frequencies are not covered
source: kind-pasteur-2026-07-18-S128 (cont.75; owner: compute the continuum limit and settle the width growth)
depends_on:
  - THM-1172    # whose (IV)/(V) growth question this settles
  - THM-1162    # the whole-safe-set counting argument being re-costed
  - THM-1142    # the linear descent
script: 04-computation/continuum_limit_kps_S128c75.py, exact_width_kps_S128c75.py (+ .out)
---

# THM-1173 — the continuum limit, exactly

## (I) The limit model

Normalise the k₁-gap to [0,1] by dividing by its length G = 6/(7k₁). For kᵢ = k₁ + dᵢ:

- a tooth of kᵢ has half-width 1/(14kᵢ), normalised **k₁/(12kᵢ) → 1/12** (full width 1/6);
- its centre is ((k₁s − j dᵢ)/(k₁kᵢ) − 1/(14k₁))·7k₁/6 → **(7/6)·frac(−dᵢu) − 1/12**, u = j/k₁;
- the threshold 1/(7k₄) normalises to **k₁/(6k₄) → 1/6**.

So the limit is: three teeth of width 1/6 inside [0,1], positions driven by the single
parameter u, and *bad* means the longest surviving piece is ≤ 1/6.

## (II) The profile is exactly linear

For u ∈ [0,1/4] the d = 3 tooth has centre 13/12 − (7/2)u, so its **left edge is
1 − (7/2)u**, and the piece [0, 1 − (7/2)u] to its left is the longest. Hence

> **F_∞(u) = 1 − (7/2)·u  on [0, 1/4].**

Verified as an exact rational identity at u = 0, 1/100, …, 25/100 — no mismatches. (At
u = 1/4 it gives 1/8, matching the measured minimum.)

## (III) The endpoints, exactly

| | u | F_∞(u) |
|---|---|---|
| entry | **5/21** = 0.2380952 | **1/6** exactly |
| exit | **2/7 = 6/21** = 0.2857143 | **1/6** exactly |

with F_∞ just inside each endpoint below 1/6 and just outside above it. So

> **per-run width = 2/7 − 5/21 = 1/21 = 0.0476190, exactly.**

## (IV) It converges — the growth question is settled

THM-1172(V) left open whether the per-run width converges or grows past 0.082. It converges:

| k₁ | 157 | 207 | 257 | 307 | 357 | 407 | **∞** |
|---|---|---|---|---|---|---|---|
| per-run fraction | 0.0382 | 0.0386 | 0.0428 | 0.0423 | 0.0448 | 0.0467 | **1/21 = 0.047619** |

The finite values rise toward 1/21 **from below** and do not pass it.

## (V) The counting argument holds

Two mirror-symmetric runs (THM-1172) give

> **total bad measure = 2/21 = 0.0952381**,

against S(P) of measure ≥ 0.164 spread over 14–26 components. Margin **0.0688**. A fine
exact grid gives 0.0952834, differing from 2/21 by 4.5·10⁻⁵ — grid resolution.

## Honest status

This settles the question asked, exactly and in the favourable direction. Two limits on the
scope remain and both matter:

- the whole analysis is for **consecutive killers** (d = 1,2,3); non-consecutive frequencies
  give different drift rates and are not covered, and THM-1172 already flagged that the
  reflection symmetry may not survive there;
- S(P) ≥ 0.164 is the atlas figure for eight-speed cores.

**Uniform r=5 remains open.** What is now exact is the consecutive-killer branch of the
argument: bad measure 2/21, safe measure ≥ 0.164, margin 0.0688.

## Named next
- Redo (I)–(III) for general (d₂,d₃,d₄). The drift rates become dᵢ/k₁, so the profile is
  F(u) with three frequencies; the linear branch should persist with slope (7/2)·(d_max/3),
  giving per-run width ∝ 1/(7 d_max) and total bad ∝ 2 d_max/21 — which would *grow* with
  d_max and could exceed 0.164 once d_max ≳ 5. That is the next thing to check, and it may
  be where the consecutive case is genuinely easiest rather than representative.
- If so, the argument splits: small d_max by this counting, large d_max by the spread cone
  of THM-1140(II) (adjacent ratios ≥ 7/3). Whether those two ranges meet is then the
  remaining question.
