---
id: HYP-3791
title: THE MULTI-FAR CORRECTION IS SELF-SIMILAR OVER THE 13-LATTICE (= n-1 = missing speed = first CF convergent of t*), AND RESONANCE => REDUNDANCY (positive correction) => no multi-far beater. Working the next lever (r-far correction) recursively/fractally with 2 as the first prime. RESULTS: (1) the resonance lattice is 13ℤ (13=n-1=the DROPPED speed = 1/first-convergent of t*=n/Phi6=14/183=[0;13,14]); it governs resonance at EVERY order r self-similarly: single-far peaks at k∈13ℤ (S66), pairwise 2-far corr2(a,b) peaks iff the comb DIFFERENCE δ=|a-b|∈13ℤ (δ=13:+0.099, 26:+0.081, 39:+0.067; non-13 δ all |.|<=0.017), r-far via signed speed-combinations ∈13ℤ. (2) the 13-spaced pair resonance does NOT decay in scale W (corr2(W,W+13)~0.08-0.11 for W=200..50000) -- because δ=13 is itself a SLOW resonant speed (13t* ~ integer) that couples to L_C at ANY scale = the OPEN-Q-108 residual locus, pinned to 13-lattice-spaced far combs. (3) BUT every strong resonance is POSITIVE = REDUNDANT (the two dangers LOCK: (a+13)t=at+13t, 13t slow ~const near t*, so D_{a+13}~D_a => overlap => cover the SAME part of L_C), while coverage-SPREADING (negative) correlations are all WEAK (<=0.016 vs +0.099). Far combs either resonate REDUNDANTLY (13-spaced, pile up) or EQUIDISTRIBUTE (non-resonant, cover 2r' each, decay) -- neither SPREADS coverage over L_C => no multi-far beater (recursive extension of S66/HYP-3787 far-element impotence). PRIME LENS (2 first): prime 2 = the ± antipode ι = the ODD sin kernel s(j)=sin(2πj r')/(πj) = why the correction is SIGNED; prime 7=n/2 = the 7-vanishing (THM-503), gates hat1 (8x suppression at k=14=2·7=n); resonance denominators 13=n-1 & Phi6=183=3·61 = the CF convergents (multi-scale/fractal). REFUTED: survival-correction sign is NOT (-1)^r (data +,+,+,+,-); the prime-2 parity lives in the KERNEL, not the survival sign
status: MIXED (grid-verified structure + heuristic mechanism). VERIFIED (grid N=6e5, n=14, core={1..12}, r=r'=1/14, L=0.0341): (1) corr2 resonance iff δ∈13ℤ (δ=13:+0.0993, 26:+0.0814, 39:+0.067; sign census top-4 positive all 13|δ; top-4 negative all |.|<=0.0161); (2) corr2(W,W+13) persists ~0.08-0.11 across W=200,1000,5000,20000,50000 (no scale decay) while corr1(W) and corr2(W,W+1) decay/stay small; (3) hat1 8x-suppressed at k=14=n (7-gating); two-atom hat1(k)~L cos(2πk t*) with envelope crossing zero by k=183 (two-scale). MECHANISM (resonance=>redundancy=>no spreading=>no beater) is HEURISTIC, not a proof. The persistent 13-resonance is a REAL residual (OPEN-Q-108), here characterized (redundant sign) not eliminated. Grid-approximate at large W (resolution edge ~N/20).
source: klein-2026-07-01-S67
depends_on:
  - HYP-3790   # S66: single-far signed correction, exact Fourier identity + O(1/w) (renumbered from 3787; = mac-mini's HYP-3787); this is the r-far recursion
  - HYP-3786   # S65: far-element impotent / equidistribution on L_C
related:
  - THM-503    # the 7-vanishing (prime 7=n/2 gating); L is a singular INTEGRAL, NOT an Euler product (no prime-multiplicativity)
  - THM-515    # L_C measure = singular series theta-form (the on-diagonal; corr = off-diagonal)
  - HYP-3762   # three-gap / continued fraction (the CF [0;13,14] source of the 13-lattice)
  - HYP-3788   # mac-mini: |H|>=2 multi-patch equidistribution (this adds the resonance/redundancy structure)
  - OPEN-Q-108 # the r>=2 multi-far uniform bound (13-spaced residual pinned here, not closed)
  - HYP-3132   # kind-pasteur multi-far crux; kps-S4 residual "near-tight cores x incommensurate combs" (complementary: this is the COMMENSURATE 13-spaced case)
  - HYP-3715   # t* = n/Phi6 hexagonal point
results:
  - 04-computation/multifar_recursion_fractal_prime_klein.py
  - 05-knowledge/results/multifar_recursion_fractal_prime_klein.out
---

# HYP-3791 [renumbered from HYP-3789 by mac-mini-S77: mac-mini-S76 committed HYP-3789 (moment relaxation) 13min earlier; kps-S6 consensus] — the multi-far correction is self-similar over the 13-lattice; resonance ⇒ redundancy

## The next lever, recursively: the r-far correction
S66 (HYP-3790, = mac-mini's independently-derived HYP-3787) wrote the single-far correction as an exact
Fourier sum on `hat1 = FT(1_{L_C})`. The r-far
correction is the r-fold version (a cumulant/convolution of single-far). Attacked here with core
`{1..n-2}`, `r=r'=1/n` (so `2r'=1/7`, the prime 7), `L = meas(L_C) = 0.0341`.

## Result 1 — the resonance lattice is `13ℤ`, self-similar at every order (the fractal)
`t* = n/Phi6 = 14/183 = [0;13,14]` (continued fraction). Its convergents are `1/13` and `14/183`. The
denominator **`13 = n-1`** — the speed the construction DROPS — is `1/first-convergent ≈ 1/t*`. This single
lattice governs resonance at **every** order `r`, self-similarly:
- **single-far** (S66): `hat1(k)` peaks at `k ∈ 13ℤ` (k=13,26,39,52,65).
- **pairwise 2-far**: `corr2(a,b) = (1/L)∫_{L_C}(D_a-2r')(D_b-2r')` peaks **iff the comb difference
  `δ=|a-b| ∈ 13ℤ`**: `δ=13 → +0.099`, `26 → +0.081`, `39 → +0.067`; every non-13 `δ` has `|corr2| ≤ 0.017`.
- **r-far**: resonance when a signed speed-combination `Σ c_i w_i ∈ 13ℤ`.
So the resonance set is a union of `13ℤ`-hyperplanes in speed-space, the SAME at every scale — the
recursive/fractal skeleton. (The second convergent `183 = Phi6 = 3·61` sets the envelope: `hat1(k)` combs
decay `~0.96, 0.94, 0.92,…` and the envelope reverses sign by `k=183` — a two-scale CF structure.)

## Result 2 — the 13-resonance PERSISTS to all scales (the OPEN-Q-108 residual, pinned)
`corr2(W, W+13)` does **not** decay in `W`: `≈ 0.083, 0.099, 0.103, 0.113, 0.074` for
`W = 200, 1000, 5000, 20000, 50000`, while `corr1(W)` and non-resonant `corr2(W,W+1)` decay/stay small.
Reason: the difference `13` is a **slow** speed (`13t ≈ integer` at `t*`), so a `13`-spaced far pair stays
phase-locked to `L_C` at **any** scale `W`. This is why single-far impotence (S66, `O(1/w)`) does **not**
trivially extend: far-far combs CAN resonate at arbitrary scale, provided they are `13`-lattice-spaced.
This is exactly the OPEN-Q-108 multi-far residual, here localized to `13ℤ`-spaced far combs (the
COMMENSURATE case; complementary to kind-pasteur-S4's incommensurate-combs residual).

## Result 3 — but resonance ⇒ REDUNDANCY (wrong sign for a beater): the impossibility mechanism
Every strong pairwise resonance is **positive** = the dangers overlap = the combs cover the **same** part
of `L_C` (redundant), NOT different parts. Structural reason: `(a+13)t = at + 13t`, and `13t` is nearly
constant across the fast `at`-cycle near `t*`, so `D_{a+13} ≈ D_a` — the two dangers **lock**. Sign
census (`δ=1..40`): the 4 most-positive `corr2` are `δ=13,26,39` (all `13ℤ`) `+ 28`; the 4 most-negative
are `δ=35,6,36,22`, all `|corr2| ≤ 0.016` (**6× weaker** than the `+0.099` redundant peak).
> To **beat**, far combs must SPREAD their danger over `L_C` (cover different parts) = strong **negative**
> correlation. But strong correlations occur only on `13ℤ` and are **positive/redundant**; the
> coverage-spreading (negative) correlations are all weak. So far combs either pile up redundantly
> (`13`-spaced) or equidistribute (non-resonant, cover `2r'` each) — **neither spreads coverage** ⇒ no
> multi-far beater. (Recursive extension of S66: single-far can't cover `L_C` by narrow arcs; multi-far
> can't cover it because its strong correlations are redundant.)

## The prime lens (2 as the first prime; the others echo it one scale finer)
Canon guard (THM-503): `L` is a singular **integral**, **NOT** an Euler product — so there is **no**
prime-multiplicative factorization of `L`. The primes enter as **structural gradings**, not factors:
- **Prime 2** (first): the **antipode** `ι` (`x ↔ 1-x`, `a ↔ -a`) built into `‖·‖` (the `±`, the two-sided
  danger band `2r'`). It makes the danger kernel `s(j)=sin(2πj r')/(πj)` **odd** in `j` — this ι-odd `ℤ/2`
  grading is *why the correction is signed*. `hat1` (FT of `1_{L_C}`) is real & **even** (`cos`) = the
  ι-even partner. Prime 2 = the sign template.
- **Prime 7** = `n/2`: the **7-vanishing** (THM-503, `s(t)=0` at `7|t`; sign period `14=2·7`). It **gates**
  `hat1`: an `8×` suppression at `k = 14 = 2·7 = n` (verified). Prime 7 = the zeros.
- **`13 = n-1`, and `Phi6 = 183 = 3·61`**: the **resonance denominators** (CF convergents of `t*`). Primes
  `3, 61` refine the resonance one scale finer (the envelope/second scale).
So the first prime `2` sets the `ℤ/2` sign-grading; each other prime plays the analogous role one level
finer (7 → zeros; 3, 61, 13 → resonances). The user's "consider how the other primes relate similarly"
is: `2` grades by antipode/sign, `7` by vanishing, `{3,13,61}` by resonance scale — a prime-indexed
hierarchy, self-similar over the CF.

## Refuted
The survival-correction sign is **not** `(-1)^r` (data: `survival_r - (6/7)^r = +,+,+,+,-` for `r=1..5`).
The prime-2 parity is a property of the **kernel** (`s` odd), not of the survival statistic.

## Net
The multi-far correction lives on a single self-similar skeleton — the `13ℤ`-lattice (`n-1` = the dropped
speed = first CF convergent of the hexagonal point `t*`) — at every order `r`. The `13`-spaced resonance
persists to all scales (the OPEN-Q-108 residual), but it is **redundant** (positive: dangers lock and
double-cover), while the coverage-spreading correlations a beater would need are weak — so no multi-far
beater emerges. The primes grade the structure (2 = sign, 7 = zeros, 3·61·13 = resonance), with `2` the
template. Heuristic mechanism + grid-verified structure; the uniform `r≥2` bound is still OPEN-Q-108.
