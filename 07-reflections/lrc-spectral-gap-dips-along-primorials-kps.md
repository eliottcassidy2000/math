---
source: kind-pasteur-2026-06-19-S9
status: reflection on a verified-computational result (THM-539, HYP-2623)
tags: [lonely-runner, spectral-gap, stern-brocot, primorial, highly-composite, basin-width, lrc14]
---

# The LRC spectral basin is arithmetic, not analytic

**Prompt (user):** the second LRC-spectrum point is `M_k={1,…,k-1,2k} -> 2/(2k+1)`, so
`1/(k+1) < sigma_2(k) <= 2/(2k+1)`, mediant-tight for `k<=5`, strict for `k>=6`; the gap is
`Theta(1/k^2)` at most — does anything dip below it unboundedly?

## The answer, and the surprise

**Yes — but the dip is arithmetic.** The values of the max-min collar just above the floor
`1/(k+1)` are the Stern-Brocot mediants `a/(a(k+1)-1)`, `a=2,3,4,...`, marching down to the
floor. The gap `g(k)·k^2 -> 1/a`. So the whole question is: *how large can the level `a` be?*

The answer is governed not by `k` but by the **factorization of `k-1`**. The witnessing family
`F(k,a) = {1,…,k-2, k, a(k-1)}` reaches level `a` exactly when `k-1` carries the first `a-1`
primes:

```
   a = 3  at  k-1 divisible by 2·3   (k = 7, 13, 19, …)
   a = 4  at  k-1 = 30 = 2·3·5       (k = 31)
   a = 5  at  k-1 = 210 = 2·3·5·7    (k = 211)
   a = 6  at  k-1 = 2310 = 2·3·5·7·11 (k = 2311)
```

So `a_max(k) ~ ω(k-1)` is **unbounded** (along primorial `k-1`), but grows only like
`log k / log log k`. The spectral gap dips below `Theta(1/k^2)` infinitely often, yet for almost
every `k` it sits exactly at `~1/(2k^2)` (the doubled apex). **`liminf_k g(k)·k^2 = 0`, but
`limsup = 1/2`.**

## Why a single large speed does so much: it is a clock-killer

The mechanism is clean. The lone large speed `a(k-1)` is divisible by **every** divisor `d` of
`k-1`. At every coarse time `t = j/d` it lands exactly on an integer — distance zero — so it
*annihilates that clock*: at `t=j/d` the collar is `0`, that time can never be the max-min
witness. The AP's loneliness lives precisely on these coarse rational times `j/d`. Replace one
mid-AP speed by `a(k-1)` and you delete, in one stroke, **all** the coarse tight times whose
denominator divides `k-1`. The more primes divide `k-1`, the more clocks die, the closer the
surviving max-min creeps to the floor.

This is the same "kill the clock at the tie-wall" instinct that codex's apex-lift sheaf
(HYP-2101) circles around — but here it is explicit and quantitative: *one* speed kills *all*
clocks dividing `k-1`, and the residual level is `ω(k-1)`.

## What it says about LRC

The recurring lesson of this project — "patterns that hold at small `n` break at a specific
arithmetic threshold" — recurs in its purest form. The old repo claim that the gap is *exactly*
`2/(2n-1)` (a clean analytic edge) was a small-box mirage. The truth is number-theoretic: the
**width of the LRC extremal basin** around the AP is `1/((a_max(k)(k+1)-1)(k+1))`, and
`a_max(k)` is read off the prime factorization of `k-1`. The hardness of LRC stays pinned at the
AP, but the *geometry of its neighborhood* is a Stern-Brocot staircase whose depth is arithmetic.

For LRC(14) (`k=13`): `k-1 = 12 = 2^2·3` has only the primes `2,3`, so `a_max(13) = 3`, and the
basin width is `g(13) = 1/574` at `3/41`. The fact that `12` is *not* divisible by `5` is exactly
why `k=13` does not dip to level 4 — and is a small, concrete instance of why `n=14` is "just
hard enough" but no harder.

## The constants

The mediant `a/(a(k+1)-1)` is the `a`-th convergent of the Stern-Brocot descent from `1/k` to
`1/(k+1)`; its denominator `a(k+1)-1` is the **lift depth** `c(k)/(k+1)` of the user's G2 lemma.
Generic `c(k) = (k+1)(2k+1)` (a=2); along primorials it is `(k+1)(a(k+1)-1)` with `a` unbounded.
The pseudo-doubling ratio `2 - 1/(k-1)` of the triangle foundation is the `a=2` shadow of this
same staircase — the higher rungs `a=3,4,…` are the deeper Cayley-Dickson-like descents that only
"unlock" when `k-1` supplies the requisite primes.
