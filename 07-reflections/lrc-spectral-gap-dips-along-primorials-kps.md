---
source: kind-pasteur-2026-06-19-S9
status: reflection on a verified-computational result (THM-539, HYP-2623)
tags: [lonely-runner, spectral-gap, stern-brocot, primorial, highly-composite, basin-width, lrc14]
---

# The LRC spectral basin is arithmetic, not analytic

**Prompt (user):** the second LRC-spectrum point is `M_k={1,…,k-1,2k} -> 2/(2k+1)`, so
`1/(k+1) < sigma_2(k) <= 2/(2k+1)`, mediant-tight for `k<=5`, strict for `k>=6`; the gap is
`Theta(1/k^2)` at most — does anything dip below it unboundedly?

## The answer, and the honest correction

**Partly — the dip is arithmetic, but bounded as far as we can prove.** The values of the max-min
collar just above the floor `1/(k+1)` are the Stern-Brocot mediants `a/(a(k+1)-1)`, `a=2,3,4,...`,
marching down to the floor. The gap `g(k)·k^2 -> 1/a`. So the whole question is: *how large can
the level `a` be?* The level is governed not by `k` but by the **factorization of `k-1`**: the
witnessing family `F(k,a) = {1,…,k-2, k, a(k-1)}` reaches level `a` when `k-1` carries the right
initial primes:

```
   a = 3  at  6  | (k-1),  k ≢ 1 mod 30   (k = 7, 13, 19, 25, …)   — exact, infinite family
   a = 4  at  30 | (k-1)                  (k = 31, 61, …, 181)     — exact, infinite family
   a = 5  at  k = 211 (k-1 = 210 = 2·3·5·7)  ->  NOT a dip: M = 1/212 = the FLOOR (tight!)
```

**The correction (MISTAKE-079).** I first read a covering test `M(F(211,5)) < 5/1059` as a
*deeper* dip and announced `a_max(k) ~ ω(k-1)` UNBOUNDED, `liminf g·k^2 = 0`. The adversarial
verification workflow computed the exact value: `M(F(211,5)) = 1/212`, the **floor** — the family
does not dip to level 5, it *collapses to a tight configuration* (`g=0`). `F(k,5)` is tight for
every `k` with `2·3·5·7 | (k-1)`. So the natural family **tops out at `a=4`**. The confirmed
picture: `g(k)·k^2 ∈ [≈1/4, 1/2]`, `limsup = 1/2`, and `liminf` is **not known to be 0** —
whether any construction reaches `a >= 5` is OPEN. This realigns with codex's S16/S17 "no
`o(1/k^2)` dip found."

**The lesson worth keeping:** `M < c` localizes but does not identify — a value below a mediant
can be a slightly-smaller special value OR the global minimum (the floor) in disguise. Compute the
exact extremal value at the extrapolation point before naming an asymptotic.

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
clocks dividing `k-1`. The residual level rises with the small-prime content of `k-1` — but only
up to `a=4` for this family. Push the multiplier one rung too far and the construction over-kills:
at `a=5` (`k=211`) it lands back ON the floor, becoming a second *tight* set rather than a deeper
dip. So the "more primes → deeper dip" intuition is real but saturates; it is not a monotone
ladder to the floor.

## What it says about LRC

The recurring lesson of this project — "patterns that hold at small `n` break at a specific
arithmetic threshold" — recurs in two ways here. First, the old repo claim that the gap is
*exactly* `2/(2n-1)` (a clean analytic edge) was a small-box mirage: below-edge dips (`a=3,4`)
exist on arithmetic subsequences of `k`. Second — and this is the cautionary half — my own
"the dip is unbounded along primorials" was itself a too-clean extrapolation, broken at the very
first untested primorial (`k=211`) where the family collapses to the floor. The honest statement:
the **width of the LRC extremal basin** around the AP is `1/((a(k)(k+1)-1)(k+1))` with the realized
`a ∈ {2,3,4}`; whether `a` can exceed 4 (whether the basin can be made arbitrarily narrow relative
to `1/k^2`) is OPEN. The hardness of LRC stays pinned at the AP; the *geometry of its neighborhood*
is a Stern-Brocot staircase whose first three rungs are arithmetic in `k-1`.

For LRC(14) (`k=13`): `k-1 = 12 = 2^2·3` has only the primes `2,3`, so `a_max(13) = 3`, and the
basin width is `g(13) = 1/574` at `3/41`. The fact that `12` is *not* divisible by `5` is exactly
why `k=13` does not dip to level 4 — and is a small, concrete instance of why `n=14` is "just
hard enough" but no harder.

## The constants

The mediant `a/(a(k+1)-1)` is the `a`-th convergent of the Stern-Brocot descent from `1/k` to
`1/(k+1)`; its denominator `a(k+1)-1` is the **lift depth** `c(k)/(k+1)` of the user's G2 lemma.
Generic `c(k) = (k+1)(2k+1)` (a=2); on the `a=3,4` dip families it is `(k+1)(a(k+1)-1)` (a up to 4 confirmed).
The pseudo-doubling ratio `2 - 1/(k-1)` of the triangle foundation is the `a=2` shadow of this
same staircase — the higher rungs `a=3,4,…` are the deeper Cayley-Dickson-like descents that only
"unlock" when `k-1` supplies the requisite primes.
