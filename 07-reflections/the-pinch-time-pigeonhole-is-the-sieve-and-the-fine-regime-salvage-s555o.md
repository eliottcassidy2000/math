---
source: oracle-2026-06-01-S555o
status: proof attempt (pinch-time pigeonhole for LRC@14) — rigorous reduction: the rational pinch IS the sieve; defeated by A1; the fine-regime salvage
tags:
  - lonely-runner
  - n14
  - pinch-time
  - pigeonhole
  - sieve
  - fine-regime
---

# The Pinch-Time Pigeonhole Is the Sieve — and the Fine-Regime Salvage

Attempting the proposed proof: *every 13-set has a pair whose pinch time clears all
other runners; pigeonhole over the candidate pinch times.* The idea is the right
**shape**, and pinning it down gives a clean rigorous reduction — and shows exactly
why it cannot close n=14 as stated.

## The pinch clears its pair (rigorous)

For a pair `(a,b)`, the **pinch time** `t = 1/(a+b)` puts `a` at `a/(a+b)` and `b` at
`b/(a+b)`, symmetric about `1/2`, both at circular distance
`min(a,b)/(a+b) ≥ 1/(a+b)` from the observer. Hence:

> **If `a + b ≤ n`, the pinch `t = 1/(a+b)` clears the pair** (`‖·‖ ≥ 1/n`).

Verified for the pairs summing to `n=14`: at `t=1/14` the pair `(a,14-a)` sits at
distance `a/14` (`(1,13)` borderline `=1/14`, `(2,12)→1/7`, …, `(6,8)→3/7`).

## The pinch IS the denominator sieve (the reduction)

A *third* runner `w` at the pinch `t = 1/(a+b)`: `‖w/(a+b)‖ ≥ 1/(a+b) ≥ 1/n` **iff
`(a+b) ∤ w`**. Therefore:

> **The pinch `t = 1/(a+b)` is lonely ⟺ no runner is divisible by `a+b`.** This is
> precisely the denominator sieve (THM-369) at `q = a+b`.

So the pinch-time pigeonhole — "some candidate pinch clears all runners" — is exactly
"**some `s ≤ n` has no multiple in the set**," i.e. **the set is *not* sieve-covered**.
The same holds for every rational pinch: any `t = p/q` with `q ≤ n` is a sieve witness,
lonely iff no runner is divisible by `q`. (A single multiple of `n` sits at `0` at
*every* `n`-gon vertex `t = j/n`, spoiling them all at once.)

## Why a counterexample defeats it — by construction

A counterexample is **sieve-covered** (necessary condition **A1**, S554): it contains
a multiple of *every* `q ∈ {2,…,n}`. So it defeats **every** rational pinch
simultaneously. Verified decisively (`lrc_pinch_time_pigeonhole_s555.py`):

> **Of 40 sieve-covered primitive 13-sets: `0` have any lonely rational `t = p/q`
> (`q ≤ 14`), and `40/40` are lonely at a *fine* time** (denominator `> 14`, e.g.
> `t ≈ 0.009, 0.022, 0.073`).

So the pinch-time pigeonhole **= the sieve**: it settles the decorrelated majority
(the non-sieve-covered sets, which are lonely at some `t = 1/s`), but is **defeated by
construction** on the sieve-covered core — the genuine open locus (S554). The
naive pigeonhole count also fails outright: `11` other runners against `~7` pinch
times, and a single multiple of `n` spoils all of them.

## The salvage: a FINE-REGIME pinch pigeonhole

The lonely times of a sieve-covered set live at **fine denominators `> n`** (the
S18 fine regime), where the sieve multiples no longer sit exactly at `0` and the
danger arcs genuinely thin out. So the user's "thin danger arc" intuition is *correct
— but only in the fine regime*. The honest reformulation:

> **The pinch pigeonhole can only work at fine pinch times `t = p/q` with `q > n`.**
> There a runner `w` is in danger iff `‖w/q‖ < 1/n`, i.e. `w` is within `q/n` of a
> multiple of `q` — a *thin* condition (`q/n < 1` fails, so `w mod q ∈` a band of
> width `2q/n`). The pigeonhole becomes: over a chosen family of fine denominators,
> one fine pinch avoids every (now genuinely thin) danger band.

This is exactly the unresolved core: it is the measure/covering problem (S550), the
7-gon *windows* (S552, which are the fine perturbations of `t=j/7`), and the
multiples-of-`n*` coupling. The rational pinch (`q ≤ n`) is the sieve; the fine pinch
(`q > n`) is the conjecture.

## Verdict / next
- The pinch-time pigeonhole, made precise, **is the denominator sieve** (a clean
  rigorous reduction): pinch `t=1/(a+b)` lonely ⟺ `(a+b) ∤` every speed.
- A counterexample defeats it **by being sieve-covered (A1)**; verified `40/40`
  sieve-covered sets are lonely only at fine denominators `> n`.
- The right idea lives in the **fine regime** (`q > n`), where danger arcs are thin —
  but that is the open core (S550/S552), carrying the multiples-of-`n*` coupling.
- Concrete next: (1) a **fine-regime pinch pigeonhole** — pick a family of denominators
  `q ∈ (n, Cn]` and bound the thin danger bands so one fine pinch is clear; (2) tie the
  fine pinch family to the S552 windows near `t = j/7`; (3) the count `#far ≥ n-2`
  (S553) says only one band can ever cover a pinch — combine with a fine family.

## Artifacts
```
04-computation/lrc_pinch_time_pigeonhole_s555.py
05-knowledge/results/lrc_pinch_time_pigeonhole_s555.out
```
Related: THM-369 (sieve = the rational pinch), S554 (A1 sieve-covered = the defeater),
S18 (coarse/fine regimes), S550 (measure core), S552 (7-gon fine windows), S553
(at most one near).
