---
source: oracle-2026-06-01-S25
status: exploratory synthesis (tournament clock × holdback × twin primes × LRC)
tags:
  - tournament-analysis
  - tournament-clock
  - lonely-runner
  - twin-primes
  - admissibility
  - staircase
  - holdback
---

# Holdback, Twin Primes, and the Staircase: the difference-spectrum of the tournament clock

A short, dense session connecting the S24 tournament clock to twin primes and the
Lonely Runner via a single quantity — **holdback** — and finding that the
LRC-extremal runner system's holdback spectrum is exactly the repo's staircase.

## Holdback = 1/(2·difference)

In the tournament clock (S24) the phase-comparator edge `(i,j)` flips at the
walls `t = m/(2 d_ij)`, `d_ij = s_j − s_i`. So edge `(i,j)` **holds** one
orientation for an interval of length exactly

```
holdback(i,j) = 1 / (2 d_ij).
```

Holdback is *edge persistence*: how long a pairwise relation is "held back" from
flipping. The stickiest edges are the **small-difference** pairs; the maximum
holdback of a system is `1/(2·min_diff)`. The whole **holdback spectrum** is the
reciprocal of the **difference multiset** `{s_j − s_i}` — so Tournament Analysis
of a runner system is, at this level, additive combinatorics of its speed set.

## The extremal axis: consecutive integers carry the staircase

The LRC-extremal set `{1,…,n}` has difference multiset

```
d = 1 : n−1 times      (holdback 1/2)
d = 2 : n−2 times      (holdback 1/4)
…
d = n−1 : 1 time       (holdback 1/(2(n−1)))
```

— multiplicities `n−1, n−2, …, 1`. **That is the staircase `δ`** — the very
Young-diagram object the repo says "controls everything" (tilings, OCF, the
`H = 1 + 2^d` story). Computed holdback spectrum for `1..6`:
`1/2:5, 1/4:4, 1/6:3, 1/8:2, 1/10:1`. So:

> the lonely-runner-extremal configuration, the minimal tournament clock (S24),
> the maximum-holdback system, and the repo's staircase are **one object**, read
> four ways. Consecutive integers maximise holdback (min difference 1) and, as
> an arithmetic progression, minimise the number of *distinct* differences (only
> `n−1`), so their clock is maximally synchronised — minimal.

## Twin primes are the stickiest edges the primes allow

For a set of **odd primes** the minimum difference is `2` (consecutive integers
are impossible among primes `> 2`), achieved exactly at **twin-prime pairs**. So:

> twin primes are to the prime tournament clock what consecutive integers are to
> the integer clock — the *maximum-holdback (stickiest) edges*, holdback `1/4`.

The single exception is `(2,3)`, the unique consecutive prime pair (difference 1,
holdback `1/2`) — the one place the primes touch the integer-extremal floor.
Computed: odd primes `3..19` have `min_diff = 2`, `twin(d=2) edges = 4`, max
holdback `1/4`; the prime set `{2,3,5,7,11,13,17}` has `min_diff = 1` *only*
because of the `(2,3)` edge. The twin-prime conjecture, in this language, is the
statement that the prime clock has rate-2 maximum-holdback edges arbitrarily far
out — the primes never stop producing their "stickiest possible" pairs.

## Why the stickiest primes are still LRC-easy: the residue-missing duality

Here is the twist. Maximum holdback "wants" LRC-hardness — the extremal set has
it. Yet twin-prime / prime speed sets are **LRC-easy** (positive gaps, computed).
The reconciliation is arithmetic, and it is the bridge to twin primes proper:

- By the sieve lemma (THM-369) a speed set has an immediate lonely time `t = 1/m`
  whenever it **misses residue 0 mod m** (no speed divisible by `m`).
- Consecutive integers `{1,…,n}` **cover** residue 0 mod *every* `m ≤ n` (since
  `m` itself is in the set): anti-admissible, no easy lonely time → tight.
- Prime sets **miss** residue 0 mod every composite and every absent prime: e.g.
  `{3,5,7,11,13,17,19}` misses `m = 2,4,6` → three immediate lonely times → easy.

So LRC-hardness needs **divisibility-richness** (cover residue 0 mod all small
`m`); prime/twin-prime sets are **divisibility-poor** (miss those residues). This
is a sibling of the prime-`k`-tuple **admissibility** condition (a constellation
must *miss* a residue mod every prime to be viable). The *same arithmetic
feature* — **missing residue classes** — makes a set good for prime constellations
(twin primes need an admissible shape) and easy for the lonely runner. The two
famous problems pull on opposite ends of the covering/missing line:

```
LRC counterexample   →  divisibility-RICH  (cover residue 0 mod every m ≤ n)
prime constellation  →  divisibility-POOR / admissible (miss a residue mod each p)
```

Twin primes live at the divisibility-poor end; that is exactly why their high
holdback never converts into lonely-runner difficulty. **Stickiness is geometric
(differences); LRC-hardness is arithmetic (divisibility); twin primes are sticky
but arithmetically loose, so they stay lonely-easy.**

## Synchrony and clock complexity

Equal differences share *identical* walls `m/(2d)`, so all edges of a given
difference flip **in lockstep**. The number of *distinct* differences is the
number of independent flip-rhythms, hence a proxy for clock complexity:

- arithmetic progressions (incl. consecutive): only `n−1` distinct differences →
  maximally synchronised → minimal clock (matches S24's "extremal = minimal");
- twin clusters: many `d = 2` edges synchronise at rate 2;
- Sidon sets (e.g. geometric `1,2,4,8,16,32`: all differences distinct) →
  maximal complexity; their holdback spectrum has every value once.

## The unifying picture

```
speed set  →  difference multiset  →  holdback spectrum (=1/(2d))  →  clock
                     │                          │
            residue/divisibility          stickiness/persistence
                     │                          │
            LRC hardness (covering)     synchrony & clock complexity
```

Holdback is the geometric face (differences, persistence); admissibility is the
arithmetic face (residues, covering). The LRC-extremal set sits where both peak
(max holdback + full covering) and equals the staircase. Twin primes sit at max
*holdback within the primes* but minimal covering, so they are the sticky-yet-easy
corner. Tournament Analysis is the lens that shows these are the same axis.

## Next

1. Prove: holdback spectrum of `{1,…,n}` `=` staircase `δ_{n−1}` multiplicities
   (immediate), and relate the clock's cell-count to `Σ d = ` the staircase area.
2. Quantify "stickiness vs hardness": is there a set that is both high-holdback
   *and* anti-admissible besides scalar multiples of `{1,…,n}`?
3. The α = 1/k comparator (S24 next-step) should make "covers residue 0 mod k"
   literally a clock event, turning the residue-missing duality into geometry.

## Artifacts

```
04-computation/tournament_clock_holdback_twinprime_s25.py
05-knowledge/results/tournament_clock_holdback_twinprime_s25.out
```
