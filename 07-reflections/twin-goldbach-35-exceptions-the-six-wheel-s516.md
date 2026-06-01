---
source: oracle-2026-06-01-S516
status: exploratory (twin-prime Goldbach exceptions; the 6k±1 wheel)
tags:
  - twin-primes
  - goldbach
  - additive-multiplicative
  - residue-wheel
  - complement-necklace
  - lonely-runner
---

# Twin-Prime Goldbach: the 35 missed evens and the six-wheel they ride

**Question.** Goldbach with twin primes — write an even number as a sum of two
*twin primes* (members of a twin pair). Which even values are missed, and what
do they point to?

## The 35 (computed, conjecturally complete = OEIS A007534)

Even `N` not a sum of two twin primes (twins = {3,5,7,11,13,17,19,29,31,…}):

```
2, 4,
94, 96, 98,     400,402,404,    514,516,518,    784,786,788,
904,906,908,    1114,1116,1118, 1144,1146,1148, 1264,1266,1268,
1354,1356,1358, 3244,3246,3248, 4204,4206,4208.
```

Exactly **35** values; largest **4208** (nothing missed in `[4210, 200000]`, so
the list is conjecturally final). Decompose as

```
35 = 2 + 11·3 :  the trivial pair {2,4}  +  ELEVEN consecutive-even TRIPLES.
```

## The structure: triples centered on multiples of 6

Every non-trivial exception is a triple `{6m−2, 6m, 6m+2}`. The **eleven centers
are all ≡ 0 (mod 6)** and, refining, all ≡ **6 or 12 (mod 30)**:

```
centers 6m:  96,402,516,786,906,1116,1146,1266,1356,3246,4206
   c/6:      16,67,86,131,151,186,191,211,226,541,701   (mostly ≡1 mod 5)
   c mod 30: {6,12}  — never 0,18,24
```

Why triples? **All twin primes `>3` live at `6k±1`** (the residues coprime to 6).
A sum of two of them lands at `(±1)+(±1) ∈ {4,0,2} (mod 6)` — *all three even
residues*, but by different "twin types":

```
N ≡ 4 (mod 6) = (6a−1)+(6b−1)   lower + lower
N ≡ 0 (mod 6) = (6a−1)+(6b+1)   lower + upper
N ≡ 2 (mod 6) = (6a+1)+(6b+1)   upper + upper
```

A triple `{6m−2, 6m, 6m+2}` is exactly one number of each type. So an exception
triple is a value of `m` at which **lower+lower, lower+upper, and upper+upper
twin-sums all fail simultaneously** — a local "twin-prime desert" wide enough to
kill all three channels at once. The `mod 30` refinement (`centers ≡ 6,12`) is the
next layer of the same sieve: the `mod 5` structure further selects which `m`
can be triple-deserts.

`{2,4}` is degenerate — below the minimal twin-sum `3+3 = 6`.

## The "complement necklace" reading

The right object is the **residue wheel** — a cyclic necklace of beads `0..5`
(mod 6), or `0..29` (mod 30), etc. (the primorial wheels `2, 2·3, 2·3·5, …`).
Twin primes occupy only the **unit beads** `{1,5} (mod 6)` (coprime to 6); the
**complement beads** `{0,2,3,4}` are forbidden to them. So:

> twin primes = an indicator supported on the unit-beads of the wheel;
> "sum of two twin primes" = the **self-convolution** of that indicator on the
> wheel; the 35 missed evens are the **holes of the convolution**.

The wheel/complement split is the additive↔multiplicative interface in miniature:
the **wheel is multiplicative** (residues coprime to `6 = 2·3`, the small-prime
sieve), the **convolution is additive** (Goldbach). This is exactly the addition/
multiplication duality threaded through the LRC work — and "complement" here is
the forbidden (non-unit) beads, the wheel's self-complementary clasp at the
multiple-of-6 center of each triple.

## What it points to (deeper)

1. **The 6k±1 (primorial) wheel** is the load-bearing structure: every feature —
   triples, the `{4,0,2} mod 6` channels, the `mod 30` refinement — is the wheel
   showing through. Twin-Goldbach is convolution on the wheel.
2. **A finite exceptional set from a sparse additive set.** Twin primes are sparse
   (density `~ C/log²`), yet their sumset misses only 35 evens — a *covering*
   phenomenon. This is the additive-combinatorics sibling of the LRC sieve /
   admissibility (S25, HYP-1851): a sparse set that nonetheless covers all but
   finitely much, with the exceptions governed by local (wheel) deserts.
3. **Triples = simultaneous failure across all residue channels** — the same
   "must cover every channel at once" logic as a Lonely Runner counterexample
   (which must hit every modulus). The exceptions are where the channels conspire.

## Caveat / next
- Finiteness (exactly 35) is conjectural (verified to 2·10^5 here; standard).
- Next: read the convolution on the `mod 210 = 2·3·5·7` wheel and check whether
  the eleven centers separate cleanly by the `mod 7` layer; relate the triple-
  desert condition to twin-prime gaps directly.

## Artifacts
```
(computation inline; OEIS A007534)
```
