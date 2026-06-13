---
source: oracle-2026-06-01-S521
status: result (computational, exact reduction verified) + conjecture
tags: [twin-primes, goldbach, additive-basis, necklace, mod-6, A007534]
---

# Twin-prime Goldbach: the 35 exceptions are 11 mod-6 triples on a necklace

**Question (user):** twin-prime version of Goldbach — every even number a sum of
two twin primes. Which 35 even values are missed, and what deeper structure do
they point to?

A "twin prime" = a prime in some twin pair; T = {3,5,7,11,13,17,19,29,31,...}.
The even numbers NOT of the form t1+t2 (t1,t2 in T) are OEIS **A007534**,
conjectured finite. Computed exhaustively to 2,000,000:

```
EXACTLY 35 exceptions, largest 4208:
  2, 4,
  94, 96, 98,        400, 402, 404,     514, 516, 518,
  784, 786, 788,     904, 906, 908,     1114, 1116, 1118,
  1144, 1146, 1148,  1264, 1266, 1268,  1354, 1356, 1358,
  3244, 3246, 3248,  4204, 4206, 4208
```

## The structure: everything is forced by mod-6 centering

The exceptions are not scattered. Apart from the trivial pair {2,4} they form
**eleven runs of exactly three consecutive even numbers** — triples
{6m-2, 6m, 6m+2}. This is not a coincidence; it is forced.

**Twin pairs are centered at multiples of 6.** For p>=5 a twin pair is
{6k-1, 6k+1}, centered at c = 6k. Only the very first pair {3,5} is anomalous
(center 4). Let C = the set of twin-pair **centers** = {4, 6, 12, 18, 30, 42,
60, 72, 102, ...}; C \ {4} ⊂ 6Z.

A sum of two twin primes picks one endpoint from each of two pairs:
a in {c1-1, c1+1}, b in {c2-1, c2+1}, so
```
a + b ∈ { c1+c2-2,  c1+c2,  c1+c2+2 }.
```
**Each pair of centers (c1,c2) covers a full triple** centered at s = c1+c2.
Because C ⊂ 6Z (ignoring the anomaly), every covering sum s is a multiple of 6,
and the three evens nearest 6m — namely 6m-2, 6m, 6m+2 — are covered **as one
unit** iff 6m ∈ C+C. So representability is constant on each mod-6 triple, and
**a triple is all-exceptional iff its center 6m is not a sum of two twin-pair
centers.** That is *why* the misses come in threes.

Verified exactly: on every even n in (8, 6000], n is a sum of two twin primes
**iff** the unique multiple of 6 in {n-2,n,n+2} lies in C+C. **Zero mismatches.**

## The necklace: the real object is a thinned Goldbach problem

Divide out the 6. Let the **twin-center necklace** be
```
K = C/6 = { k : 6k-1 and 6k+1 are both prime } = {1, 2, 3, 5, 7, 10, 12, 17,
            18, 23, 25, 30, 32, 33, 38, 40, 45, 47, 52, 58, 70, 72, ... }
```
Then twin-prime Goldbach is **exactly** ordinary binary Goldbach on K:

> even n (>8) is a sum of two twin primes  ⇔  round(n/6) = k1 + k2 for some
> k1,k2 ∈ K.

The 11 exception-triples correspond to the 11 integers that are **not** a sum of
two necklace elements:
```
m = 16, 67, 86, 131, 151, 186, 191, 211, 226, 541, 701
( ×6 = 96, 402, 516, 786, 906, 1116, 1146, 1266, 1356, 3246, 4206 )
```
The largest is **m = 701**. The twin-Goldbach conjecture is therefore equivalent
to the single statement: **every integer m > 701 is a sum of two elements of the
twin-center necklace K** — an additive-basis-of-order-2 claim for a set K of
density ~ (twin-prime constant)·n/log^2 n, i.e. *sparser than the primes*.

## What it points to (the "deeper structure")

1. **The mod-6 quotient is the natural stage, not Z.** The "complement necklace"
   the user has in mind is `K`: the twin-Goldbach problem is a binary additive
   problem one level *down*, on the indices of twin-pair centers. The fattening
   k ↦ {6k-2,6k,6k+2} is the only thing that makes the surface look 3-wide.
2. **Triples = the +-2 thickening of a sumset.** A007534's "triples" are an
   artifact of representing each abstract center-sum 6m by a 3-element residue
   window. Strip the window and the data is clean: 11 holes of K+K.
3. **It is harder than ordinary Goldbach, in the expected way.** K is sparser
   than P, and binary Goldbach on a sparse set can genuinely fail finitely often
   before the density wins — which is exactly the observed 11-hole comb ending
   at 701. The tail (last hole 701, then a 3000-wide gap to nothing) mirrors the
   ordinary-Goldbach "the exceptions die once the count of representations
   r(m) >= 1 for good") heuristic: r_K(m) grows and the holes stop.
4. **Connection to the repo's pair-lens.** This is the same "work on PAIRS, not
   points" move as `pair-first-twin-prime-lens-s502` / HYP-1965/1966: twin primes
   ARE the pair objects; here the pair-center coordinate 6k is exactly the lens
   coordinate, and Goldbach-on-pairs becomes Goldbach-on-K.

## Open

- **HYP-1994:** m=701 is the last hole of K+K (twin-Goldbach has no exception
  above 4208). Needs: r_K(m) (number of necklace representations) is positive
  for all m>701 — verify r_K(m)>0 far out, and bound r_K below via twin-prime
  density (Hardy–Littlewood). Plausibly provable conditionally on a quantitative
  twin-prime conjecture even though twin-prime infinitude itself is open.
- Why these particular 11 holes? Are the m-values (16,67,86,131,151,186,191,
  211,226,541,701) themselves structured (residues, near-K-gaps)? First pass:
  they sit in stretches where K is locally thin.

## Convergence note (two sessions, same structure)

A concurrent session, `oracle-2026-06-01-S516`
(`07-reflections/twin-goldbach-35-exceptions-the-six-wheel-s516.md`), reached the
same picture independently from the "6k±1 residue wheel" angle and the same
"complement necklace" language the user used. It also flagged that the eleven
necklace holes m = 16, 67, 86, 131, 151, 186, 191, 211, 226, 541, 701 are
**mostly ≡ 1 (mod 5)** — checks out: 10 of 11 are ≡1 mod 5 (only 67 ≡ 2). So the
holes are not just sparse, they sit in a single residue lane mod 5 — a second
modular constraint stacked on the mod-6 one. Worth chasing: is the ≡1 mod 5
bias a real feature of K+K holes or small-sample? Two independent derivations
converging on the same reduction is itself the signal that the mod-6 necklace
is the right object, not an artifact of one viewpoint.

Scripts: `04-computation/twin_goldbach_exceptions_s521.py` (the 35 + run/residue
structure), `04-computation/twin_goldbach_necklace_reduction_s521.py` (the exact
reduction + the 11 K+K holes). Results in `05-knowledge/results/`. See also
HYP-1994 and the S516 six-wheel reflection.
