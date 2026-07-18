---
id: THM-983
title: Prime common-scale Hamming-six scalar obstruction above seventeen
status: PROVED STRUCTURAL + REFEREED FINITE-EXACT — the exact floor formula and B6/B5 table leave only p=23,29; p=23 is excluded by an exact Cayley-spectrum bound, p=29 by multiplicative closure, and a C++ primary plus two independently structured Python referees replay both complete 924-support banks
source: codex-2026-07-17-S66 prime-scale residue-capacity synthesis and exact referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-981]
related: [THM-962, THM-974, THM-980, THM-982, HYP-6820]
verification:
  - 04-computation/lrc13_prime_scale_hamming_six_scalar_obstruction_codex_prime.cpp
  - 05-knowledge/results/lrc13_prime_scale_hamming_six_scalar_obstruction_codex_prime.out
  - 04-computation/lrc13_prime_scale_hamming_six_scalar_referee_codex_prime.py
  - 05-knowledge/results/lrc13_prime_scale_hamming_six_scalar_referee_codex_prime.out
  - 04-computation/lrc13_prime_common_scale_hamming_six_scalar_referee_codex_S66.py
  - 05-knowledge/results/lrc13_prime_common_scale_hamming_six_scalar_referee_codex_S66.out
---

# THM-983 — every prime common scale above seventeen is scalar-impossible

Let `p>=19` be prime.  The primitive proper AP-centred common-scale-`p`
Hamming-six sheet bank is empty already at unit-independent scalar owner
capacity.

The only effective orders are `1,p`.  For an order-`p` provider/owner ratio
`r in F_13^*`, unit choice does not change the number of covered sheets.  The
cardinality has the exact residue-class description

```text
a_p(r)=#{x in Z : -p < x <= p and x == p*r (mod 13)}.
```

Write `p=13q+s`, `1<=s<=12`.  Each `a_p(r)` is `2q` plus a bonus depending
only on `(s,r)`.  Explicitly, if `a in {1,...,12}` represents `sr mod 13`,
then the half-open endpoint convention gives the exact all-`q` identity

```text
a_p(r)=floor((p-a)/13)-floor((-p-a)/13)
      =2q+floor((s-a)/13)-floor((-s-a)/13).
```

If `B_6(s)` and `B_5(s)` are the sums of the six and five largest bonuses,
their complete twelve-row table is

```text
s        1  2  3  4  5  6  7  8  9 10 11 12
B_6(s)   1  3  5  6  6  6  7  9 11 12 12 12
B_5(s)   1  3  5  5  5  5  6  8 10 10 10 10.
```

An order-one provider covers all `p` sheets only at its matching owner and
zero sheets elsewhere.  If a hereditary row is mixed, choose an owner whose
provider has order `p` (there are at least two such providers).  Every
order-one provider vanishes at that owner, while the at most five order-`p`
providers cover at most

```text
10q+B_5(s)<13q+s=p.
```

Indeed the deficit is `3q+s-B_5(s)>=1`.  Hence a scalar row must be
all-order-`p`.

For an all-order-`p` row, any owner has capacity at most
`12q+B_6(s)`.  Since `B_6(s)-s<=2`, this is strictly below `p` whenever
`q>=3`.  The primes below 39 are `17,19,23,29,31,37`.  THM-981 classifies the
exceptional scale seventeen.  The same table immediately excludes
`19,31,37`.  Only `23` and `29` reach the numerical threshold.  Both
exceptions have short structural obstructions.

### The `p=23` Cayley-spectrum obstruction

Here the seven high ratios

```text
H_23={1,2,3,6,7,10,11}
```

have bonus two and the five low ratios have bonus one.  An owner reaches 23
only if it sees at least five high ratios, or equivalently at most one low
provider.  Identify `F_13^*` with `Z/12` using generator two.  The low ratios
have exponent set

```text
D={2,3,6,8,9}.
```

Let `A` be the directed Cayley adjacency for `D` and
`M=(A+A^T)/2`.  Its character-`k` eigenvalue is

```text
cos(pi*k/3)+2cos(pi*k/2)+cos(2pi*k/3)+(-1)^k,
```

so the exact Fourier spectrum of `M` is

```text
5,-1,-2,-1,2,-1,1,-1,2,-1,-2,-1.
```

For a six-support indicator `x=(1/2)1+y`, one has `y` orthogonal to `1` and
`||y||^2=3`.  Therefore its number of internal directed low edges satisfies

```text
x^T A x = x^T M x >= 5*6^2/12 - 2*3 = 9.
```

If all six owners were feasible, each would receive at most one such edge,
so the same count would be at most six.  This contradiction excludes `p=23`
without a support census.

### The `p=29` closure obstruction

Here the five high ratios are

```text
H_29={1,4,5,8,9}.
```

Every owner must see all five of them.  Thus a hypothetical support `S` obeys
`yH_29 subseteq S` for every `y in S`.  Iteration gives
`y generated(H_29) subseteq S`, but `H_29` generates all twelve elements of
`F_13^*`.  This cannot fit inside a six-support.

Thus no prime `p>=19` reaches even the scalar owner gate.  The proof is uniform
in `p`, not a cutoff scan: adding thirteen to the scale adds two complete
residue blocks, and `B_6(s)-s<=2<q` closes every `q>=3`.  Within THM-860's
finite primitive range this simultaneously removes all 139 primes from 19
through 839.  The primary certificate additionally scans all 129,360 supports
at those primes plus the `p=17` anchor and finds only the two expected
quadratic rows at seventeen.

## Independent exact replay

A C++ primary and a first independently written Python referee derive the
same bonus table and separately scan the two exceptional 924-support banks.
The standard-library structural referee in this packet independently derives
the integer-floor formula, reconstructs the literal CRT sheet masks, and
checks unit-independence at both exceptional primes.  It performs 9,312 direct
formula checks through the full THM-860 range `c<=840` and 7,488 literal
exceptional mask checks.  It then enumerates all 924 six-supports in both
integer-capacity and binary Cayley coordinates.  The two predicates agree
owner by owner.  Its stronger feasible-owner histograms are

```text
p=23: 0:204, 1:528, 2:138, 3:48, 4:6;
p=29: 0:846, 1:72,  2:6.
```

Thus the finite replay strengthens the structural result: no `p=23` support
has more than four feasible owners, and no `p=29` support has more than two.
Multiplication gives 80 support orbits, with size histogram
`2:1, 4:1, 6:3, 12:75`; the quotient is telemetry only and is not used to
skip any support.

Frozen SHA-256 values are

```text
C++ primary source          86b34f9d404e10831ef9a86003dbf9f6fbc062054a9dfa884e20e853947911e8
C++ primary output          87eb95da5e554ab7ffa960ed8a485ca97ae5ffc75b78df7444119810b01620cc
Python census source        4ffd04a51c9c4584ec4c2e43806cd22703765f7921af42436ed41718b57c9883
Python census output        2570e3b5b285fb2b901e89113734ffdcbe8b7586865976c2e5b69c9d0bf19112
Python structural source    c699c773e20bec740a0f132b29ed93231c855e26bc44720d2c8d2b15290d7073
Python structural output    02a2a34271ff2fecd005b64433996f09d874fce6cb9243fe4e32d0158e11f758
```

Normal and `python -O` runs of both Python referees reproduce their respective
outputs byte-for-byte.  The C++ and first Python implementations landed
independently of the spectral/closure argument; agreement on the cardinality
table and exceptional emptiness is therefore not a shared-census assumption.

## Tournament and carrier audit

The proof-facing pair observable is the provider/owner ratio bonus.  At the
two exceptional primes its high/low switch is faithful because each bonus row
has only two heights.  It produces a directed Cayley graph, not a tournament:
reciprocal ratios can tie or point differently.  Both Cayley graphs are one
SCC; the `p=23` graph has degree six and 64 directed triangles, while the
`p=29` graph has degree four and 16.

For telemetry, compare the two directed bonuses on each unordered pair and
break ties along label order.  At both primes the completed tournament has
score histogram

```text
2:1, 3:3, 5:2, 6:2, 8:3, 9:1,
```

40 directed triangles, one SCC, 42 tie edges, 12 flips from the tie order,
and 124,961 Hamiltonian paths.  This tournament cannot prove the theorem: it
forgets the absolute capacity threshold and all symmetric high/low ties.

The challenged vertex assumption matters here.  Runners/providers alone are
faithful only when the complete directed ratio kernel and embedded support are
retained.  Ratio residues are the smallest proof carrier.  Unlabelled sheets,
units, gaps, divisor states, or projective pairs destroy owner incidence.
Individual sheets and units may be discarded only at this scalar gate because
unit-independence was proved first; they would be essential after a scalar
survivor.

There is a second, inequivalent tournament shadow: ordering the twelve ratio
weights and breaking ties by label gives a transitive tournament for every
residue class.  It checks the bonus order statistics but likewise forgets the
embedded six-support and absolute owner sums.

## Certificate cross-checks

The primary C++ implementation proves the complete-residue-block lemma,
reconstructs the bonus rows by literal integer counts, maximizes all five- and
six-subsets, scans the exceptional banks, and audits every THM-860 prime.  The
first independently written Python referee derives the same recurrence and
tables, independently scans `p=17,23,29`, and supplies the reciprocal-capacity
tournament fingerprint.  C++ `-O3`, `-O0`, and ASan/UBSan outputs agree
byte-for-byte; normal and optimized runs of both Python referees agree.

The theorem concerns only prime common-scale Hamming-six faces.  Composite
scales, H5 ramification, non-AP/deep sheets, and global sporadic emptiness
remain open.  ∎
