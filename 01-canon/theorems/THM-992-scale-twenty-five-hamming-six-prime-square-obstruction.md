---
id: THM-992
title: Scale-twenty-five Hamming-six prime-square obstruction
status: PROVED STRUCTURAL + REFEREED FINITE-EXACT — exactly 36 scalar rows survive in three multiplication orbits; a five-coset overlap proof makes every owner miss at least three sheets, and independent Python-set and literal-CRT C++ sorted-vector implementations agree on every proof-bearing digest
source: codex-2026-07-17-S66 scale-twenty-five continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_five_hamming_six_prime_square_obstruction_codex_c25.py
  - 05-knowledge/results/lrc13_scale_twenty_five_hamming_six_prime_square_obstruction_codex_c25.out
  - 04-computation/lrc13_scale_twenty_five_hamming_six_referee_codex_c25.cpp
  - 05-knowledge/results/lrc13_scale_twenty_five_hamming_six_referee_codex_c25.out
---

# THM-992 — the scale-twenty-five prime-square face is empty

For `c=25`, the effective orders are `1,5,25`.  Hereditary leave-one-out
lcm is equivalent to having at least two order-25 coordinates.  This gives
473 order words and 243,750,000 literal state words per support, hence

```text
924*243,750,000 = 225,225,000,000
```

unquotiented labelled state contexts.

## Structural scalar collapse

Let `B` be the order-five providers and write `c25` for the number of
order-25 providers.  An order-one coordinate already forces capacity at most
23 at an order-25 owner.  Once only orders five and twenty-five remain, the
capacity at owner `y` has the form

```text
30 - c25 - 5*z_y - delta_y,
```

where `z_y` counts order-five provider/owner ratios in `{4,9,12}`, and
`delta_y` records an antipodal order-25 provider.  Requiring capacity at all
six owners gives `z_y=0` and `c25 in {2,3,4,5}`.  Applying the same forbidden
ratio condition in both orientations puts at most one order-five provider in
each quadratic class: the symmetric closure of the forbidden-ratio Cayley
relation is exactly `K6 disjoint-union K6` on the quadratic and nonquadratic
classes.  This eliminates `c25=2,3`.  If `c25=5`, the simultaneous
conditions `z_y=delta_y=0` demand an antipodal-free support even though the
single order-five provider forbids `{-b,±3b}`, removing a whole other
antipodal class.  Hence `c25=4`.

The two order-five labels lie in opposite quadratic classes and force the
support to be the complement of their `{3,10,12}` multiples.  Exactly 36
labelled rows remain.  Multiplication by `F_13^*` splits them into three
orbits of size twelve, represented by

```text
B={1,2}, support={1,2,4,5,8,9};
B={1,5}, support={1,4,5,6,7,9};
B={1,6}, support={1,2,4,6,9,11}.
```

## Five-coset owner obstruction

Partition the 25 sheets through the quotient `pi: Z/25 -> Z/5`, with fibres
`Q_0,...,Q_4`.  An order-five mask is one complete five-sheet `Q_j`.  At self
ratio it is `Q_3`; at any other allowed ratio its four unit choices occupy
`Q_0,Q_1,Q_2,Q_4`.  An order-25 mask has the following invariant signature:

- ratios `4,9`: one point in each non-`Q_3` class;
- ratio `12`: three points in three non-`Q_3` classes;
- every other ratio: one `Q_3` point and one point in three of the other four
  classes.

At an order-five owner, the two ratios `4,9` among the order-25 providers hit
the moving order-five coset and the two nonresidue providers hit `Q_3`.
Therefore at least four of their sixteen points overlap the ten order-five
sheets, and the union has size at most `10+(16-4)=22`.

At an order-25 owner, the order-25 ratios include `1,12` and two nonresidues.
If the two order-five cosets differ, all four order-25 masks meet their union,
so the total is at most `10+(15-4)=21`.  If the cosets coincide, the thick
part has size only five, so even the overlap-free relaxation gives the
stronger bound `5+15=20`.  Thus every owner misses at least three of the 25
sheets.
No global unit word exists.

This is a finite Kakeya/toothpick incidence proof: order-five masks are thick
parallel fibres and order-25 masks are rigid three/four-point needles.  Every
surviving owner either pays a four-incidence overlap tax against two distinct
fibres or collapses the thick starting mass from ten to five when those fibres
coincide.  The quotient discards exact needle positions but preserves exactly
this disjunction.

The primary's separate complete literal DP sharpens the maxima to 72
order-five owner rows of size 22 and 144 order-25 owner rows of size 21; all
36 contexts have zero feasible owners.  The structural coset proof, not that
DP, is the logical obstruction.

## Frozen primary and carrier audit

The standard-library primary reconstructs every CRT base both algebraically
and by literal search, proves the closed cardinality formulas, checks all 729
order words against the prime-square grammar, and compares the structural
scalar predicate against all `924*473=437,052` labelled order contexts.  Its
exact scalar feasible-owner histogram is

```text
0:1156, 1:149868, 2:171636, 3:90884, 4:21864, 5:1608, 6:36.
```

The 36 survivors represent 92,160,000 literal unit words.  The 216 full union
banks contain 23,338,080 reachable masks and have size histogram

```text
45200:24, 48380:24, 48540:24,
133390:48, 140330:48, 141430:48.
```

Normal and `python -O` runs reproduce the frozen 38-line output byte-for-byte.
The independent C++ referee reconstructs every CRT base by literal congruence
search, discovers the scalar rows from cardinalities rather than calling the
structural predicate, and constructs each reachable bank as a sorted vector.
Forward and reverse provider traversals agree exactly in all 216 owner rows.
It derives the forbidden Cayley graph and audits the distinct/coincident fibre
case split directly from literal masks.  `-O2`, `-O0`, and ASan+UBSan
(`detect_leaks=0`) outputs are byte-identical.  Frozen SHA-256 values are

```text
Python primary source  2b40edafa026dbc94f75a21ae2fbc764588d64e923779996760e746006b57870
Python primary output  f0c3cb63459c6cfccd024f05cc76d8d7819ed53217d469aaaf9682cd91e2ccdb
C++ referee source     a133eb9624d6c8c1b1118b89c220a52397e0203bc9685ffaca7c178a8530cc85
C++ referee output     b927a01c6439ea0db9b1e4d8a366b8c9ebdd141abc4fd6a5ab51a7b0bafbaaec
```

The proof-facing pair relation is the directed forbidden-ratio Cayley graph:
36 arcs, 6 reciprocal antipodal pairs, and symmetric closure
`K6 disjoint-union K6`.  It is not a tournament, and its independent-set
bound is the scalar obstruction.  For owner telemetry, compare
`(locally feasible, maximum union, scalar capacity)` lexicographically and
break ties in coordinate order.  Every one of the 36 completions is
transitive, with score multiset `{0,1,2,3,4,5}`, seven tie edges, no directed
triangle, six singleton SCCs, and one Hamiltonian path.  This completion
forgets the absolute threshold and therefore cannot prove the theorem.

Quadratic classes plus the two order-five obligations are faithful at the
scalar layer; quotient fibres plus owner labels are faithful for the overlap
tax.  Providers alone, isolated sheets, gaps, wall events, unlabelled residues,
Fano points, or chi-seven colours destroy shared-unit, owner, or fibre
incidence.

This theorem does not cover H5 ramification, non-AP/deep sheets, or global
sporadic emptiness.  Since THM-860 excludes every multiple of thirteen, scale
26 is primitive-impossible; the next legal untreated common scale is `c=27`.
