---
id: THM-992
title: Scale-twenty-five Hamming-six prime-square obstruction
status: PROVED STRUCTURAL + TRIPLE-CERTIFICATE FINITE-EXACT — exactly 36 scalar rows survive in three multiplication orbits; the owner-normalized five-coset proof makes every owner miss at least three sheets, and one primary plus independently written Python/literal-search and C++/sorted-vector referees agree on every decisive census
source: codex-2026-07-17-S66 scale-twenty-five continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_five_hamming_six_prime_square_obstruction_codex_c25.py
  - 05-knowledge/results/lrc13_scale_twenty_five_hamming_six_prime_square_obstruction_codex_c25.out
  - 04-computation/lrc13_scale_twenty_five_hamming_six_structural_referee_codex_c25.py
  - 05-knowledge/results/lrc13_scale_twenty_five_hamming_six_structural_referee_codex_c25.out
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

Fix an owner and partition the 25 sheets into their five residue classes
modulo five.  The raw residue of the owner's self order-five mask depends on
the owner.  Relabel the five classes **for this owner** so that this self class
is called `Q_3`.  This normalization is essential: `Q_3` is not one global
sheet residue shared by all six owner projections.  An order-five mask is one
complete five-sheet `Q_j`.  At self ratio every unit gives `Q_3`; at any other
allowed ratio its four unit choices occupy `Q_0,Q_1,Q_2,Q_4`.  An order-25
mask has the following signature in the same owner-local names:

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
so the total is at most `10+(15-4)=21`; if the cosets coincide, the stronger
bound is twenty.  Thus every owner misses at least three of the 25 sheets.
No global unit word exists.

This is a finite Kakeya/toothpick incidence proof: order-five masks are thick
parallel fibres, order-25 masks are rigid three/four-point needles, and every
surviving packet pays a four-incidence overlap tax.  The quotient discards
exact needle positions but preserves precisely the intersection lower bound
needed for the theorem.

## Frozen primary and two independent referees

The standard-library primary reconstructs every CRT base both algebraically
and by literal search, audits all 729 order words against the prime-square
grammar, derives the forbidden-ratio Cayley graph `K6 disjoint-union K6`, and
compares the structural scalar predicate against all 437,052 labelled
support/order contexts.  Its separate exact union DP produces the same sharp
owner maxima.  A fresh execution directly from the frozen git object matches
its stored output byte-for-byte.

The primary's exact scalar feasible-owner histogram over all 437,052 rows is

```text
0:1156, 1:149868, 2:171636, 3:90884, 4:21864, 5:1608, 6:36.
```

The 36 scalar survivors represent 92,160,000 literal unit words.  This full
histogram is additional primary telemetry; the structural classification and
the 36-row scalar digest are what the independent referee rederives.

A standard-library-only referee, written without reading a scale-25 primary,
reconstructs each CRT representative by literal search and checks all 3,600
label/order/unit/owner masks against the independent period-cardinality
formula.  It verifies 316,008 instances of the scalar formula above and
reconstructs the 36 rows directly from all `924*473=437,052` labelled
support/order contexts.  Their owner capacities are

```text
25:144, 26:72.
```

The referee then performs 2,880 order-25 owner-normalized coset-signature
checks.  As a guard against accidentally treating the normalization as
global, it records the raw self-coset residues for owners `1,...,12`:

```text
(3,1,2,0,4,3,1,0,4,2,3,1).
```

Finally it builds every owner-local projected union bank without quotienting
supports, owners, or unit masks.  The terminal banks contain 23,338,080
reachable-mask incidences (counting a mask once in every owner bank where it
occurs), with exact histograms

```text
maximum union:       21:144, 22:72
reachable-bank size: 45200:24, 48380:24, 48540:24,
                     133390:48, 140330:48, 141430:48.
```

Thus the elementary five-coset bounds are sharp at both owner types, and all
216 owner rows are infeasible.  Normal and `python -O` executions agree
byte-for-byte.  The two programs independently agree on the grammar and scalar
bank hashes as well as the obstruction-bearing counts:

```text
grammar SHA-256  7ae50439ddbd7e09d37516d067fe20f35d1f36b7830c25dc7156c900c6fde62f
scalar SHA-256   ad266f55f820615eb2f7b4e323b6599024842c56352ff65c1b6dcdd117c250f9
```

Their mask and owner-bank hashes use deliberately different serialization
orders, so equality is not asserted for those digests.  Frozen artifact hashes
are

```text
primary source             2b40edafa026dbc94f75a21ae2fbc764588d64e923779996760e746006b57870
primary output             f0c3cb63459c6cfccd024f05cc76d8d7819ed53217d469aaaf9682cd91e2ccdb
structural referee source  bdd745a20ff17116707e72bfb6af24cfda961ab913a8572b130a5925a6d073cd
structural referee output  2f9aab211505a41fd3442dcd5f88e9b94d9481851792b2e9aafe4f80a1b3cd16
```

The separately written C++ referee gives a third implementation.  It finds
every CRT base by literal congruence search, discovers the scalar rows from
cardinalities rather than calling the structural predicate, and represents
each reachable bank as a sorted vector.  Forward and reverse provider
traversals agree on all 216 owner rows.  It independently derives the
forbidden Cayley graph and checks the distinct/coincident fibre case split
from literal masks.  Optimized, unoptimized, GCC, and address/undefined-
sanitized runs are byte-identical.  Its frozen hashes are

```text
C++ referee source  a133eb9624d6c8c1b1118b89c220a52397e0203bc9685ffaca7c178a8530cc85
C++ referee output  b927a01c6439ea0db9b1e4d8a366b8c9ebdd141abc4fd6a5ab51a7b0bafbaaec
```

The structural coset proof, not the projected DP, is the logical obstruction.
The primary checks the raw five-coset table at owner one; the Python referee
closes the covariance/presentation gap by checking all twelve owners after
explicitly renaming each owner's self coset; the C++ referee changes both the
CRT construction and the reachable-bank representation.  Together they
promote the result from a scratch claim to a triple-certificate theorem.

## Tournament and carrier audit

The proof-facing pair relation is the directed forbidden-ratio Cayley graph:
36 arcs, 30 edges in its symmetric closure, six reciprocal antipodal pairs,
and two components of size six.  It is not a tournament; its `K6+K6`
independent-set bound is the scalar obstruction.

On owner vertices, the primary and C++ referee orient
`(locally feasible, maximum union, scalar capacity)` and the structural Python referee orients
`(maximum union, scalar capacity, reachable-bank size)`, in each case
lexicographically with coordinate-order tie breaks.  Both gauges give a
transitive tournament in every one of the 36 rows: score sequence
`(0,1,2,3,4,5)`, zero directed triangles, six singleton SCCs, and one
Hamiltonian path.  The primary has seven tie edges in every row; the richer
referee gauge has tie histogram `2:24,7:12`.  Its flip histogram is

```text
1:1, 2:1, 3:2, 4:10, 5:2, 6:7, 7:7, 8:1, 9:4, 10:1.
```

These total-order shadows are diagnostic only.  They forget both the absolute
25-sheet threshold and the incidence statement that an order-25 mask occupies
three or four of the four nonself five-cosets, which is exactly what forces
overlap.

The challenged vertex choice sharpens the proof object.  The scalar collapse
is carried by ratio residues, their two quadratic classes, and antipodal
pairs.  After that collapse, the faithful small carrier is the five owner-local
sheet cosets with provider masks as hyperedges.  Runners/owners alone lose
coset incidence; individual sheets retain unnecessary affine offset data;
units alone lose which coset is the owner's self class.  This is a colored
five-vertex incidence system, not naturally a tournament.

A compact Lean terminal quotient is therefore feasible: formalize an abstract
five-class partition and the two intersection lemmas giving bounds 22 and 21,
then discharge the finite ratio/signature tables separately.  It should not
pretend to formalize the 225,225,000,000 raw contexts, and no such unverified
Lean module is claimed here.

This theorem does not cover H5 ramification, non-AP/deep sheets, or global
sporadic emptiness.  THM-860 already rules out primitive `c=26`, because a
primitive packet has `13` not dividing its common scale.  THM-993 contains the
frozen primary attack on the next numerically possible common scale, 27.
