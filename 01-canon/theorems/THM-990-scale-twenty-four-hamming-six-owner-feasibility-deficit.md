---
id: THM-990
title: Scale-twenty-four Hamming-six owner-feasibility deficit
status: PROVED STRUCTURAL + REFEREED FINITE-EXACT — independent algebraic-CRT Python/set and literal-CRT C++/sorted-vector certificates exhaust all 154,461,339,648 labelled state contexts, prove that every scalar survivor has at least two impossible owners, and classify the thirty extremal rows by a cubic-coset three-sheet-class nerve
source: codex-2026-07-17-S66 scale-twenty-four continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989]
related: [THM-978, THM-980, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_four_hamming_six_frontier_scout_codex_c24.py
  - 05-knowledge/results/lrc13_scale_twenty_four_hamming_six_frontier_scout_codex_c24.out
  - 04-computation/lrc13_scale_twenty_four_hamming_six_referee_codex_c24.cpp
  - 05-knowledge/results/lrc13_scale_twenty_four_hamming_six_referee_codex_c24.out
---

# THM-990 — scale twenty-four has at least two impossible owners

This theorem closes the next legal primitive proper AP-centred common-scale
Hamming-six face after scale 22; scale 23 is uniformly prime-excluded by
THM-983.

For `c=24`, the effective orders are
`1,2,3,4,6,8,12,24`, with twenty-four literal `(D,e)` states.  Hereditary
leave-one-out lcm enumeration gives 108,813 divisor words and 167,165,952
literal state words per support, hence

```text
924*167,165,952 = 154,461,339,648
```

unquotiented labelled state contexts.  A complete algebraic-CRT scratch
reconstruction leaves 66,984 scalar contexts on 854 supports across 202 order
multiplicity profiles.  Exact owner-local set-union reachability gives

```text
0 feasible owners: 64,962 contexts,
1 feasible owner :  1,800 contexts,
2 feasible owners:    192 contexts,
4 feasible owners:     30 contexts.
```

Thus every scalar row has at least two empty owner projections.  The 401,904
owner rows have maximum-union histogram

```text
12:72, 14:2136, 15:1644, 16:15876, 17:24420, 18:76296,
19:94872, 20:104592, 21:53040, 22:24948, 23:1704, 24:2304.
```

There are 20,302 distinct owner maximum vectors, and no reachable-mask bank
exceeded 7,728 states.  A global unit word would project to a feasible word at
all six owners, so the local deficit is terminal for this common-scale face.

## Cubic structure of the thirty extremal rows

The thirty rows with four feasible owners are not an unstructured tail.  Let
`g=2` generate `F_13^*`, and write

```text
C0 = {g^e : e=0 mod 3} = {1,5,8,12},
C1 = {g^e : e=1 mod 3} = {2,3,10,11},
C2 = {g^e : e=2 mod 3} = {4,6,7,9}.
```

The independent referee proves the following exact classification.  After a
common dilation by some `a in F_13^*`, every extremal row has

```text
support                 a*C2 union E,
four orders on a*C2     3,3,3,3,
two equal orders on E   D,D,
absent cubic coset      a*C1,
```

where `E` is a two-subset of `a*C0`.  There are precisely three cases:

```text
E adjacent in the cyclic order-four coset (ratio g^(+/-3)), D=8  : 12 rows,
E adjacent in the cyclic order-four coset (ratio g^(+/-3)), D=24 : 12 rows,
E antipodal in that coset (ratio g^6=-1),              D=8       :  6 rows.
```

Thus the three multiplication-orbit sizes are `12,12,6`.  The feasible owners
are exactly the four labels in `a*C2`; the two impossible owners are exactly
`E`.  For `D=8`, the scalar capacities are `30` on `a*C2` and `25` on `E`;
for `D=24`, they are respectively `32` and `24`.  In both cases the true
maximum union is `24` on `a*C2` and only `19` on `E`.

This five-sheet deficit has a small three-class nerve proof.  Partition the
twenty-four sheets as

```text
S_j = {s in Z/24 : s=j mod 3},  j=0,1,2.
```

At an owner in `a*C2`, the four order-three provider ratios are `C0`.  Their
exact distinct mask fibres are

```text
ratio 1       : {S2},
ratios 5,8    : {S0,S1},
ratio 12      : {empty}.
```

The ratio-`1`, ratio-`5`, and ratio-`8` providers therefore cover
`S2,S0,S1`, so this owner is feasible without using either high-order
provider.

At an owner in `E`, the order-three ratios are `C2`, with fibres

```text
ratios 4,9    : {S0,S1},
ratios 6,7    : {empty}.
```

If the two active order-three providers choose different classes, they cover
sixteen sheets.  The two high-order ratios are `{1,h}`, with
`h in {5,8,12}` for `D=8` and `h in {5,8}` for `D=24`.  A literal
`4 x 4` (`D=8`) or `8 x 8` (`D=24`) fibre lemma proves, for every `j`,

```text
max |(M_1 union M_h) intersect S_j| = 3.
```

Hence a distinct-class choice covers at most `16+3=19`.  A same-class choice
covers at most `8+(6+3)=17` for `D=8` and `8+(4+4)=16` for `D=24`.
The bound nineteen is attained.  The obstruction is therefore exact overlap
debt in the mod-three sheet nerve, not scalar slack.

The thirty missing-owner pairs form a useful alternate tournament carrier.
As a multiplicative Cayley multigraph on `F_13^*`, exponent steps `+/-3`
have multiplicity two (orders eight and twenty-four) and the antipodal step
`6` has multiplicity one (order eight).  Every vertex has weighted degree
five; the underlying simple graph is three disjoint `K_4`s, one on each cubic
coset.  This graph retains the character class and the exact pair of failed
proof obligations.  The transitive summary tournament below forgets both.

## Independent certificates and carrier audit

The exact primary checks every algebraic CRT representative against literal
search, checks every local mask cardinality against an independent period
formula, audits the prime-power hereditary grammar against all six
leave-one-out lcms, and enumerates all support/order capacity rows without a
symmetry quotient.  Python normal and `-O` runs are byte-identical.  The
primary artifact hashes are

```text
primary source  c3d20203dea9c36396db4a9975759b3d65b101c0aaf5a2c053d1370669861db1
primary output  3bb20b04576bcdb293ef474c4a083b43e8ab2e228ce87107b7bae36c02311ee4
```

The referee was designed and its invariants were frozen before its author read
the primary source or output.  It generates divisors and exact-order residues
from gcd/lcm definitions, finds every CRT base by bounded literal search,
proves the exact owner-to-ratio cyclic sheet gauge exhaustively, and joins
distinct unit-mask images with sorted-vector union DP in both provider orders.
It audits multiplication covariance but does not delete any labelled row.

The referee independently reproduces the `66,984 / 854 / 202` scalar census,
the complete feasible-owner and maximum-union histograms, all 20,302 owner
maximum vectors, 101,961,528 total reachable masks in 674 bank-size bins, and
the maximum bank size 7,728.  It additionally serializes all thirty extremal
rows, proves the cubic classification and three-class nerve lemma above, and
finds multiplication-orbit sizes `12,12,6`.  Its internal exact fingerprints
include

```text
literal CRT bases FNV64       bef281a507266cdb
literal owner masks FNV64     b700276da083bd75
cubic three-class nerve FNV64 e20047236f3ef6ab
scalar bank FNV64             c68102afc4f81501
canonical owner bank FNV64    cf86c2c61edbacb9
thirty extremal rows FNV64    d40b8aa5ee61cb18
```

Optimized, unoptimized, and address/undefined-sanitized C++20 builds complete
with byte-identical stdout.  The frozen referee hashes are

```text
referee source  e8e1bc41e451957180d61967af634ba0584a2969f622fea3fab08f773a901446
referee output  2702e91969a6aaf2387de8941004479094f330066c98313585dd5911ec34f95e
```

Owner obligations are the terminal-faithful tournament vertices.  The primary
uses ordered observable `(feasible,max-union,capacity)`; the referee strengthens
it to `(feasible,max-union,capacity,reachable-count,maximum-mask-count)`.
Lexicographic switching with the coordinate order as tie Hamiltonian path
produces a transitive tournament in all 66,984 rows, with score word
`0,...,5`, no directed triangle, six singleton SCCs, and one Hamiltonian path.
This is diagnostic telemetry: either total-order gauge forgets the absolute
24-sheet threshold, exact masks, and unit-witness incidence.  The cubic
missing-owner multigraph retains substantially more structure.  Provider,
divisor, residue, isolated-sheet, and wall-event vertices lose shared-unit
incidence even earlier unless equipped with an incidence sidecar.

This theorem concerns only the primitive proper AP-centred common-scale-
twenty-four Hamming-six face.  It does not cover scale 25 or higher, H5
ramification, non-AP/deep sheets, or global sporadic emptiness.
