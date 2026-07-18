---
id: THM-993
title: Scale-twenty-seven Hamming-six owner-feasibility deficit
status: PROVED FINITE-EXACT + INDEPENDENT STRUCTURAL REPLAY — all 450 scalar survivors have at most four owner-local feasible projections; a standard-library primary and an independently developed NumPy-batched nested-fibre certificate agree on the scalar bank and every exact headline census
source: codex-2026-07-17-S66 scale-twenty-seven continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990, THM-992]
related: [THM-978, THM-980, THM-981, THM-994, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_seven_hamming_six_owner_deficit_codex_c27.py
  - 05-knowledge/results/lrc13_scale_twenty_seven_hamming_six_owner_deficit_codex_c27.out
  - 04-computation/lrc13_scale_twenty_seven_hamming_six_nested_fibre_obstruction_codex_c27.py
  - 05-knowledge/results/lrc13_scale_twenty_seven_hamming_six_nested_fibre_obstruction_codex_c27.out
---

# THM-993 — scale twenty-seven has at least two impossible owners

Assume the primitive proper AP-centred common-scale Hamming-six reduction of
THM-860, now with common scale `c=27`.  The effective orders are

```text
1, 3, 9, 27.
```

Hereditary leave-one-out lcm is equivalent to having at least two order-27
coordinates.  Hence there are 1,909 labelled order words.  With the literal
unit fibres retained, those words carry 380,511,756 state words per support,
or

```text
924 * 380,511,756 = 351,592,862,544
```

unquotiented labelled state contexts.

## Scalar collapse

For a provider/owner ratio `r in F_13^*`, the exact mask cardinalities are

```text
D=1  : 27 if r=1, else 0;
D=3  :  9 if r in {1,4,5,8,9}, else 0;
D=9  :  6 if r in {1,2,5,8,11}, else 3;
D=27 :  5 if r=1, else 4.
```

The frozen primary checks these formulas against every literal CRT mask.  It
then tests all `924*1,909=1,763,916` support/order contexts without quotienting
by multiplication.  Exactly 450 rows on 84 supports have scalar capacity at
least 27 at all six owners.  In multiplicity coordinates
`(#D1,#D3,#D9,#D27)`, their complete histogram is

```text
(0,0,3,3):  12    (0,0,4,2):  18
(0,2,1,3):  60    (0,2,2,2): 294
(0,3,0,3):  12    (0,3,1,2):  36
(0,4,0,2):  18.
```

Thus order one is already absent from every scalar survivor.  The 450 rows
form 36 free multiplication orbits of size twelve and three stabilized orbits
of size six; this orbit census is telemetry and is not used as a quotient.

## Exact owner-local deficit

For every scalar row and each of its six owners, the primary reconstructs the
literal CRT masks and forms the complete reachable union bank, deduplicating
only equal bit masks.  Across all 2,700 owner projections the exact maximum
union histogram is

```text
20:120, 21:336, 22:192, 23:336,
24:528, 25:432, 26:588, 27:168.
```

The resulting context census is

```text
0 feasible owners: 336 rows,
1 feasible owner :  96 rows,
4 feasible owners:  18 rows.
```

In particular, no scalar row has six feasible owner projections.  A global
unit word would project to a covering word in every owner bank, so this local
deficit is terminal for the scale-27 common-sheet H6 face.  The computation
visits 13,598,160 distinct reachable masks counted in their owner banks; there
are 51 bank-size bins and the largest bank has 128,880 masks.  There are 319
distinct maximum-union vectors.  The frozen digests are

```text
primary source          e1a08766308e8c98bce38630310fc597b75643e394dfda2ae8daf461dd464089
primary output          c1516ca0021c2cb95a3b4049fb41f6c4d3ba6507d74eff96a6e4d7601c7a0b5e
literal mask table      1f6007a1f21d0d3f0c4382634ad27c005139e175e67afbb7f7cafe130533c249
hereditary order bank   0eec286f7d0032bcbbed2b4ece0b62638562301484a907940fab374109c46df4
scalar bank             141ec6d8c551c2e0ebac31dc102d3f38ad257be268dfe24403835afa7613dc05
owner summaries         6471aa37d9a0d6a630a830d36865d1449c7d7361d068041a60e7482425bed085
mod-3/mod-9 signatures  e7e460b659cc589e88a47bf40149d3ffbe8e50cfdc004187ae5cc32e2fb63073
```

The primary source is standard-library only.  It audits algebraic CRT representatives
against literal search, the prime-power hereditary grammar against all six
leave-one-out lcms, mask sizes against an independent period count, and full
coverage against maximum union.  THM-994 supplies the independently developed
replay: it batches the scalar bank through NumPy rather than nested Python
loops, hashes every sorted reachable mask rather than only owner summaries,
and separately proves a sound nested-fibre upper relaxation.  The two programs
agree on the exact scalar digest
`141ec6d8c551c2e0ebac31dc102d3f38ad257be268dfe24403835afa7613dc05`,
all seven multiplicity counts, all 450 labelled rows, the 2,700 owner maxima,
the feasible-owner census `0:336,1:96,4:18`, all 51 bank-size bins, and the
13,598,160 reachable-mask incidence total.  This promotes the finite-exact
claim; a short conceptual classification of the 18 four-feasible rows remains
desirable but is not needed for terminal emptiness.

## Tournament and alternate carriers

On owner vertices, orient the exact ordered observable

```text
(feasible, maximum union, scalar capacity, reachable-bank size),
```

breaking an exact tie by owner coordinate.  Every one of the 450 completed
tournaments is transitive: score sequence `(0,1,2,3,4,5)`, no directed
triangles, six singleton SCCs, and one Hamiltonian path.  The tournament is a
useful ranking of failed obligations, but it destroys the absolute
27-sheet threshold, overlap multiplicity, and the colours of the ramified
sheet classes.

The challenged vertex choice points to a more faithful object.  At scale 27,
the natural small vertices are the residue classes of the 27 sheets modulo
three and modulo nine, with provider masks as coloured hyperedges.  This
quotient retains the prime-power ramification signatures while discarding
irrelevant affine offsets.  The primary records 96 exact combined signature
types.  Owners retain the terminal deficit vector but forget why a sheet is
missed; individual sheets retain too much offset data; units alone forget
which owner-local ratio produced the mask.  A structural referee should work
on this coloured coset-cover system rather than infer meaning from the
transitive owner tournament.

This claim does not cover H5 ramification, non-AP-centred or deep-sheet
packets, lift-dependent shell-five certificates, or global sporadic
emptiness.  The next numerically possible primitive common scale is 28.
