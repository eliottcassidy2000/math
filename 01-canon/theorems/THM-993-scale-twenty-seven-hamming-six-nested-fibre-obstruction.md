---
id: THM-993
title: Scale-twenty-seven Hamming-six nested-fibre obstruction
status: CLAIMED + FROZEN PRIMARY — a complete deterministic primary checks all 351,592,862,544 labelled state contexts through exact weighted grammar, reduces to 450 scalar rows, and gives a terminal two-stage Z/3 then Z/9 fibre relaxation; normal/optimized Python is byte-stable, while independent replay remains pending
source: codex-2026-07-17-S66 scale-twenty-seven continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990, THM-992]
related: [THM-963, THM-969, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_seven_hamming_six_nested_fibre_obstruction_codex_c27.py
  - 05-knowledge/results/lrc13_scale_twenty_seven_hamming_six_nested_fibre_obstruction_codex_c27.out
---

# THM-993 — scale twenty-seven has a terminal nested-fibre deficit

The primitive proper AP-centred common-scale-27 Hamming-six face is empty,
subject to independent replay of the frozen primary.

For `c=27`, the effective orders are `1,3,9,27`, with literal unit counts
`1,2,6,18`.  Hereditary leave-one-out lcm is equivalent to at least two
order-27 coordinates.  Hence there are 1,909 hereditary order words and

```text
27^6 - 9^6 - 6*18*9^5 = 380,511,756
```

literal state words per support, or

```text
924*380,511,756 = 351,592,862,544
```

unquotiented labelled state contexts.

## Exact scalar layer

For provider/owner ratio `r in F_13^*`, unit choice does not change mask
cardinality.  Literal interval counting gives

```text
D=1 : 27 at r=1, otherwise 0;
D=3 :  9 at r in {1,4,5,8,9}, otherwise 0;
D=9 :  6 at r in {1,2,5,8,11}, otherwise 3;
D=27:  5 at r=1, otherwise 4.
```

The primary checks those formulas against every literal mask.  Its complete
`924*1909=1,763,916` scalar order-context census has feasible-owner histogram

```text
0:11452, 1:391056, 2:849558, 3:418988,
4:83472, 5:8940, 6:450.
```

Exactly 450 rows on 84 supports reach capacity at all six owners.  None uses
order one.  Their multiplicity histogram in order `(1,3,9,27)` is

```text
(0,2,2,2):294, (0,2,1,3):60, (0,3,1,2):36,
(0,4,0,2):18,  (0,0,4,2):18, (0,0,3,3):12,
(0,3,0,3):12.
```

They represent 46,002,816 literal unit words.  Multiplication by
`F_13^*` splits the labelled contexts into three orbits of size six and 36
orbits of size twelve; the 84 supports split into two size-six and six
size-twelve orbits.  These orbit counts are telemetry only: the primary checks
all 450 rows and every owner.

## Nested-fibre relaxation lemma

The structural consumer is the following general finite inequality.  Let
`d|c`, let `A` be a set of anchor providers whose local masks are periodic
modulo `d`, and fix an owner.  For an anchor unit tuple `u_A`, write

```text
Q(u_A) = union_(i in A) M_i(u_i).
```

For any choices of all remaining units,

```text
|union_i M_i(u_i)|
 <= |Q(u_A)| + sum_(i notin A) |M_i(u_i) minus Q(u_A)|
 <= |Q(u_A)| + sum_(i notin A) max_e |M_i(e) minus Q(u_A)|.       (1)
```

Therefore

```text
U_d = max_(u_A) [ |Q(u_A)|
                  + sum_(i notin A) max_e |M_i(e) minus Q(u_A)| ]
```

is a sound upper bound on every literal union at that owner.  The relaxation
deliberately lets each nonanchor provider choose its best unit independently;
it can only enlarge the attainable union.  Thus `U_d<c` is terminal.

For an effective-order-`d` provider, shifting sheet `t` to `t+d` changes the
CRT expression by a multiple of `13d`.  Its mask is therefore the pullback of
a subset of `Z/dZ`, proving the periodicity premise rather than assuming it.

## The exact two-stage certificate

Among the 450 scalar rows, 420 contain an order-three provider.  Take all
order-three coordinates as anchors in (1), retaining the thick fibres of
`Z/27 -> Z/3`.  Across their 2,520 owner obligations, the exact anchor-union
banks have only one, two, or three elements:

```text
1:36, 2:2112, 3:372.
```

The remaining 30 rows have only orders nine and twenty-seven.  Take every
order-nine coordinate as an anchor and pass to `Z/27 -> Z/9`.  Their 180
owner obligations have anchor-union bank-size histogram

```text
45:24, 46:36, 68:24, 77:24, 134:72,
```

and exact upper-bound histogram

```text
21:12, 22:48, 25:96, 26:24.
```

Combining the two carriers, the number of owners whose relaxed bound still
reaches 27 has context histogram

```text
0:48, 1:96, 3:96, 4:210.
```

In particular every one of the 450 rows has at least two impossible owner
projections.  A global unit word would have to project feasibly at all six
owners, so the common-scale-27 face is empty.  This is an exact finite
structural certificate over all labelled scalar rows, not a sample and not an
orbit-quotiented search.

## Literal sharpness sidecar

The same primary separately performs immutable-set reachability at all 2,700
owner rows.  It constructs 13,598,160 reachable masks, with largest bank
128,880, and obtains

```text
feasible owners/context: 0:336, 1:96, 4:18;
maximum union: 20:120, 21:336, 22:192, 23:336,
               24:528, 25:432, 26:588, 27:168.
```

Some individual owners are literally feasible, but no context has more than
four feasible owners.  This DP checks sharpness and the loss ledger; inequality
(1), not the DP, is the theorem-bearing obstruction.

## Cayley, tournament, and carrier-loss audit

Switching the order-three ratio cardinality at its high value gives a directed
Cayley relation with 48 arcs, 36 symmetric edges, 12 reciprocal pairs, 16
directed triangles, and one SCC.  The corresponding order-nine relation has
the same first three counts, no directed triangle, and one SCC.  Its symmetric
switch is exactly the nonquadratic class, hence its undirected shadow is
`K6,6` on quadratic versus nonquadratic residues.  This is the bipartite dual
of THM-992's `K6 disjoint-union K6` prime-square shadow.  Neither Cayley
relation is a tournament.

For telemetry, use owner obligations as vertices and compare the exact keys

```text
(locally feasible, maximum union, scalar capacity, fibre bound)
```

lexicographically, breaking equality along coordinate order.  All 450
completed tournaments are transitive, with score multiset `{0,1,2,3,4,5}`,
no directed cycle, six singleton SCCs, and one Hamiltonian path.  The
tournament forgets the absolute threshold and cannot prove the theorem.

The `Z/3` quotient retains order-three thick fibres but forgets positions
inside each fibre, pair overlaps among the relaxed order-nine/order-27 needles,
and their shared unit constraints.  The `Z/9` quotient repairs the rows with no
order-three anchor while still forgetting point-level needle incidence.  This
loss is safe only because (1) is an upper relaxation.  The literal DP is the
sidecar that measures the discarded overlap.  Owners with their absolute
bounds are faithful; runners, isolated residues, raw sheets, Fano points,
chi-seven colours, or completed tournament orientations are not.

This nested `Z/27 -> Z/9 -> Z/3` carrier is the exact 3-adic
toothpick/Kakeya self-similarity: order-three masks are nine-sheet fibres,
order-nine masks are three-sheet fibres, and order-27 masks are point needles.

## Frozen primary

The primary reconstructs CRT bases algebraically and by literal search,
checks mask periodicity and the closed scalar formulas, scans every labelled
order context, evaluates every structural owner bound, and hashes every sorted
reachable-mask bank.  Normal and `python -O` executions reproduce the frozen
45-line output byte-for-byte.  Frozen SHA-256 values are

```text
Python primary source  d6bc42ac3f2d9609e37df0e11bfbc1eba7ffd6eabf1dd931c801c1a94fd959e1
Python primary output  0b1fc1b19dbbf8ea2e8295a1d3fb19f50144ee9cdf4b1fb6f980230bdfcfe294
```

Independent reconstruction is still required before promotion from claimed
to proved.  This theorem does not address H5 ramification, non-AP/deep sheets,
or global sporadic emptiness.  The next legal untreated AP-centred common scale
is `c=28`; scale 29 is already prime-excluded by THM-983.
