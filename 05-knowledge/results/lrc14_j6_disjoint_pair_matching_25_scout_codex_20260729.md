# Disjoint-pair matching scout on 25 marked suffixes

**Status: FINITE-EXACT SCOPED SCOUT.  NO UNIFORM SEVEN-BODY THEOREM.**

This scout refines THM-2897's parent certificate on the `25` scalar-hard
marked suffixes from the four-root THM-2895 battery.  On a carrier of mass
`h`, let `q_5` be the fifth allowed singleton order statistic, let `B_2` be
a global pair-union cap, and put

```text
T=h-q_5,                    L=T-B_2.
```

Give each pair edge `{u,v}` its intrinsic undirected weight

```text
U({u,v})=|C intersect (D_u union D_v)|.
```

If two vertex-disjoint edges have weights summing to at least `T`, each
edge has weight at least `L`: the other edge weighs at most `B_2`.  Write
`q_1` for the global singleton maximum and
`gamma=(99/70)r/7`.  When

```text
L>q_1+h/7,
```

every pair containing a label

```text
w >= ceil(gamma/(L-q_1-h/7))
```

has weight strictly below `L`.  Hence the `L`-heavy pair graph is finite.
Testing its maximum-weight two-edge matching is stronger than replacing
both disjoint edges by the same reusable global maximum `B_2`.

The exact census is

```text
scalar-hard branches                              25
finite-heavy-graph eligible                       13
q_5+2B_2 closed                                    3
disjoint-pair matching closed                      9
new closures beyond q_5+2B_2                       6
exact L-heavy edges                               59
base pair-cap evaluations                        242
heavy-edge evaluations                           434
literal/direct reconstruction controls            13
```

By root, the fields are
`(hard, eligible, q_5+2B_2, matching, new)`:

```text
(2,8,9,10,11,13,14):       (6,4,1,2,1)
(1,3,9,10,11,12,14):       (5,3,2,3,1)
(2,5,9,11,12,13,14):       (8,2,0,1,1)
(2,3,4,5,6,7,8):           (6,4,0,3,3)
```

The smallest positive finite-tail gap is

```text
574489/206158680
```

at `E=(2,8,9,10,11,13,14)`, rank `6`, apex `29`; this row needs the
largest exact tail entry `2611`.  The closest eligible matching failure
exceeds its target by

```text
241/360360
```

at `E=(2,3,4,5,6,7,8)`, rank `5`, apex `36`, with disjoint edges
`(20,44)` and `(33,78)`.

The complete output preserves all `25` marked prefixes, `q_1`, `q_5`,
`B_2`, finite levels, tail cutoffs, heavy-edge digests, and maximum
matchings.  Ordinary and optimized replays are byte-identical:

```text
04-computation/lrc14_j6_disjoint_pair_matching_25_scout_codex_20260729.py
SHA-256 ce0465fb1c17912e0dec76a96f18d6259aa35bb2cae55842c873a5eec1c2069b

05-knowledge/results/lrc14_j6_disjoint_pair_matching_25_scout_codex_20260729.out
SHA-256 45426624b4779abf2153838e2d8effb83b49d258299fec5fce921ad5c084ad1a

complete ledger
SHA-256 ccd5ac755132dc574925345018f25833aae54535fcae2b47a204203d663df129
```

The relation is an intrinsic symmetric weighted graph.  Its useful
operation is disjoint matching, so forcing orientations would discard the
theorem-bearing incompatibility.  The matching cap still forgets
cross-block overlap and is only sufficient.  The `12` ineligible rows and
the `4` eligible matching failures still require a tail-count split, a
child rank-selective cap, or THM-2895's literal H4 descent.
