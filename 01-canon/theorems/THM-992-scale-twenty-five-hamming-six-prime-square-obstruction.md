---
id: THM-992
title: Scale-twenty-five Hamming-six prime-square obstruction
status: CLAIMED STRUCTURAL FROM A COMPLETE SCRATCH RECONSTRUCTION — exactly 36 scalar rows survive, in three multiplication orbits, and a five-coset overlap argument proves every owner misses at least three sheets; a frozen primary and independent replay are in progress
source: codex-2026-07-17-S66 scale-twenty-five continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990]
related: [THM-978, THM-980, THM-981, HYP-6820]
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
each quadratic class, eliminating `c25=2,3`.  If `c25=5`, the simultaneous
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

Partition the 25 sheets into their five residue classes `Q_0,...,Q_4`
modulo five.  An order-five mask is one complete five-sheet `Q_j`.  At self
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
so the total is at most `10+(15-4)=21`; if the cosets coincide, the stronger
bound is twenty.  Thus every owner misses at least three of the 25 sheets.
No global unit word exists.

The complete literal DP independently sharpens the maxima to 72 order-five
owner rows of size 22 and 144 order-25 owner rows of size 21; all 36 contexts
have zero feasible owners.  The structural coset proof, not that DP, is the
logical obstruction.  Promotion requires a frozen self-auditing primary and
an independently structured replay.  This theorem does not cover `c>=26`,
H5 ramification, non-AP/deep sheets, or global sporadic emptiness.
