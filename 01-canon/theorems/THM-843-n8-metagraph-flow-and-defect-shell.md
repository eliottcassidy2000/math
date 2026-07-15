---
id: THM-843
title: The n=8 two-phase flow reaches every balanced node, and the n=9 defect support has an almost-global black radius-one shell
status: PROVED FINITE-EXACT FULL n=8 MERGED-FLOW + DEFECT-SHELL CLASSIFICATION
source: codex-2026-07-15-S13e
depends_on: [THM-781, THM-785, THM-790, THM-834]
related: [THM-811, THM-814, THM-830, THM-842, HYP-6825, HYP-6880]
verification:
  - 04-computation/n8_metagraph_flow_defect_shell_codex_S13e.py
  - 05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.out
  - 05-knowledge/results/n8_metagraph_flow_defect_shell_codex_S13e.json
---

# THM-843 — n=8 flow completeness and the defect shell

In the exact `2^21` size-eight atlas, the converse-merged metagraph has

```text
3,528 nodes = 3 pure blue + 173 mixed + 3,352 pure black.       (1)
```

Its complement-line layer contains `2,048` blue and `1,046,528` black
literal lines.  After projection these occupy respectively `1,794` and
`436,580` node-pair supports, including loops; exactly `158` projected
supports carry both colours.  Literal lines and projected supports are kept
separate throughout.

## 1. The n=8 transitivity-flow question

Orient a non-tied support toward increasing cyclic-triangle count `C3`, or
equivalently decreasing score energy

```text
x=sum_i(2d_i-7)^2=168-8 C3.                              (2)
```

Starting at the unique transitive node, unrestricted nondecreasing-`C3`
reach and reach by a two-phase colour word

```text
blue* black*                                                   (3)
```

are the identical set of `3,069` nodes.  Thus the equality observed by
THM-785 through `n=7` survives at `n=8`, although monotone reach is not the
whole metagraph.

The decisive endpoint statement is stronger:

```text
C3-maximal nodes                         50
phase split                              15 mixed + 35 pure black
unrestricted monotone reach              50/50
blue* black* monotone reach               50/50.          (4)
```

Indeed every node with `C3>=15` is reached.  This settles THM-785's finite
`n=8` balanced-locus question affirmatively.  It does not prove the all-size
claim.

The same calculation also gives an exact nonmonotonicity boundary for the
THM-834 support.  Of its 155 nodes, 153 are reached monotonically.  The two
exceptions are

```text
rank 1008: pure black, C3=11, x=80,
rank 3512: mixed,      C3= 6, x=120.                     (5)
```

They are not disconnected; each requires a final rebound in `x`.  Thus even
the most-transitive defect sector contains marked presentations that no
globally monotone axis order can place correctly.

## 2. Blue symmetry and black imbalance at n=8

All three literal and all three projected blue boundary supports point from
pure blue to mixed as `C3` increases.  At the mixed--pure-black boundary the
black census is

```text
                         outward       reverse       C3 tie
literal lines              16,118        24,370        12,584
projected supports           7,298        10,634         4,836. (6)
```

So raw flow again has a left/right reversal.  It is a fibre-volume effect,
not a reversal of source-normalized tendency.  Dividing by the available
black endpoint mass in each source phase gives

```text
mixed -> pure black     8059/27308  = .2951149846...,
pure black -> mixed     2437/203844 = .0119552207....    (7)
```

This is consistent with, and aggregated from, THM-814's sharper source-`q`
disintegration.  Equation (6) concerns quotient direction; grid reflection
still cancels signed endpoint current exactly.

## 3. Radius one around the 155 defect nodes

THM-834's selected bank itself consists of `153` pure-black and two mixed
nodes.  Its selected 87 complement-line supports are exclusively black, but
the full radius-one shell is much larger:

```text
colour  literal incident lines  projected supports  external nodes  all nodes
blue                         30                  28              25         27
black                   128,652              47,872           3,153      3,308
union                         -                   -            3,155      3,310. (8)
```

The black projected supports split into `1,956` nonloop internal supports,
`45,889` boundary supports, and `27` loops.  One incident projected support
has both colours.  Blue touches exactly the two mixed seed nodes and no
pure-black seed.  They are the only blue portals:

```text
rank 3411: C3=12, x=72,  23 literal / 21 projected blue neighbours,
rank 3512: C3= 6, x=120,  7 literal /  7 projected blue neighbours. (9)
```

All 28 external blue incidences land in the mixed phase, on 25 distinct
nodes.  Thus global blue symmetry re-enters through two sharply localized
gateways rather than through a diffuse rail.

Only 218 of 3,528 nodes lie outside the combined radius-one shell:

```text
outside phase split       3 pure blue + 15 mixed + 200 pure black,
combined distance profile 0:155, 1:3155, 2:215, 3:3.            (10)
```

Every omitted node has `x>=56`, the n=8 Ornstein--Uhlenbeck centre from
THM-833.  Consequently (8) contains the entire distributed tail `x<56` and
all 50 balanced nodes.  The collision bank is black-only internally but is
not a small or closed black island: one black step sees 93.8% of all merged
nodes, and three coloured steps see every node.

## 4. Curvature refines an edge but does not identify it

Mark every literal black endpoint presentation whose node lies in the seed
support and measure the change in `C3` to its complement-line mate.  The exact
curvature disintegration is

```text
 q   marked ends     lower     tie    higher   sum Delta C3       mean
 0        13,828      4,120   1,210     8,498          3,642    1821/6914
 1        43,480     21,380   8,594    13,506        -20,670   -2067/4348
 2        46,976     30,618  13,214     3,144        -39,292  -9823/11744
 3        24,954     10,946   9,424     4,584        -10,420  -5210/12477
 4         5,910        436     826     4,648          6,220      622/591
 5           696          0       0       696          1,536        64/29
 6            70          0       0        70            280            4. (11)
```

The aggregate black drift is `-29352/67957`, whereas the 30 marked blue
endpoints have mean drift exactly `+2`.  The black sign is not monotone in
curvature: it is positive at `q=0`, negative at `q=1,2,3`, and positive again
at `q=4,5,6`.  This is a local disintegration, not a contradiction to
THM-814's source-normalized phase-boundary theorem; the conditioning events
are different.

For each incident black literal line, also attach THM-811's endpoint curvature
`q` to its projected endpoints and retain `|epsilon|`.  On the `128,652`
lines the bare node pair has `47,872` cells.  The attached carrier

```text
(projected node pair, endpoint q values, |epsilon|)      (12)
```

has `58,937` cells with fibre profile

```text
54,247 x 2, 4,124 x 4, 452 x 6, 97 x 8, 15 x 10, 2 x 12. (13)
```

Thus curvature resolves 11,065 additional projected collisions, but no
literal line: every fibre remains even, as reflection predicts.  This is the
local form of THM-811/814's preservation boundary.  Curvature is a genuine
black-edge coordinate, not a substitute for mirror position or the full
positional `B2/B3` address.

The verifier audits the curvature formula independently on all `2^11`
endpoint-leg states by counting cyclic triples containing both fixed-path
endpoints.

## 5. Convention correction inherited by THM-834

The full-shell audit exposed one root-convention error in THM-834's depth
sidecar.  That verifier uses the legacy bit convention in which `FULL`, not
zero, is the transitive fixed-path tiling.  Re-rooting at
`atlas[FULL]=3527` changes the local depth of 114/155 selected nodes and the
global maximum from seven to eight.  It does **not** change any node fibre,
selected line, category, score axis, pairwise bridge, defect-cube statement,
or the sector order, whose primary axis means are all distinct.  The THM-834
artifact and theorem have been corrected and replayed.

## Tournament Analysis and preservation boundary

The primary tournament vertices are the eight marked path positions before
quotienting and the 3,528 merged classes after quotienting.  A secondary
five-vertex information tournament compares phase pair, `C3` pair, node pair,
node pair plus curvature, and literal line.  Raw separated-pair retention has
the transitive path

```text
literal > node+curvature > node > C3 > phase,            (14)
```

while retention per log-cell moves `C3` to the front and creates five edge
flips.  Both gauges are transitive; their disagreement quantifies economy,
not a cycle in the underlying merged metagraph.

The quotient preserves tournament isomorphism, converse, line colour and
multiplicity, `C3`, endpoint curvature/current, and membership in the marked
THM-834 support.  It destroys unrecorded size-nine gluing and all LRC metric
gaps, owners, wall chronology, loneliness predicates, and continued-fraction
carry.  No statement here resolves LRC(14).  The challenged assumption is
that a conspicuously black selected bank should be studied as a closed
subgraph; its exact shell proves the opposite. ∎
