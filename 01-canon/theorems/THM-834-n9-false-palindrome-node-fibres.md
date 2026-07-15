---
id: THM-834
title: The 58 n=9 false-palindrome pairs occupy 53 converse-merged tournament nodes
status: PROVED FINITE-EXACT FORWARD/INVERSE-FIBRE MAP + n8 BLACK-FLOW PLACEMENT
source: codex-2026-07-15-S13b
depends_on: [THM-810, THM-828, THM-832]
related: [THM-796, THM-809, THM-833, THM-838, HYP-3809, HYP-6880]
verification:
  - 04-computation/n9_false_palindrome_node_forward_map_codex_S13b.py
  - 05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.out
  - 05-knowledge/results/n9_false_palindrome_node_forward_map_codex_S13b.json
  - 04-computation/n9_false_palindrome_metagraph_flow_codex_S13b.py
  - 05-knowledge/results/n9_false_palindrome_metagraph_flow_codex_S13b.out
  - 05-knowledge/results/n9_false_palindrome_metagraph_flow_codex_S13b.json
---

# THM-834 — exact node fibres and black-flow placement

For a size-nine tiling `u`, let `T(u)` be its tournament in the fixed
Hamiltonian-path chart.  For any tournament `T`, let `c(T)` be its least
upper-triangle adjacency integer among every vertex relabelling.  Write

```text
nu(u) = {c(T(u)),c(T(u)^op)}                             (1)
```

as a sorted pair, with a repeated entry for a self-converse class.  On the
116 endpoints of THM-828's 58 false-palindrome pairs, `nu` has exactly 53
values, with fibre profile

```text
48 nodes x one pair + 5 nodes x two pairs.               (2)
```

All five doubled nodes are self-converse and all ten of their pairs lie in
the dominant defect sector `D=4c41818`.  Forty-nine of the 53 nodes are
self-converse.  On 54 pairs the two reflected tilings present the same
self-converse ordinary class; on four pairs they present two distinct
ordinary classes exchanged by converse.

The bare projection `nu` is therefore not injective on collision pairs.  It
becomes injective after retaining the complement-partner node.  Distinguish
the ordered and unordered versions

```text
P^->(u) = (nu(u),nu(u xor FULL)),
bar P(u) = {nu(u),nu(u xor FULL)}.                        (3)
```

The first coordinate of (3) has 53 values and the second has 54, but both the
ordered pair and its unordered endpoint version have 58 singleton values.
This is the exact finite bridge between THM-828's literal tiling fibres and
the merged-metagraph nodes.

## The two maps, including the inverse fibre

The forward function is completely canonical:

```text
tiling mask
  -> fixed-path labelled tournament
  -> least code under S9
  -> sorted pair with the least code of the converse.     (4)
```

The verifier uses a subset dynamic program.  If `v` is the first vertex of a
canonical order, its arcs to the remaining vertices occupy the low `k-1`
bits and the recursively canonical induced tournament occupies the higher
bits.  All minimizing tail orders are retained.  Thus this is exact
canonicalization, not a hash or graph-library heuristic.

Conversely, for an ordinary unlabelled class `C`, every directed Hamiltonian
path gives a presentation in the fixed-path tiling chart, and two paths give
the same tiling exactly when an automorphism identifies them.  Hence

```text
F(C) = {u : c(T(u))=C},    |F(C)| = H(C)/|Aut(C)|.        (5)
```

For a merged node `m={C,C^op}` its full tiling fibre is

```text
F(m)=F(C) union F(C^op),
|F(m)| = H/|Aut|                 if C=C^op,
         2H/|Aut|                otherwise.              (6)
```

Equations (4)--(6) are the requested bidirectional dictionary.  In practice
one obtains the set in (5) either by enumerating Hamiltonian paths modulo
automorphisms or by scanning the exact tiling atlas.  The JSON witness table
also gives the restricted inverse

```text
m -> F(m) intersect {the 116 THM-828 endpoints}.          (7)
```

The five nontrivial fibres of (7) have merged keys

```text
0000a68d4, 0008840b2, 000884195, 0008848a3, 0008c51f3.   (8)
```

Each contains four endpoints, hence two distinct false-palindrome pairs.
Their ordinary-class invariants `(H,|Aut|,H/|Aut|)` are respectively

```text
(1581,1,1581), (1903,1,1903), (1711,1,1711),
(1931,1,1931), (1639,1,1639).                            (9)
```

This proves concretely that a node is an unmarked tournament, whereas a
collision pair remembers two marked Hamiltonian-path presentations.

## A sufficient two-marginal repair

On one chosen presentation from each of the 58 pairs, the two coordinate
fibres in (3) are

```text
nu(u):            53 values, multiplicity {1:48,2:5},
nu(u xor FULL):   54 values, multiplicity {1:51,2:2,3:1}.
```

Neither marginal is sufficient.  The ordered `P^->` coordinate has 58
values; forgetting which marginal is first still leaves 58 values, and none
is a loop.  Reflection swaps the literal endpoints coherently and creates no
failure.  THM-838 further proves that all 58 coupled cells remain distinct
under the centered-CF copy `X_9 -> X_10`.

This does not say that `P^->` or `bar P` is a global invariant of every
metagraph line or that either always descends under continued-fraction
actions.  It is an exact
two-marginal repair on this 58-pair bank, with no claim of minimality among
all possible repairs.  The 53-versus-54 counts also depend on
which complement chart is called direct; the coupled coordinate does not.

## Placement through the six size-eight face presentations

The six face endpoints of every collision pair give 348 distinct size-eight
tilings.  In the exact `2^21` size-eight atlas they occupy

```text
155 merged nodes = 153 pure-black + 2 mixed nodes.        (10)
```

Project each face tiling to the merged node of its complement partner.  The
174 marked-face complement-line selections (three per witness) collapse to
87 distinct projected complement-line node-pair codes.  Every code occurs in
the black code set and not in the blue code set.  These are projected codes,
not 87 asserted literal lines.  Their occurrence multiplicities are

```text
58 x 1 + 29 x 4 = 174,                                   (11)
```

and their endpoint categories are 172 pure-black--pure-black occurrences
and two mixed--pure-black occurrences.  Along THM-810's score axis `x`, the
absolute jumps are

```text
|Delta x| : 0  8  16  24  32  48
count     :32 52  72   4  13   1.                        (12)
```

Thus the selected node-pair codes have no internal blue left-right rail:
after projection, their complete incidence is exclusively black.  The
imbalance is also quantitative—142 of the 174 marked-face occurrences change
the transitivity axis, while 32 are level.  This is a statement about the
selected face bank, not a claim that the global size-eight metagraph lacks
blue symmetry.

## An objective sector order and its obstruction

For each realized `D` sector, average over its supporting face-node
presentations.  Order first by decreasing mean score-axis coordinate, then
by increasing mean local metagraph depth, increasing mean `H/|Aut|`, and the
four-bit cube coordinate as a final exact tie-breaker.  The resulting order
is

```text
11b4600, 5df5e18, 4c41818, 08c2c0c, 4483414, 4511092,
18e4e8a, 1976a0c, 095088a, 0192486, 5c67a9e.             (13)
```

It runs from more transitive to more distributed support in the precise
THM-810 coordinate, but it is an order of decorated sectors, not a total
order of all isomorphism classes.  The JSON gives for every sector its exact
axis mean, local-depth mean, `H/|Aut|` mean, node count, black-line count,
category histogram, and jump histogram.  A compatible proposed total order
on individual nodes is

```text
(-x, score sequence, colour stratum, local depth,
 H/|Aut|, canonical merged key).                         (14)
```

Incoming THM-833 supplies the dynamical centre of this order.  For the
size-eight single-arc-flip walk the exact OU equilibrium is

```text
x_* = n(n-1) = 56.
```

Exactly the first two sector barycentres in (13), `220/3` and `68`, lie above
56; the other nine, beginning with `1660/39`, lie below it, and none is level.
Thus “toward distributed” is not one global drift direction.  Above 56 the
mean flip drift decreases `x`; below 56 it increases `x`, back toward the
random-tournament centre.  The maximally distributed end is an entropic tail,
not the stationary endpoint of the flow.

The defect cube does **not** supply the missing metagraph geometry.  Among
the 14 cube-adjacent sector pairs, all 14 have a black bridge and none a blue
bridge; their local metagraph distances are `0:1, 1:11, 2:2`.  Among the 41
non-cube pairs, 38 also have black bridges, none has a blue bridge, and the
distances are `0:7, 1:21, 2:12, 3:1`.  The mean absolute axis-barycentre gaps
are exactly

```text
cube:    (1724/91)   = 18.945054...,
noncube: (38720/1599)= 24.215134....                     (15)
```

Orienting cube edges down the axis produces three local maxima and three
local minima; the top sector reaches only two of eleven sectors, including
itself, by directed cube paths.  Hence a monotone transitive-to-balanced flow cannot be recovered
from defect-coordinate adjacency alone.  The black metagraph bridges are
the missing non-cubical incidence.

## Tournament Analysis and preservation boundary

For every marked size-eight face tournament, the tournament vertices are its
eight path positions, the pairwise observable is arc orientation, converse
reverses every edge, and the fixed path `8->...->1` is a Hamiltonian path.
Aggregated over the 348 selected endpoint occurrences, the audit records
4,643 climbs, 2,548 level flips, and 2,553 descents for full-arc flips; in the
marked-wiggly relation it records 3,133 climbs, 2,027 levels, and 2,148
descents.  The 87 projected codes are then disintegrated by sector and by
(12).  The directed defect cube is not a tournament and is therefore
reported by its maxima, minima, reachability, and bridge distances rather
than assigned spurious tournament fingerprints.

The challenged assumption is that tournament vertices must be runners—or
even the tilings themselves.  Here the useful vertices alternate among
unmarked tournament classes, marked Hamiltonian paths, complement-line node
pairs, and defect sectors.  The quotient (4) preserves tournament
isomorphism and converse, while destroying the observer path and collision
pairing; (3) restores exactly the latter information on this bank.  No map
in this theorem preserves LRC gaps, owners, walls, or loneliness, and no
statement here resolves the fourteen-runner case. ∎
