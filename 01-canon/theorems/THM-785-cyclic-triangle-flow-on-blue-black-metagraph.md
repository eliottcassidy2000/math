---
id: THM-785
title: Cyclic-triangle transitivity flow on the blue/black merged metagraph
status: PROVED general score, flux, blue-symmetry, line-distribution, and category-incidence laws; VERIFIED FINITE oriented quotient flow n=3..7
source: codex-2026-07-14-S9
depends_on:
  - HYP-6825
  - THM-781
related: [THM-550, THM-646, THM-773, THM-778, THM-784, THM-787]
artifacts:
  - 04-computation/merged_metagraph_transitivity_flow_codex_S9.py
  - 05-knowledge/results/merged_metagraph_transitivity_flow_codex_S9.out
  - 05-knowledge/results/merged_metagraph_transitivity_flow_codex_S9.json
---

# THM-785 — Cyclic-triangle flow on the blue/black merged metagraph

The merged tournament nodes have a natural horizontal coordinate that is not
an arbitrary drawing convention.  It is the number `C3(T)` of cyclic
triangles.  It runs from the unique transitive class to exactly the
regular/near-regular score locus, survives tournament isomorphism and
converse, and has a closed flux law on every blue or black complement line.

The resulting picture has two different kinds of symmetry:

1. blue grid symmetry kills the endpoint defect and leaves a one-coordinate,
   centrally symmetric binomial law;
2. the whole black labelled ensemble is still centrally symmetric, but it has
   a second endpoint-defect coordinate.  After quotienting to iso-class
   categories and orienting toward larger `C3`, this hidden coordinate produces
   a measurable categorical left/right imbalance.

The second point is essential.  Black asymmetry is not failure of the
complement involution; it is asymmetry of an oriented quotient flow.

## 1. The exact transitivity spectrum

Let `d_i` be the outdegrees of an `n`-vertex tournament `T` and put
`mu=(n-1)/2`.  Then

```text
sum_i (d_i-mu)^2 = n(n^2-1)/12 - 2 C3(T).                 (1)
```

Equivalently, for the integral energy

```text
E4(T) := sum_i (2d_i-(n-1))^2,
E4(T) = n(n^2-1)/3 - 8 C3(T).                            (2)
```

Consequences:

- `C3=0` if and only if the tournament is transitive, so there is one root
  merged node.
- Moving right in `C3` is exactly moving down in score variance.
- The maximum locus is precisely the most evenly distributed score locus:

```text
n odd:  ((n-1)/2,...,(n-1)/2),       C3_max=n(n^2-1)/24;
n even: (n/2 repeated n/2 times,
         n/2-1 repeated n/2 times),   C3_max=n(n^2-4)/24. (3)
```

Thus the examples `(2,2,2,2,2)`, `(3,3,3,3,4,4,4,4)`, and
`(1,1,2,2)` are literally the right endpoint of this spectrum.

Proof of (1): every transitive triple has a unique source, while a cyclic
triple has none, so the number of transitive triples is
`sum_i binom(d_i,2)`.  Subtract from `binom(n,3)` and use
`sum_i d_i=n(n-1)/2`.  The integer vectors with this fixed sum minimize the
quadratic energy exactly at (3); regular and near-regular tournaments realize
those vectors.

Converse can merge two distinct score sequences.  The correct node datum is
therefore the unordered score-sequence orbit under `d -> n-1-d`, not an
arbitrarily selected representative.  Formula (1) is constant on this orbit.

## 2. Complement-line flux and the missing black coordinate

Use the explorer convention with fixed directed path
`0 -> 1 -> ... -> n-1`.  A tiling mask `t` chooses the remaining
`m=binom(n-1,2)` arcs; `bar(t)` flips all of those arcs and preserves the path.
The labelled degrees transform as

```text
bar(d_0)     = n-d_0,
bar(d_(n-1)) = n-2-d_(n-1),
bar(d_i)     = n-1-d_i                 (0<i<n-1).
```

Substitution in (1) gives the exact line flux

```text
C3(bar(t))-C3(t) = d_0(t)-d_(n-1)(t)-1.                 (4)
```

Introduce the two endpoint coordinates

```text
a       = d_0-n/2,
epsilon = d_0+d_(n-1)-(n-1).
```

Then

```text
Delta C3 = 2a-epsilon,                                  (5)
(a,epsilon)(bar(t)) = (-a,-epsilon).                    (6)
```

Equation (6) is the exact underlying left/right symmetry.  The black sector
does not break it.  What black lines add is the transverse coordinate
`epsilon`, which permits the `C3` orientation and the pure/mixed/black phase
orientation to disagree after projection.

## 3. Blue reflection law and its closed binomial distribution

A blue tiling is fixed by reflection in the grid anti-diagonal.  Vertex
reversal is then a tournament anti-isomorphism, hence

```text
d_i+d_(n-1-i)=n-1 for every i,
epsilon=0,
Delta C3=2d_0-n.                                        (7)
```

In particular every blue step has parity `n`.  Odd `n` has no blue `C3` tie;
even `n` can have one.

More is true: the entire blue step distribution is closed form.  Let

```text
m = binom(n-1,2),  q=n-2,
f = floor((n-1)/2),  r=(m+f)/2.
```

There are `r` reflection-orbits of tiles, hence `2^r` blue masks and
`2^(r-1)` blue complement lines.  The `q` top-row orbit bits are independent.
If exactly `z` of them are wins for vertex zero, then `d_0=1+z` and the line
step is `|2z+2-n|`.  Taking one member of each pair `z < q-z`, the number of
such lines is

```text
binom(q,z) 2^(r-q).                                      (8)
```

When `z=q-z` (only the even-`n`, zero-step case), divide (8) by two.  This is
the precise blue left/right symmetry suggested by the drawings.

Two further consequences reconcile the concurrent THM-787 formulation.  The
largest blue step is exactly

```text
max |Delta C3|=n-2,
max |Delta E4|=8(n-2).                                  (8a)
```

Hence blue `|Delta E4|` is `8 (mod 16)` for odd `n` and `0 (mod 16)` for even
`n`; the parity and maximum laws conjectured in THM-787 are general theorems.
By THM-781 the transitive node has fibre size one, so it is incident to exactly
one complement line.  That line is blue and attains (8a), landing at
`C3=n-2`.  Its landing node is mixed in the exact `n=4..7` atlases; that last
phase assertion is not promoted here beyond the audited range.

## 4. Closed laws for all steps and for black defect

Put `p=n-3`.  Separate the endpoint bits into:

- `R in {0,1}`, the arc between vertices `0` and `n-1`;
- `X`, the number of the `p` left-endpoint wins against middle vertices;
- `Y`, the number of the `p` right-endpoint wins against middle vertices.

Then `X,Y` are independent binomials at the level of mask counts and

```text
Delta C3 = 2R-1+X-Y.                                    (9)
```

Define the number of masks with signed flux `s` by

```text
A_n(s) = 2^(m-(2n-5))
         sum_[2r-1+x-y=s] binom(p,x)binom(p,y).          (10)
```

The number of unoriented complement lines of absolute step `k` is `A_n(k)`
for `k>0` and `A_n(0)/2` for `k=0`.  Subtracting the blue counts in (8) gives
the black step distribution for every `n`.

The defect itself is even simpler: `epsilon=X+Y-p`.  For `e>0`, the number of
black lines with `|epsilon|=e` is

```text
D_n(e)=2^(m-2p) binom(2p,p+e).                           (11)
```

For zero defect it is

```text
D_n(0)=2^(m-2p-1) binom(2p,p) - 2^(r-1),                (12)
```

where the subtracted term is the complete blue locus.  Equations (8),
(10)--(12) were asserted against every line at `n=3..7` with zero failures.

At `n=7` they specialize to

```text
blue |Delta C3|:  1:160, 3:80, 5:16;
black |Delta C3|: 0:3584, 1:6112, 2:4096, 3:1776, 4:512, 5:48;
black |epsilon|:  0:4224, 1:7168, 2:3584, 3:1024, 4:128.
```

## 5. Exact categorical topology

Call a merged node `pure_blue`, `mixed`, or `pure_black` according as its
THM-781 tiling fibre contains only blue masks, both colours, or only black
masks.  Complement preserves grid symmetry.  Therefore:

```text
a blue line has no pure_black endpoint;
a black line has no pure_blue endpoint.                 (13)
```

Consequently every cross-category line lies in the path

```text
pure_blue --blue-- mixed --black-- pure_black,           (14)
```

with no pure-blue/pure-black shortcut.  Lines within one category and
self-lines are allowed.  Statement (14) is a theorem for every size; it is
topological, not a claim that increasing `C3` always points from left to right.

## 6. Exact oriented quotient census through n=7

Orient a non-tied projected line instance toward increasing `C3`.  The finite
atlas gives:

| n | nodes `(PB,M,PK)` | `C3_max` | balanced-node phases | monotone reach |
|---:|---:|---:|---|---:|
| 3 | `(2,0,0)` | 1 | `PB` | 2/2 |
| 4 | `(1,1,1)` | 2 | `M` | 2/3 |
| 5 | `(3,5,2)` | 5 | `PB` | 7/10 |
| 6 | `(2,10,22)` | 8 | `M,M,M,PK` | 22/34 |
| 7 | `(4,84,184)` | 14 | `M,M,M` | 212/272 |

Here “monotone” allows movement across a `C3` tie.  Requiring the colour word
to have the two-phase form `B* K*` reaches exactly the same node set at every
audited size.  Every balanced node is reached.  At `n=6`, three chosen paths
are `BBB` and the fourth is `BKK`; at `n=7`, all three are `BBBBBB`.  These are
finite facts, not yet a theorem for arbitrary `n`.

The oriented line-instance phase counts `(outward,inward,same,tie)` are

| n | blue | black |
|---:|---:|---:|
| 3 | `(0,0,1,0)` | `(0,0,0,0)` |
| 4 | `(1,0,0,1)` | `(0,2,0,0)` |
| 5 | `(2,3,3,0)` | `(0,12,4,8)` |
| 6 | `(2,0,18,12)` | `(38,78,256,108)` |
| 7 | `(6,0,250,0)` | `(1254,2798,8492,3584)` |

At the `n=7` blue boundary all six instances and all six projected supports
point `pure_blue -> mixed`.  At the black boundary,

```text
mixed -> pure_black: 1254 instances / 522 supports;
pure_black -> mixed: 2798 instances / 1075 supports;
boundary C3 ties:      992 instances / 297 supports.
```

Thus the reverse black drift is exactly `1544` line instances and `553`
supports; the reverse/forward ratios are about `2.231` and `2.059`.  This is
the requested quotient-level left/right imbalance.  It coexists with the
exact labelled central symmetry (6).

The black current across the `C3` cut `k | k+1` at `n=7`, for
`k=1,...,13`, is

```text
4, 26, 86, 294, 626, 1520, 2440, 3512, 4256, 4290,
3430, 1366, 70.                                         (14a)
```

It rises by three orders of magnitude from the transitive end, peaks at the
late-middle cuts `9|10` and `10|11`, and falls again at the balanced endpoint.
The two largest individual depth-pair packets are `C3 10->11` with `1504`
lines and `C3 9->10` with `1206`, agreeing with THM-787's energy-axis census.

Local one-tile depth and blue/black root distance should not replace `C3`.
At `n=7`, their Pearson correlations with `C3` are respectively `0.898915843`
and `0.465133288`; the latter has 4,587 discordant unordered node pairs even
after removing coordinate ties.  The spectrum and the metagraph position are
genuinely different axes.

## 7. An objective flow address and ordering

The result JSON gives every node the lexicographic address

```text
(C3,
 phase rank pure_blue < mixed < pure_black,
 rooted blue/black line distance,
 least rooted blue/black path word,
 Landau-slack orbit,
 score-sequence orbit,
 stable weighted blue/black WL colour,
exact HYP-6825 rank).                                  (15)
```

Sorting (15) produces `flow_rank`.  The first coordinate is the exact
transitivity spectrum; the next three place a node in the blue-to-black
geometry; the score/slack orbit resolves majorization shape inside a `C3`
layer; weighted line refinement records recursive neighbourhood position;
and the old canonical address is an admitted exact final tie-breaker.  Landau
slack is a readable transform of the score orbit, not falsely claimed as an
independent invariant.

Together with THM-781, this supplies both directions requested by the
explorer programme:

```text
tiling -> merged node -> flow address/rank,
flow node -> exact HP(T)/Aut(T) tiling fibre.
```

## Preservation and challenged assumption

For the organizational question, the vertices are merged iso-class nodes.
For the flux theorem, however, the more faithful vertices are complement-line
instances or endpoint-bit packets `(R,X,Y)`: quotient nodes alone erase line
multiplicity, the sign of `epsilon`, and which fibre element realizes a line.
This explicitly challenges the assumption that tournament vertices must be
runners, arcs, or even iso classes.

The quotient preserves tournament transitivity depth and fibre membership.
It does **not** preserve the LRC lonely-time predicate, endpoint owners,
metric gaps, threshold side, wall chronology, Euclidean transport, or sheet
carry.  THM-773 and THM-778 still require those data as stalks over the node.

The carrier Tournament Analysis in the computation treats position carriers,
not runners, as vertices.  Its pairwise observable is the number of unordered
merged-node pairs separated by a carrier; retention and economy gauges have
25 edge flips at `n=7`.  This confirms that category, `C3`, score, local depth,
line distance, and exact address are competing controlled-forgetting
quotients, not interchangeable names for one order.

## Open frontier

The next structural questions are now sharp:

1. Does every balanced node at every `n` admit a nondecreasing `C3` path with
   colour word `B* K*`?
2. Does two-phase reach always equal unrestricted nondecreasing reach?
3. What invariant predicts the sign and size of black categorical drift after
   the exact symmetric defect law (11) is conditioned on iso-class phase?
4. At `n=8`, do the new balanced score nodes remain reachable, and does the
   blue boundary stay outward while the black boundary retains reverse drift?
5. Can signed defect and line-orbit data replace the canonical-code fallback
   in (15), or exhibit the first unavoidable twins?

These are no longer questions about how to draw the metagraph.  They ask how
an exact symmetric line measure disintegrates over asymmetric quotient fibres.

## Reconciliation with the concurrent THM-787 census

THM-787 independently found the same energy axis `E4`, blue parity spectra,
unique transitive pipe, and black-sea concentration while this theorem was in
progress.  The normalizations are exactly `Delta E4=-8 Delta C3`.  THM-785
adds the general flux identity, proves THM-787's conjectural blue parity/max
law, gives closed distributions for all sizes, works on converse-merged nodes,
and resolves directional category flow and monotone reach.

One bookkeeping distinction prevents a false discrepancy: THM-787's phase
histograms count **unmerged tournament classes**, whereas the table in §6
counts **converse-merged nodes**.  At `n=7`, for example, its 368 pure-black
classes form the 184 pure-black nodes reported here.  Line-instance totals and
energy/`C3` spectra agree exactly.
