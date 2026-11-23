# A square-sum cycle of 72 attained tournament counts

**FINITE-EXACT.** For every \(v\in[36]\), one can select a tournament with exactly \(100v-1\) directed Hamilton paths and reverse one specified arc to obtain exactly \(100v+1\). All 36 selections fit on ten vertices, and ten is the smallest common order for these prescribed heights. On the resulting set of 72 distinct heights, join two heights precisely when their sum is a square. This graph has 118 edges and **exactly 389 Hamilton cycles invariant under swapping the two heights at every midpoint simultaneously**. Consequently the total number of its Hamilton cycles is odd.

## A general lifting theorem

**PROVED.** Fix an even positive integer \(c\) and positive \(n\). Give \(v\in[n]\) the two odd labels
\[
                    f_+(v)=c^2v+1,\qquad f_-(v)=c^2v-1 .
\]
The labels are all distinct. Two labels on the same side have sum 2 modulo 4. Opposite labels at base vertices \(v,w\) have sum \(c^2(v+w)\). For distinct \(v,w\), an edge of the base square-sum graph therefore gives precisely its two lifts, and all cross-vertex edges arise this way. The two copies of a single vertex are joined precisely at \(v=2r^2\). These are the only edges fixed as sets by the involution exchanging the two copies.

When \(n\ge2\) is even, the sheet-invariant Hamilton cycles on the \(2n\) labels correspond bijectively to unoriented Hamilton paths on \([n]\) whose two endpoints are of the form \(2r^2\). Given such a path \(v_0,\ldots,v_{n-1}\), travel forward on \(c^2v_i+(-1)^i\), cross at the end, return on \(c^2v_i-(-1)^i\) in reverse order, and cross at the start. The two crossings are squares because the endpoints qualify.

Conversely, on an invariant cycle the sheet swap is a fixed-point-free involution of the cycle, reversing the bipartition. The half-turn would move by \(n\) positions and preserve the bipartition since \(n\) is even. The involution must therefore be a reflection through two opposite edges. Those edges are the two crossings at equal base labels. One half of the cycle contains one representative of every vertex orbit; projecting it recovers the unique base path. Exchanging the sheets gives the same unoriented cycle. All other Hamilton cycles pair off under the swap, proving the parity assertion whenever the base path count is odd. On a nonbipartite base, a single sign per vertex at this same scale cannot preserve every square edge. For example, the triangle on 6, 19, 30 already requires both signs at a vertex.

The parallel signed chart gives a second exact representation. Its map on odds is \(\rho(x)=(3x+\chi_4(x))/4\). At scale \(c=10\) it sends \(100v\pm1\) to \(75v\pm1\). Here the types move from the two odd classes modulo four to the two unit classes modulo three. Edges become pairs summing to three times a square, since opposite types cancel their signs. Same-type sums are units modulo three. Thus all 389 symmetric cycles below also have a second realization as integer-label cycles in that chart. [Parallel four-down chart and the 19th-vertex decoder](../chanpak_signed_square36_20261123/README.md)

## The 36-vertex instance

Here \(c=10\). There are 57 base edges and four crossings, at 2, 8, 18 and 32. Enumerating simple paths from the smaller endpoint to the larger endpoint, with the terminal used only after all others, gives the full count in each case.

| Endpoints | Base Hamilton paths |
|---|---:|
| 2, 8 | 20 |
| 2, 18 | 18 |
| 2, 32 | 42 |
| 8, 18 | 81 |
| 8, 32 | 66 |
| 18, 32 | 162 |

One base path, in order, is

```text
2 34 30 19 17 32 4 12 24 25 11 5 20 16 9 27 22 14
35 29 7 18 31 33 3 6 10 26 23 13 36 28 21 15 1 8
```

Its two terminal crossings add to \(20^2\) and \(40^2\). The vertical edges at these midpoints are actual one-arc tournament edits. All 72 labels in this and every other lifted cycle are actual counts. The horizontal square-sum condition compares the numerical counts. A single coherent arc walk through specified tournament representatives requires further data.

## Which square symmetries survive an arc reversal?

The [all-orders arc census](PHYSICAL.md) keeps 101 of the 118 numerical edges. All four same-label crossings survive. Among the 57 base edges, 47 have both lifts realized by arc reversals, seven have neither, and three have exactly one. The three asymmetric pairs of lifts expose the data that the arithmetic sheet swap loses.

| Realized by one reversal | Its sheet image, excluded at all orders |
|---|---|
| 299 and 3301 | 301 and 3299 |
| 399 and 2101 | 401 and 2099 |
| 501 and 3099 | 499 and 3101 |

Despite this loss of symmetry in the whole graph, its **unique** Hamilton cycle is invariant. It is one of the 389 arithmetic cycles, and the only one that survives the one-arc test edge by edge. Its fixed edges are the crossings at 2 and 8. Projecting one half gives this path through the 47 base edges with two realized lifts:

```text
2 7 18 31 33 16 9 27 22 14 35 29 20 5 11 25 24 12
4 32 17 19 30 34 15 21 28 36 13 23 26 10 6 3 1 8
```

In a realization of each edge separately, all 72 edits fit on 13 vertices, and 13 is necessary. At height 199 the two incident edges force different strong-component orders, so one cannot follow through that vertex by keeping a common representative. This obstruction survives relabeling and a sequence of edits that stay at count 199. [Component grammar and finite certificates](PHYSICAL.md)

Each [certificate row](midpoints.rows) gives the midpoint, the tail and head of the flipped arc (zero based), then all outgoing-neighbor bitmasks of the lower tournament. Midpoints 300 through 1500 in steps of 100 come from the [previous target table](../chanpak_ternary_buckets_20261122/target_flips.rows). The two smaller cases came from order-eight-or-smaller tournaments. For the remaining centers we varied every apex mask on sampled nine-vertex deletions. Enumeration completeness of the search is unnecessary for these positive witnesses. The replay checks the tournament axioms and counts full paths directly at both endpoints of all 36 rows. A new source has to be first on every full path, so sources pad every certificate to order ten. Count 3601 exceeds the order-nine maximum 3357, which also bounds all smaller orders by the same padding argument. [THM-329, maximum path-count table](../../../01-canon/theorems/THM-329-a038375-new-terms.md)

There are also unbounded versions. For any tournament of odd count \(h\), put two fixed copies before each selected tournament in an ordinal sum (all edges forward between blocks). Every full path concatenates full paths of its blocks. The 72 counts are all multiplied by \(h^2\), preserving the entire square-sum graph in both directions, its two mod-four types and each marked single-arc edit. The resulting differences within pairs are \(2h^2\). Taking ordinal sums of directed triangles makes this scale unbounded. Unit differences are certified at the original scale.

Run `python3 -B lift.py` from this directory to check all midpoint witnesses, enumerate the six complete endpoint fibers, and verify all 389 lifted cycles edge by edge.
