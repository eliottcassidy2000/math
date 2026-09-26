# From 233 to 223 by two local flips, and what a ten-vertex lift preserves

**Status: FINITE-EXACT complete distance-two search around one eight-vertex
tournament; PROVED elementary conditional path-contraction identity and
ordered-join transport. No arithmetic or Collatz closure.**

2026-09-26. The difference `233-223=10` prompted the search; it is not an
assumption or a number-theoretic explanation of its outcome.

## 1. Inheritance and the board

* Closest proved mechanism: [THM-082, Hamiltonian path deletion/contraction](../../01-canon/theorems/THM-082-deletion-contraction-ham-paths.md).
  A prescribed contiguous directed block contracts by retaining incoming
  edges from its first vertex and outgoing edges from its last vertex.
* Starting object: [THM-316, staircase antipalindrome](../../01-canon/theorems/THM-316-staircase-antipalindrome.md),
  independently replayed in [the 233 graph lane](crossroads233_20260926_graph.md).
  Its eight-vertex staircase has 233 Hamiltonian paths, not eight or 233
  vertices. The literal ten-vertex staircase has 2489 paths.
* Corrected near miss: MISTAKE-009, `sympy_proof_n8.py Used Simplified n<=7
  Formula`, in `01-canon/MISTAKES.md`. At eight vertices, disjoint 3-cycle/
  5-cycle pairs contribute to the odd-cycle formula. Omitting that term
  gives the wrong explanation of the difference below.
* Canonical hostile: the complete 28-neighbor search shows no one-edge flip
  of this staircase has 223 paths. The closest absolute change is twelve.
* Least-used sidecar: retain both orientations of a **contracted contiguous
  three-vertex block**. Their sum is the mixed response; the two single-flip
  changes by themselves do not determine the simultaneous change.

Board: the staircase; local edge-flip coordinates; contiguous blocks;
disjoint odd cycles; ordered strong components; arithmetic interpretation.
The anchor is the exact 223/233 neighborhood, the niche is conditional
two-edge response, and the wildcard is whether the ten-vertex count can
carry information beyond numerical coincidence.

## 2. Exact search and the two surviving neighbors

On vertices `0,...,7`, pair `(2p,2p+1)` has arc `2p+1 -> 2p`. Define
`rank(2p)=p` and `rank(2p+1)=4+p`; all other arcs point from lower rank to
higher rank. The resulting adjacency bitmasks are

```
T = [252,169,242,164,202,144,42,64],       H(T)=233.
```

Enumerating all 28 one-edge flips finds no 223. Among all 378 unordered
two-edge flips, **exactly two** give 223:

```
flip {0,3} and {2,3}: [244,169,250,161,202,144,42,64];
flip {4,5} and {4,7}: [252,169,242,164,106,128,42,80].
```

The antiautomorphism `i -> 7-i` exchanges the two descriptions. This proves
that the distance from this marked staircase to the level set `H=223` in
the edge-flip cube is exactly two. It does not classify all tournaments
with either count.

For the first pair let `e={0,3}` and `f={2,3}`. The complete response square is

| Reversed pairs | Hamiltonian paths |
|---|---:|
| none | 233 |
| e only | 291 |
| f only | 123 |
| both | 223 |

The separate changes sum to `58-110=-52`, but the simultaneous change is
`-10`. The missing interaction is exactly

```
223-291-123+233 = 42.
```

Every one of these four counts was checked by subset DP and independently
by literal enumeration of all `8!` vertex permutations.

## 3. A conditional two-edge path identity

**Proposition.** In any tournament T containing `u -> v -> w`, let e and f
be the two displayed arcs, and write `H_ij` for the Hamiltonian-path count
after flipping e when `i=1` and f when `j=1`. Then

```
H_11-H_10-H_01+H_00
 = H(T/[u,v,w]) + H(T/[w,v,u]) >= 0.                 (1)
```

Here the notation on the right denotes the **boundary contraction**: remove
the three vertices and insert one block vertex, with incoming edges from
the block's first vertex and outgoing edges from its last. The second
contraction corresponds to the reversed block in the doubly flipped
tournament; outside the block that tournament agrees with T. The contracted
digraph need not be a tournament.

**Proof.** Expand each count over all vertex permutations. Any permutation
not using both underlying adjacencies cancels in the alternating sum. If
both occur in a Hamiltonian path, the shared vertex v must be internal in
the contiguous block, so the only possibilities are `[u,v,w]` and `[w,v,u]`.
The first is legal at corner `00`, the second at `11`, both with positive
sign. Contracting either block is a bijection with Hamiltonian paths of the
corresponding boundary contraction, including when the block occurs first
or last. This proves (1). QED.

This is an elementary conditional form of the inherited block-contraction
mechanism, not a priority claim for an unfamiliar graph invariant. It has
two useful consequences:

* The mixed response is unchanged when any outside edge incident with the
  middle vertex v is reversed: that vertex is internal to the contracted
  block and none of its outside edges is used.
* The directed-path hypothesis is essential for the sign. Taking either
  one-flip corner as the new base gives two arcs both entering or both
  leaving v, and reverses the sign of the same alternating difference.
  Thus there is no universal positive mixed response for arbitrary pairs
  of incident edges.

For the staircase, the original block is `[0,3,2]`. Exactly 31 original
Hamiltonian paths contain it contiguously, and exactly 11 doubly flipped
paths contain `[2,3,0]`. The two six-vertex contractions independently give
the same counts, proving `42=31+11`. All 32 assignments to the five other
edges incident with vertex three preserve this response.

The proposition is also checked exhaustively on every labelled tournament
on three, four and five vertices: 12, 384 and 15,360 directed-wedge tests.
The proof supplies the all-order quantifier; the finite checks do not.

## 4. Why the difference is ten: full odd-cycle information

The inherited [THM-002 odd-cycle collection formula](../../01-canon/theorems/THM-002-ocf.md)
expresses H as the sum of `2^(number of cycles)` over vertex-disjoint
collections of directed odd cycles, including the empty collection. The
primary Irving--Omar [current rendered paper, Corollary 19](https://arxiv.org/html/2412.10572v1#S3.SS4)
states the equivalent permutation formula and attributes it to
Grinberg--Stanley. The current rendering numbers it 19; older repository
citations call it Corollary 20. Reversing Hamiltonian paths identifies the
tournament and its opposite. Cycles here are distinct up to cyclic rotation,
not reversal, and all nontrivial odd lengths are included.

At eight vertices at most two nontrivial odd cycles can be disjoint. Direct
cycle enumeration gives:

| Object | c3 | c5 | c7 | disjoint 3+3 | disjoint 3+5 | H |
|---|---:|---:|---:|---:|---:|---:|
| staircase | 12 | 28 | 28 | 16 | 8 | 233 |
| e reversed | 14 | 36 | 37 | 19 | 10 | 291 |
| f reversed | 9 | 15 | 11 | 9 | 4 | 123 |
| both reversed | 12 | 25 | 24 | 16 | 9 | 223 |

Thus the net change is

```
2*((25-28)+(24-28)) + 4*(9-8) = -14+4 = -10.        (2)
```

Both the triangle count and the number of disjoint triangle pairs stay
unchanged. The difference comes from longer odd cycles together with their
compatibility. This is a precise instance where local triangle statistics
lose the desired global count. Omitting the 3+5 term would incorrectly
predict minus fourteen.

The same two labelled flips applied to staircases of other sizes give:

| vertices | original H | modified H | original minus modified |
|---|---:|---:|---:|
| 4 | 5 | 5 | 0 |
| 6 | 29 | 29 | 0 |
| 8 | 233 | 223 | 10 |
| 10 | 2489 | 2265 | 224 |
| 12 | 33773 | 29367 | 4406 |
| 14 | 562685 | 468673 | 94012 |

This is FINITE-EXACT evidence, not a fitted recurrence. In particular, the
local move does not have a constant numerical cost of ten.

## 5. Literal ten-vertex transport and the strong-tournament boundary

Adjoin a source singleton before the core and a sink singleton after it,
orienting all cross-block edges forward. Every Hamiltonian path must start
at the source, traverse one Hamiltonian path of the core, and end at the
sink. This proves an exact count-preserving bijection for every core.

The ten-vertex adjacency arrays are

```
H=233: [1022,1016,850,996,840,916,800,596,640,0],
H=223: [1022,1000,850,1012,834,916,800,596,640,0].
```

Both have strongly connected component sizes `[1,8,1]`. Padding therefore
preserves H, the core's entire odd-cycle collection polynomial and the
response square, while deliberately adding two trivial components. It does
not supply a strong ten-vertex example. It works at every larger order too,
so the bare condition of having ten vertices gives no special arithmetic
meaning to the difference of ten.

The nearby canon has a different, carefully scoped ten-vertex phenomenon:

* [THM-4137, complete order-ten strong centrality](../../01-canon/theorems/THM-4137-strong-tournament-centrality-complete-order-ten.md)
  checks all 9,355,949 strong isomorphism classes. The central layer uniquely
  optimizes two **certified support floors**, but 3,146,972 classes have only
  noncentral actual response maxima. Floors and maxima cannot be identified.
* [THM-4163, order-eleven homogeneous-pair centrality](../../01-canon/theorems/THM-4163-order-eleven-homogeneous-pair-johnson-centrality.md)
  transports a typed rooted quotient, preserving the full layer lattice.
  Its 93,559,490 rooted presentations are a covering census, not a count of
  distinct eleven-vertex tournaments.
* [THM-4133, cyclic-substitution counterexample](../../01-canon/theorems/THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample.md)
  gives a strong order-twelve failure. Small positive cases do not justify
  all-order centrality; the theorem's exact normalized exposure data, not
  the vertex count, measures the surviving obstruction.

The padded examples lie outside the strongness hypotheses of these results.

## 6. Connection contract, Collatz boundary, and reproduction

| Source and operation | Preserved target | Lost information / needed sidecar |
|---|---|---|
| Eight-vertex staircase, two named flips | Exact response square and block-contraction counts | The scalar H forgets which edges and cycle collections caused it |
| Ordered singleton padding | Hamiltonian paths bijectively and all core odd-cycle collections | Strongness fails; retain component decomposition |
| Odd-cycle collection formula | Exact global H from a finite compatibility object | Triangle totals omit long cycles and disjointness |
| A hypothetical arithmetic interpretation | None established by equality of integers alone | Requires an explicit map preserving a number-theoretic predicate and its exponent/height data |

For Collatz, the earlier graph lane retains valuation words, affine carries,
source height and labelled path cuts. Replacing that full object by a small
tournament's Hamiltonian-path count discards all four. Neither `H=223`,
`H=233` nor the interaction value 42 proves descent, orbit periodicity or
multiplicative independence. The genuine reusable move here is to expose
the **conditional interaction** through its contracted boundary data before
interpreting a scalar numerical match.

Reproduce from the repository root:

```
python -B 04-computation/experiments/crossroads10_20260926_graph.py
python -B -O 04-computation/experiments/crossroads10_20260926_graph.py
```

The [output](crossroads10_20260926_graph.out) retains all 28 one-flip counts,
the complete two-flip hits, all four cycle packets, literal block counts,
32 middle-neighborhood controls, ten-vertex lifts, six family rows, and the
complete small-order mixed-response audit. Integer checks remain active
under optimization; no floating-point calculation or imported counting
engine is used.
