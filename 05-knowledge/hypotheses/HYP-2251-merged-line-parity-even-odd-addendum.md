---
id: HYP-2251
status: PARTIALLY-TRUE line-lift theorem; false for simple merged-edge projection
source: codex-2026-06-05-S675b
related:
  - HYP-2245
  - HYP-2244
  - HYP-2240
  - HYP-2187
  - HYP-1764
---

# HYP-2251: Black-Even And Blue-Odd Live In The Complement-Line Lift

## Claim

This addendum sharpens HYP-2250.  The user's parity intuition is correct after
the missing line/address coordinate is retained:

```text
black complement-line lift = even graph / Eulerian multigraph
blue complement-line lift  = odd-degree graph on its active SC support
```

But it is false if "black" and "blue" mean the older simple merged-metagraph
edge coloring by SC/non-SC endpoint type.

Equivalently, parity is not a property of the naked simple graph
`G_n/Z_2`.  It is a property of the lifted complement-pair line graph over the
fixed-base-path tiling cube.

## Evidence

S675b enumerates the fixed-base-path tiling cube exactly for `n=3..7`,
canonicalizes tournament isomorphism classes, merges tournament complements,
and checks against the known A000568 and merged-edge counts:

```text
n  m  classes  merged  E(G_n/Z_2)  blue lines  black lines
3  1        2       2            1           1            0
4  3        4       3            3           2            2
5  6       12      10           21           8           24
6 10       56      34          143          32          480
7 15      456     272         2123         256        16128
```

For every checked `n`, the complement-line lift has:

```text
active blue degrees all odd = True
black degrees all even      = True
```

The same script audits the simple SC/NS-colored merged edge projection.  That
version fails immediately:

```text
n=4: black odd witnesses [0, 2]
n=5: black odd witnesses [1, 2, 3, 8]
n=6: black odd witnesses [2, 3, 4, 6, 7, 11, 17, 23]
n=7: black odd witnesses [3, 5, 7, 8, 35, 43, 56, 58]
```

So the proof-relevant object is not just a color on metagraph edges; it is the
line-lift with multiplicity and loop half-edges.

## Proof Mechanism

Let `Q_m` be the fixed-base-path tiling cube, where `m=C(n-1,2)`, and let
`C(x)` be all-free-tile complement.  The blue/black line graph has one line for
each unordered pair `{x,C(x)}`.

Blue lines are those whose tiling endpoint is grid-reflection fixed.  Such
tilings live over SC merged nodes.  In an SC class, non-grid-fixed tilings pair
under grid reflection, while the fixed-base tiling count is odd.  Therefore
the grid-fixed endpoint count over each active SC node is odd.  Counting a
self-loop twice, the blue line multigraph has odd degree at every active node.

Black lines are the complement.  In an SC node, black count equals odd total
minus odd blue count, hence even.  In an NS merged pair, black count equals
odd class A plus odd class B, hence even.  Thus the black line lift is
Eulerian.

## Even/Odd Number Interpretation

Even graphs are kernel objects: their boundary is zero over `F_2`.  The black
line lift is literally in that cycle-space kernel.

The blue line lift is an affine odd coset: its boundary is `1` on the SC
support and `0` elsewhere.  Since the number of active SC nodes is even in the
checked rows, that all-ones support is a legal boundary.  Adding any matching
or other odd-boundary repair on the SC support translates the blue lift back
into the even cycle space.

This is the more useful reading of "tournaments are equinumerous with even
graphs": the labelled fixed-base tiling cube and labelled degree-even graphs
share the vector-space size `2^C(n-1,2)`, but quotient-level parity statements
survive only when the line/address coordinate is retained.

## Assumption Challenge

This session explicitly considered alternate vertices: merged tournament
nodes, tiling-complement pairs, single-flip edges, even-graph cycle
coordinates, and proof-obligation addresses.  The predicate preserved by the
successful quotient is degree parity of complement-pair endpoints.  The simple
merged graph preserves adjacency but destroys multiplicity and loop half-edge
data, so it destroys the parity theorem.

## Next Moves

1. Promote the proof sketch into a theorem file, probably next to THM-345/346
   because it uses the same odd-SC / `2 mod 4` NS bucket arithmetic.
2. Test whether the black even graph has stable cycle-space invariants:
   component cycle rank, Euler circuits, and projection to the degree-even
   graph quotient `E_n`.
3. Use the odd blue boundary as a selector for SC support.  If a later quotient
   forgets which SC nodes are odd-boundary active, attach the smallest address
   that restores it.
4. For LRC14, search for owner/carry line lifts whose black/even carrier is an
   Eulerian certificate while blue/odd support marks the exact floor atoms.
