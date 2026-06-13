# The Deep Symmetry: G_n/Z_2 as a Tournament on (n-2) Effective Vertices

**Session**: kind-pasteur-2026-03-22-S20cc

## The Observation

G_4/Z_2 = K_3 = the complete graph on 3 vertices = the simplest tournament.

The merged meta-graph at n=4 IS a tournament at n=3. The meta-theory reduces by 2 dimensions.

## The Pattern

| n | G_n/Z_2 vertices | n-2 | Match? |
|---|-----------------|-----|--------|
| 3 | 2 | 1 | K_2 = edge = "tournament on 2" |
| 4 | 3 | 2 | K_3 = "tournament on 3" |
| 5 | 10 | 3 | 10 >> 3, but C(5,2) = 10... |
| 6 | 34 | 4 | 34 >> 4, but close to C(6,2)+... |

At n=3: 2 merged vertices = a line segment = the "tournament" on 2 items.
At n=4: 3 merged vertices = K_3 = the tournament on 3 items.

But at n=5 it breaks: 10 merged vertices, not 3. However, 10 = C(5,2) = the number of ARCS in a 5-tournament. And C(n,2) is the dimension of the tournament space.

So G_n/Z_2 has vertices equal to:
- n=3: C(2,1) = 2. Or just 2.
- n=4: C(3,2) = 3. The number of arcs in K_3.
- n=5: C(4,2) + C(4,1) = 6 + 4 = 10? Or (A000568(5) + 8) / 2 = 10.
- n=6: (A000568(6) + 12) / 2 = 34.

The merged vertex count is (A000568(n) + SC_count(n)) / 2. The SC count sequence: 2, 2, 8, 12... This doesn't give a clean formula for the merged count.

## The Deeper Meaning

What IS G_n/Z_2 structurally? It's the space of **unoriented tournament types** -- where we don't distinguish between a tournament and its complement.

In physics, this is the **projective** version: not distinguishing between a state and its "opposite." Just as projective space RP^n is the quotient of the sphere S^n by the antipodal map, G_n/Z_2 is the quotient of G_n by the complement map.

The complement map sends T to T^op (reverse all arcs). On the iso class graph, it maps class C to its complement class. SC classes are fixed points. NS pairs are free orbits.

## Why K_3 at n=4

At n=4 there are 4 iso classes: transitive (H=1), two H=3 classes, and H=5. The complement map:
- Transitive (H=1, SC) -> itself
- H=3 class A (NS) -> H=3 class B (NS)
- H=5 (SC) -> itself

After merging: 3 nodes. Every pair is connected because:
- Transitive connects to both H=3 classes (black edges) -> in merged graph, connects to the merged H=3 pair
- Transitive connects to H=5 (blue edge)
- H=3 merged connects to H=5 (black edges)

So: K_3 = complete graph = every merged node connects to every other.

K_3 with 3 nodes is also the Cayley graph of Z/3Z with generator 1. And the cyclic group Z/3Z is the automorphism group of the 3-cycle tournament at n=3. So G_4/Z_2 has the structure of the automorphism group of G_3.

## The Recursive Principle

G_4/Z_2 = K_3 = the "tournament graph" at the n=3 level.
G_3/Z_2 = K_2 = the "tournament graph" at the n=2 level.

If this continues: G_n/Z_2 should have a structure related to the (n-2)-level theory.

At n=5: G_5/Z_2 has 10 vertices and 21 edges. Compare with the tournament space at n=3: 8 labeled tournaments, 2 iso classes. The connection isn't a simple "G_n/Z_2 = G_{n-2}" but something more nuanced.

10 = C(5,2) = number of arcs in a 5-tournament. 21 = C(7,2) = edges in K_7. Or 21 = C(7,3) = blocks in STS(7) = lines of the Fano plane. These numerological coincidences may or may not be meaningful.

## The Real Structure

What the complement-merge reveals is the **underlying symmetry class of tournaments**. By factoring out both vertex permutation (S_n) AND complement (Z/2Z), we get the irreducible building blocks.

Each merged node represents a pair {T, T^op} of complementary tournament types. These are the objects that are genuinely structurally different -- not just "the same tournament viewed from the other side."

The blue edges in the merged graph connect types that are related by SC-preserving flips. The black edges connect types related by SC-breaking flips. The structure is:
- Blue subgraph: the "orientation-preserving" connections
- Black subgraph: the "orientation-reversing" connections

This is exactly the structure of a **Klein geometry**: a space with a group of motions (flips) that splits into orientation-preserving (blue) and orientation-reversing (black).

## What This Means for the Edge Formula

The user's insight: count blue and black edges separately, then add.

In the merged graph:
- Blue edges connect "same-type" nodes (both SC, or both merged-NS)
- Black edges connect "cross-type" nodes

The blue edge count may be expressible in terms of the SC class structure alone (an internal property of the blue skeleton).
The black edge count may be expressible as a bipartite matching between SC classes and merged-NS pairs.

If both have clean formulas, the total edge count = blue + black + collapsed + twin, giving the exact E(G_n).

## The Grand Vision

Tournament theory at level n has a meta-theory at level n-2 (via the complement-merged quotient). This meta-theory itself has a meta-meta-theory at level n-4. And so on, descending by 2 each time, until hitting the trivial level.

This DESCENT is the Cayley-Dickson tower in reverse: from S (n=17) to O (n=9) to H (n=5) to C (n=3) to R (n=2). Each step removes one layer of complexity. The complement merge is the mechanism of descent.

The full analytical understanding of G_n requires understanding this descent -- how the structure at level n is built from the structure at level n-2, with the complement merge as the connecting map.
