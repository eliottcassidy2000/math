# Endpoint Collisions As Incidence Not Adjacency S95

The support-3 SC collision looked at first like another triangle.

That was the natural guess.  The repo already has many honest triangles:
directed 3-cycles, Petersen/Kneser root orthogonality, beta-3 generators,
Krawtchouk `B_3`, disjoint 3-cycle packings.  A non-private endpoint column
with three odd parent owners seemed likely to be a literal triangle somewhere
in the parent merged metagraph.

The geometry probe says no.

At `5->6`, the single collision owner triple is a clique.  But at `6->7`, the
14 owner triples have induced parent-metagraph edge counts:

```text
{0: 1, 1: 6, 2: 5, 3: 2}
```

One support-3 collision has no parent-metagraph adjacency among its owners at
all.  So the ternary object is not the clique complex speaking through endpoint
transfer.  It is endpoint-transfer incidence itself.

## What This Refines

The cubic-obstruction slogan survives, but in a stricter form:

```text
three does not always mean triangle;
sometimes it means three preimages of a boundary state.
```

That distinction matters.  Petersen/root orthogonality is adjacency geometry:
which root pairs are disjoint.  Directed 3-cycles are orientation geometry:
which three arcs fail to rank.  Endpoint SC collisions are quotient-incidence
geometry: which three parent nodes map oddly to the same SC child boundary.

They rhyme because all three are first failures of a binary simplification, but
they are not the same object.

## Feedback Into The Proof Search

The merged endpoint-rank conjecture should not be reduced too early to the
ordinary merged metagraph `G_n/Z_2`.  The parent metagraph is a shadow of the
same tiling cube, but it forgets endpoint-transfer multiplicities.  The
collision block lives in the transfer matrix.

The second pass gives a positive replacement for the failed triangle picture.
The support-3 collision hypergraph leaf-peels completely through `6->7`:

```text
peel_removed = [0, 0, 0, 1, 14]
peel_core    = [[], [], [], [], []]
```

So the collision block may still be triangular, but only after moving into the
right category.  A parent owner that is not private globally can become private
after earlier collision columns are removed.

So HYP-1791 is better stated as a 3-uniform hypergraph incidence problem:

```text
rows    = parent merged nodes
columns = non-private SC child boundary nodes
entry   = odd endpoint-transfer incidence
```

The metagraph can still provide features: edge count among owners, SC/NSC
types, good-cut/SCC-defect height, H-values.  But it is not the carrier of the
relation.

## Connection Back To Older Threads

This is another example of the boundary-state principle:

- scalar `H` failed to recurse until bridge exposure `Q_T` was added;
- fixed `t3` failed to determine `H` until the lifted object `Omega(T)` was
  inspected;
- projection defects became meaningful only after comparing two quotient maps
  along structured moves;
- endpoint collisions are invisible if we only ask for adjacency in the parent
  metagraph.

The useful object is the smallest interface that remembers how the operation
acts.  For endpoint insertion, that interface is the transfer-incidence matrix,
not the first graph shadow cast by its rows.
