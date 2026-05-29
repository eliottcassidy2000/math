# Recursive Boundary S95

The fixed-path explorer has a literal product recursion:

```text
Q_{m(n+1)} = Q_{m(n)} x Q_{n-1}
m(n) = C(n-1,2).
```

Add the new endpoint to the Hamiltonian path. The old tiling stays fixed; the
new free choices are the `n-1` arcs from the new endpoint to all old non-last
vertices.

This turns the vague "recursive metagraph" idea into an exact transfer matrix.

## Class Transfer

Let `A_n(C,D)` count endpoint children from parent class `C` at `n` to child
class `D` at `n+1`.

```text
row sum:     Σ_D A_n(C,D) = 2^(n-1) F_n(C)
column sum:  Σ_C A_n(C,D) = F_{n+1}(D)
```

The first line says every parent tiling gets all endpoint signatures. The
second line says every child tiling has a unique parent by deleting the final
endpoint.

For merged nodes the same statement holds with `M` replacing `F`.

Measured support:

```text
class support:  [2, 5, 22, 181, 2338]
merged support: [2, 5, 19, 112, 1312]
```

The transfer is sparse, but not thin: by `6->7`, each child class has on
average about five parent classes, while each parent class reaches about 42
child classes.

## Parity Boundary

The row sums are even; the unmerged child column sums are odd. Therefore:

```text
XOR(parent parity rows) = all child classes.
```

In the merged quotient, child columns are odd exactly for self-complementary
nodes:

```text
XOR(parent parity rows) = SC child indicator.
```

This is the same boundary phenomenon as the edge law

```text
2 lambda_u + Σ tau_uv = m M_u,
```

but one categorical level higher. Tile-flip edges give a boundary inside one
level; endpoint transfer gives a boundary between adjacent levels.

## Full Row Rank

The surprise is rank:

```text
unmerged GF2 rank: [1, 2, 4, 12, 56]
parent classes:    [1, 2, 4, 12, 56]

merged GF2 rank:   [1, 2, 3, 10, 34]
parent merged:     [1, 2, 3, 10, 34]
```

Through `n=7`, endpoint insertion is parity-injective on the quotient tower.
No parent-level parity signal vanishes in the children.

This suggests a new recursive invariant: the quotient tower may be a chain of
full-row-rank `GF(2)` transfer maps.

## Cross-Thread Interpretation

1. **Bucket parity.** The previous bucket theorem says odd merged masses are
   exactly SC nodes. Transfer parity recovers the same SC set at the next level.

2. **Two reductions.** Mode A is now a transfer operator. Mode B is a twisted
   overlap bundle. The difference is clear: Mode A has a unique parent per
   child tiling; Mode B forgets a labeled overlap and creates twist.

3. **Insertion exposure.** Arbitrary vertex insertion needs start counts, end
   counts, and bridge exposure `Q_T`. Endpoint transfer is the fixed-path
   aggregate where those local cut sums become a class-level matrix.

4. **Even graphs.** The same transfer was built for the even-graph quotient,
   and it fails full row rank:

   ```text
   even parent counts: [1, 2, 3, 7, 16]
   even GF2 ranks:     [1, 1, 2, 6, 8]
   ```

   The cycle-space lens keeps the row/column transfer theorem but loses
   tournament parity-injectivity.

5. **Betti/deletion threads.** Vertex deletion has repeatedly been the tool for
   proving path-homology constraints. This transfer matrix is the tiling-side
   adjoint of deletion: every child has a unique endpoint-deletion parent, but
   classes have many parents and children.

## New Questions

1. Prove or refute full row rank of endpoint transfer over `GF(2)`.

2. Explain the even-graph rank defect. The child boundary is the set of
   classes with odd orbit size `n!/|Aut(G)|`, so the answer should involve
   Sylow-2 structure of graph automorphism groups.

3. Look for a Smith normal form of the integer transfer matrix. The mod-2 rank
   may be the first shadow of a stronger integral divisibility pattern.

4. Relate transfer support to automorphism tax. Rigid parent classes should
   spread more freely; symmetric parents may produce structured collisions.

5. Compose transfers:

   ```text
   A_{n+1} A_n : level n -> level n+2.
   ```

   Compare this fast-time composition with the Mode B source-sink descent.

The pattern is recursive but not scalar. The right object is not a sequence; it
is a tower of sparse integer transfer matrices whose mod-2 boundary remembers
self-complementarity.
