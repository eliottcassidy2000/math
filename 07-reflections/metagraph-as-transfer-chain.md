# The Metagraph as a Transfer Chain

**Filed by:** opus-2026-04-04-S7

## The Vertex-Addition Construction

The metagraph G_n is built by a CHAIN of vertex additions:

```
1 class → 2 classes → 4 classes → 12 classes → 56 classes → ...
 (k=2)     (k=3)       (k=4)        (k=5)         (k=6)
```

At each step k → k+1, we add vertex k+1 to the tournament on {1,...,k}.
The new vertex connects to all existing vertices via:
- One base-path arc (k+1 → k, always present)
- k-1 tile arcs (to vertices 1,...,k-1, each independently forward/backward)

Each k-class extends to multiple (k+1)-classes depending on the k-1 binary
choices for the new vertex's connections.

## The Transfer Matrix

T_k: an |A000568(k+1)| × |A000568(k)| matrix where
  T_k[C', C] = number of (k-1)-bit extensions mapping class C to class C'

The full metagraph G_n is the PRODUCT of these transfer matrices:
  G_n = T_{n-1} · T_{n-2} · ... · T_3

This product builds the metagraph INCREMENTALLY, adding one vertex at a time.

## Key Properties

1. **Fan-out grows**: each class at step k maps to O(2^{k-1}/symmetry) classes at k+1.
   At k=5→6: 6-16 successor classes per predecessor.

2. **H-range expands monotonically**: the range [H_min, H_max] grows at each step.
   The transitive class (H=1) always generates all H values up to some threshold.

3. **The H-maximizer propagates**: the class achieving max H at step k always
   generates the max-H class at step k+1. This is the HEREDITARY MAXIMIZER property.

4. **The strip tiles form a GROWING COLUMN**: strip k contains k-1 tiles connecting
   the new vertex to ALL existing internal vertices.

## Explicit Representation at n=5

The n=5 metagraph has 12 vertices, 30 edges, and the following weighted adjacency:

```
       H: 1  3  3  5  3  5  9  9 11 13 15 15
H= 1:    0  1  1  1  1  1  1  0  0  0  0  0  (transitive → 6 neighbors)
H=15:    0  0  0  0  0  0  0  0  1  1  0  0  (regular → 2 neighbors)
```

The transitive is the HUB (max degree = m = 6, connected to all reachable classes).
The regular is the ENDPOINT (min degree = 2, connected only to near-maximizers).

The metagraph DIAMETER is 3: transitive → H≤9 (d=1) → H≤15 (d=2) → regular (d=3).

## Tiling-Fibration Identity

For every isomorphism class C:
  **size(C) × |Aut(T_C)| = H(T_C)**

This connects the three fundamental quantities: class size in the tiling model,
automorphism group order, and Hamiltonian path count. The product is always H.

## Connection to the Multilinear Polynomial

The multilinear polynomial H(t₁,...,t_m) can be computed via the transfer chain:
  H(t) = ⟨1| T_{n-1}(t_{strip_{n-1}}) · ... · T_3(t_{strip_3}) |1⟩

where each transfer matrix depends on the tile variables of its strip.
This gives a FACTORED representation of H(t) as a chain of small linear maps.

The DEGREE of the multilinear polynomial (2⌊(n-1)/2⌋) corresponds to the
maximum degree contributed by the transfer chain: each transfer T_k contributes
at most ⌊(k-1)/2⌋ tile arcs that participate in the max-degree cycle
(the reversal cancellation limits the effective degree).
