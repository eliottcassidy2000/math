# THM-303: Grid Reflection Symmetry of the Multilinear Polynomial

**Status:** PROVED (opus-2026-04-04-S3)
**Verified computationally:** n=5,6,7 (100% coefficient matching)

## Statement

The multilinear polynomial H(t₁,...,tₘ) is invariant under the grid reflection (x,y) → (n+1-y, n+1-x), which is the anti-diagonal reflection of the staircase Young diagram δ_{n-2}.

For every tile subset S, letting σ(S) denote the reflected tile subset:
$$c_S = c_{\sigma(S)}$$

Equivalently: H(t) = H(σ(t)) for all tilings t.

## Proof

The grid reflection (x,y) → (n+1-y, n+1-x) corresponds to the **path reversal** involution on the base Hamiltonian path. The base path n→n-1→...→1 becomes 1→2→...→n under reversal, which is the base path for the relabeled tournament π(T) where vertex v maps to n+1-v.

Under this relabeling:
1. The base path is preserved (same path, reversed direction = same path with relabeled vertices)
2. Tile (x,y) maps to tile (n+1-y, n+1-x) — the grid reflection
3. The forward/backward direction of each tile is preserved (both endpoints swap high↔low)
4. H(T) = H(π(T)) because the vertex relabeling v→n+1-v is a graph isomorphism that maps Hamiltonian paths to Hamiltonian paths bijectively

Therefore H is invariant under grid reflection, and the multilinear coefficients inherit this symmetry.

## Verification

| n | Total nonzero coefficients | Reflection-invariant |
|---|---------------------------|---------------------|
| 5 | 34 | 34/34 |
| 6 | 199 | 199/199 |
| 7 | 1781 | 1781/1781 |

All non-constant coefficients are perfectly matched by their reflections.
