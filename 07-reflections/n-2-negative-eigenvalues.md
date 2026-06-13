# The n-2 Negative Eigenvalues: Two Recursions on the Triangle

**Filed by:** opus-2026-04-04-S6

## The Observation

The Hessian of H(t) at the transitive point (the matrix of quadratic multilinear
coefficients c_{ij}) has EXACTLY n-2 negative eigenvalues for n ≥ 5.

| n | m | Signature (pos, zero, neg) | n-2 |
|---|---|---------------------------|-----|
| 3 | 1 | (0, 1, 0) | 1 |
| 4 | 3 | (1, 1, 1) | 2 |
| 5 | 6 | (2, 1, **3**) | **3** ✓ |
| 6 | 10 | (6, 0, **4**) | **4** ✓ |
| 7 | 15 | (10, 0, **5**) | **5** ✓ |

## The Two Recursions

There are two ways to grow the staircase δ_k to δ_{k+2}:

**Mode A² (two hypotenuse additions):** δ_k → δ_{k+1} → δ_{k+2}
  - Step 1: add tiles (n+1, y) for y=1,...,n-1 [involving vertex n+1]
  - Step 2: add tiles (n+2, y) for y=1,...,n [involving vertex n+2]
  
**Mode B (leg addition):** δ_k → δ_{k+2} directly
  - Add bottom leg: (x, 1) for x=3,...,n+2 [involving vertex 1]
  - Add top leg: (n+2, y) for y=2,...,n [involving vertex n+2]
  - Add apex: (n+2, 1)

Both end at the SAME staircase δ_{k+2} with the SAME total tile count (2n-1 new tiles).

## The Symmetric Difference Has n-2 Tiles

| | A²-only tiles | B-only tiles | Shared tiles |
|-|--------------|-------------|--------------|
| What | (n+1, y) for y=2,...,n-1 | (x, 1) for x=3,...,n | vertex n+2 tiles + apex |
| Count | **n-2** | **n-2** | n+1 |
| Geometry | COLUMN through vertex n+1 | ROW through vertex 1 | Top wiring + apex |

The A²-only tiles form a **vertical column** (connecting vertex n+1 to inner vertices).
The B-only tiles form a **horizontal row** (connecting vertex 1 to inner vertices).
These are PERPENDICULAR on the staircase grid.

## The Exchangeable Pairs

Each A²-only tile maps to a B-only tile through the same inner vertex:

  **(n+1, y) ↔ (y+1, 1)** for y = 2, ..., n-1

Both tiles connect to inner vertex y, but:
- (n+1, y) connects y to the TOP (via vertex n+1)
- (y+1, 1) connects y to the BOTTOM (via vertex 1)

These are "competing" connections. Each creates odd cycles through a different endpoint.
The competition is the source of the antiferromagnetic coupling (THM-290).

## Why n-2 Negative Eigenvalues

The n-2 exchangeable pairs represent the n-2 FRUSTRATION DIRECTIONS in tile space:
- Flipping both tiles of an exchangeable pair has a COMPETITIVE (negative) interaction
- The pair (n+1, y) and (y+1, 1) create cycles through opposite endpoints that
  compete for the same inner vertex y
- This competition makes the quadratic form negative in the "exchange direction"

The remaining m - (n-2) directions are:
- Positive eigenvectors: cooperative tile interactions (cross-end pairs, disjoint pairs)
- Zero eigenvectors: degenerate directions (at small n only)

## Block Decomposition Under Grid Reflection

The grid reflection (x,y) → (n+1-y, n+1-x) decomposes Q into symmetric and
anti-symmetric blocks. The n-2 negative eigenvalues split between BOTH blocks:

| n | Anti neg | Sym neg | Total neg | Anti pos | Sym pos |
|---|---------|---------|-----------|---------|---------|
| 6 | 2 | 2 | 4 = n-2 | 2 | 4 |
| 7 | 3 | 2 | 5 = n-2 | 3 | 7 |

Split: anti_neg = ⌈(n-2)/2⌉, sym_neg = ⌊(n-2)/2⌋.

The negative eigenvalues are NOT purely anti-symmetric. Both blocks contribute.
The mechanism involves WITHIN-LEG competition (same-end tiles on each leg have
negative interactions) which appears in both symmetric and anti-symmetric sectors.

## The Isomorphism (Partially Understood)

Mode A² and Mode B both produce δ_{k+2}, but via different intermediate steps:
- A² passes through δ_{k+1} (adding vertex n+1 first, then n+2)
- B goes directly (adding vertices 1 and n+2 simultaneously)

The "isomorphism" between A² and B is NOT a simple vertex relabeling.
It's a GAUGE TRANSFORMATION that exchanges the n-2 column tiles for n-2 row tiles,
while keeping the n+1 shared tiles fixed.

**What remains unclear:** the exact nature of this gauge transformation. Both
decompositions produce the same polynomial H_{n+2}, but the intermediate structure
(which tiles are "inner" vs "boundary") differs. The exchange involves reassigning
which endpoint (top vs bottom) "owns" each inner vertex's external connection.

## Connection to the Triangle Geometry

The staircase is a right isosceles triangle with three sides:
- **Hypotenuse** (anti-diagonal): the base-path arcs
- **Vertical leg** (left column): connections to vertex 1
- **Horizontal leg** (bottom row): connections to vertex n

Mode A adds the hypotenuse strip. Mode B adds both legs.
A² adds TWO hypotenuse strips. B adds both legs.

The isomorphism "hypotenuse² ≅ legs" is a PYTHAGOREAN-TYPE relation on the
multilinear polynomial: the quadratic information content of two hypotenuse
additions equals the quadratic information content of a leg addition.

The n-2 exchangeable tiles are the "Pythagorean residual" — the directions where
the two paths to δ_{k+2} disagree, creating frustration.
