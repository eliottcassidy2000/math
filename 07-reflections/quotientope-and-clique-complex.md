# The Quotientope Connection and the Clique Complex

*opus-2026-03-23-S240/S241*

## G_n = Q_m / S_n

The arc-reversal graph on labeled tournaments is exactly the m-cube Q_m where m = C(n,2). Each tournament is a vertex of Q_m (binary string of arc orientations), each arc reversal flips one bit. The isomorphism class graph G_n is the orbit graph of Q_m under S_n acting on pair coordinates.

This is NOT a Pilaud-Santos quotientope (which arises from lattice congruences of the weak order on S_n). Our construction quotients the CUBE, not the permutohedron. The analogy:

| Quotientope (Pilaud-Santos) | Our construction |
|---|---|
| Weak order on S_n (n! elements) | m-cube Q_m (2^m elements) |
| Lattice congruence ≡ | S_n × Z_2 group action |
| Edge contraction | Orbit identification |
| Quotient lattice | Orbit graph G_n/Z_2 |
| ALWAYS a polytope | NOT a polytope |

## The Clique Complex: Topology of Tournament Space

The clique complex Δ(G_n/Z_2) has faces = cliques in the merged meta-graph.

| n | f-vector | dim | χ | β_1 | β_2 | genus ≥ |
|---|----------|-----|---|-----|-----|---------|
| 5 | (10, 21, 12, 2) | 3 | -1 | 2 | 0 | 0 |
| 6 | (34, 143, 139, 38, 1) | 4 | -7 | 15 | 7 | 8 |

The Betti numbers explode: β_1 grows from 2 to 15, and β_2 appears for the first time at n=6 with 7 independent 2-cycles (cavities).

## Why It's Not a Polytope

Three independent reasons:
1. **h-vector has negative entries**: (1, 6, -3, -4, 2) at n=5. Polytope boundaries have non-negative h-vectors.
2. **Nontrivial homology**: β_1 > 0 at n=5. Polytope boundaries are spheres (β_1 = 0 for d ≥ 3).
3. **Low connectivity**: vertex connectivity = 2. A d-polytope skeleton has connectivity d (Balinski).

## What the Betti Numbers Mean

β_1 = 2 at n=5 means the clique complex has TWO independent 1-cycles that are not boundaries of any 2-chain. These are "essential loops" in tournament space that cannot be filled. They represent fundamental topological obstructions.

β_1 = 15 at n=6 means FIFTEEN such independent loops. The growth 2 → 15 (factor 7.5) suggests exponential blowup.

β_2 = 7 at n=6 means SEVEN independent 2-cavities — "hollow bubbles" in tournament space that cannot be filled by any 3-chain. These are higher-dimensional topological features that first appear at n=6.

## The Deep Connection

The clique complex Δ(G_n/Z_2) is the **quotient space** Q_m / (S_n × Z_2) in a topological sense. The topology of this quotient space encodes the "shape" of tournament isomorphism — the fundamental question of when two tournaments are structurally equivalent.

The Betti number explosion suggests that the topology of this quotient space becomes increasingly complex, reflecting the increasing difficulty of the tournament isomorphism problem at larger n.

## New Sequences

| Invariant | Values | OEIS |
|-----------|--------|------|
| β_1(Δ(G_n/Z_2)) | 2, 15, ... | NEW |
| β_2(Δ(G_n/Z_2)) | 0, 7, ... | NEW |
| χ(Δ(G_n/Z_2)) | -1, -7, ... | NEW |
| ω(G_n/Z_2) (clique number) | 4, 5, ... | NEW |
| genus(G_n/Z_2) | 0, ≥8, ... | NEW |
