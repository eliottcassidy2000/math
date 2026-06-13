# THM-205: Cone H-Invariance — H(Cone(T)) = H(T)

**Status:** PROVED + VERIFIED (exhaustive n <= 6)
**Source:** opus-2026-03-14-S89

## Statement

For any tournament T on n vertices:
- H(Cone_top(T)) = H(T)
- H(Cone_bot(T)) = H(T)

where Cone_top(T) adds a new vertex v* beating all vertices of T, and Cone_bot(T) adds v* losing to all.

## Proof

**Top cone:** v* beats all vertices => no vertex has an edge to v*. In any Hamiltonian path of Cone_top(T), v* cannot appear at any position except the first (it has no predecessor). Since v* beats all, the edge v*→w exists for every w ∈ V(T). After v*, the remaining path visits all vertices of T — this is exactly a Hamiltonian path of T. So:

H(Cone_top(T)) = |{v* → P : P is a Ham path of T}| = H(T)

**Bottom cone:** Symmetric argument. v* loses to all => v* has no outgoing edges within Cone_bot(T). So v* must be the last vertex in any Hamiltonian path. Every vertex w beats v*, so any Ham path of T can be extended by v*. ∎

## Verification

Exhaustively verified for all tournaments at n = 3, 4, 5, 6 (up to 2^15 = 32768 tournaments at n=6).

## Corollary

H is a **stable invariant** under vertex suspension (coning). Adding a universal source or universal sink to a tournament preserves the Hamiltonian path count.

## Consequence for Fourier Analysis

The Walsh-Fourier structure of H on the (m+n)-edge tournament cube, with n edges fixed (the cone edges), projects down to the same structure on the m-edge cube. The cone edges are "spectators" in the Walsh decomposition.
