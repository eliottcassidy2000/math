# Dark Tournaments and the Sign Representation

**Session:** kind-pasteur-2026-03-25-S20gb
**Date:** 2026-03-25

## The Royle-Praeger Definition (2022)

An **even graph** (Royle-Praeger sense) is one where no automorphism reverses an odd number of oriented edges. An **odd graph** has at least one such sign-reversing automorphism.

This is NOT about vertex degree parity. It's about the **sign representation** of the automorphism group acting on the edge space.

The sign: sgn_X(g) = product over {u,v} in E(X) of (u^g - v^g)/(u - v) = ±1.
Even: sgn = +1 for ALL automorphisms. Odd: some automorphism has sgn = -1.

## The Equinumerosity (PROVED, 2022-2023)

|even graphs on n vertices| = |tournaments on n vertices| = A000568(n)

Proved by Devillers, Freedman, Glasby, Praeger, Royle via Burnside counting.
**The bijection is an OPEN PROBLEM.** No natural bijection is known.

## The Dark Tournament Sequence (NEW, not in OEIS)

|odd graphs on n vertices| = A000088(n) - A000568(n) = 0, 1, 2, 7, 22, 100, 588, 5466

This sequence is **absent from OEIS** and appears unstudied.

## The Four-Object System

| Object | Count | Burnside sum over | Even/Odd |
|--------|-------|-------------------|----------|
| Tournaments | A000568 | All-odd-cycle σ only | Even (light) |
| Dark tournaments | A000088 - A000568 | Even-cycle σ only | Odd (dark) |
| Even graphs | A000568 | Same Burnside count! | Even (light) |
| Odd graphs | A000088 - A000568 | Complementary count | Odd (dark) |

The key: tournaments are counted by the all-odd-cycle part of the Burnside sum. The remaining terms (even-cycle σ) count the "dark" objects. The graph Burnside sum uses ALL σ, while the tournament sum excludes even-cycle σ (which fix 0 tournaments).

## The Bipartite Metagraph Structure

In the graph edge-flip metagraph:
- **Even classes have ZERO internal edges** (edge flip always changes degree parity)
- **Odd classes have rich internal structure** (the dark metagraph)
- **Cross edges** connect even to odd only

The tournament world (wiggly metagraph on arc-flips) and the graph world (edge-flip metagraph) are ORTHOGONAL structures. Tournaments communicate through arc-direction flips; graphs communicate through edge-presence flips. These are different operations on different spaces.

## Open Questions

1. **What is the natural bijection** between even graphs and tournaments? (Royle-Praeger open problem)
2. **Does the dark metagraph have a "Hamiltonian path count"** analog? (Some function D(G) playing the role of H(T))
3. **Are there "dark tilings"** — a staircase model for odd graphs?
4. **Does the fractal codec structure** extend to dark tournaments?
5. **Is the dark sequence 0,1,2,7,22,100,588,5466 Burnside-computable?** (Yes, via the even-cycle terms)

## The Deep Symmetry

Tournaments and dark tournaments are the YIN and YANG of graph theory:
- Tournaments = even graphs = all-odd-cycle Burnside = the "light" side
- Dark tournaments = odd graphs = even-cycle Burnside = the "dark" side
- Together they make all graphs = the complete Burnside sum
- They don't interact directly (different metagraph structures) but are bound by the constraint that their counts sum to A000088.

The sign representation of the automorphism group is the BRIDGE between them: even graphs have trivial sign representation, odd graphs have non-trivial sign. This is the simplest possible representation-theoretic distinction — the binary split of graphs into two classes based on the most basic invariant of the automorphism action.
