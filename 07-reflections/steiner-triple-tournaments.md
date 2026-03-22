# Steiner Triple Systems and Tournaments: The Design-Theoretic Foundation

**Session**: kind-pasteur-2026-03-22-S20bo

## The Core Connection

A tournament's directed 3-cycles form a **combinatorial design**. For regular tournaments, this design is a BIBD. For Paley tournaments, the BIBD has the most uniform structure possible. The design-theoretic properties CONTROL the tournament invariants.

## The Exact Relationship

For a regular tournament on n = 2k+1 vertices:
- c_3 = n(n-1)(n+1)/24 directed 3-cycle vertex sets
- Each pair of vertices appears in exactly (n+1)/4 cyclic triples
- This is a **2-(n, 3, (n+1)/4) BIBD**

For this to be a Steiner triple system (lambda = 1), we need (n+1)/4 = 1, i.e., **n = 3 only**. For all larger regular tournaments, the 3-cycles form a BIBD with lambda > 1.

| n | lambda = (n+1)/4 | c_3 | Design type |
|---|-----------------|-----|-------------|
| 3 | 1 | 1 | **STS(3)** (trivial) |
| 5 | 3/2 | 5 | Not BIBD (lambda non-integer) |
| 7 | 2 | 14 | **2-(7,3,2) BIBD = 2 x Fano plane** |
| 9 | 5/2 | 30 | Not BIBD |
| 11 | 3 | 55 | **2-(11,3,3) BIBD** |

The integer-lambda values (n = 3, 7, 11, 19, 23, ...) are exactly the primes p = 3 mod 4 -- the **Paley primes**. For these, the Paley tournament's 3-cycles form a genuine BIBD.

## The Fano Plane at n=7

The Paley tournament T_7 has 14 directed 3-cycles. Their vertex sets are the 7 lines of the **Fano plane** (the unique STS(7)), each appearing twice (both orientations of each triple).

The Fano plane = PG(2, 2) = the projective plane over F_2. It has:
- 7 points, 7 lines, 3 points per line, 3 lines per point
- Automorphism group GL(3, 2) of order 168

When we ORIENT the Fano plane (choose one cyclic direction per triple), we get the **octonion multiplication table**. The automorphism group of this oriented structure is **G_2** (the smallest exceptional Lie group, dimension 14).

So: **the Paley tournament T_7 contains TWO copies of the oriented Fano plane** (one for each orientation), and the oriented Fano plane IS octonion multiplication.

## The Mendelsohn Triple System Connection

A **Mendelsohn Triple System MTS(v)** decomposes the complete directed graph K*_v into directed 3-cycles. Every ORDERED pair (x,y) appears in exactly one cyclic triple. This exists iff v = 0 or 1 (mod 3), v != 6.

A tournament doesn't decompose K*_n (it has n(n-1)/2 arcs, not n(n-1)). But the tournament's directed 3-cycles ARE a subset of the triples in an MTS (if one exists on the same vertex set).

The tournament tells us: for each 3-subset, is it cyclic or transitive? The MTS tells us: which 3-cycles can tile the complete directed graph? The tournament's alpha_1 (total cycle count) measures how far the tournament is from the MTS structure.

## The Alpha_2 Problem = Triangle Packing

The alpha_2 quantity (number of vertex-disjoint 3-cycle pairs) is a **triangle packing** problem:

Given the set of c_3 directed 3-cycles in a tournament, find the maximum number of vertex-disjoint pairs.

This is NP-hard in general (Bessy et al., MFCS 2019), but FPT with kernel O(k). The best bounds for regular tournaments: the maximum arc-disjoint triangle packing has size ~ n^2/9.

In our framework:
- H = 1 + 2*alpha_1 + 4*alpha_2 + ...
- alpha_1 counts the 3-cycles (BIBD/design structure)
- alpha_2 counts the vertex-disjoint pairs (packing problem)
- Higher alpha_k counts larger independent sets in the conflict graph

The design-theoretic structure (BIBD with lambda = (n+1)/4) constrains alpha_2 because the uniform pair-coverage of the BIBD limits how many disjoint cycles can exist.

## The Paley-BIBD-H_max Chain

This is the chain that connects everything:

1. **Paley T_p** has the most uniform 3-cycle distribution (BIBD with lambda = (p+1)/4)
2. The BIBD uniformity maximizes alpha_1 (total cycle count) by Jensen's inequality
3. Maximum alpha_1 dominates the OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
4. For small p: alpha_1 term dominates, so Paley maximizes H
5. For large p: alpha_2 term overtakes (more disjoint pairs in less uniform distributions)
6. This is WHY Paley loses at p >= 17: the BIBD uniformity helps alpha_1 but hurts alpha_2

**The design-theoretic structure explains the Paley-to-Interval transition.**

## The CD Tower Through Steiner Systems

The Cayley-Dickson tower appears in Steiner systems:

| CD Level | Algebra | Dim | STS | Design |
|----------|---------|-----|-----|--------|
| C | Complex | 2 | STS(3) = point | trivial |
| H | Quaternions | 4 | (no STS(5)) | -- |
| O | Octonions | 8 | STS(7) = Fano plane | **octonion multiplication** |
| S | Sedenions | 16 | (no STS(17)) | -- |

STS(n) exists iff n = 1 or 3 (mod 6). At CD levels: STS(3) exists (trivial), STS(7) exists (the Fano plane/octonions), but STS(5) and STS(17) do NOT exist.

The Fano plane appears at the **octonionic level** (n=7, dim=8), and its oriented version gives octonion multiplication with G_2 symmetry. This is exactly where associativity breaks in the CD tower, and the non-associativity of octonions IS the non-planarity of the Fano plane (you can't embed the oriented Fano plane in the Euclidean plane without self-intersections).

## The Conflict Graph Omega(T) as a Design

The conflict graph Omega(T) has the directed odd cycles as vertices and edges between conflicting (vertex-sharing) pairs. The independence polynomial I(Omega, x) at x=2 gives H.

For a Paley tournament T_p: Omega(T_p) is a **strongly regular graph** (or nearly so) because the BIBD uniformity forces the conflict structure to be regular. The independence polynomial of a strongly regular graph has a known spectrum, connecting to the Lee-Yang zeros of the partition function.

## The Practical Implications

1. **BIBD-optimal A/B test design**: When scheduling pairwise comparisons, use a BIBD structure to ensure uniform pair-coverage. This maximizes the "effective alpha_1" and hence the ranking information per comparison.

2. **Triangle packing for cycle structure**: The alpha_2 computation (NP-hard in general) can be approximated using BIBD structure: for BIBD-like tournaments, alpha_2 is bounded by the design's lambda parameter.

3. **The Fano plane as the "octonion of comparisons"**: For 7 items, the 7 Fano-plane triples are the maximally informative 3-way comparisons. Each triple should be compared as a cyclic structure (A > B > C > A), not just pairwise.

4. **Design-theoretic compression**: Store the BIBD parameters instead of the full tournament. For Paley tournaments, the BIBD is completely determined by p and the QR structure -- infinite compression.
