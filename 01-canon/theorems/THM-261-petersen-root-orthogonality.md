# THM-261: Petersen Graph = A_4 Root Orthogonality Graph

**Status:** PROVED (verified computationally, follows from standard definitions)
**Filed by:** kind-pasteur-2026-03-21-S18

## Statement

The Petersen graph K(5,2) is isomorphic to the orthogonality graph of the A_4 positive root system. More generally, the Kneser graph K(n,2) is isomorphic to the orthogonality graph of the A_{n-1} positive root system for all n >= 3.

Specifically: the positive roots of A_{n-1} are {e_i - e_j : 0 <= i < j <= n-1}, indexed by 2-subsets of [n]. Two roots alpha = e_i - e_j and beta = e_k - e_l are orthogonal (inner product 0) if and only if {i,j} and {k,l} are disjoint, which is exactly the adjacency condition in K(n,2).

## Proof

Direct computation. The inner product (e_i - e_j, e_k - e_l) = delta_{ik} - delta_{il} - delta_{jk} + delta_{jl}. This equals 0 iff {i,j} ∩ {k,l} = empty (all four Kronecker deltas vanish). The Kneser graph K(n,2) has edges between disjoint 2-subsets. QED.

## Significance

This places the Petersen graph (and the entire K(n,2) family) firmly within Lie theory:
- K(5,2) = Petersen: orthogonality graph of 10 positive roots of A_4 = sl(5)
- K(n,2): orthogonality graph of C(n,2) positive roots of A_{n-1} = sl(n)

Since dim so(n) = C(n,2) = # positive roots of A_{n-1}, and the standard basis of so(n) IS the set of root vectors {E_{ij} - E_{ji}}, tournaments live at the intersection of so(n) and A_{n-1}.

## Computational verification

- Adjacency matrices of Petersen and A_4 orthogonality graph are IDENTICAL (10x10)
- Both are 3-regular with 15 edges
- Petersen spectrum: 3^1, 1^5, (-2)^4 matches standard result
- Verified that K(n,2) matches A_{n-1} orthogonality for n=3,4,5,6,7,8

## Connection to tournament theory

A tournament T on [n] assigns eps_{ij} in {+1,-1} to each positive root alpha_{ij} of A_{n-1}, defining a skew-adjacency matrix B = sum eps_{ij} (E_{ij} - E_{ji}) in so(n). This is simultaneously:
- A vertex of the hypercube {+/-1}^{C(n,2)} in so(n)
- A signed function on the positive roots of A_{n-1}

The Petersen graph structure on these roots determines which root pairs are "orthogonal" (disjoint arcs) vs "conflicting" (overlapping arcs), directly connecting to the tournament conflict graph via the Kneser/Johnson duality.

## Source
petersen_lie_bridge_s18.py, Part 1
