# THM-262: Dual Lie Algebra Embedding of Tournaments

**Status:** PROVED (structural theorem)
**Filed by:** kind-pasteur-2026-03-21-S18

## Statement

A tournament T on [n] simultaneously embeds in two Lie algebras:

**(A) In so(n):** The skew-adjacency matrix B_T = sum_{i<j} eps_{ij}(E_{ij} - E_{ji}) is an element of so(n). The root system of so(n) is B_{floor(n/2)} (n odd) or D_{n/2} (n even).

**(B) In sl(n):** The weight w(T) = sum_{i<j} eps_{ij}(e_i - e_j) is an element of the root lattice of A_{n-1} = sl(n). In coordinates, w(T) = (d_0, ..., d_{n-1}) where d_k = 2*score(k) - (n-1) is the score deviation.

The bridge: the standard basis of so(n), namely {E_{ij} - E_{ji} : i<j}, coincides with the root vectors of sl(n). Thus dim so(n) = C(n,2) = # positive roots of A_{n-1} for all n.

## Key properties

1. **Weight = score deviation:** w(T) = (2s_0-(n-1), ..., 2s_{n-1}-(n-1)) where s_k is the score of vertex k. Sum = 0 (traceless).

2. **Zero weight = regular tournament:** w(T) = 0 iff T is regular. At n=5, ALL 24 regular tournaments have H=15 (Paley max).

3. **Weight norm = score spread:** ||w||^2 in the A_{n-1} Killing form measures how far T is from regularity. At n=5:
   - ||w||^2 = 0: H = 15 (regular only)
   - ||w||^2 = 8: H in {11, 13, 15}
   - ||w||^2 = 16: H = 9
   - ||w||^2 = 24: H = 5
   - ||w||^2 = 32: H = 3
   - ||w||^2 = 40: H = 1 (transitive)

   The weight norm ALMOST determines H at n=5.

4. **Casimir invariants:** Tr(B^2) = -n(n-1) is CONSTANT for all tournaments. Tr(B^4) varies and gives 2 classes at n=5. The joint (Tr(B^4), ||w||^2) gives 8 classes at n=5 (still doesn't fully determine H).

## Significance

This "fundamental triangle" — A_{n-1} positive roots, K(n,2) orthogonality, so(n) basis — unifies three perspectives on tournament structure:
- Combinatorial: arcs as signed roots
- Algebraic: skew-adjacency in so(n)
- Graph-theoretic: Kneser/Johnson duality

The weight norm anticorrelation with H (higher norm = lower H, with equality for transitive and regular extremes at n=5) suggests that Lie-algebraic distance from zero weight controls Hamiltonian path count.

## Source
petersen_lie_bridge_s18.py, Parts 2, 3, 5
