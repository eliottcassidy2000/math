# THM-262: Dual Lie-Algebra Carriers of Tournaments

**Status:** PROVED formulas; carrier terminology corrected by MISTAKE-507
**Filed by:** kind-pasteur-2026-03-21-S18

## Statement

A tournament `T` on `[n]` has two related but inequivalent Lie-theoretic
carriers:

**(A) Faithful `so(n)` encoding:** The skew-adjacency matrix
`B_T=sum_{i<j} eps_{ij}(E_{ij}-E_{ji})` is an element of `so(n)` and uniquely
recovers every arc sign.  For `n>=3`, the root system of `so(n)` is type
`B_floor(n/2)` for odd `n` and `D_(n/2)` for even `n`.

**(B) Nonfaithful `A_(n-1)` score projection:** The weight
`w(T)=sum_{i<j} eps_{ij}(e_i-e_j)` is an element of the `A_(n-1)` root
lattice.  In coordinates, `w(T)=(d_0,...,d_(n-1))`, where
`d_k=2*score(k)-(n-1)`.  Different tournaments can have the same score
weight, so this is not an embedding of the tournament set.

**The bridge:** both the skew basis `E_{ij}-E_{ji}` and the positive
`A_(n-1)` roots are canonically indexed by the `C(n,2)` unordered pairs.
They do not coincide as vectors or as root systems.  The shared pair index is
the exact carrier; passing from `B_T` to `w(T)` destroys all information not
visible in the score sequence.

## Key properties

1. **Weight = score deviation:** w(T) = (2s_0-(n-1), ..., 2s_{n-1}-(n-1)) where s_k is the score of vertex k. Sum = 0 (traceless).

2. **Zero weight = regular tournament:** w(T) = 0 iff T is regular. At n=5, ALL 24 regular tournaments have H=15 (Paley max).

3. **Weight norm = score spread:** `||w||^2` for the standard normalized
   `A_(n-1)` root inner product (roots have squared length two) measures how
   far `T` is from regularity. At `n=5`:
   - ||w||^2 = 0: H = 15 (regular only)
   - ||w||^2 = 8: H in {11, 13, 15}
   - ||w||^2 = 16: H = 9
   - ||w||^2 = 24: H = 5
   - ||w||^2 = 32: H = 3
   - ||w||^2 = 40: H = 1 (transitive)

   The weight norm ALMOST determines H at n=5.

4. **Casimir invariants:** Tr(B^2) = -n(n-1) is CONSTANT for all tournaments. Tr(B^4) varies and gives 2 classes at n=5. The joint (Tr(B^4), ||w||^2) gives 8 classes at n=5 (still doesn't fully determine H).

## Significance

This typed triangle -- `A_(n-1)` positive-root supports, `K(n,2)`
orthogonality, and the pair-indexed `so(n)` basis -- relates three
perspectives on tournament structure without identifying their ambient
vectors:
- Combinatorial: arcs as signed roots
- Algebraic: skew-adjacency in so(n)
- Graph-theoretic: Kneser/Johnson duality

The weight norm anticorrelation with H (higher norm = lower H, with equality for transitive and regular extremes at n=5) suggests that Lie-algebraic distance from zero weight controls Hamiltonian path count.

## Source
petersen_lie_bridge_s18.py, Parts 2, 3, 5
