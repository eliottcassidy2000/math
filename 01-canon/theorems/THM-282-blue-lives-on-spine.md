# THM-282: Blue Lives on the Spine; Ribs Are Always Black

**Status:** PROVED  
**Session:** opus-2026-04-03-S27  
**Depends on:** THM-280 (grid reflection = complement), THM-281 (SC sizes odd)  

## Statement

For complement-flip edges (d=m) in the tournament metagraph:

1. **B ⊆ SC-SC**: Every pure-blue edge connects two self-complementary classes.
2. **M ⊆ SC-SC**: Every mixed edge connects two self-complementary classes.
3. **SC-NS ⊆ K**: Every edge between an SC class and an NS class is pure-black.
4. **NS-NS ⊆ K**: Every edge between two NS classes is pure-black.

Equivalently: the blue/mixed color diversity exists ONLY on the SC backbone (spine). The ribs (SC-NS) and sea (NS-NS) are uniformly black.

## Proof

### Lemma: Grid-symmetric tilings are closed under complement flip.

If tiling b is grid-symmetric (R(b) = b where R is the grid reflection), then the complement tiling b' (flip all m bits) is also grid-symmetric.

Proof: isGridSym(b) checks that b_i = b_{R(i)} for all paired tile positions. Flipping all bits preserves this equality: (1-b_i) = (1-b_{R(i)}) iff b_i = b_{R(i)}. □

### Lemma: Grid-symmetric tilings only exist in SC classes.

If b is grid-symmetric, then T(b) ≅ T(b)^op (by THM-280: the grid reflection sends b to a tiling of T^op, but if b is a fixed point of R, then T(b) = T(R(b)) = σ(T(b)^op), hence T(b) ≅ T(b)^op). So the isomorphism class of T(b) is self-complementary. □

### Corollary: Grid-symmetric tilings map to SC classes under complement.

If b is grid-symmetric, its complement b' is grid-symmetric (Lemma 1), and T(b') is in an SC class (Lemma 2). So the complement of a gs tiling always lands in an SC class. □

### Proof of (1) B ⊆ SC-SC:

A pure-blue edge between classes A and B means there exist grid-symmetric tilings in A whose complement lands in B, and NO non-gs tilings in A map to B. By Lemma 2, A must contain gs tilings, so A is SC. By the Corollary, the complement of those gs tilings lands in an SC class, so B is SC. Therefore both A and B are SC. □

### Proof of (3) SC-NS ⊆ K:

An edge between SC class A and NS class B has lines from tilings in A whose complement lands in B. If any such tiling were grid-symmetric, its complement would land in an SC class (Corollary), but B is NS — contradiction. So all such tilings are non-gs, meaning all lines are black. Therefore the edge is pure-black. □

### Proof of (4) NS-NS ⊆ K:

NS classes contain zero grid-symmetric tilings (verified: NS = pure-black purity). So all complement lines from NS tilings are black. □

### Proof of (2) M ⊆ SC-SC:

A mixed edge has at least one blue line (from a gs tiling) and one black line. The blue line requires both endpoints to be SC (by the same argument as (1)). □

## Verified

Computationally confirmed at n = 3, 4, 5, 6, 7. Zero exceptions.

## Consequence

The three-part decomposition of the merged metagraph:
- **Spine (SC-SC)**: carries B, M, and some K edges — the only place with color diversity
- **Ribs (SC-NS)**: uniformly K (black) — the bridge between spine and sea  
- **Sea (NS-NS)**: uniformly K (black) — the dominant bulk at large n

Filtering to show only blue (B) edges reveals a subgraph of the SC backbone — the "colored skeleton" of the metagraph. This skeleton vanishes as n grows (SC fraction → 0).
