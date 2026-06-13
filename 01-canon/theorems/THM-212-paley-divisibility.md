# THM-212: p divides H(T) for Circulant Tournaments on Prime Vertices

**Status:** PROVED
**Found by:** opus-2026-03-14-S89c
**Verified in:** `04-computation/pi_maxH_89c.py`

## Statement

For every prime p and every circulant tournament T on Z/pZ (including the Paley tournament P_p when p ≡ 3 mod 4):

**p | H(T)**

More generally: if a tournament T on n vertices has an automorphism of prime order p dividing n, then p | H(T).

## Verified Cases

| p | H(P_p) | H(P_p)/p |
|---|--------|----------|
| 3 | 3 | 1 |
| 7 | 189 = 3³×7 | 27 |
| 11 | 95095 = 5×7×11×13×19 | 8645 |
| 19 | 1172695746915 | 61720828785 |

## Connection to THM-209

By THM-209, H(P_p) = IP(G(P_p), 2) where G(P_p) is the odd-cycle disjointness graph.

The Paley tournament P_p has a transitive Z/pZ symmetry group (cyclic shift i ↦ i+1 mod p). This symmetry acts on the odd cycles and hence on independent sets. If the action has no fixed points (other than the empty set), then each orbit has size p, and H(P_p) - 1 = 2·(sum of cycle counts) + 4·(disjoint pair counts) + ... would be divisible by p if all cycle counts and packing counts are divisible by p.

Since t_k(P_p) ≡ 0 (mod p) for all k (because the cyclic shift permutes cycles in orbits of size p), the level-1 contribution 2·Σt_k is divisible by 2p. Similarly for higher levels.

## Proof Sketch (needs verification)

1. The cyclic shift σ: i ↦ i+1 mod p is an automorphism of P_p.
2. σ acts freely on the set of directed odd cycles (no cycle is fixed by σ).
3. Therefore each orbit has size p, so t_k ≡ 0 (mod p) for all k.
4. σ also acts freely on disjoint k-packings, so d_{...} ≡ 0 (mod p).
5. Therefore H = 1 + 2·(Σt_k) + 4·(Σd) + 8·(Σtriplets) + ...
6. Each term after the initial 1 is divisible by 2p.
7. So H ≡ 1 (mod 2p)... but we need H ≡ 0 (mod p), not H ≡ 1 (mod p).

**ISSUE**: Step 7 gives H ≡ 1 (mod p), but we OBSERVE H ≡ 0 (mod p). This means either:
- The proof sketch has an error (likely: the action on Hamiltonian paths, not on the IP decomposition, gives the divisibility)
- OR: there's a subtler mechanism.

**Better approach**: σ acts on Hamiltonian paths directly. Each orbit of paths under σ has size p (since σ has order p and p is prime). If NO path is fixed by σ, then p | H(P_p). A path π is fixed by σ iff π(i+1) = π(i) + 1 mod p for all i, which gives only one path (the identity). But actually σ acts by relabeling: (σπ)(i) = π(i)+1 mod p. A path v₁→v₂→...→vₚ is mapped to (v₁+1)→(v₂+1)→...→(vₚ+1). This is still a Hamiltonian path (since P_p is circulant). The path is fixed iff vᵢ+1 = vᵢ₊₁... hmm, this isn't quite right.

The correct statement: σ maps the path (v₁,...,vₚ) to (v₁+1,...,vₚ+1) mod p. This is a DIFFERENT path (different vertex sequence) unless all vᵢ shift to give the same path. Since p is prime and the shift has order p, either the orbit has size p or the path is a "periodic" sequence, which is impossible for a permutation when p is prime.

**Therefore p | H(P_p).** The proof is: the cyclic automorphism group Z/pZ acts freely on Hamiltonian paths, partitioning them into orbits of size p. □

## Additional Observations

- H(P_p) achieves the MAXIMUM among all tournaments on p vertices for p = 3, 7, 11 (OEIS A038375)
- But H(P_19) < H(C_19) (cyclic beats Paley at p=19!)
- H(P_p) / Mean(H) ≈ 2.4-2.5 for large p
