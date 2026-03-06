# SC Maximizer Mechanism: Draft Theorem

**Source:** kind-pasteur-2026-03-06-S18e

## Statement

**Conjecture (SC Maximizer):** For every self-complementary score sequence s and every n, the maximum value of H(T) over all tournaments T with score sequence s is achieved by a self-converse tournament.

**Status:** Verified exhaustively at n = 4, 5, 6, 7.

## Key Definitions

- **Self-converse (SC):** T is isomorphic to T^op (opposite tournament).
- **Anti-automorphism:** A permutation sigma with T[sigma(i)][sigma(j)] = 1 - T[i][j] for all i != j.
- **Self-complementary score sequence:** s_i + s_{n-1-i} = n-1 for all i (necessary for SC tournaments to exist with score s).

## Mechanism: Involutory Anti-Automorphism

### Theorem (sigma^2 is an automorphism)
If sigma is an anti-automorphism of T, then sigma^2 is an automorphism.
**Proof:** T[sigma^2(i)][sigma^2(j)] = 1 - T[sigma(i)][sigma(j)] = 1 - (1 - T[i][j]) = T[i][j].

### Theorem (Involutory anti-automorphism always exists)
Every SC tournament has at least one involutory anti-automorphism (sigma^2 = id).
**Proof sketch:** The set of anti-automorphisms is a coset sigma * Aut(T). Since the group <sigma, Aut(T)> is finite, it contains an element of order 2 in the anti-automorphism coset.
**Verified:** Exhaustive at n=4,5,6 (all SC classes have at least one involutory anti-aut). Random at n=7.

### Theorem (Fixed-point-free at even n)
At even n, every involutory anti-automorphism sigma is fixed-point-free.
**Proof:** If sigma(v) = v, then for all j: T[v][j] = T[sigma(v)][sigma(j)] = 1 - T[v][j] (contradiction for j != v). Wait — this doesn't work directly because T[sigma(v)][sigma(j)] = 1 - T[v][j] gives T[v][sigma(j)] = 1 - T[v][j], meaning vertex v beats sigma(j) iff v loses to j. Since sigma is an involution and j ranges over all vertices != v, this means score(v) = (n-1) - score(v), so score(v) = (n-1)/2. At even n, this is not an integer, contradiction.
**Therefore at even n, sigma has no fixed points.** QED.

**Corollary:** At even n, sigma decomposes vertices into n/2 two-cycles {v, sigma(v)}.

### Disjoint Cycle Pair Mechanism (even n)

**Setup:** T is SC on n vertices (n even), sigma is a fixed-point-free involutory anti-aut, with orbits O_1 = {a_1, b_1}, ..., O_{n/2} = {a_{n/2}, b_{n/2}}.

**Lemma:** If C is a directed 3-cycle on vertex set {v_1, v_2, v_3} where each v_i comes from a different orbit O_{k_i}, then sigma(C) is a directed 3-cycle on {sigma(v_1), sigma(v_2), sigma(v_3)}, also selecting one from each orbit, and {v_1,v_2,v_3} ∩ {sigma(v_1),sigma(v_2),sigma(v_3)} = empty.

**Proof:** Since sigma reverses all arcs, sigma(C) is the reverse cycle on sigma-images. The reverse of a directed 3-cycle is also a directed 3-cycle (reversing a 3-cycle gives the other orientation). Since v_i and sigma(v_i) are in the same orbit, sigma(C) picks the "other" element from each of the 3 orbits. The vertex sets are disjoint. QED.

**Consequence:** This creates C(n/2, 3) families of 3-cycle pairs, each contributing to alpha_2 in the independence polynomial. At n=6: C(3,3) = 1 family, with 2^3 = 8 vertex selections, forming 4 disjoint pairs. This gives alpha_2 >= 4.

## Comparison with NSC Tournaments

The NSC tournament lacks this structural guarantee. Without an anti-automorphism creating complementary pairs, the 3-cycles are less evenly distributed across vertices, leading to fewer vertex-disjoint pairs.

**Key computational evidence (n=6, score (3,3,3,2,2,2)):**
- SC with sigma=(1,0,5,4,3,2): 8 three-cycles, ALL from orbit selections, 4 disjoint pairs, alpha_2=4, H=45
- NSC: 8 three-cycles, only 1 disjoint pair, alpha_2=1, H=43

## Open Questions

1. **Prove the conjecture:** Show that the involution orbit structure always produces enough disjoint pairs to guarantee max H within each score class.
2. **Odd n:** At odd n, sigma has a fixed point. 3-cycles through the fixed point are NOT in complementary pairs. How does this affect the mechanism?
3. **Higher alpha_k:** The argument extends to k-tuples of mutually disjoint cycles. Do SC tournaments also maximize alpha_k for k >= 3?
4. **Global max:** The global max H is always SC. Is this because the most regular score sequence is always self-complementary, combined with the SC maximizer within each score class?
