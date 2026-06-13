---
theorem: THM-332
name: Dominant Vertex Position Theorem
status: PROVED
session: opus-2026-05-27-S1
verified: computationally n=3..5 (all deg=n-1 cases)
---

## Statement

Let T be a tournament on n vertices and Q a vertex with d⁺(Q) = n−1 (Q beats all other vertices — Q is the **absolute source**). Then:

1. Q appears at **position 0** (first position) in every directed Hamiltonian path of T.
2. H(T) = H(T−Q), i.e., deleting Q does not change the Hamiltonian path count.
3. Q is NOT contained in any directed odd cycle.

## Proof

(1): Every HP is a permutation (v₀, v₁, ..., v_{n-1}) where v_i → v_{i+1}. Since Q beats all others, the only valid position for Q is position 0 (any later position would require some other vertex v to beat Q, contradicting d⁺(Q) = n−1).

(2): Every HP starts Q → (remaining n−1 vertices in some HP-order of T−Q). The map HP of T ↦ its suffix gives a bijection with HP of T−Q. Hence H(T) = H(T−Q).

(3): Any odd cycle through Q would require some vertex to beat Q. None exist. □

## Duality (Absolute Sink)

Let P be the **absolute sink**: d⁺(P) = 0 (P beaten by all others). Then:
- P appears at position n−1 (last) in every HP of T.
- H(T) = H(T−P).
- P is not in any odd cycle.

## Connection to King Theorem

If d⁺(Q) = n−1, then Q is trivially a king (reaches everything in 1 step, not 2). The King theorem's "2 steps" is vacuous for Q. More interesting: the unique vertex P with d⁺(P) = 0 is the vertex that the king theorem applies to most weakly — P is not a king at all (it can't reach anyone).

## Computation

Verified at n=3..5: for all tilings where the source vertex (n-1 in 0-indexed) beats all others (all "vertical leg" tiles in the default direction), Q always appears first in 100% of Hamiltonian paths (frac_first = 1.000, frac_last = 0.000). □
