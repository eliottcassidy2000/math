# Even-graph bridge and limitations

**Source.** G. F. Royle, C. E. Praeger, S. P. Glasby, S. D. Freedman and A. Devillers, “Tournaments and even graphs are equinumerous,” *Journal of Algebraic Combinatorics* 57 (2023), 515–524, https://doi.org/10.1007/s10801-022-01197-0 ; https://arxiv.org/abs/2204.01947 .

## Precise definition and result

Fix an arbitrary orientation of the edges of a finite simple graph X. For an automorphism g, count edges e whose image orientation opposes the chosen orientation of g(e). Its parity is independent of the chosen orientation. X is **even** if every automorphism has even count, and **odd** otherwise. This is not the standard notion of a graph with every vertex degree even.

Theorem 1.1 of Royle et al. states that, for n ≥ 2, the number of isomorphism classes of n-vertex even graphs equals the number of isomorphism classes of n-vertex tournaments. Their proof uses the Cauchy–Frobenius theorem and double counting to obtain #graphs = #tournaments + #odd graphs, then subtracts. They explicitly leave a natural bijection as an open problem in their paper. This is a statement about what that 2023 paper proves, not a claim that no later bijection exists. The identity itself gives no natural free 2-group action or orbit decomposition.

## Rigorous obstructions to simple candidate symmetries

1. **Graph complement is not closed on even graphs.** The empty graph on n ≥ 2 vertices is even, because it has no edges. Its complement K_n is odd: a transposition of two vertices reverses exactly their joining edge modulo 2. For example, choose the transposed pair's edges to every other vertex compatibly; then those edges contribute zero reversals, while the joining edge contributes one. Thus complement cannot define an involution of even-graph isomorphism classes.
2. **Reversing all tournament arcs is not fixed-point-free on isomorphism classes.** The transitive tournament, with vertices ordered 1 < ⋯ < n and edges directed forward, is isomorphic to its reverse by i ↦ n+1−i. This holds for every n.
3. **The n = 2 endpoint is an absolute divisibility obstruction.** There is exactly one tournament isomorphism class, so a(2) = 1 is not divisible by 2^floor(2/2) = 2; an orbit decomposition into blocks of that size is impossible there. Any such conjecture must exclude n = 2.
4. **No label-preserving bijection of the complete labeled sets at n = 3.** There are 2^3 = 8 labeled tournaments. The even graphs are exactly the empty graph and the three labeled copies of P3, hence four. For P3, the only nonidentity automorphism swaps its two leaves; orient both edges toward the center and it reverses neither. A one-edge graph and K3 are odd, witnessed by swapping the endpoints of an edge (orient the remaining K3 edges compatibly). Consequently a bijection between unlabeled classes cannot arise as a bijection between these labeled sets.

These counterexamples rule out only the stated simple constructions; they do not rule out a different natural symmetry or prove any general divisibility result for n ≥ 3.
