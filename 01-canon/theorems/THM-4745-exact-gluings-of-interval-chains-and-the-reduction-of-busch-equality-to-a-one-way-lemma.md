---
id: THM-4745
title: "Exact gluings of interval chains, and the reduction of Busch's equality case to a one-way lemma. Let A, B be strong tournaments on vertex sets meeting in one vertex v, T a tournament on their union inducing A and B, A' = A - v, B' = B - v. (1) H(T) >= H(A) H(B), with equality iff every Hamiltonian path of T restricts to Hamiltonian paths of A and of B and every pair of such paths has exactly one merge. (2) If A and B are interval chains, every cross arc goes from A' to B', and H(T) = H(A) H(B), then v is the last spine vertex of A and the first spine vertex of B (any vertex of a 3-cycle factor), so T is the concatenated chain; symmetrically when every cross arc goes from B' to A'. (3) FINITE-EXACT one-way lemma: for all chains with at most 6 vertices and (|A|-1)(|B|-1) <= 16, and for C_3, R_4 with the optimal 7-chain IC(3,3), every exact gluing has all its cross arcs in one direction, and every exact gluing is a concatenation. (4) Reduction: if every exact gluing of two Busch-extremal chains A, B with f(|A|) f(|B|) = f(|A| + |B| - 1) has one-way cross arcs, then for every m the strong m-tournaments with f(m) Hamiltonian paths are exactly the optimal interval chains (HYP-9356)"
status: >
  PROVED: statements 1, 2, 4 (all orders; proofs below; statement 1 is Busch's merge argument, Lemma 2 of his note,
  with its equality case made explicit). FINITE-EXACT: statement 3 (exact_gluings.py: every orientation of the cross
  arcs, every gluing vertex). The hypothesis of statement 4 (the one-way lemma for extremal chains) is OPEN beyond
  the finite range; together with THM-4744 (4) (orders <= 12) it is the remaining step of HYP-9356.
  Session opus-2026-10-21-S50; pre-audit.
session: opus-2026-10-21-S50
source: 05-knowledge/results/h_spectrum_multiplicities_20261021.md
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/exact_gluings.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/gluing_failure_modes.py (+ .out)
related:
  - THM-4744 (interval chains: product law, recognition, the extremal classes through order 12)
  - HYP-9356 (Busch's minimum is attained exactly by the optimal interval chains)
  - 05-knowledge/results/spectrum_basins_20261020e/spectrum.md (Busch's theorem, CITED)
---

# THM-4745 — exact gluings of interval chains

## Setting

* `A`, `B` are strong tournaments with `|A|, |B| ≥ 3` and `V(A) ∩ V(B) = {v}`. `T` is any tournament on `V(A) ∪ V(B)` with `T[V(A)] = A` and `T[V(B)] = B`.
* The **cross arcs** are the `(|A|−1)(|B|−1)` arcs between `A' = A − v` and `B' = B − v`.
* The gluing is **exact** if `H(T) = H(A)·H(B)`.
* **Merges.** For Hamiltonian paths `P = P_1 v P_2` of `A` and `Q = Q_1 v Q_2` of `B`, a merge is a Hamiltonian path `R_1 v R_2` of `T` where:
  * `R_1` interleaves `P_1` with `Q_1` and `R_2` interleaves `P_2` with `Q_2`;
  * each interleaving keeps the orders of both sequences.
* **Chains.** Interval chains, spines `x_0, …, x_N` and blocks are as in THM-4744. Block `i` has its start `x_(l_i)`, its end `x_(l_(i+1))`, and the vertices strictly between, its interior.
* Busch's minimum is `f`: Busch, *EJC* 13 (2006) N3 (CITED).

## Statements

1. **Merge inequality (PROVED).**
   * `H(T) ≥ H(A) H(B)`.
   * Equality holds iff both of the following hold:
     * **(R)** every Hamiltonian path of `T` restricts to Hamiltonian paths of `A` and of `B`;
     * **(U)** every pair of Hamiltonian paths of `A` and `B` has exactly one merge.
2. **One-way exact gluings of chains are concatenations (PROVED).**
   * Let `A = IC(σ)` and `B = IC(τ)`, and suppose every cross arc goes from `A'` to `B'` and the gluing is exact.
   * Then `v = x_N` is the last spine vertex of `A` and `v = y_0` is the first spine vertex of `B`, so `T ≅ IC(σ τ)`. If a factor is `C_3 = IC(2)`, every vertex counts as both, by rotation.
   * The case where every cross arc goes from `B'` to `A'` is the converse picture, and gives `IC(τ σ)`.
3. **One-way lemma, finite range (FINITE-EXACT).**
   * Take chains with parts in `{2, 3, 4, 5}`, at most 6 vertices each and `(|A|−1)(|B|−1) ≤ 16`: 28 ordered pairs. Add the four ordered pairs of `C_3` or `R_4` with the optimal 7-chain `IC(3,3)`, for 32 pairs in all.
   * For every gluing vertex and every cross orientation, Busch's inequality holds.
   * Every exact gluing has one-way cross arcs and is a concatenation.
   * Both concatenations of every pair are exact.
   * There are 136 exact gluings in all (`exact_gluings.py`, 443 s).
4. **Reduction (PROVED).** Suppose that for all `a, b ≥ 3` with `f(a)f(b) = f(a+b−1)`, every exact gluing of a Busch-extremal chain of order `a` with one of order `b` has one-way cross arcs. Then for every `m ≥ 3` the strong `m`-tournaments with `H = f(m)` are exactly the optimal interval chains. That is HYP-9356, and with THM-4744 (2) it gives `E(3q) = q`, `E(3q+1) = 1`, `E(3q+2) = q(q+3)/2`.

## Proofs

**1.**
* **The inequality.** Lemma 2 of Busch's note merges any `P_1` with `Q_1` and any `P_2` with `Q_2`. Since `last(R_1)` is the last vertex of `P_1` or of `Q_1`, both of which precede `v` in their paths, `last(R_1) → v`. Likewise `v → first(R_2)`. So every pair has at least one merge.
* **Restriction recovers the pair.** A merge restricts to `A` and to `B` as the original pair. So the merges of different pairs are different paths.
* **Counting.** Hence `H(T)` equals:
  * the number of pairs, which is `H(A)H(B)`;
  * plus the surplus merges of pairs with more than one merge;
  * plus the Hamiltonian paths of `T` that do not restrict, that is, those violating (R).
* So equality holds iff both surpluses are zero. ∎

**2.** Assume every cross arc goes from `A'` to `B'` and the gluing is exact.

* **Shape of a path.** No vertex of `B'` is followed by a vertex of `A'` on a path of `T`. So every Hamiltonian path of `T` is `X_a X_b v Y_a Y_b`, where each `X` and `Y` is a block of vertices of `A'` or of `B'`.
* **A bad shape.** If `X_a` and `X_b` are nonempty and `v → last(X_a)`, the path violates (R), because its restriction to `A` contains the step `last(X_a), v`.
* **The B-side ingredient.** Since `B` is strong it has a Hamiltonian cycle. Deleting `v` from it gives a Hamiltonian path `X_b` of `B'` that ends at an in-neighbour of `v`.
* **Bad splits.** Hence (R) forbids every **bad split** of `A'`: a pair of disjoint paths `X_a ≠ ∅` and `Y_a` (possibly empty) covering `A'`, with `v → last(X_a)` and (`Y_a = ∅` or `v → first(Y_a)`). Given one, `X_a X_b v Y_a` is a Hamiltonian path of `T`:
  * `last(X_a) → first(X_b)` is a cross arc;
  * `last(X_b) → v`;
  * `v → first(Y_a)`.

  But its restriction to `A` is not a path.

Let `A = IC(p_1, …, p_k)` with `N ≥ 3`, and let `v = x_t`. A bad split exists unless `t = N`:

* **(a) `A − v` is strong.** Then `A'` has a Hamiltonian cycle, so it has a Hamiltonian path ending at any out-neighbour of `v`. Take `Y_a = ∅`. `A − v` is strong in three cases:
  * `v` is interior to a block with `p_j ≥ 3`: `A − v = IC(…, p_j − 1, …)`.
  * `v` is the single interior vertex of a block `j < k` with `p_j = 2`: the block start `x_(l_j)` reaches the block end through block `j+1` and its back arc.
  * `v` is the single interior vertex of the last block, `k ≥ 2`: `x_N → x_(l_k)`, and the back arcs lead from there to `x_0`.
* **(b) `v = x_(l_j)` is the start of a block `j < k`**, including `x_0` when `k ≥ 2`.
  * The vertices `x_(l_(j+1)), …, x_N` form a strong chain, and every arc between them and the rest of `A'` points into them. So they form the terminal strong component of `A − v`, and they contain `x_N`.
  * `v → x_N`: this arc is forward, and it is not the back pair of the last block, since `j < k`.
  * In a tournament every vertex of the terminal component ends some Hamiltonian path. Take that path as `X_a`, with `Y_a = ∅`.
* **(c) `v = x_(l_k)` is the start of the last block.**
  * If `p_k ≥ 3`: take `X_a = x_0 ⋯ x_(l_k − 1) x_(l_k + 1)` (just `(x_1)` if `l_k = 0`) and `Y_a = x_(l_k + 2) ⋯ x_N`. The arcs used are spine arcs and forward arcs into interior vertices. `v` beats both `x_(l_k + 1)` and `x_(l_k + 2)`.
  * If `p_k = 2` and `k ≥ 2`: with `c` the interior vertex, take `X_a = x_0 ⋯ x_(l_(k−1) − 1) c` (just `(c)` if `l_(k−1) = 0`) and `Y_a = x_(l_(k−1)) ⋯ x_(l_k − 1) x_N`. Here `v → c` is forward, and `v → x_(l_(k−1))` is the back arc of block `k − 1`.

So `v = x_N`.

* **The `B` side.** The converse of `T` reverses every cross arc. By THM-4744 (2) it is the gluing of the converses `IC(σ reversed)` and `IC(τ reversed)`, whose spines are the reversed spines, and it is still exact. Applying the argument with the roles exchanged shows that `v` is the last spine vertex of the converse of `B`, which is `y_0`.
* **Conclusion.** With every cross arc forward, `v = x_N = y_0` and the back arcs of both chains, `T` is `IC(σ τ)` by definition. ∎

**4.** Induction on `m`.
* **Base and the nearly transitive case.**
  * For `m ≤ 5` the census agrees, and THM-4744 (4) covers `m ≤ 12`.
  * For `m ≥ 6`, the nearly transitive tournament has `2^(m−2) + 1 > f(m)` paths, so it is not extremal.
* **The decomposition.**
  * By Busch's Lemma 1, an extremal strong `T` therefore splits as `A ∪ B` with `A`, `B` strong, `|A ∩ B| = 1` and `|A|, |B| ≥ 3`.
  * Then `f(m) = H(T) ≥ H(A)H(B) ≥ f(a)f(b) ≥ f(m)`. The last step holds because `f(a)f(b)` is a product `3^x 5^y` with `2x + 3y = m − 1`.
  * So `A` and `B` are extremal, `f(a)f(b) = f(m)`, and the gluing is exact.
* **Conclusion.**
  * By induction `A` and `B` are optimal chains.
  * By hypothesis the cross arcs are one-way, so by statement 2, `T` is a concatenated chain. It has `f(m)` paths, so its composition is optimal.
  * The converse direction is THM-4744 (1). ∎

## What this does and does not say

* It does not prove the one-way lemma.
  * Mixed cross arcs can still give unique merges for a single pair: `(q_1, q_2)` merged with `(p_1)` has the unique merge `q_1 p_1 q_2` when `q_1 → p_1 → q_2`. So pairwise uniqueness alone cannot exclude mixed arcs.
  * Conversely, (R) alone does not suffice either. In `gluing_failure_modes.py`, at every gluing position of `C_3 ∪ C_3`, four of the fourteen mixed patterns violate only (U), with no path breaking (R).
  * With an `R_4` factor, at the end-to-start positions, every mixed pattern has a path breaking (R). The least surplus is 4 for `C_3` with `R_4`, and 6 for `R_4` with `R_4`.
  * A proof must therefore use both conditions.
* The finite census covers all chain pairs whose sizes enter the induction up to `m = 10`:
  * `(3,3) → 5`, `(3,4) → 6`, `(4,4) → 7`;
  * `(3,6), (4,5) → 8`;
  * `(4,6), (3,7) → 9`;
  * `(4,7) → 10`.

  THM-4744 (4) then verifies the classification directly through `m = 12`.
