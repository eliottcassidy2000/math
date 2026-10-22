---
id: THM-4745
title: "Exact gluings of interval chains, and an independent computer-assisted reproof of the Moon-Yang uniqueness theorem for Busch's equality case. Let A, B be tournaments on vertex sets meeting in one vertex v, T a tournament on their union inducing A and B, A' = A - v, B' = B - v. (1) H(T) >= H(A) H(B), with equality iff every Hamiltonian path of T restricts to Hamiltonian paths of A and of B and every pair of such paths has exactly one merge. (2) If A and B are interval chains, every cross arc goes from A' to B', and H(T) = H(A) H(B), then v is the last spine vertex of A and the first spine vertex of B (any vertex of a 3-cycle factor), so T is the concatenated chain; symmetrically when every cross arc goes from B' to A'. (3) FINITE-EXACT: for all chains with at most 6 vertices and (|A|-1)(|B|-1) <= 16, and for C_3, R_4 with the optimal 7-chain IC(3,3), every exact gluing is one-way and a concatenation. (4) Reduction: if exact gluings of Busch-extremal chains of compatible orders are one-way, the Busch-extremal strong tournaments are exactly the optimal interval chains. (5) The one-way lemma holds for all chains whose parts lie in {2, 3, 4}: every exact gluing of two such chains is one-way and a concatenation (PROVED, computer-assisted at one finite step). (6) Hence (KNOWN: Moon-Yang 2007) for every m the strong m-tournaments with exactly f(m) Hamiltonian paths are exactly the optimal interval chains"
status: >
  PROVED: statements 1, 2, 4 (all orders; statement 1 is Busch's merge argument, Lemma 2 of his note, with its equality
  case made explicit). PROVED, computer-assisted: statement 5. Its only finite input is Step 3, the single-block end
  positions for blocks of 3, 4, 5 vertices, which exact_gluings.py checks exhaustively; oneway_proof_checks.py confirms
  the excluded configurations by exact counts. FINITE-EXACT: statement 3.
  Statement 6 is KNOWN: J. W. Moon and L. L. M. Yang, Utilitas Math. 72 (2007) 193-198, CITED via El Sahili and Ghazo
  Hanna (arXiv:2101.00713); the primary paper was not accessed. The proof here is an independent reproof, not a new
  result.
  The S50 audit confirmed statements 1-4: every bad-split construction was checked at all 745 positions for N = 3..10,
  and the census was reproduced. Statements 5 and 6 were added after that audit and have had a separate audit pass
  (see below). The novelty first claimed for the reduction was retracted (MISTAKE-652).
  Session opus-2026-10-21-S50.
session: opus-2026-10-21-S50
source: 05-knowledge/results/h_spectrum_multiplicities_20261021.md
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/exact_gluings.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/gluing_failure_modes.py (+ .out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/oneway_proof_checks.py (+ .out)
related:
  - THM-4744 (Moon's interval chains: product law, recognition, the extremal classes through order 13)
  - HYP-9356 (the Moon-Yang characterization; KNOWN)
  - 05-knowledge/results/spectrum_basins_20261020e/spectrum.md (Busch's theorem, CITED)
---

# THM-4745 — exact gluings of interval chains, and a reproof of Moon–Yang

## Setting

* `A`, `B` are tournaments with `|A|, |B| ≥ 3` and `V(A) ∩ V(B) = {v}`. `T` is any tournament on `V(A) ∪ V(B)` with `T[V(A)] = A` and `T[V(B)] = B`.
* The **cross arcs** are the `(|A|−1)(|B|−1)` arcs between `A' = A − v` and `B' = B − v`. "Forward" means `A' → B'`.
* The gluing is **exact** if `H(T) = H(A)·H(B)`.
* **Merges.** For Hamiltonian paths `P = P_1 v P_2` of `A` and `Q = Q_1 v Q_2` of `B`, a merge is a Hamiltonian path `R_1 v R_2` of `T` where:
  * `R_1` interleaves `P_1` with `Q_1` and `R_2` interleaves `P_2` with `Q_2`;
  * each interleaving keeps the orders of both sequences.
* **Chains.** Interval chains, spines `x_0, …, x_N` (for `B`: `y_0, …, y_M`) and blocks are as in THM-4744. Block `i` has its start `x_(l_i)`, its end `x_(l_(i+1))`, and the vertices strictly between, its interior. `C_3 = IC(2)` and `R_4 = IC(3)`.
* Busch's minimum is `f`: Busch, *EJC* 13 (2006) N3 (CITED).
* For a tournament `X`:
  * `Init(X)` is its initial strong component, the possible first vertices of its Hamiltonian paths;
  * `Term(X)` is its terminal strong component, the possible last vertices.

## Statements

1. **Merge inequality (PROVED).** `H(T) ≥ H(A) H(B)`, with equality iff both of the following hold:
   * **(R)** every Hamiltonian path of `T` restricts to Hamiltonian paths of `A` and of `B`;
   * **(U)** every pair of Hamiltonian paths of `A` and `B` has exactly one merge.

   No strongness is needed.
2. **One-way exact gluings of chains are concatenations (PROVED).**
   * Let `A = IC(σ)` and `B = IC(τ)`, and suppose every cross arc goes from `A'` to `B'` and the gluing is exact.
   * Then `v = x_N` and `v = y_0`, so `T ≅ IC(σ τ)`. If a factor is `C_3`, every vertex counts as both, by rotation.
   * The backward case is the converse picture, and gives `IC(τ σ)`.
3. **Finite range (FINITE-EXACT).**
   * Take chains with parts in `{2, 3, 4, 5}`, at most 6 vertices each and `(|A|−1)(|B|−1) ≤ 16`: 28 ordered pairs. Add the four ordered pairs of `C_3` or `R_4` with `IC(3,3)`, for 32 pairs in all.
   * For every gluing vertex and every cross orientation, every exact gluing is one-way and a concatenation.
   * There are 136 exact gluings in all (`exact_gluings.py`). The audit reproduced this, plus 8 further pairs.
4. **Reduction (PROVED).** If every exact gluing of two Busch-extremal chains of orders `a, b` with `f(a)f(b) = f(a+b−1)` is one-way, then for every `m` the strong `m`-tournaments with `H = f(m)` are exactly the optimal chains.
5. **One-way lemma (PROVED, computer-assisted).** Let `A`, `B` be chains whose parts lie in `{2, 3, 4}`, so every block has 3, 4 or 5 vertices. If `T = A ∪_v B` is exact, then either every cross arc is forward with `v = x_N = y_0`, or every cross arc is backward with `v = x_0 = y_M`.
6. **The Moon–Yang theorem (KNOWN; reproved).** For every `m ≥ 3`, the strong `m`-tournaments with exactly `f(m)` Hamiltonian paths are exactly the optimal interval chains:
   * the parts lie in `{2, 3, 4}`, with product `f(m)`;
   * `E(3q) = q`, `E(3q+1) = 1`, `E(3q+2) = q(q+3)/2`;
   * self-converse iff palindromic.

## Proofs

**1.**
* **The inequality.** Lemma 2 of Busch's note merges any `P_1` with `Q_1` and any `P_2` with `Q_2`. Since `last(R_1)` is the last vertex of `P_1` or of `Q_1`, both of which precede `v` in their paths, `last(R_1) → v`. Likewise `v → first(R_2)`.
* **Restriction recovers the pair.** A merge restricts to `A` and to `B` as the original pair.
* **Counting.** So `H(T)` equals:
  * the number of pairs, which is `H(A)H(B)`;
  * plus the surplus merges of pairs with more than one merge;
  * plus the Hamiltonian paths of `T` that violate (R).
* Hence the inequality, with equality iff (R) and (U). ∎

**2.** Assume every cross arc goes from `A'` to `B'` and the gluing is exact.

* **Shape of a path.** No vertex of `B'` is followed by a vertex of `A'` on a path of `T`. So every Hamiltonian path of `T` is `X_a X_b v Y_a Y_b`, where each `X` and `Y` is a block of vertices of `A'` or of `B'`.
* **The B-side ingredient.** Since `B` is strong, deleting `v` from a Hamiltonian cycle of `B` gives a Hamiltonian path `X_b` of `B'` that ends at an in-neighbour of `v`.
* **Bad splits.** Hence (R) forbids every **bad split** of `A'`: disjoint paths `X_a ≠ ∅` and `Y_a` (possibly empty) covering `A'`, with `v → last(X_a)` and (`Y_a = ∅` or `v → first(Y_a)`).
  * Given one, `X_a X_b v Y_a` is a Hamiltonian path of `T`: `last(X_a) → first(X_b)` is a cross arc, `last(X_b) → v`, and `v → first(Y_a)`.
  * Its restriction to `A` is not a path.

Let `A = IC(p_1, …, p_k)` with `N ≥ 3`, and let `v = x_t`. A bad split exists unless `t = N`:

* **(a) `A − v` is strong.** Then `A'` has a Hamiltonian path ending at any out-neighbour of `v`; take `Y_a = ∅`. `A − v` is strong in three cases:
  * `v` is interior to a block with `p_j ≥ 3`: `A − v = IC(…, p_j − 1, …)`;
  * `v` is the single interior vertex of a block `j < k` with `p_j = 2`;
  * `v` is the single interior vertex of the last block, `k ≥ 2`.
* **(b) `v = x_(l_j)` is the start of a block `j < k`**, including `x_0` when `k ≥ 2`.
  * The vertices `x_(l_(j+1)), …, x_N` form the terminal strong component of `A − v`, and they contain `x_N`.
  * `v → x_N`.
  * So take `X_a` to be a Hamiltonian path of `A'` ending at `x_N`, with `Y_a = ∅`.
* **(c) `v = x_(l_k)` is the start of the last block.**
  * If `p_k ≥ 3`: take `X_a = x_0 ⋯ x_(l_k − 1) x_(l_k + 1)` (just `(x_1)` if `l_k = 0`) and `Y_a = x_(l_k + 2) ⋯ x_N`.
  * If `p_k = 2` and `k ≥ 2`: with `c` the interior vertex, take `X_a = x_0 ⋯ x_(l_(k−1) − 1) c` and `Y_a = x_(l_(k−1)) ⋯ x_(l_k − 1) x_N`. Here `v → x_(l_(k−1))` is the back arc of block `k − 1`.

So `v = x_N`. The converse of `T` reverses every cross arc and every spine (THM-4744 (2)), which gives `v = y_0`. With every cross arc forward, `T` is `IC(σ τ)` by definition. ∎

**4.** Induction on `m`.
* **Base.** For `m ≤ 5` the census agrees.
* **Not nearly transitive.** For `m ≥ 6` the nearly transitive tournament has `2^(m−2) + 1 > f(m)` paths.
* **The decomposition.** By Busch's Lemma 1, an extremal `T` splits as `A ∪ B` with `A`, `B` strong, `|A ∩ B| = 1` and `|A|, |B| ≥ 3`. Then:
  * `f(m) = H(T) ≥ H(A)H(B) ≥ f(a)f(b) ≥ f(m)`;
  * the last step holds because `f(a)f(b)` is a product `3^x 5^y` with `2x + 3y = m − 1`;
  * so `A` and `B` are extremal, and therefore optimal chains by induction, and the gluing is exact.
* **Conclusion.** By hypothesis the gluing is one-way, so by statement 2, `T` is a chain. It has `f(m)` paths, so its composition is optimal.
* **The converse direction** is THM-4744 (1). ∎

**5.** Induction on `|T|`. Let `T = A ∪_v B` be exact, with `A`, `B` chains whose parts lie in `{2, 3, 4}`.

*Step 0: tools.*

* **(T1)** Statement 1, for any two subtournaments meeting in exactly one vertex.
* **(T2)** A chain splits at a block boundary `w` as `A_1 ∪_w A_2`, and `H(A) = H(A_1)H(A_2)` (THM-4744 (1)).
* **(T3) Single-path constraints.** Each is witnessed by a Hamiltonian path of `T` that violates (R).
  * Strong factors have Hamiltonian paths starting and ending anywhere (Camion).
  * A Hamiltonian path of `A'` ending at `a` exists iff `a ∈ Term(A − v)`; one starting at `a` exists iff `a ∈ Init(A − v)`. The same holds for `B'`.

  | constraint | forbidden configuration | witness | broken restriction |
  |---|---|---|---|
  | (F1) | a forward arc `a → b` with `v → a` and `a ∈ Term(A − v)` | `P Q`, with `P` a path of `A'` ending at `a` and `Q` a path of `B` starting at `b` | `H|A = P v` |
  | (F2) | a forward arc `a → b` with `b → v` and `b ∈ Init(B − v)` | `P Q`, with `P` a path of `A` ending at `a` and `Q` a path of `B'` starting at `b` | `H|B = v Q` |
  | (B1) | a backward arc `b → a` with `a → v` and `a ∈ Init(A − v)` | `Q P`, with `Q` a path of `B` ending at `b` and `P` a path of `A'` starting at `a` | `H|A = v P` |
  | (B2) | a backward arc `b → a` with `v → b` and `b ∈ Term(B − v)` | `Q P`, with `Q` a path of `B'` ending at `b` and `P` a path of `A` starting at `a` | `H|B = Q v` |

* **(T4) Restriction.** Suppose `A_1` is the first block of `A` and `v ∉ A_1`. Then `T = A_1 ∪_w T_2` with `T_2 = A_2 ∪_v B`, and by (T1), (T2):
  `H(A)H(B) = H(T) ≥ H(A_1)H(T_2) ≥ H(A_1)H(A_2)H(B) = H(A)H(B)`.
  So `T_2` is exact. The same holds for the last block of `A` and for the outer blocks of `B`. Iterating, `T_loc = A_loc ∪_v B_loc` is exact, where `A_loc` and `B_loc` consist of the blocks containing `v`: one block, or two when `v` is a block boundary.

*Step 1: junctions.*

* **The split.** Suppose `A_loc = A_left ∪_v A_right`, where `v` is the end of the block `A_left` and the start of the block `A_right`.
* **Two exact sub-gluings.** Grouping two ways with (T1), (T2):
  * `H(T_loc) ≥ H(A_left) H(A_right ∪_v B_loc) ≥ H(A_loc) H(B_loc)`;
  * `H(T_loc) ≥ H(A_right) H(A_left ∪_v B_loc) ≥ H(A_loc) H(B_loc)`.

  So `A_right ∪_v B_loc` and `A_left ∪_v B_loc` are exact and smaller.
* **What induction gives.** Both are concatenations. Writing `B'` for `B_loc'`:
  * if `A_right ≠ C_3`, then `v` is the end of `B_loc` and `B' → A_right'`;
  * if `A_left ≠ C_3`, then `v` is the start of `B_loc` and `A_left' → B'`.

  A `C_3` factor also allows the other order.
* **The cases.** Every combination is excluded:
  * **Both one-way.** If every cross arc of `T_loc` is forward, or every one backward, then `v`, a junction of `A_loc`, is excluded by statement 2.
  * **The configuration `F`.** Otherwise `v` is both the start and the end of `B_loc`, so `B_loc = C_3`. The arcs are `A_left' → B' → A_right'`, or, when both blocks are `C_3`, `F'`: `A_right' → B' → A_left'`.
    * Take a back-arc path `L_X v L_Y` of `A_left`, whose `L_Y` starts at `x_0`, the out-neighbour of `v`.
    * Take a path `R_X v R_Y` of `A_right` whose `R_X` ends at the in-neighbour of `v`.
    * Then `L_X b_1 R_X v L_Y b_2 R_Y` is a Hamiltonian path of `T_loc` whose restriction to `B` is `b_1 v b_2`, where `v → b_1`. So (R) fails.
  * **The configuration `F'`.** The path `b_1 x_0 x_1 v z_1 z_2 b_2` does the same.
  * **`v` a junction of `B_loc`.** Then it is neither the start nor the end of `B_loc`, which contradicts what induction gives.
  * **Symmetry.** Junctions of `B` are excluded by the same argument.

*Step 2: interior positions.* Now `A_loc = NT_s` and `B_loc = NT_r` are single blocks.

* **Both interior.** Suppose `v` is interior to both, so `s, r ≥ 4` and `A − v`, `B − v` are strong.
  * (F1) and (B1) make every in-neighbour of `v` in `A` beat all of `B'`.
  * (F2) and (B2) make every in-neighbour of `v` in `B` beat all of `A'`.
  * This is a contradiction.
* **Interior in `A`, an end of `B`.** (A 3-cycle counts as an end.) Let `v = x_t`, `0 < t < N`.
  * (F1) and (B1) force `a → b` iff `a ∈ {x_0, …, x_(t−1)}`.
  * Take the path `x_1 ⋯ x_t x_N x_0 x_(t+1) ⋯ x_(N−1)` of `NT_s`, and a path `Q_1 v b` of `B`.
  * The vertex `b` can follow `v` directly (since `b → x_N`), or follow `x_0` (since `x_0 → b` and `b → x_(t+1)`, or `b` comes last). That gives two merges, so (U) fails.
* **The mirror case**, an end of `A` with an interior position of `B`, follows by symmetry.

*Step 3: end positions of single blocks (FINITE-EXACT).* `exact_gluings.py` checks every ordered pair of blocks among `NT_3`, `NT_4`, `NT_5`, every gluing position and every cross orientation. Exactness occurs only at `(x_N, y_0)` with all arcs forward, or at `(x_0, y_M)` with all arcs backward.

*Step 4: globalization.* Take `v = x_N = y_0` with the local arcs forward; the backward case is the converse picture.
* **The A side.** If `A` has `k ≥ 2` blocks, then:
  * `Init(A − v) = {x_0, …, x_(l_k)}`;
  * `N^-_A(v) = A' − {x_(l_k)}`;
  * so (B1) makes every arc between `{x_0, …, x_(l_k − 1)}` and `B'` forward.
* **The B side.** If `B` has `k' ≥ 2` blocks, then:
  * `Term(B − v) = {y_(m_2), …, y_M}`, where `y_(m_2)` is the end of `B`'s first block;
  * `N^+_B(v) = B' − {y_(m_2)}`;
  * so (B2) makes every arc between `A'` and `{y_(m_2 + 1), …, y_M}` forward.
* **The rest.** The remaining cross arcs join the last block of `A` to the first block of `B`. These are the local arcs, forward by Step 3.

So every cross arc is forward, and `T = IC(σ τ)` by statement 2. ∎

**6.** Optimal chains have parts in `{2, 3, 4}`, so statement 5 supplies the hypothesis of statement 4. The count and the palindromes follow from THM-4744 (2). ∎

## Audit trail

* The S50 audit (`audit_S50/`) confirmed statements 1–4 with independent code.
  * It checked every bad-split construction at all 745 positions for `N = 3..10`.
  * It reproduced the gluing census, and the failure modes: four of the fourteen mixed patterns of `C_3 ∪ C_3` violate only (U), with least surpluses 4 and 6.
* It also found the prior art. The product law is Moon 1972, and the characterization is Moon–Yang 2007, so the novelty claims were retracted (MISTAKE-652).
* Statements 5 and 6 were written after that audit. `oneway_proof_checks.py` confirms the excluded configurations by exact counts:
  * `F(s,r)` for `s, r = 3..6`;
  * `F'(3,3)`;
  * every interior configuration with `s ≤ 6`, `r ≤ 5`;
  * the witness paths.
