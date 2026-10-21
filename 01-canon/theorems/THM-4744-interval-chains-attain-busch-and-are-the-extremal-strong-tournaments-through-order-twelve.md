---
id: THM-4744
title: "Interval chains and the extremal strong tournaments of Busch's theorem. (1) For every composition N = p_1 + ... + p_k with parts p_i >= 2, the interval chain IC(p) - vertices x_0 < ... < x_N, every arc forward except one back arc x_(l_(i+1)) -> x_(l_i) per block, with l_1 = 0 and l_(i+1) = l_i + p_i - has exactly prod (2^(p_i - 1) + 1) Hamiltonian paths. (2) For N >= 3, IC(p) is rigid, distinct compositions give non-isomorphic chains, and the converse of IC(p) is IC(reversed p), so IC(p) is self-converse iff p is a palindrome; all chains of order m = N + 1 share the score sequence (1, 1, 2, ..., m-2, m-2) and have exactly m - 2 cyclic triples. (3) FINITE-EXACT, complete in the count: every strong tournament with fewer than 375 Hamiltonian paths has order at most 11, and the strong raw-count polynomials P_d are determined for every odd d < 375; P_d = 0 exactly for d = 7, 21. (4) FINITE-EXACT: for every m = 3..12 the strong m-tournaments with exactly Busch's minimum f(m) paths are exactly the optimal chains (parts in {2, 3, 4}, product f(m)), so E(m) = 1, 1, 2, 2, 1, 5, 3, 1, 9, 4"
status: >
  PROVED: statements 1 and 2, all orders (proofs below). FINITE-EXACT, complete for the stated ranges by Moon's
  strong-deletion theorem, the insertion lemma (both in THM-1370) and Busch's minimum (CITED): statements 3 and 4.
  The inputs are the exhaustive n <= 9 class census (h_multiplicities.py) and exact one-vertex ears to orders 10, 11, 12
  (strong_ears.py 225 and 375). E(12) = 4 was committed as a prediction in fad6a40581 (HYP-9356) before it was
  computed. Statement 3 extends the concurrent Session24 inventory, complete for H <= 73
  (actual_hamilton_multiplicities_20261021), and the class series of spectrum_basins_20261020e then gives the exact
  multiplicity polynomials N(n, h) for every odd h < 375. The all-order version of statement 4 is HYP-9356.
  Session opus-2026-10-21-S50; pre-audit.
session: opus-2026-10-21-S50
source: 05-knowledge/results/h_spectrum_multiplicities_20261021.md
scripts:
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears.py (+ strong_ears_225.out, strong_ears_375.out)
  - 04-computation/experiments/h_spectrum_multiplicities_20261021/h_multiplicities.py (+ .out)
related:
  - THM-1370 (the omissions 7, 21; Moon's monotonicity; the insertion lemma)
  - HYP-9356 (all-order characterization of Busch's extremal tournaments)
  - 05-knowledge/results/spectrum_basins_20261020e/spectrum.md (Busch, CITED; the class series F_N)
  - 05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md (concurrent)
---

# THM-4744 — interval chains and the extremal strong tournaments

## Setting

* `H(T)` is the number of directed Hamiltonian paths.
* `f(m)` is Busch's minimum over strong `m`-tournaments: `f(3q) = 3·5^(q−1)`, `f(3q+1) = 5^q`, `f(3q+2) = 9·5^(q−1)`.
  * A. H. Busch, *Electron. J. Combin.* 13 (2006) N3, Corollary 1 (CITED).
  * The paper proves the minimum and that Moon's 1972 construction attains it. It does not describe or count the tournaments that attain it.
* `NT_s` is the nearly transitive tournament on `s` vertices: the transitive order with its first and last vertex arc reversed. It has `2^(s−2) + 1` paths (Moon; Busch's note, p. 2): the unique path avoiding the back arc, and one path for each split of the `s − 2` interior vertices.
* **Interval chain** `IC(p_1, …, p_k)`:
  * vertices `x_0, …, x_N` with `N = Σ p_i`;
  * block `i` is the interval `[l_i, l_(i+1)]`, with `l_1 = 0` and `l_(i+1) = l_i + p_i`;
  * every arc goes forward except the `k` back arcs `x_(l_(i+1)) → x_(l_i)`.

  Each block alone is `NT_(p_i + 1)`, and consecutive blocks share an end vertex.

## Statements

1. **Product law (PROVED).** `H(IC(p)) = ∏_i (2^(p_i − 1) + 1)`.
2. **Recognition (PROVED, `N ≥ 3`).**
   * The out-degree of `x_t` is `N − 1` for `t = 0` and `N − t` for `1 ≤ t ≤ N − 1`, and `x_N` has out-degree 1.
   * So the score sequence is `(1, 1, 2, …, N−1, N−1)` for every composition, and there are `m − 2` cyclic triples.
   * The spine order is determined by the tournament. Hence:
     * `IC(p)` is rigid;
     * `IC(p) ≅ IC(p')` iff `p = p'`;
     * the converse of `IC(p)` is `IC(p reversed)`.
3. **Complete strong inventories below 375 (FINITE-EXACT).**
   * The strong classes with `H < 375` have order at most 11.
   * The ears give 1739 strong classes of order 10 and 411 of order 11 with `H ≤ 375`.
   * For every odd `d < 375` the polynomial `P_d(z) = Σ z^|T|` over strong classes with `H(T) = d` is determined, and `P_d = 0` exactly for `d ∈ {7, 21}`. In particular, every odd `d ∈ [3, 373]` other than 7 and 21 is the count of a **strong** tournament.
   * Converse pairing holds on these inventories: `P_d ≡ SC_d (mod 2)` coefficientwise at orders 10 and 11, where `SC_d` counts the self-converse strong classes. Every self-converse strong 10-class with `H ≤ 375` has `H ≡ 1 (mod 4)`, as THM-791 requires of every self-converse tournament of even order.
4. **The extremal classes through order 12 (FINITE-EXACT).** For `m = 3, …, 12` the strong `m`-tournaments with `H = f(m)` are exactly the chains `IC(p)` whose parts lie in `{2, 3, 4}` and whose product is `f(m)`.

| `m` | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| `f(m)` | 3 | 5 | 9 | 15 | 25 | 45 | 75 | 125 | 225 | 375 |
| optimal compositions | (2) | (3) | (2,2), (4) | (2,3), (3,2) | (3,3) | 2,2,3 in 3 orders; (4,3), (3,4) | 2,3,3 in 3 orders | (3,3,3) | 2,2,3,3 in 6 orders; 4,3,3 in 3 orders | 2,3,3,3 in 4 orders |
| `E(m)` | 1 | 1 | 2 | 2 | 1 | 5 | 3 | 1 | 9 | 4 |
| self-converse | 1 | 1 | 2 | 0 | 1 | 1 | 1 | 1 | 3 | 0 |

## Proofs

**1.** Induction on `k`. For `k = 1` the chain is `NT_(p_1+1)`. For `k ≥ 2` put:
* `A = {x_0, …, x_(p_1)}`, the first block;
* `v = x_(p_1)`, the shared vertex;
* `B = {x_(p_1), …, x_N}`, the chain `IC(p_2, …, p_k)`;
* `A' = A − v` and `B' = B − v`.

Every arc between `A'` and `B'` goes from `A'` to `B'`, since back arcs stay inside blocks.

* **Shape of a path.** In a Hamiltonian path `X v Y`, no vertex of `B'` is followed by a vertex of `A'`. So `X` and `Y` each have the form (vertices of `A'`)(vertices of `B'`).
* **Merges are unique.** For Hamiltonian paths `P_1 v P_2` of `A` and `Q_1 v Q_2` of `B`, the only order-preserving merge forming a path is `P_1 Q_1 v P_2 Q_2`. It is a path because every cross arc is available, `last(Q_1) → v` and `v → first(P_2)`.
* **Every path restricts.** Conversely, every Hamiltonian path `X_a X_b v Y_a Y_b` restricts to the paths `X_a v Y_a` of `A` and `X_b v Y_b` of `B`, unless one of two shapes occurs. Both are impossible:
  * **(i)** `X_a` and `X_b` are nonempty and `v → last(X_a)`. The only out-neighbour of `v` in `A` is `x_0`, so `last(X_a) = x_0`. The only in-neighbour of `x_0` in `T` is `v`, so `x_0` is the first vertex of the path and `X_a = (x_0)`. Then `Y_a = A' − x_0` is nonempty, since `p_1 ≥ 2`. But `Y_a` follows `v`, and `x_0` is the only out-neighbour of `v` in `A'`. Contradiction.
  * **(ii)** `Y_a` and `Y_b` are nonempty and `first(Y_b) → v`. The only in-neighbour of `v` in `B` is `x_(l_3)`, the end of block 2, so `first(Y_b) = x_(l_3)`. If `X_b` were nonempty, its last vertex would precede `v` and be another in-neighbour of `v`. So `X_b` is empty, and the interior `M` of block 2 (nonempty, since `p_2 ≥ 2`) lies in `Y_b`, after `x_(l_3)`. A vertex of `M` can be entered only from `v`, from `A'`, or from interior vertices of block 2 below it. So the smallest-index vertex of `M` in `Y_b` has no possible predecessor. Contradiction.

Thus `H(T) = H(A)·H(B)`, and induction gives the product. ∎

**2.**
* **Scores.** A forward arc leaves `x_t` toward each of the `N − t` later vertices. A block start loses one out-arc (its arc to the block end is reversed), and a block end gains one (its back arc). A shared vertex is both, so it nets zero. This gives the stated scores.
* **The spine is recognizable.** For `N ≥ 3` the scores `N − t`, `2 ≤ t ≤ N − 2`, occur once each. The two vertices of score `N − 1` are `x_0, x_1`, with `x_0 → x_1`: the back arc of block 1 ends at `x_(p_1)` and `p_1 ≥ 2`. The two vertices of score 1 are `x_(N−1), x_N`, with `x_(N−1) → x_N`.
* **Consequences.**
  * Every isomorphism preserves the spine order, so `IC(p)` is rigid.
  * The arcs that point backward in the spine order are exactly the back arcs, and these spell out `p`. So different compositions give non-isomorphic chains.
  * Reversing every arc reverses the spine order and maps the back arcs to back arcs, so the converse of `IC(p)` is `IC(p reversed)`.
* **Triangles.** The cyclic-triple count `C(m,3) − Σ C(s_v, 2)` depends only on the scores. It equals `m − 2` (checked for `m ≤ 12`, and an elementary sum in general). This is Moon's minimum for strong tournaments. ∎

**3, 4.**
* **Completeness.** Let `T` be strong of order `m ≥ 4`. By Moon (1966), `T` has a vertex `v` with `T − v` strong, and `H(T − v) ≤ H(T)` by the insertion lemma. So every strong class with `H ≤ X` arises as a nonconstant one-vertex ear of a strong class of one order lower with `H ≤ X`.
* **Computation.**
  * Starting from the census, the ears of the 3,492 strong 9-classes with `H ≤ 375` give every strong 10-class with `H ≤ 375`.
  * Their ears give every strong 11-class with `H ≤ 375`.
  * Their ears give every strong 12-class with `H = 375`, and show that no strong 12-tournament has `H < 375`.
  * Canonical forms: the S36 refined-invariant code for `n ≤ 11`, and an exact individualisation-refinement code at `n = 12`.
* **Order bound.** `f(12) = 375` and `f` is nondecreasing (THM-1370), so no strong tournament of order 12 or more has fewer than 375 paths. This gives statement 3.
* **Identification.** The extremal classes were compared with the canonical forms of the optimal chains. ∎

## What this does and does not say

* Statement 4 is proved for `m ≤ 12` only. The all-order characterization is HYP-9356.
  * Statement 1 gives one direction for all `m`: every optimal chain is extremal.
  * Statement 2 shows that distinct chains are distinct classes.
  * So `E(m) ≥ #optimal compositions` holds for all `m`: `E(3q) ≥ q`, `E(3q+1) ≥ 1`, `E(3q+2) ≥ q(q+3)/2`.
* Statement 1 holds for any parts `p_i ≥ 2`, but only parts 2, 3, 4 can be optimal. Per added vertex, block `NT_(p+1)` costs `log(2^(p−1)+1)/p`, which is 0.549, 0.536, 0.549, 0.567, … for `p = 2, 3, 4, 5`. So `p = 3` is cheapest, one 4 costs the same as two 2s, and nothing else is optimal.
