---
id: THM-4652
title: "Paley pattern character sums: reversal signs, odd and isolated-edge vanishing, odd-reversal symmetry (complete patterns K_m vanish for q = 3 mod 4), and series-parallel polynomiality"
status: PROVED (elementary character-sum identities); FINITE-EXACT checks in the scripts listed below
source: opus-2026-10-14-S34
depends_on: [orthogonality of the quadratic character; Jacobi sums of the quadratic character (sum_u chi(u(1-u)) = -chi(-1))]
verification: 04-computation/experiments/paley_rate_c1_ingredients_20261014.py, paley_A6_polynomiality_20261014.py, paley_K4_charsum_cm_20261014.py, paley_K5_charsum_20261014.py, paley_W4_charsum_20261014.py
note: 05-knowledge/results/paley_h31_exact_rate_and_apery_wheel_20261014.md
related: [THM-438 (Paley cluster integrals are Catalan; R(p) -> e), THM-582 (palindromic involution), HYP-9274, HYP-9273]
---

# THM-4652 — Paley pattern sums

## Setting

- `q` is an odd prime power with `q ≡ 3 (mod 4)`, and `χ` is the quadratic character of `F_q` (`χ(0) = 0`, `χ(−1) = −1`).
- A **pattern** is an oriented multigraph `G` on the vertex set `[m]`; each edge `e` has a tail `t(e)` and a head `h(e)`.
- The **pattern sum** is `Σ*(G) = Σ_v Π_e χ(v_(h(e)) − v_(t(e)))` over injective maps `v : [m] → F_q`. The unrestricted sum `Σ(G)` is the same sum over all maps.

These sums are the building blocks of every count of patterns in the Paley tournament `P_q` (arc `x → y` iff `χ(y − x) = 1`). Expanding `[x → y] = (1 + χ(y − x))/2` over the arcs of any pattern expresses the count through the `Σ*` of its sub-patterns. In particular the Hamiltonian-path ratio `R(q) = H(P_q) 2^(q−1)/q!` (THM-438) is a sum of normalised `Σ*` of unions of paths.

## Statement

1. **(Reversal.)** Reversing one edge multiplies `Σ*(G)` by `−1`. Up to sign, `Σ*` depends only on the undirected multigraph.
2. **(Odd size.)** If `G` has an odd number of edges, `Σ*(G) = 0`.
3. **(Isolated edges.)** If `G` has a connected component consisting of one edge, `Σ*(G) = 0`. Taking that edge with multiplicity `k` gives the same conclusion for odd `k`. For even `k` the component contributes a factor `χ² = 1`, i.e. it is a pair of distinct free points.
4. **(Odd-reversal symmetry.)** Suppose a permutation `π` of `[m]` carries the underlying multigraph of `G` to itself and reverses an odd number of edges. Then `Σ*(G) = 0`. In particular, all of the following give `Σ*(G) = 0` for every `m >= 2`:
   - the complete pattern `K_m` (every pair joined once, any orientation): its summand is `±χ(Vandermonde)`;
   - `K_5` minus one edge;
   - `K_5` minus two edges sharing a vertex.
5. **(Series-parallel polynomiality.)** Let `G` be a two-terminal series-parallel multigraph, built from single edges by series and parallel composition. Then `Σ(G) ∈ Z[q]`: a polynomial in `q` whose integer coefficients depend only on `G`, not on the arithmetic of `q` beyond `q ≡ 3 (mod 4)`. The same holds for every graph whose blocks are series-parallel.

## Proof

1. `χ(x − y) = χ(−1)χ(y − x) = −χ(y − x)`.
2. `v ↦ −v` is a bijection of injective maps, and it multiplies the summand by `(−1)^|E|`.
3. Fix the values of all other vertices, a set `D` of distinct values. Then
   `Σ_(a ≠ b; a, b ∉ D) χ(b − a) = Σ_(a ≠ b) χ(b − a) − Σ_(a or b ∈ D) χ(b − a) = 0 − (0 + 0 − Σ_(a ≠ b ∈ D) χ(b − a))`.
   The last sum vanishes because the pairs `(a, b)` and `(b, a)` cancel (`χ(−1) = −1`). For multiplicity `k`, `χ(x)^k = χ(x)` for odd `k` and `= 1` for even `k` when `x ≠ 0`.
4. `Σ*` is a sum over all injective maps, so it is invariant under relabelling the vertices. Relabelling by `π` gives the same multigraph with `r` edges reversed, so `Σ*(G) = (−1)^r Σ*(G)` by part 1.
   - **`K_m`.** The transposition of two consecutive labels in the low-to-high orientation reverses exactly one edge, the edge between them. For any other orientation, use part 1 first.
   - **`K_5 − e`, `e = {i, j}`.** Take the transposition of two labels outside `{i, j}`. The reversal count is odd: one edge between the swapped labels, plus two for each label strictly between them. The Vandermonde formulation below is the clean argument.
   - **`K_5` minus `{i, j}, {i, k}`.** Take the transposition `(j k)`. It fixes the missing pair `{i, j}, {i, k}` as a set and reverses the edge `jk`.

   In each case the sign bookkeeping reduces to: the summand is `±χ(Vandermonde)·Π_(missing edges) χ(·)`, and the transposition negates the Vandermonde while fixing the product over the missing edges.
5. For two terminals `s, t`, put `T_G(d) = Σ` over the internal vertices (unrestricted) with `v_t − v_s = d`. The class
   `𝓕 = {α χ(d) + γ + β [d = 0] : α, β, γ ∈ Z[q]}`
   contains a single edge (`T = χ`), and it is closed under:
   - **series composition, which is convolution:**
     - `χ * χ = 1 − q[d = 0]`, by the Jacobi sum `Σ_u χ(u(1 − u)) = −χ(−1) = 1` for `d ≠ 0`, and `Σ_x χ(x)χ(−x) = −(q − 1)` at `d = 0`;
     - `χ * 1 = 0`, `1 * 1 = q`, `δ * f = f`;
   - **parallel composition, which is pointwise product:**
     - `χ·χ = 1 − [d = 0]`;
     - `χ·δ = 0`, `1·f = f`.

   So `T_G ∈ 𝓕`, and `Σ(G) = q·Σ_d T_G(d) ∈ Z[q]`. A cut vertex multiplies the block sums (divided by `q`), which gives the final sentence. ∎

## Checks (FINITE-EXACT)

- Part 2: `A_3 = A_5 = 0` for runs of three and five steps, at every prime tested.
- Part 3: the cross terms of `paley_rate_c1_ingredients_20261014.py`.
- Part 4: `S(q) = Σ_(x,y) χ(x(x−1)y(y−1)(y−x)) = 0` (`K_4`) for all 24 primes `q ≡ 3 (mod 4)` below 200. The same `K_4` sum equals `4a² − 2q` (`q = a² + b²`, `a` odd) for all 21 primes `q ≡ 1 (mod 4)` below 200, where `χ(−1) = +1` and part 4 does not apply. `T(q) = 0` (`K_5`) for all 19 primes `q ≡ 3 (mod 4)` below 160.
- Part 5 is consistent with the exact polynomial run sums of HYP-9274: `A_2 = q(q−1)`, `A_4 = 2q(q−1)(q−3)`, `A_6 = q(q−1)(q−3)(5q−29)`. (Those are injective sums, so part 5 does not prove them by itself; Möbius inversion over set partitions can create non-series-parallel quotients.)

## Reading

- The first 5-vertex pattern that survives parts 2–4 and is not series-parallel is `K_5` minus two disjoint edges, the wheel `W_4`. There a genuine Hasse–Weil quantity appears.
- HYP-9274 (B) holds for every odd prime (PROVED + CITED: Ahlgren–Ono Lemma 2.2 and Theorem 6, plus the collision correction; codex's wheel packet and the S34 audit). The statement is `Σ*(W_4) = −q(q−1)(a_q + 4)` for `f = η(2z)⁴η(4z)⁴ = Σ a_n q^n`, the weight-4 level-8 newform of Beukers' Apéry-number supercongruence. Deligne's theorem (the Riemann hypothesis for varieties) bounds it by `|a_q| <= 2q^(3/2)`.
- The same core enters single runs from 10 steps on: `A_10 = poly + 34q(q−1)a_q`.
- Parts 4–5 explain why every pattern sum computed below `W_4` is a polynomial for Paley tournaments.
