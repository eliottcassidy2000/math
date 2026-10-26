---
id: THM-4780
title: "The Paley census to six vertices. Write T_G for the injective character sum of an oriented graph G over F_q, q = 3 (mod 4). After the degree-at-most-two elimination calculus of the 2026-10-14 run packets and THM-4652's vanishing laws, the only dense graph sums (minimum degree >= 3) that survive in the transitive k-census of the Paley tournament P_q, k <= 6, are those without an odd-reversal automorphism: the wheel W4 = K5 - 2K2 (five vertices), and the octahedron K_(2,2,2) and the Hamiltonian-path complement K6 - P6 (six vertices). The octahedron is (q - 5) wheels by PGL(2) invariance, which proves HYP-9275's octahedron identity at every odd prime. With signed multiplicities +1 (wheel) and -8 = (1/2) sum_(S6) sgn(s)(-1)^des(s) (K6 - P6), for every prime power q = 3 (mod 4): TT5(P_q) = q(q-1)[(q-3)(q^2-21q+105) + S_W4(q)]/2^10 and TT6(P_q) = q(q-1)[(q-3)(q^3-46q^2+682q-3239) + 7(q-15)S_W4(q) - 8 S_(K6-P6)(q)]/2^15, where S = T/(q(q-1)). With the wheel identity (HYP-9274 (B)), TT5(P_p) = p(p-1)(p^3-24p^2+168p-319-a_p)/2^10 for primes p = 3 (mod 4), a_p the coefficient of eta(2z)^4 eta(4z)^4. The four-vertex level: a regular tournament has TT4 >= n(n-1)(n-3)(n-7)/64, with equality iff it is doubly regular, which for P_q is the Gauss identity g^2 = -q."
status: >
  PROVED. Parts 1-4 and 7 are elementary. Part 5 is the reduction calculus with a priori degree bounds, plus exact
  interpolation from finitely many computed values: the linear systems have full column rank, so their unique solution
  is the true one. Part 6 uses HYP-9274 (B), which is PROVED + CITED (Ahlgren-Ono) for odd primes. FINITE-EXACT
  consistency checks:
  TT5 at 127 values of q (3, 7, 11 and the 119 primes 19..1499);
  TT6 at 75 values of q (3, 7, 11 and the 67 primes 19..751). Both counts include the prime powers 27, 243, 343, 1331, 2187,
  computed by direct census;
  the octahedron identity at all 27 odd primes 5..109, in both classes mod 4.
  Audited 2026-10-25 (04-computation/experiments/pythagorean_hourglass_20261025/audit/REPORT.md): CORRECT WITH FIXES
  (attributions and wording only; MISTAKE-670). Part 5 re-derived by exact symbolic elimination, with no interpolation; TT5, TT6 and the
  graph sums recomputed by different algorithms; part 4 tested exhaustively on P^1(F_7) and on further even-degree graphs.
session: mac-mini-2026-10-25-pyth
source: 05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md
depends_on:
  - THM-4652   # reversal signs, odd and isolated-edge vanishing, odd-reversal symmetry
  # degree <= 2 elimination calculus: 05-knowledge/results/wheel_run_interface_20261014h.md (codex-tiling) and
  # 05-knowledge/results/paley_next_cores_20261014i.md (tiling-modular-atoms)
  - THM-4653   # K6 - P6 terminal = Legendre bilinear form M_q (used in the M_q form of part 5)
  - HYP-9274   # (B) wheel W4 = Apery form at odd primes: PROVED + CITED (Ahlgren-Ono 2000); used in part 6 only
scripts:
  - 04-computation/experiments/pythagorean_hourglass_20261025/census_interpolation.py (+ .out): part 5 by exact interpolation
  - 04-computation/experiments/pythagorean_hourglass_20261025/census_level_forecast.py (+ .out): part 3 for k = 4..7, signed multiplicities
  - 04-computation/experiments/pythagorean_hourglass_20261025/graph_sums6.py (+ .out): the nine dense six-vertex sums at 23 primes
  - 04-computation/experiments/pythagorean_hourglass_20261025/octahedron_pgl2_check.py (+ .out): part 4 and the HYP-9275 normalisation
  - 04-computation/experiments/pythagorean_hourglass_20261025/paley_prime_powers.py (+ .out): census at q = 3, 7, 11, 27, 243, 343, 1331, 2187
  - 04-computation/experiments/pythagorean_hourglass_20261025/paley_tt5_arith.py (-> paley_tt5_data.out), paley_tt6.py (-> paley_tt6_data.out)
  - 04-computation/experiments/pythagorean_hourglass_20261025/paley_graph_sums.py (+ .out): K4, K5, K5 - e, W4 sums at all odd primes <= 419
  - 04-computation/experiments/pythagorean_hourglass_20261025/qr_four_census.py (+ .out): part 7
related:
  - 05-knowledge/hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md
  - 05-knowledge/hypotheses/HYP-9275-paley-runs-one-modular-form-octahedron-k6p6.md
  - 05-knowledge/hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md   # S_(K6-P6) = |b(q)|^2 - q^2 - 7
  - 05-knowledge/hypotheses/HYP-9371-four-loop-law-doubly-regular-tournaments-see-the-five-vertex-census.md
  - 01-canon/theorems/THM-1820-intransitivity-is-schur-convexity-not-covariant-vanishing.md   # the Schur-convexity argument one level down
---

# THM-4780 — The Paley census to six vertices

## Setting

- `q ≡ 3 (mod 4)` is an odd prime power and `χ` is the quadratic character of `F_q`, so `χ(−1) = −1`.
- `P_q` is the Paley tournament: `x → y` iff `χ(y − x) = 1`.
- `TT_k(T)` is the number of transitive `k`-subsets of a tournament `T`.
- For an oriented graph `G` on `[k]`, `T_G(q) = Σ ∏_{(i→j)∈G} χ(x_j − x_i)`, summed over **injective** `x : [k] → F_q`.
- When `|E(G)|` is even, `T_G` is invariant under `AGL(1, q)`, so `T_G = q(q−1)·S_G` with `S_G` the sum over `x_1 = 0`, `x_2 = 1`.
- Unless stated otherwise, graphs carry the **order orientation** (`i → j` for `i < j`).
- The **dense** graphs are those with minimum degree `≥ 3`.
- `W4 = K5 − 2K2` is the wheel (a hub joined to a 4-cycle).
- `K_{2,2,2} = K6 − 3K2` is the octahedron. `K6 − P6` is `K6` minus the edges of a Hamiltonian path.

## Statement

1. **Expansion.** `2^{C(k,2)} TT_k(P_q) = Σ_{G ⊆ E(K_k)} T_G(q)`, with the order orientation.
2. **Reduction (the degree-≤2 elimination calculus of `wheel_run_interface_20261014h.md`, with degree bounds).**
   - Every `T_G` is a combination `Σ_H c_{G,H}(q) T_H(q)`.
   - For `G` not dense, `H` runs over dense oriented graphs on fewer vertices, together with the empty graph on no vertices (`T = 1`); a dense `G` is its own reduction.
   - Each `c_{G,H}` is in `Z[q]`, of degree `≤ |V(G)| − |V(H)|`, and the same polynomial for every `q ≡ 3 (mod 4)`.
3. **Survivors.** `T_G = 0` if `|E(G)|` is odd, or if the underlying graph has an automorphism that reverses an odd number of arcs. Among dense graphs with an even number of edges:

   | `k` | dense graphs | even `\|E\|` | killed by an odd automorphism | survivors |
   |---|---|---|---|---|
   | 4 | 1 | 1 | 1 (`K4`) | none |
   | 5 | 3 | 2 | 1 (`K5`) | `W4` |
   | 6 | 19 | 9 | 7 | `K6 − P6`, `K_{2,2,2}` |
   | 7 | 150 | 75 | 37 | 38 (one with all degrees even) |
4. **PGL(2) lemma.**
   - Suppose `G` is connected, every degree is even and `|E|` is even.
   - Then `T_G = ε·(q + 1 − k)·T_{G−u}` for every vertex `u`, with a sign `ε = ±1` fixed by the orientations.
   - In particular `T_{K_{2,2,2}} = (q − 5)·T_{W4}`. Here the octahedron has antipodal pairs `{1,2}, {3,4}, {5,6}` and order orientation; `W4` is labelled as in `paley_graph_sums.py`, and `ε = +1`.
   - The same identity holds for `q ≡ 1 (mod 4)`.
   - Consequence: HYP-9275's octahedron identity `U(K_{2,2,2}) = q(q−1)((5 − q)a_q + 3q)` holds at every odd prime. Here `U` is the all-maps sum and `a_q` is the coefficient of `η(2z)⁴η(4z)⁴`.
5. **Census formulas.** For every prime power `q ≡ 3 (mod 4)`:
   ```
   TT5(P_q) = q(q−1) [ (q−3)(q² − 21q + 105) + S_W4(q) ] / 2^10
   TT6(P_q) = q(q−1) [ (q−3)(q³ − 46q² + 682q − 3239) + 7(q − 15) S_W4(q) − 8 S_(K6−P6)(q) ] / 2^15
   ```
   - The coefficients `+1` and `−8` are the signed multiplicities of `W4` in `K5` and of `K6 − P6` in `K6`.
   - Explicitly, `−8 = ½ Σ_{σ∈S_6} sgn(σ)(−1)^{des σ} = ½·2³A_3(−1)`, where `A_3(t) = 1 + 4t + t²` is the Eulerian polynomial.
   - By THM-4653, `S_(K6−P6) = M_q + 2q − 6`, so equivalently `TT6 = q(q−1)[(q−3)(q³ − 46q² + 682q − 3255) + 7(q−15)S_W4 − 8M_q]/2^15`.
6. **Modular form at five vertices.** For primes `p ≡ 3 (mod 4)`, `S_W4(p) = −4 − a_p` (HYP-9274 (B)). Hence
   `TT5(P_p) = p(p−1)(p³ − 24p² + 168p − 319 − a_p)/2^10`, with `a_p` the coefficient of `η(2z)⁴η(4z)⁴`.
7. **Four vertices.**
   - Let `T` be a regular tournament on `n` vertices. Then `TT4(T) ≥ n(n−1)(n−3)(n−7)/64`, with equality iff `T` is doubly regular.
   - `TT4(T)` determines all four 4-vertex class counts. With `k = (n−1)/2`: `#(source + C3) = #(sink + C3) = n·C(k,3) − TT4`, and `#(strong) = C(n,4) − TT4 − 2(n·C(k,3) − TT4)`.
   - For `P_q`: equality ⟺ `S² = J − qI` (`S = A − Aᵀ`) ⟺ `g² = −q` for the quadratic Gauss sum.
   - So the Paley value `TT4(P_q) = q(q−1)(q−3)(q−7)/64` (McCarthy–Springfield, Cor. 2.3) is the Gauss identity.

## Proof

**Part 1.**
- A transitive `k`-set has exactly one ordering in which every pair is an arc `x_i → x_j` with `i < j`.
- So `TT_k = Σ_x ∏_{i<j} (1 + χ(x_j − x_i))/2` over injective `x`.
- Expand the product.

**Part 3, vanishing (THM-4652).**
- `x ↦ −x` multiplies every factor by `χ(−1) = −1`. So `T_G = (−1)^{|E|} T_G`.
- If `σ` is an automorphism of the underlying graph, relabelling the variables gives `T_G = (−1)^{rev(σ)} T_G`, where `rev(σ)` is the number of arcs that `σ` reverses.
- Reorienting a set `F` of arcs changes `rev(σ)` by `|F| + |σF| ≡ 0` (mod 2). So the parity, and hence the vanishing, does not depend on the orientation.
- The table is an exhaustive computation over `geng -d3 k` (`census_level_forecast.py`).
- For `k = 6` it agrees with a direct evaluation: the seven killed sums are `0` at all 23 primes `7..199` (`graph_sums6.py`).

**Part 2, reduction.** Let `v` have degree `d ≤ 2`, and let `X` be the set of the other `k − 1` values.
- `d = 0`: the sum over `x_v ∉ X` gives `(q − k + 1)·T_{G−v}`.
- `d = 1`, neighbour `u`: `Σ_{t∉X} χ(t − x_u) = −Σ_{z≠u,v} χ(x_z − x_u)`. Each term is `T` of `G − v` with the edge `{z,u}` toggled; a doubled edge gives `χ² = 1` on injective maps, or `χ(−1) = −1` when the two factors are oppositely oriented. The coefficient is `±1`.
- `d = 2`, neighbours `u, w`: by the Jacobi sum `Σ_t χ(t − a)χ(t − b) = −1` (`a ≠ b`), `Σ_{t∉X} χ(t − x_u)χ(t − x_w) = −1 − Σ_{z≠u,w,v} χ(x_z − x_u)χ(x_z − x_w)`. This gives `T_{G−v}` and the toggled graphs, with coefficients `±1`.
- The signs depend only on the orientations and on `χ(−1)`.
- Each step removes one vertex and multiplies by a polynomial of degree `≤ 1`. Iterating gives the degree bound.
- Dense graphs on `k` vertices are never produced from graphs on `k` vertices.

**Part 4, PGL(2).**
- Extend `T_G` to injective maps `x : V → P¹(F_q)`, dropping the factors of arcs at the vertex placed at `∞`. Call the summand `F(x)`.
- Choose an orientation with every in-degree even. A connected graph with `|E|` even has one, by the parity-orientation lemma; reorienting changes `T_G` by a sign.
- Affine maps `x ↦ ax + b` multiply `F` by `χ(a)^{|E'|}`, where `E'` is the set of arcs not at `∞`. Then `|E'| = |E|` or `|E| − deg(u)`, which is even in both cases.
- For the inversion `ι(x) = −1/x`, the identity `−1/x_j + 1/x_i = (x_j − x_i)/(x_i x_j)` gives `F(ιx) = F(x)·∏_v χ(x_v)^{deg v} = F(x)` when no vertex sits at `0` or `∞`.
- When a vertex `u` sits at `0`, the factors at `u` contribute `(−1)^{in(u)}`. When vertices `u` and `w` sit at `0` and `∞`, they contribute `(−1)^{deg'(w) + in'(w) + in'(u)}`, primes denoting degrees into the other `k − 2` vertices. Both exponents are even: if `u ~ w`, then `deg'(w)` and exactly one of `in'(u)`, `in'(w)` drop by one, so the second exponent has the parity of `deg(w) + in(w) + in(u)`. So `F` is `PGL(2)`-invariant.
- `PGL(2)` is transitive on `P¹`. So, for each vertex `u`, the configurations with `x_u = P` have the same total for every `P ∈ P¹`. At `P = ∞` that total is `T_{G−u}`.
- Hence the full `P¹` sum is `(q+1)·T_{G−u}` for every `u`. Subtracting the `k` classes with a vertex at `∞` gives `T_G = (q + 1 − k)·T_{G−u}`.
- The octahedron in order orientation has in-degrees `0, 0, 2, 2, 4, 4`, so `ε = +1`. Deleting `u = 6` leaves `W4` with hub `5`; relabelling the hub from the largest to the smallest label reverses 4 arcs, so the sign stays `+1`.
- For `q ≡ 1 (mod 4)` the in-degree condition is not needed.
- Check: `S_oct = (q − 5)·S_W4` at all 27 odd primes `5 ≤ q < 110` (13 in the class 1 mod 4, 14 in the class 3 mod 4).

**HYP-9275's normalisation.**
- `U(G)` sums over all maps, so non-adjacent vertices may collide. Decompose `U` by kernel partitions into independent blocks.
- For the octahedron, merged antipodal pairs carry same-direction double arcs. This gives `U(K_{2,2,2}) = T_oct + 3(q−4)T_{C4} + 3(q)_4 + (q)_3`.
- Here `T_{C4} = q(q−1)(1 − (q−2)) = −q(q−1)(q−3)` for the 4-cycle with two sources and two sinks.
- For the wheel the same decomposition gives `U(W4) = T_W4 − q(q−1)(q−4)`.
- HYP-9274 (B) states the injective identity `T_W4 = −q(q−1)(a_q + 4)`; with the line above it is `U(W4) = −q(q−1)(q + a_q)` (as recorded in HYP-9275 (A); `U(W4)` is `q(q−1)` times Ahlgren–Ono's moment). Together with `T_oct = (q−5)T_W4`, this yields exactly `U(K_{2,2,2}) = q(q−1)((5−q)a_q + 3q)`.
- `octahedron_pgl2_check.py` confirms every identity of this paragraph by brute force at `q = 5, 7, 11` (and the wheel ones at `q = 13`).

**Part 5.**
- By parts 1–3, `2^{10} TT5 = A5(q) + m_{W4}·T_{W4}`, with `deg A5 ≤ 5`.
- `m_{W4}` is the signed number of labelled copies of `W4` in `K5`: `m_H = |Aut H|^{-1} Σ_{σ∈S_k} (−1)^{rev_H(σ)}`. It equals `+1`.
- Similarly `2^{15} TT6 = A6(q) + B6(q)·T_{W4} + m_oct·T_oct + m_{K6−P6}·T_{K6−P6}`, with `deg A6 ≤ 6` and `deg B6 ≤ 1`.
- Here `m_oct = +1`, and `T_oct = (q−5)T_W4` is absorbed into `B6`.
- For `K6 − P6` with the path `1−2−…−6`, a relabelling `σ` reverses `inv(σ) − des(σ)` complement arcs. So `m = ½ Σ sgn(σ)(−1)^{des σ}`, and by the Désarménien–Foata identity `Σ_{σ∈S_n} sgn(σ) t^{des σ} = (1 − t)^{⌊n/2⌋} A_{⌈n/2⌉}(t)` (checked for `n ≤ 8`), this is `½·8·(−2) = −8`.
- The unknown coefficients (7 for `TT5`, 10 for `TT6`) were then solved exactly from computed values of `TT_k(P_q)`, `S_W4(q)` and `S_{K6−P6}(q)`: 41 primes for `TT5` and 24 for `TT6`. In the `TT6` run the multiplicity `D` was left free.
- Both systems have full column rank, so the solution is unique. The true coefficients satisfy the system, so the solution is the truth.
- The solved multiplicities come out as `+1` and `−8`, as computed independently.
- The polynomials are the same for all `q ≡ 3 (mod 4)`, so the formulas hold for prime powers too.
- The `M_q` form follows from THM-4653: `V_q = q(q−1)(−εM_q + 7 + 13ε − 2εq)` with `ε = −1`.

**Part 6.**
- HYP-9274 (B) states `S_W4 = −4 − a_p` at odd primes. Its proof: `U(W4) = q(q−1)·Σ_b χ(b)(Σ_c χ(c(c−1)(c−b)))²` exactly, Ahlgren–Ono evaluate the moment as `−a_p − p`, and the collision correction `U(W4) = T_W4 − q(q−1)(q−4)` converts this to the injective sum.
- Equivalently, `S_W4(q) = Σ_{x≠0,1} χ(x)a_x² + (q − 3)`, where `a_x = −Σ_y χ(y(y−1)(y−x))` is the Legendre-family trace.
  With the hub at `0`, one rim vertex at `1` and its non-neighbour at `x`, the two remaining rim vertices each contribute `χ(y(y−1)(y−x))`.
  Their distinct-pair sum is `a_x² − (q − 3)`, and `Σ_{x≠0,1} χ(x) = −1`. This is checked numerically in `census_interpolation.py`.
- Substitute into part 5.

**Part 7.**
- `TT4 = Σ_v t3(N⁺(v))` and `t3(Y) = Σ_{u∈Y} C(d_Y(u), 2)`.
- Inside `N⁺(v)`, a tournament on `k = (n−1)/2` vertices, the out-degrees sum to `C(k, 2)`. Convexity gives the bound, with equality iff every arc `v → u` has `(n−3)/4` common out-neighbours, i.e. `T` is doubly regular.
- The class counts follow from `Σ_v c3(N⁺(v)) = #(source + C3)`, `Σ_v c3(N⁻(v)) = #(sink + C3)` and `c3 + t3 = C(k,3)`.
- For `P_q`, `S` acts on additive characters with eigenvalues `χ(a)g` (`a ≠ 0`) and `0`. Doubly regular means `SSᵀ = qI − J`, i.e. `|g|² = q`, and since `S` is skew this is `S² = J − qI`, i.e. `g² = −q`.
- Checks: exhaustive over all regular tournaments at `n = 7, 11`; random regular tournaments at `n = 15, 19, 23`; `P_q` for `q ≤ 43` (`qr_four_census.out`).

## Data and consequences (FINITE-EXACT)

- **Zeros of the census.**
  - `TT5(P_q) = 0` at `q = 3, 7, 11`.
  - `TT6(P_q) = 0` at `q = 3, 7, 11, 19, 23, 27`. `P_27` is ST27, the unique largest `TT6`-free tournament (Sanchez-Flores; Neiman–Mackey–Heule).
  - In part 5 these zeros are arithmetic identities. For example, `TT5(P_7) = 0` ⟺ `S_W4(7) = −28` ⟺ `a_7 = 24`. Similarly `TT5(P_11) = 0` ⟺ `a_11 = −44`.
- **Prime powers.** With `a_q` read as the Frobenius trace `α^r + β^r`, the formula of part 6 holds at `q = 27, 243, 343, 1331, 2187` (`paley_prime_powers.out`). There the wheel identity itself is FINITE-EXACT, not cited; in this Frobenius form it was already FINITE-EXACT, by direct evaluation, at 9, 25, 27, 49, 81, 121, 125, 169 (`wheel_apery_extension_20261014h.md` §4).
- **Two-step reading.** Every surviving dense sum lacks an odd-reversal symmetry. The five-vertex survivor is the `₄F₃(1)` wheel. The six-vertex level adds exactly one new sum, the `₂F₁·₃F₂` terminal of THM-4653. This is the run-expansion ladder of the 2026-10-14 notes, seen through subtournament counts.

## Prior art and credits

- **TT3 and TT4.** `TT3(P_q) = q(q−1)(q−3)/8` and `TT4(P_q) = q(q−1)(q−3)(q−7)/64` are McCarthy–Springfield, arXiv:2311.02135, Cor. 2.7 and 2.3. They count orders 3 and 4 in `k`-th power Paley digraphs. Both values are older facts for the classes containing `P_q`: `TT3 = n(n−1)(n−3)/8` holds for every regular tournament, and the `TT4` value for every doubly regular one (part 7). We found no published order-5 or order-6 count.
- **Earlier repo work (2026-10-14).** The vanishing laws are THM-4652 (opus-2026-10-14-S34). The degree-≤2 injective elimination calculus is codex-tiling's `wheel_run_interface_20261014h.md`; the ten-edge terminal classification, including the survival of `K6 − P6` and the vanishing of the other ten-edge six-vertex terminals, is tiling-modular-atoms' `paley_next_cores_20261014i.md`. The wheel and octahedron are cores of the Paley run expansion: HYP-9274 (B), PROVED + CITED, and HYP-9275 (opus-2026-10-14-S35), where the octahedron identity was FINITE-EXACT and open. The `K6 − P6` reduction is THM-4653 and HYP-9276 (mac-mini-2026-10-14-wilson).
- **New here.**
  - The PGL(2) counting lemma and the resulting proof of HYP-9275's octahedron identity. The invariance behind it is classical: with an even in-degree orientation, `F(x) = χ(∏_{(i→j)} det(x_i, x_j))` in homogeneous coordinates, which is well defined on `P¹(F_q)` because every degree is even and `GL(2)`-invariant because `|E|` is even (the triangle case is the oriented Paley two-graph on `P¹(F_q)`; see e.g. Gunderson–Semeraro, arXiv:2204.10775, and Babai–Cameron on switching classes of tournaments).
  - The census formulas of part 5, with the signed multiplicities and the Désarménien–Foata value `−8`.
  - The survivor census for `k = 7`.
- **Part 7.** It is elementary and classical in substance: for regular `T`, `TT4 = n·C(k,3) − Σ_v c3(N⁺(v))` and `c4 = n(n²−1)(n−3)/48 − Σ_v c3(N⁺(v))`, so the bound is the known lower bound `c4(T) ≥ c4(DR_n)` on 4-cycles of regular tournaments, with the same equality case (see Savchenko, arXiv:2403.07629, §1, and the references there). The class-count relations are standard for regular tournaments.
