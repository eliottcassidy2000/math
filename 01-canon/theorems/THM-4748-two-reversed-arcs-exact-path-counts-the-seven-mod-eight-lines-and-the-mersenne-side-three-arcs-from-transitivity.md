---
id: THM-4748
title: "Two reversed arcs: exact Hamiltonian path counts, the seven-mod-eight lines, and the Mersenne side three arcs from transitivity. For the transitive order with two arcs reversed (spans >= 2), H has an exact formula in each of the five interval configurations: disjoint, touching, crossing (with a 3^r term), strictly nested, shared endpoint. Consequently every tournament with feedback arc number fa <= 2 and H = 7 (mod 8) has H in {15, 23} or on one of the two lines 3*2^j - 1, 5*2^j - 1; so for k >= 5 every tournament with H = 2^k - 1 has fa >= 3. The bound is attained for every k: the Fibonacci fan (an apex losing exactly to floor(k/2), k + 1 - floor(k/2) and k + 1 of a transitive chain) has H = 2^k - 1 and three back arcs, while 2^k + 1 = H(NT_(k+2)) needs one"
status: >
  PROVED: statements 1-5 (all orders). Statement 1 is a run-decomposition count; statements 2 and 3 are its 2-adic case
  analysis; statement 4 uses the one-apex formula, which is KNOWN in the repo (05-knowledge/variables/single-core-cycle-count.md,
  Session16/25/26 one-apex packets). FINITE-EXACT: statement 6. Prior art: two-arc counts may appear in M. Abi Aad's
  formulas for special reversed-arc families (Discrete Appl. Math. 2020, doi 10.1016/j.dam.2019.09.005, and DAM 2026,
  doi 10.1016/j.dam.2026.07.017; neither paper accessed); no literature priority is claimed for statement 1. The
  observation that within one-apex tournaments 2^k - 1 needs three apex in-arcs, against one for 2^k + 1 and two for
  3*2^k - 1 and 5*2^k - 1, was made the same day by mac-mini's wide-net survey (MSG-2400, brute force n <= 8); statements
  3 and 4 prove it for all tournaments and all k. Session opus-2026-10-21-S51.
session: opus-2026-10-21-S51
source: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md
scripts:
  - 04-computation/experiments/repunits_hyp9357_20261021/two_arc_formulas.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/three_arc_census.py (+ .out)
related:
  - THM-1370 (H never 7 or 21)
  - THM-344 (the two n = 8 classes with H = 63 are one-apex, Omega = K_31)
  - THM-466 (2-adic digits of H are odd-cycle collection counts)
  - THM-4744 (Moon's interval chains; NT_s has 2^(s-2) + 1 paths)
  - THM-4749 (the repunit degree law)
  - HYP-9359 (the arithmetic-support law of strong counts)
---

# THM-4748 — two reversed arcs, the seven-mod-eight lines, and the Mersenne side

## Setting

* The **transitive order** on `{0, 1, …, m−1}` has `i → j` for `i < j`.
* A **back arc** `y → x` (`y > x`) is a reversed pair. Its interval is `[x, y]` and its span is `y − x`.
* Every tournament `T` is the transitive order of some ordering with `fa(T)` back arcs, where `fa(T)` is the minimum feedback arc number.
  * In such a minimum ordering every back arc has span at least 2: swapping two adjacent vertices removes a span-1 back arc and changes nothing else.
* `NT_s` is the transitive order on `s` vertices with the single back arc `(s−1) → 0`. Its count `H(NT_s) = 2^(s−2) + 1` is classical.

## Statements

1. **Two-arc formulas (PROVED).** Reverse two arcs of spans at least 2.
   * **Disjoint or touching intervals.** The tournament has two strong blocks, and `H = (2^(a) + 1)(2^(b) + 1)`, where `a + 1` and `b + 1` are the spans. Touching intervals give Moon's chain.
   * **Crossing intervals**, `x_1 < x_2 < y_1 < y_2`. Put `p = x_2 − x_1 − 1`, `r = y_1 − x_2 − 1` and `q = y_2 − y_1 − 1`. Then
     `H = 1 + 2^(p+r+1) + 2^(q+r+1) − [q = 0]·2^p − [p = 0]·2^q + 2^(p+q+1)·3^r`.
   * **Strictly nested intervals**, `x_1 < x_2 < y_2 < y_1`. Put `s = y_2 − x_2`, `x = x_2 − x_1 − 1` and `y = y_1 − y_2 − 1`. Then
     `H = 1 + 2^(s−1) + 2^(x+y+1)·(2^s − 1 + 3^(s−1))`.
   * **Intervals sharing one endpoint**, with outer span `L` and inner span `s`, `2 ≤ s ≤ L − 1`:
     `H = 1 + 2^(L−1) − 2^(L−1−s) + 2^(s−1) − [s = L − 1]`.
2. **The seven-mod-eight lines (PROVED).** Suppose `fa(T) ≤ 2` and `H(T) ≡ 7 (mod 8)`. Then `H(T) ∈ {15, 23} ∪ {3·2^j − 1 : j ≥ 3} ∪ {5·2^j − 1 : j ≥ 3}`, and every listed value occurs.
   * The two lines come from shared endpoints with `s = L − 1` (giving `3·2^(L−2) − 1`) and `s = L − 2` (giving `5·2^(L−3) − 1`).
   * `3·2^j − 1` is the Thābit sequence.
3. **The Mersenne side needs three arcs (PROVED).** For `k ≥ 5`, every tournament with `H = 2^k − 1` has `fa ≥ 3`.
   * For `k = 3`, `H = 7` never occurs (THM-1370).
   * For `k = 4`, `15` is the chain `IC(2,3)` with `fa = 2`.
   * For `k = 2`, `3` is `C_3` with `fa = 1`.
4. **The Fibonacci fan (PROVED).** For `k ≥ 4`, put `a = ⌊k/2⌋` and `b = k + 1 − a`, so that `b − a ∈ {1, 2}`.
   * Let `F_k` be the tournament on `{0, 1, …, k+1}` in which:
     * `1, …, k+1` form a transitive chain;
     * vertex `0` loses exactly to `a`, `b` and `a + b = k + 1`, and beats every other vertex.
   * Then `F_k` is strong and has exactly three back arcs in the order `0 < 1 < … < k+1`, and `H(F_k) = 2^k − 1`.
   * With statement 3: for every `k ≥ 5` the least feedback arc number of a tournament with `H = 2^k − 1` is exactly 3.
5. **The Fermat side (KNOWN).** `H(NT_(k+2)) = 2^k + 1` with one back arc.
   * Both sides of `2^(2k) − 1 = (2^k − 1)(2^k + 1)` are therefore realized at the same order `k + 2`: the Fermat side one arc from transitivity, the Mersenne side three.
6. **Finite record (FINITE-EXACT).**
   * The class census for `n ≤ 9`: the classes with `H = 31, 63, 127, 255` have least feedback arc number 3 (`two_arc_formulas.py`).
   * Three reversed arcs (all configurations, `m ≤ 13`) reach `2^k − 1` for every `4 ≤ k ≤ 11` and for `k = 15`. Their values `≡ 7 (mod 8)` are no longer on finitely many lines: the largest odd part of `H + 1` is 13, 37, 97, 343, 829, 3737, 7267 at `m = 7, …, 13` (`three_arc_census.py`).

## Proofs

**1. The run decomposition.**
* **Runs and descents.** Write a Hamiltonian path as a permutation, cut into maximal increasing runs. Every descent of the path is a back arc, and each back arc is used at most once. Inside a run, consecutive vertices `u < v` need `u → v`, so a run may not contain the two ends of a back arc consecutively.
* **The count.** Fix the sequence of back arcs used as descents. Assigning the vertices to runs determines the path: each run is increasing. A vertex inside the intervals of `c` used arcs has `c + 1` admissible runs. The assignments that place a reversed pair consecutively inside one run are subtracted.
* **No descent.** One path, `0 1 ⋯ m−1`. A span of at least 2 keeps every reversed pair non-adjacent.
* **A single back arc `y → x`.** One descent: the first run ends at `y` and the second starts at `x`. Every vertex strictly inside `[x, y]` chooses a run, and every other vertex is forced. That gives `2^(span−1)` paths, so `H(NT_s) = 1 + 2^(s−2)`.

Vertices outside the union of the two intervals are forced singletons and change nothing. So take the union to be `[0, m−1]`.

* **Crossing** (`0 = x_1 < x_2 < y_1 < y_2 = m−1`).
  * One descent at `y_1 → 0` gives `2^(y_1 − 1)` assignments. Of these, `2^(x_2 − 1)` put `x_2` immediately before `m−1` in one run, and that happens only if `y_1 = m − 2`.
  * One descent at `(m−1) → x_2` gives `2^(m − x_2 − 2)` assignments, less `2^(m − 2 − y_1)` if `x_2 = 1`, where `0` would precede `y_1` directly.
  * For two descents in either order, the middle run joins the two arcs. Vertices left of `x_2` have 2 runs, vertices in `(x_2, y_1)` have 3, and vertices right of `y_1` have 2. No reversed pair can become consecutive, so each order gives `2^(p+q)·3^r`.
  * Summing gives the crossing formula.
* **Strictly nested** (outer arc `(m−1) → 0`, inner arc `y_2 → x_2`, `0 < x_2 < y_2 < m−1`).
  * One descent at the outer arc gives `2^(m−2)` assignments. Of these, `2·2^(m−3−s)` put `x_2, y_2` in one run with every vertex between them in the other run.
  * One descent at the inner arc gives `2^(s−1)`. The vertices `0` and `m−1` are then in different runs.
  * Two descents in either order give `2^(x+y)·3^(s−1)` each.
  * With `m − 1 = L` and `L − 1 − s = x + y + 1`, the total is the nested formula.
* **Shared endpoint** (outer `(m−1) → 0`, inner `s → 0`; the other case is the converse picture).
  * One descent at the outer arc gives `2^(L−1)` assignments, less `2^(L−1−s)` that put `0, s` consecutively in the second run.
  * One descent at the inner arc gives `2^(s−1)`, less 1 when `s = L − 1`, because then `0` and `m−1` can be adjacent in the second run.
  * Two descents would make `0` start two runs, so there are none.
* **Disjoint and touching intervals.** They multiply (Moon).

The check is in `two_arc_formulas.py`: every formula equals an exact Hamiltonian-path DP in all 9,009 two-arc configurations with `m ≤ 14`. ∎

**2.** Since `H(T) ≡ 7 (mod 8)`, `T` is not transitive. If `fa = 1`, `H = 2^a + 1 ≡ 3, 5, 1 (mod 8)`, so `fa = 2`. Reduce each formula of statement 1 modulo 8.
* **Products** `(2^a + 1)(2^b + 1)`. The factors are `≡ 3, 5, 1` for `a = 1, 2, ≥ 3`, so the product is `≡ 7` only for `{a, b} = {1, 2}`: `H = 15`.
* **Crossing.**
  * If `p, q ≥ 1`, every power of two in `H − 1` has exponent at least 2, so `H ≡ 1 (mod 4)`.
  * If `p = q = 0`, `H − 1 = 2^(r+2) − 2 + 2·3^r`, which is `≡ 4` or `0 (mod 8)`.
  * If `p = 0 < q`, `H − 1 = 2^(r+1) + 2^(q+r+1) − 2^q + 2^(q+1)·3^r`. This is `≡ 6 (mod 8)` exactly at `(q, r) = (2, 0)`, giving 15, and at `(q, r) = (1, 1)`, giving 23.
  * The case `q = 0 < p` is symmetric.
* **Strictly nested.** For `s ≥ 3`, `H − 1 ≡ 0 (mod 4)`. For `s = 2`, `H − 1 = 2 + 3·2^(x+y+2)`, which is `≡ 6 (mod 8)` only when `x = y = 0`: `H = 15`.
* **Shared endpoint.** Put `t = L − 1 − s`.
  * `t = 0`: `H = 3·2^(L−2) − 1`, which is `≡ 7` for `L ≥ 5`.
  * `t = 1`: `H = 5·2^(L−3) − 1`, which is `≡ 7` for `L ≥ 6`, while `L = 4, 5` give 9 and 19.
  * `t ≥ 2`: for `s ≥ 3`, `H ≡ 1 (mod 4)`; for `s = 2`, `H = 3·2^(L−3) + 3`, which is `≡ 7` only at `L = 5`, giving 15.

A symbolic scan in `two_arc_formulas.py` covers every formula with parameters up to 40–80. ∎

**3.** For `k ≥ 3`, `2^k − 1 ≡ 7 (mod 8)`. Suppose `fa(T) ≤ 2`. By statement 2, `2^k − 1` would be 15, 23, `3·2^j − 1` or `5·2^j − 1`. The last two would force `3·2^j = 2^k` or `5·2^j = 2^k`, which is impossible, and 23 is not a repunit. So `k = 4`. ∎

**4.** `T − 0` is transitive, so every cycle passes through vertex 0. The odd-cycle conflict graph is therefore complete, and `H = 1 + 2α_1` (OCF; THM-466).
* **Counting the odd cycles.** Let 0 lose exactly to the set `A ⊆ {2, …, m−1}` of chain positions. An odd cycle is `0 → i → (an even-size subset of the chain strictly between i and t) → t → 0`, with `i ∉ A` and `t ∈ A`. There are `w(t − i)` of them, where `w(1) = 1` and `w(d) = 2^(d−2)`.
* **The one-apex formula.** Since `Σ_{i=1}^{t−1} w(t − i) = 2^(t−2)`,
  `H − 1 = Σ_{t∈A} 2^(t−1) − Σ_{s<t in A} 2^max(t−s−1, 1)`
  (KNOWN: the single-core formula; re-checked here against a DP for every `A` with `m ≤ 12`).
* **The fan.** For `A = {a, b, a + b}` with `a ≥ 2` and `b − a ∈ {1, 2}`, the pair terms are `2^(b−1)` (from `a, a+b`), `2^(a−1)` (from `b, a+b`) and `2` (from `a, b`). The first two cancel the singleton terms of `a` and `b`, leaving `H − 1 = 2^(a+b−1) − 2`.
* **Strong, with three back arcs.** `0 → 1 → ⋯ → k+1 → 0` is a Hamiltonian cycle, since `1 ∉ A`. The back arcs in the order `0 < 1 < ⋯ < k+1` are exactly `a → 0`, `b → 0` and `(k+1) → 0`.
* **Checks.** `two_arc_formulas.py` checks `H(F_k) = 2^k − 1` by DP for `k ≤ 14`, and by the formula for `k ≤ 400`. ∎

**5.** This is the single-back-arc count in the proof of statement 1. ∎

## Remarks

* **Why the asymmetry.** One back arc of span `s` contributes `+2^(s−1)` paths: each vertex under it has two choices. A sum of positive powers of two plus 1 lands on the Fermat side `2^k + 1` immediately. The Mersenne side `2^k − 1 = 2^k + 1 − 2` needs a subtraction. Subtractions come only from arcs sharing an endpoint (the negative pair coefficients `−2^max(gap, 1)`). Removing everything but `2^(k−1)` and the final `−2` takes exactly the Fibonacci relation `a + b = c` among three apex spans.
* **Tiling-cube reading** (mac-mini's lane: the Möbius coefficients are all even). For one-apex tournaments the cube expansion stops at pairs. A singleton coefficient is `+2^(s−1)`, and an endpoint-sharing pair coefficient is `−2^max(t−s−1, 1)`. The fan is the smallest exact cancellation.
* **What this does not say.** It does not say that `fa ≥ 3` alone makes a value rare. Multiplicities are measured in THM-4749 and HYP-9359. 15 and 255 are common because they are products of Fermat numbers (chains).
