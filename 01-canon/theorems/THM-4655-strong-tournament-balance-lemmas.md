---
id: THM-4655
title: "Uniform Hamiltonian paths of strong tournaments: the 1/3-2/3 property holds for every strong 'transitive-plus-one' tournament, with the explicit balanced pair (z, u_2); there the extra vertex z is balanced against every interior vertex, every vertex has a balanced partner and the unique Hamiltonian cycle carries at least two balanced arcs; the exact transfer chain at a cut vertex, and with it: a cut vertex z is balanced against every vertex of each interior strong component of T - z with at most 7 vertices, whatever the rest of T (computer-assisted); the vertex-wise constant is at most 1/3 for every n >= 5 (= 1/3 for 5 <= n <= 10); a rotation lemma for Hamiltonian cycles; Stanley-type log-concavity fails; exhaustive verification for all strong tournaments on n <= 10 vertices"
status: >
  PROVED: statements 1, 2, 5, 6 (elementary; statement 1(f) uses an exact computer check for m <= 14 together with an
  analytic bound for m >= 9, and is reproved by statement 7). PROVED (computer-assisted): statement 7 (a general reduction to a
  local certificate on one block, checked exactly for all blocks with at most 7 vertices; audited by audit E part 2, which
  reproduced the certificate independently). FINITE-EXACT: statements 3-4 (exhaustive, gentourng + exact subset DP in C). The tournament 1/3-2/3
  conjecture HYP-9277 (opus S35) remains OPEN; this theorem proves it, and its vertex-wise strengthening HYP-9279(B), for
  a natural subclass of strong tournaments, and extends the exhaustive verification from n <= 8 to n <= 10.
source: mac-mini-2026-10-14-tourn
depends_on: [HYP-9277 and its proved reductions (opus S35, one_twenty_eight_atlas_balance_gold_20261014.md: delta(T) = max over strong components; with no balanced pair the >2/3-majority relation is a linear order)]
verification:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out): exhaustive strong n <= 9
  - 04-computation/experiments/tournament_balance_20261014/tbal2.c (+ tbal2_strong9.out, tbal2_strong10.out): exhaustive strong n <= 10, with per-vertex partners
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py, tpo_proof_steps.py (+ .out), tpo_u2_check.c (+ .out)
  - 04-computation/experiments/tournament_balance_20261014/tpo_cycle_arcs.py (+ .out), tpo_vertexwise_proof_check.py (+ .out): statements 1(e)-(g)
  - 04-computation/experiments/tournament_balance_20261014/cutvertex_chain.py (+ .out): statement 5, validated against the subset DP
  - 04-computation/experiments/tournament_balance_20261014/vw_tight_family.py (+ .out): statement 6
  - 04-computation/experiments/tournament_balance_20261014/z1_cone_certificate.py (+ z1_cert_upto5.out, z1_cert_6.out, z1_cert_7.out), z1_local_slack.py (+ .out): statement 7
  - 04-computation/experiments/tournament_balance_20261014/rotation_lemma_check.py, logconcave_test.py, analyse_minimisers.py
note: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
related: [HYP-9277, HYP-9279, THM-438 (Paley Hamiltonian paths), THM-4652 (Paley pattern sums)]
---

# THM-4655: balance of uniform Hamiltonian paths in strong tournaments

## Setting

- `T` is a tournament on `n` vertices, and `P` is a uniformly random Hamiltonian path of `T` (every tournament has one).
- `p_T(x, y)` is the probability that `x` precedes `y` in `P`. The balance is `δ(T) = max_(x≠y) min(p_T(x,y), 1 − p_T(x,y))`.
- HYP-9277 (opus S35) is the tournament 1/3–2/3 conjecture: `δ(T) >= 1/3` for every non-transitive `T`, with equality exactly at transitive compositions of singletons and cyclic triangles.
- Its proved reductions:
  - `δ(T)` is the maximum over the strong components.
  - If no pair is balanced, the relation `p > 2/3` is a linear order `σ*`.

  So the conjecture is equivalent to `δ >= 1/3` for strong tournaments.

## Statement

1. **(Transitive-plus-one tournaments.)** Let `T` consist of a transitive tournament `u_1 → u_2 → … → u_m` (that is, `u_i → u_j` for `i < j`) together with a vertex `z`, where `z → u_i` iff `w_i = 1` for a word `w ∈ {0,1}^m`.
   - **(a)** `T` is strong iff `w_1 = 1` and `w_m = 0`.
   - **(b) Paths.** The Hamiltonian paths of a strong `T` correspond bijectively to the sets `A` of `u`'s that precede `z`, each run increasing, subject to: (`A = ∅` or `w_(max A) = 0`) and (`B = ∅` or `w_(min B) = 1`), where `B` is the complement of `A`. Hence

     `H(T) = |K| + Σ_(β<α, w_β = 1, w_α = 0) 2^(α−β−1)`,   `K = {0, m} ∪ {k : w_k = 0, w_(k+1) = 1}`.

     For `w = 1^(m−1) 0` this is the nearly transitive tournament, with `H = 2^(n−2) + 1` (KNOWN: El Sahili–Ghazo Hanna, arXiv:2207.11510).

   - **(c) The balanced pair.** If `m >= 3` (`n >= 4`), then

     `1/3 < p_T(u_2 before z) <= 2/3`,

     with equality iff `w = 1 0 1^(m−3) 0` (`m >= 4`). Hence `δ(T) >= 1/3`: **the 1/3–2/3 property holds for every strong transitive-plus-one tournament, with the explicit balanced pair `(z, u_2)`.** For `m = 2`, `T = C_3` and `δ = 1/3`.
   - **(c′) Strictness (audit E).** For `m >= 3`, `δ(T) > 1/3`. So HYP-9277's equality characterisation also holds on this family.
     - If `w ≠ 1 0 1^(m−3) 0`, the pair `(z, u_2)` is strictly balanced by (c).
     - At the equality word with `m >= 5`, `p(u_3 before z) = 2^(m−4)/(1 + 2^(m−3)) ∈ (1/3, 1/2)`.
     - At `w = 1010`, `p(u_1 before u_4) = 5/9`.
   - **(d) Corollary (with S35's component reduction).** Let `T` be a non-transitive tournament in which some strong component on at least 3 vertices becomes transitive after deleting one vertex. Then `δ(T) >= 1/3`.
   - **(e) Balanced arcs on the Hamiltonian cycle.** The unique Hamiltonian cycle of a strong `T` is `z u_1 u_2 … u_m z`. Let `m >= 3`.
     - If `w_2 = 1`, then `1/2 <= p_T(z before u_1) < 2/3`.
     - If `w_2 = 0`, then `1/3 <= p_T(u_2 before u_1) <= 1/2`, with `1/3` iff `w = 1 0 1^(m−3) 0`.

     Applied to the reverse tournament, the same holds at the other end: `(u_m, z)` or `(u_(m−1), u_m)` is balanced. So **the Hamiltonian cycle carries at least two balanced arcs**.

     Two is attained for every `m >= 3` by `w = 1 0^(m−1)`, where `H = 2^(m−1) + 1`.
     - Its balanced cycle arcs are exactly:
       - `(u_1, u_2)`, with `p(u_2 before u_1) = 2^(m−2)/H`;
       - `(u_m, z)`, with `p(u_m before z) = (2^(m−2) + 1)/H`.
     - The other cycle arcs have `p(z before u_1) = 2^(m−1)/H` and `p(u_(i+1) before u_i) = 2^(m−3)/H`.
     - Audit E checked this exactly for `m = 3..16`.
   - **(f) `z` is a universal partner.** If `m >= 3`, then `1/3 <= p_T(u_i before z) <= 2/3` for every `2 <= i <= m − 1`.
   - **(g) Corollary: vertex-wise balance and connectivity.** In every strong transitive-plus-one tournament:
     - every vertex has a balanced partner: `z` and `u_2` by (c), `u_2, …, u_(m−1)` and `z` by (f), `u_1` with `z` or `u_2` and `u_m` with `z` or `u_(m−1)` by (e);
     - the graph of balanced pairs is connected, since it contains the star from `z` to `u_2, …, u_(m−1)` with `u_1` and `u_m` attached.

     So HYP-9279(B) and (B⁺) hold on this family.

2. **(Rotation lemma.)** Let `C` be a Hamiltonian cycle and `σ` a linear order of `V`. Cutting `C` before each of its `n` vertices gives `n` Hamiltonian paths (the closable paths of `C`).
   - Averaged over these paths and over all pairs, the fraction of pairs ordered as in `σ` is at most `(2n − 1)/(3n) < 2/3`, with equality iff `C` visits `V` in `σ`-order.
   - **Corollary.** If a strong tournament has no balanced pair, fewer than a fraction `n/(n+1)` of its Hamiltonian paths are closable (`last → first`).

3. **(No Stanley inequality; FINITE-EXACT.)** For uniform Hamiltonian paths, neither of the following distributions is log-concave in general:
   - the position distribution of a vertex;
   - the gap distribution of `pos(y) − pos(x)`.

   Failures appear at `n = 4`. The counts of strong tournaments failing are:

   | `n` | 4 | 5 | 6 | 7 |
   |---|---|---|---|---|
   | fail for positions | 1 of 1 | 4 of 6 | 35 of 35 | 351 of 353 |
   | fail for gaps | 1 of 1 | 6 of 6 | 35 of 35 | 353 of 353 |

   This is audit E's independent recount. So the Kahn–Saks route to `δ >= 3/11` for posets, which runs through Stanley's inequality, does not transfer verbatim.

4. **(Exhaustive; FINITE-EXACT.)** Every strong tournament on `n <= 10` vertices (1, 1, 6, 35, 353, 6008, 178133, 9355949 for `n = 3..10`) has `δ >= 1/3`. The minima are:

   | `n` | minimum `δ` |
   |---|---|
   | 3 | `1/3` |
   | 4 | `2/5` |
   | 5, 6, 7 | `4/11` |
   | 8, 9, 10 | `62/157` |

   - The minimisers are near-transitive (minimum feedback arc sets 1, 2, 3, 5, 4, 6 for `n = 4..9`).
   - With S35's reduction, **HYP-9277's inequality holds for every tournament on `n <= 10` vertices** (S35 had `n <= 8`).
   - Audit E replicated `n <= 9` with independent code: DFS on gentourng `-z` output, its own decoder, and a completeness check by `Σ n!/|Aut|`. The `n = 10` sweep rests on `tbal2.c` alone; audit E rechecked only its minimiser (`δ = 62/157`, `H = 1413`).
   - In every strong tournament with `n <= 10`, every vertex has a balanced partner (HYP-9279(B)).

5. **(The cut-vertex transfer chain and the fair-coin lemma.)** Let `T` be strong and `z` a vertex with `T − z` not strong, and let `D_1 ⇒ D_2 ⇒ … ⇒ D_k` (`k >= 2`) be the strong components of `T − z`. For `X ⊆ D_j` let `h(X)`, `h⁺(X)`, `h⁻(X)` count the Hamiltonian paths of `T[X]`, of those starting in `N⁺(z)`, and of those ending in `N⁻(z)`, all three equal to 1 at `X = ∅`. Write `D = D_j` and `X^c = D ∖ X`.
   - **(a)** `H(T) = u^T M_1 M_2 ⋯ M_k f`. The states are `pre, gap, mid, post`, with `u = e_pre + e_gap` and `f = e_gap + e_post`. The only nonzero entries of `M_j` are:

     | transition | entry |
     |---|---|
     | `pre → pre` | `h(D)` |
     | `pre → gap` | `h⁻(D)` |
     | `pre → mid` | `Σ_(X ≠ D) h(X) h⁺(X^c)` |
     | `pre → post` | `Σ_(X ≠ ∅, D) h⁻(X) h⁺(X^c)` |
     | `gap → post` | `h⁺(D)` |
     | `mid → mid` | `Σ_X h(X) h(X^c)` |
     | `mid → post` | `Σ_(X ≠ ∅) h⁻(X) h(X^c)` |
     | `post → post` | `h(D)` |

   - **(b)** For `v ∈ D_i`, the number of Hamiltonian paths with `v` before `z` is `u^T M_1 ⋯ M_i^(v) ⋯ M_k f`. Here `M_i^(v)` keeps `pre → pre` and `pre → gap`, restricts every sum over `X` to `X ∋ v`, and sets `gap → post` and `post → post` to 0.
   - **(c) (Fair coin.)** The `mid → mid` entry is invariant under `X ↔ D ∖ X`. So its `v`-part is exactly half of it. Conditionally on block `D_i` lying strictly inside the crossing window (`b < i < a` below), every vertex of `D_i` precedes `z` with probability exactly `1/2`.
   - **(d)** When `T − z` is transitive (all blocks singletons), `mid → mid = 2`, and (a) is the crossing decomposition of 1(b), with weights `2^(α−β−1)`.

6. **(The vertex-wise constant is at most `1/3` for every `n >= 5`.)** Let `V_n` (`n >= 5`, `t = n − 3`) be the transitive tournament `0 → 1 → … → t−1 → {a, b, c}`, with the cyclic triangle `a → b → c → a`, and with one arc reversed to `b → 0`. Then:
   - `V_n` is strong and `H(V_n) = 3(2^(n−4) + 1)`;
   - in a uniform Hamiltonian path, the order induced on `{a, b, c}` is uniform over the three rotations `abc, bca, cab`;
   - `p(c before a) = 2/3`, `p(c before b) = 1/3`, `p(c before 0) = 2^(t−1)/H`, and `p(c before i) = 2^(t−2)/H` for `1 <= i <= t−1`.

   So `max_w min(p(c, w), 1 − p(c, w)) = 1/3` exactly. The vertex-wise constant `min_T min_v max_w min(p, 1 − p)` is therefore at most `1/3` for every `n >= 5`. It equals `1/3` wherever HYP-9279(B) holds, in particular for `5 <= n <= 10` (exhaustive; corrected by audit E part 2). It is `1/3` at `n = 3` (`C_3`) and `2/5` at `n = 4`.
   - `V_n − b` is transitive, in the order `0, …, t−1, c, a`. So `V_n` is the transitive-plus-one tournament with `z = b`, and it is the reverse of the one with the equality word `1 0 1^(m−3) 0` of (c), where `m = n − 1`.
   - Audit E found the same sharpness from that side: at the equality word, `u_2` has `p = 2/3, 1/3` against `z, u_1` and is unbalanced against everything else (S2 in its report).
   - Within the transitive-plus-one family, only this word and its reverse attain `1/3`, for each `m >= 5`.

7. **(Cut vertices with a small interior component; computer-assisted.)** Let `T` be strong and `z` a vertex such that `T − z` has strong components `D_1 ⇒ … ⇒ D_k` with `k >= 3`. Let `2 <= i <= k − 1` with `|D_i| <= 7`. Then

   `1/3 <= p_T(v before z) <= 2/3` for every `v ∈ D_i`.

   - The other components may be arbitrary, of any size, and so may `n`.
   - Corollary: if some interior component of `T − z` has at most 7 vertices, then `z` has a balanced partner and `δ(T) >= 1/3`.
   - This proves HYP-9279 (Z1) for interior components with at most 7 vertices, and reproves 1(f), which is the case of singleton blocks.

## Proof

**1(a).**
- If `w_1 = 0`, then `u_1` beats everything, so it is a source and `T` is not strong.
- If `w_m = 1`, then `u_m` is a sink.
- Otherwise `u_1 → … → u_m → z → u_1` is a Hamiltonian cycle.

**1(b).**
- Consecutive `u`'s in a Hamiltonian path are increasing, so the `u`'s before and after `z` form increasing runs `A` and `B`. The junctions require `w_(max A) = 0` and `w_(min B) = 1`.
- Either `A = {u_1, …, u_k}` is an initial segment, which is valid iff `k ∈ K`; or `β := min B < α := max A`. In the second case:
  - every index below `β` lies in `A`;
  - every index above `α` lies in `B`;
  - the `α − β − 1` indices strictly between are free.

  Such a split is valid iff `w_β = 1` and `w_α = 0`.

**1(c).** Write:
- `C₁ = Σ_(w_α = 0) 2^(α−2)`, the crossed weight with `β = 1`;
- `C₃ = Σ_(β>=3, w_β=1) Σ_(α>β, w_α=0) 2^(α−β−1)`;
- `N₀ = #{α >= 4 : w_α = 0}`;
- `in2` for the number of valid `A` containing `u_2`.

Since `w_1 = 1`, the initial segments containing `u_2` are all of `K` except `k = 0`, and `1 ∉ K`. The crossed sets with `β = 2` exist only if `w_2 = 1`, and then their weight is `C₁/2`. With `β = 1`, `u_2` is free (`α >= 3`) or equal to `α` (`α = 2`, possible only if `w_2 = 0`, with weight 1). Hence:
- if `w_2 = 1`: `in2 = |K| − 1 + C₃ + C₁/2` and `H = |K| + 3C₁/2 + C₃`;
- if `w_2 = 0`: `in2 = |K| + C₃ + (C₁ − 1)/2` and `H = |K| + C₁ + C₃`.

Two lemmas:
- **Lemma A: `|K| <= 2 + N₀`.** Send each 01-boundary `(k, k+1)` (necessarily `k >= 2`) to the first zero after `k + 1`. That zero lies at position `>= 4`, and the map is injective.
- **Lemma B: `2C₃ <= Σ_(α>=4, w_α=0) (2^(α−2) − 2)`.** For a fixed zero `α`, `Σ_(3<=β<α) 2^(α−β−1) = 2^(α−3) − 1`.

The lower bound:
- `w_2 = 1`: `3·in2 − H = 2|K| + 2C₃ − 3 >= 1`.
- `w_2 = 0`: `3·in2 − H = 2|K| + 2C₃ + C₁/2 − 3/2 > 0`.
- In both cases `|K| >= 2`.

The upper bound:
- `w_2 = 1`: Lemma B gives `C₃ <= C₁/2 − N₀`, so `2H − 3·in2 = 3 − |K| + 3C₁/2 − C₃ >= 3 − |K| + C₁ + N₀ >= 1 + C₁ > 0` by Lemma A.
- `w_2 = 0`: `C₁ = 1 + 2[w_3 = 0] + Σ_(α>=4, w_α=0) 2^(α−2)`. Lemmas A and B then give `2H − 3·in2 = 3/2 − |K| + C₁/2 − C₃ >= 2 + N₀ − |K| + [w_3 = 0] >= 0`.
- Equality requires `w_3 = 1`, equality in Lemma B (every position between 3 and each zero is a 1), and equality in Lemma A. That forces `w = 1 0 1^(m−3) 0`, which does attain `2/3`.

`tpo_proof_steps.py` checks every identity and inequality above on all 8,190 strong words with `m <= 14`.

`tpo_u2_check.c` enumerates all valid `A` for every word with `m <= 20`: 524,286 words, no violation, and the extremal words are as stated. This check is independent of the closed-form decomposition but rests on the bijection 1(b). Audit E checked 1(b) against DFS enumeration of the actual Hamiltonian paths for all 16,382 words with `m <= 13`. ∎

**1(d).** By S35's reduction, `δ` of a tournament is the maximum of `δ` over its strong components. Apply 1(c). ∎

**1(e).** We keep the notation of 1(c).
- `u_1 ∈ B` iff `A = ∅` or the split is crossed with `β = 1`. So `p(z before u_1) = (1 + C₁)/H`.
- `u_2 ∈ A` and `u_1 ∈ B` iff `β = 1` and either `α = 2` (weight 1, only if `w_2 = 0`) or `α >= 3` with the free `u_2` on the `A` side. So `p(u_2 before u_1) = (C₁ + [w_2 = 0])/(2H)`.

Case `w_2 = 1`, where `H = |K| + 3C₁/2 + C₃`:
- `2H − 3(1 + C₁) = 2|K| + 2C₃ − 3 > 0`, since `|K| >= 2`.
- `2(1 + C₁) − H = 2 + C₁/2 − |K| − C₃ >= 0`. This uses Lemma A and `2C₃ <= C₁ − 2N₀`; the latter is Lemma B, because every zero `α >= 4` contributes `2^(α−2)` to `C₁`.

Case `w_2 = 0`, where `H = |K| + C₁ + C₃`:
- `H − (C₁ + 1) = |K| + C₃ − 1 > 0`.
- `3(C₁ + 1) − 2H = C₁ + 3 − 2|K| − 2C₃ >= 2[w_3 = 0] >= 0`, by the computation in the upper bound of 1(c). Equality forces `w = 1 0 1^(m−3) 0`.
- In this case `p(u_2 before u_1) = 1 − p(u_2 before z)`. The bijection pairs `A = ∅` with the crossing `(1, 2)`, and flips the coin of `u_2` in the crossings `(1, α)` with `α >= 3`.

For the reverse tournament `T^R`, take the word `w'_j = 1 − w_(m+1−j)`, so that `u_i ↦ u_(m+1−i)`. Every probability is replaced by its complement, and balanced pairs are the same. This gives the statement at the other end. For `m >= 3` the arc pairs `{(z,u_1), (u_1,u_2)}` and `{(u_(m−1),u_m), (u_m,z)}` are disjoint. ∎

**1(f).** Fix `2 <= i <= m − 1` and sort the valid sets by where `u_i` lies:
- an initial segment `{u_1, …, u_k}`, `k ∈ K`, has `u_i ∈ A` iff `k >= i`;
- a crossing `(β, α)` has `u_i ∈ A` if `α = i` or `β > i`, has `u_i ∈ B` if `α < i` or `β = i`, and leaves `u_i` a fair coin if `β < i < α`.

Let:
- `X = #{k ∈ K : k >= i} + W(α = i) + W(β > i)`;
- `Y = #{k ∈ K : k < i} + W(α < i) + W(β = i)`;
- `W_M = W(β < i < α)`,

where `W(·)` is the total crossing weight of the given type. Then `H = X + Y + W_M` and `p(u_i before z) = (X + W_M/2)/H`. So the claim is `X <= 2Y + W_M/2` and `Y <= 2X + W_M/2`. Both follow from `2X <= W_M` and `2Y <= W_M`.
- For `i = 2` this is 1(c). By reversal we may assume `i <= (m+1)/2`.
- For `m <= 14` the claim is checked exactly for all 8,190 strong words and all `i` (`tpo_vertexwise_proof_check.py`). The check also confirms `2X, 2Y <= W_M` in all 38,144 cases with `8 <= m <= 14` and `3 <= i <= (m+1)/2`.
- For `m >= 9` and `3 <= i <= (m+1)/2`:
  - The crossing weight factorises: `2^(α−β−1) = 2^(α−1) · 2^(−β)`. Put `Z = Σ_(α>i, w_α=0) 2^(α−1)`, which is `>= 2^(m−1)` since `w_m = 0`.
  - Then `W_M = (Σ_(β<i, w_β=1) 2^(−β)) · Z >= Z/2`, from `β = 1`.
  - `W(β > i) <= 2^(−i) Z` and `W(β = i) <= 2^(−i) Z`, since the inner sum over `β > i` is below `2^(−i)`.
  - `W(α = i) < 2^(i−1)` and `W(α < i) < 2^(i−1)`.
  - The `K`-terms are at most `m`.
  - Hence `2X, 2Y < 2m + 2^i + 2^(1−i) Z <= 2m + 2^((m+1)/2) + Z/4`. This is `<= Z/2 <= W_M` as soon as `2m + 2^((m+1)/2) <= 2^(m−3)`, which holds for `m >= 9`. ∎

**1(g).** Immediate from 1(c), 1(e) and 1(f). ∎

**2.**
- For a pair `i < j` (in `σ`), `i` precedes `j` in exactly `(c_i − c_j) mod n` of the `n` cuts, where `c` gives the positions on `C`. So the concordance count is `n·C(n,2) − Σ_(i<j) y_ij`, with `y_ij = (c_j − c_i) mod n`.
- From a fixed `i`, the forward distances to the other vertices are exactly `1, …, n − 1`. So the `n − 1 − i` larger labels have distance sum at least `1 + … + (n − 1 − i)`.
- Summing over `i` gives `Σ_(i<j) y_ij >= (n³ − n)/6`, with equality iff each `i` is followed by `i+1, …, n−1` at distances `1, …, n−1−i`, i.e. the identity cycle. The average concordance is therefore at most `(2n − 1)/(3n)`.
- For the corollary: if no pair is balanced, every pair is concordant with `σ*` with probability `> 2/3`. Non-closable paths contribute at most 1, so `2/3 < f·(2n−1)/(3n) + (1 − f)`, where `f` is the closable fraction. ∎

**3, 4.** These are the exhaustive computations listed in `verification`.

**5.** Write a Hamiltonian path as `Q_A z Q_B`, where `A` is the set of vertices before `z` and `B` its complement in `V ∖ z`.

All arcs between blocks go forward, `D_j → D_l` for `j < l`. So a Hamiltonian path of `T[A]` visits `A ∩ D_1, A ∩ D_2, …` in this order, each contiguously, and any Hamiltonian paths of the pieces concatenate. The same holds for `B`. Hence the number of paths with a given `A` is
`Π_(j<a) h(A_j) · h⁻(A_a) · h⁺(B_b) · Π_(j>b) h(B_j)`,
with `A_j = A ∩ D_j`, `B_j = B ∩ D_j`, `a = max{j : A_j ≠ ∅}` and `b = min{j : B_j ≠ ∅}`.

Blocks below `b` lie in `A` and blocks above `a` lie in `B`. Two cases:
- `a < b`: then `b = a + 1`, and `A` is an initial union of blocks.
- `b <= a`: blocks strictly between `b` and `a` are split arbitrarily, with weight `h(A_j) h(B_j)`.

Scan the blocks and record the state after block `j` (`a = 0` if `A = ∅`, `b = k + 1` if `B = ∅`):
- `pre` if `j < min(a, b)`;
- `gap` if `a = j < b`, which forces `b = j + 1`;
- `mid` if `b <= j < a`;
- `post` if `j >= max(a, b)`.

The initial state is `gap` when `A = ∅` and `pre` otherwise. The final state is `gap` when `B = ∅` and `post` otherwise. The entries of `M_j` are exactly the block factors of the transitions, which proves (a).

Restricting block `i` to the splits with `v ∈ A_i` proves (b). The involution `X ↦ D ∖ X` fixes `h(X) h(D ∖ X)` and exchanges `v ∈ X` with `v ∉ X`, which proves (c). Part (d) is direct substitution.

`cutvertex_chain.py` checks (a) and (b) against the full subset DP on 300 random cut-vertex tournaments with `n <= 12`: no mismatch. ∎

**6.** The cycle `0 → 1 → … → t−1 → a → b → 0`, together with `b → c → a`, shows that `V_n` is strong.

The only in-neighbour of `0` is `b`, so a Hamiltonian path either starts at `0` or contains `b 0`.
- **Paths starting at `0`.** `V_n − 0` is `TT_(t−1) ⇒ C_3`. Its Hamiltonian paths are `1 2 … (t−1)` followed by a rotation of `abc`, and all three start at `1`, an out-neighbour of `0`. This gives 3 paths, one per rotation.
- **Paths containing `b 0`.** Write the path as `X b 0 Y`. The tournament on `{1, …, t−1, a, c}` is transitive, in the order `1, …, t−1, c, a`. So `X` and `Y` are the unique increasing paths on their vertex sets `S` and `S'`.
  - `Y` may start anywhere, since `0` beats all of `S'`.
  - `X` must end at an in-neighbour of `b`, that is `a` or some `i >= 1`. This rules out exactly `c ∈ S` with `a ∉ S`.
  - So these paths correspond to `S ∩ {1, …, t−1}` (`2^(t−1)` choices) times one of the three patterns:
    - `a, c ∈ S`, giving the order `c a b`;
    - `a ∈ S` only, giving `a b c`;
    - neither in `S`, giving `b c a`.

Hence `H = 3 + 3·2^(t−1)`, and each rotation of `abc` occurs in exactly `1 + 2^(t−1)` paths.

The vertex `c` precedes `0` only in the pattern `c a b`. It precedes `i >= 1` only when moreover `i ∉ S`. This gives the stated values, and both are below `1/3`. `vw_tight_family.py` confirms the values exactly for `n = 5, …, 14`. ∎

**7.** Normalise the chain of statement 5 by dividing `M_j` by `h(D_j)`. For a block `D` with `O = N⁺(z) ∩ D` and `I = N⁻(z) ∩ D`, write:
- `e = h⁻(D)/h(D)`, `s = h⁺(D)/h(D)`;
- `β = Σ_(X≠D) h(X) h⁺(X^c)/h(D)`, `γ = Σ_(X≠∅,D) h⁻(X) h⁺(X^c)/h(D)`;
- `σ = Σ_X h(X) h(X^c)/h(D)`, `α = Σ_(X≠∅) h⁻(X) h(X^c)/h(D)`.

The normalised entries are `pre → pre = post → post = 1`, `pre → gap = e`, `pre → mid = β`, `pre → post = γ`, `gap → post = s`, `mid → mid = σ` and `mid → post = α`. For `v ∈ D`, `β(v)`, `γ(v)`, `α(v)` are the parts with `v ∈ X`.

*(i) Elementary block inequalities.* `γ + s <= β <= σ − 1` and `e + γ <= α <= σ − 1`.
- `s` is the `X = ∅` term of `β`, and `h⁻ <= h`.
- `σ` contains the terms `X = D` and `X = ∅`, both equal to 1.

*(ii) End blocks.* If `z` has an out-neighbour in `D`, then `β >= 1` and `γ + s >= 1`.
- Cut each Hamiltonian path of `D` just before its first vertex in `O`. This is a `β`-split.
- If the cut is at the start, it is an `s`-term. Otherwise the vertex before the cut lies in `I`, and the cut is a `γ`-split.

Dually, if `D` has an in-neighbour of `z`, then `α >= 1` and `e + γ >= 1`, by cutting just after the last vertex in `I`. All these maps are injective, since concatenation inverts them.

*(iii) Invariant region.* Since `T` is strong and `D_1` receives no arc from the other blocks, `z` has an out-neighbour in `D_1`; dually, some vertex of `D_k` beats `z`. After blocks `1, …, i−1` the left vector is `(1, g, P, Lp)`, with `0 <= g <= 1`, `P >= 1` and `1 <= Lp <= P`. The base case is the first block: by (ii), `P = β_1 >= 1` and `Lp = γ_1 + s_1 >= 1`, and `Lp <= P` by (i). Each block maps
- `g ↦ e`,
- `P ↦ β + σP >= P`,
- `Lp ↦ γ + g s + αP + Lp`,

and the new `Lp` lies between the old `Lp` and `γ + s + (α + 1)P <= β + σP`. Symmetrically, the right vector from block `i+1` on is `(Rp, s', Q, 1)`, with `0 <= s' <= 1`, `Q >= 1` and `1 <= Rp <= Q`.

*(iv) The two inequalities.* For `v ∈ D_i`, statement 5 gives
- `N = Rp + e s' + β(v) Q + γ(v) + P(σQ/2 + α(v))`,
- `H = Rp + e s' + β Q + γ + g s + P(σQ + α) + Lp`.

Both `3N − H` and `2H − 3N` are affine in each of `g, s', Lp, Rp`. So on the region (iii) their minima are attained at `g, s' ∈ {0, 1}`, `Lp ∈ {1, P}` and `Rp ∈ {1, Q}`. Each such corner is `aPQ + bP + cQ + d` with `a = σ/2`, and this is `>= 0` on `[1, ∞)²` iff `a >= 0`, `a + b >= 0`, `a + c >= 0` and `a + b + c + d >= 0`.

*(v) Finite check.* `z1_cone_certificate.py` verifies all these conditions exactly for:
- every strong tournament `D` with `|D| <= 7`: the singleton, and the 1, 1, 6, 35 and 353 classes of sizes 3–7;
- every partition `D = O ⊔ I`;
- every `v ∈ D`.

That is 47,642 block types and every `v`, with no failure (`z1_cert_upto5.out`, `z1_cert_6.out`, `z1_cert_7.out`). Audit E reproduced the certificate independently, enumerating permutations and cut points in exact integers (`audit_E_20261014/part2/z1_cert_audit.c`, 0 failures). It also brute-forced 3,000 random cut-vertex tournaments up to `n = 16` (11,630 interior vertices, no violation, minimum exactly `1/3`). ∎

*Remark (what is left in general).*
- Of the four bilinear conditions, three hold for every block by short comparisons with the complement-symmetric `σ`; for instance `β − β(v) <= σ/2` and `α − α(v) <= σ/2 − 1`.
- Only the corner `a + b + c + d >= 0` (`L4` for the lower bound, `U4` for the upper) needs the computer.
- Its minimum slack is 0 exactly for singleton blocks (realised in the equality-word tournament), and 5/3, 2, 26/9, 10/3, 298/69 for block sizes 3, 4, 5, 6, 7 (`z1_local_slack.out`; confirmed by audit E).
- A proof of `L4` and `U4` for all strong tournaments would prove HYP-9279 (Z1) in full, and hence HYP-9277 for every strong tournament with a cut vertex whose deletion leaves at least three strong components.

## Reading

- **The balanced pair is special.** In the transitive-plus-one family, pairs among the `u`'s are typically unbalanced: for `w = 1^a 0^b` with `a, b` large, most pairs `i < j` have `p(u_i before u_j) ≈ 3/4`, since a wide crossing makes both free. The balance lives on the pairs `(z, ·)`. So averaging over all pairs cannot prove the conjecture. Statement 2 quantifies this for closable paths: rotations alone cannot sustain a `> 2/3` majority order, but non-closable paths can.
- **Why `u_2`?** The crossed splits starting at `β = 1` make `u_2` a fair coin, unless `u_2` is the very first `0`. The geometric weights `2^(α−β−1)` concentrate the mass on wide crossings, where `u_2` is free. Statement 1(f) is the same mechanism for every interior `u_i`: the widest crossings `(1, α)` with `α > i` carry at least half of all weight, and inside them `u_i` is a fair coin against `z`.
- **Towards the general strong case.** The minimisers for `n <= 9` all have small feedback arc sets. A natural induction is "transitive plus `k`", tracking how the balanced pair moves as vertices are added. The `n = 8` minimiser is transitive plus 2 (deleting two vertices leaves a transitive tournament); the `n = 9` minimisers are transitive plus at most 3.
- **What does generalise.** Statement 5 carries the transitive-plus-one mechanism to every cut vertex. Interior blocks are fair coins inside the crossing window, and wide windows carry geometric weight, since `mid → mid >= 2 h(D)`.
  - The natural general statement is HYP-9279 (Z): a cut vertex is balanced against every vertex of `D_2 ∪ … ∪ D_(k−1)`, and always has a balanced partner.
  - 1(f) is its case of singleton blocks.
  - (Z) would prove HYP-9277 for every strong tournament with a cut vertex. At `n = 10` these are 2,820,361 of the 9,355,949 strong tournaments.
- **What does not generalise.** Statement 1(e) suggests that every Hamiltonian cycle of a strong tournament carries a balanced arc. That holds for all 133,632 Hamiltonian cycles of strong 8-tournaments, but it is false at `n = 9`: see HYP-9279, refuted variant (H). The vertex-wise form 1(g) survives exhaustively to `n = 10`.
