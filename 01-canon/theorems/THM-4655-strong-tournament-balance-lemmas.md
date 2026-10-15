---
id: THM-4655
title: "Uniform Hamiltonian paths of strong tournaments: the 1/3-2/3 property holds for every strong 'transitive-plus-one' tournament, with the explicit balanced pair (z, u_2); a rotation lemma for Hamiltonian cycles; Stanley-type log-concavity fails; exhaustive verification for all strong tournaments on n <= 9 vertices"
status: >
  PROVED: statements 1-2 (elementary). FINITE-EXACT: statements 3-4 (exhaustive, gentourng + exact subset DP in C).
  The tournament 1/3-2/3 conjecture HYP-9277 (opus S35) remains OPEN; this theorem proves it for a natural subclass of
  strong tournaments and extends its exhaustive verification from n <= 8 to n <= 9.
source: mac-mini-2026-10-14-tourn
depends_on: [HYP-9277 and its proved reductions (opus S35, one_twenty_eight_atlas_balance_gold_20261014.md: delta(T) = max over strong components; with no balanced pair the >2/3-majority relation is a linear order)]
verification:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out): exhaustive strong n <= 9
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py, tpo_proof_steps.py (+ .out), tpo_u2_check.c (+ .out)
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

   - **(c) The balanced pair.** If `m >= 3` (`n >= 4`), then

     `1/3 < p_T(u_2 before z) <= 2/3`,

     with equality iff `w = 1 0 1^(m−3) 0` (`m >= 4`). Hence `δ(T) >= 1/3`: **the 1/3–2/3 property holds for every strong transitive-plus-one tournament, with the explicit balanced pair `(z, u_2)`.** For `m = 2`, `T = C_3` and `δ = 1/3`.
   - **(d) Corollary (with S35's component reduction).** Suppose each strong component of a tournament becomes transitive after deleting one vertex. Then the tournament satisfies the 1/3–2/3 inequality.

2. **(Rotation lemma.)** Let `C` be a Hamiltonian cycle and `σ` a linear order of `V`. Cutting `C` before each of its `n` vertices gives `n` Hamiltonian paths (the closable paths of `C`).
   - Averaged over these paths and over all pairs, the fraction of pairs ordered as in `σ` is at most `(2n − 1)/(3n) < 2/3`, with equality iff `C` visits `V` in `σ`-order.
   - **Corollary.** If a strong tournament has no balanced pair, fewer than a fraction `n/(n+1)` of its Hamiltonian paths are closable (`last → first`).

3. **(No Stanley inequality; FINITE-EXACT.)** For uniform Hamiltonian paths, neither of the following distributions is log-concave in general:
   - the position distribution of a vertex;
   - the gap distribution of `pos(y) − pos(x)`.

   Failures appear at `n = 4`. All 35 strong tournaments at `n = 6` fail for positions, and every strong tournament at `n = 6, 7` fails for gaps. So the Kahn–Saks route to `δ >= 3/11` for posets, which runs through Stanley's inequality, does not transfer verbatim.

4. **(Exhaustive; FINITE-EXACT.)** Every strong tournament on `n <= 9` vertices (1, 1, 6, 35, 353, 6008, 178133 for `n = 3..9`) has `δ >= 1/3`. The minima are:

   | `n` | minimum `δ` |
   |---|---|
   | 3 | `1/3` |
   | 4 | `2/5` |
   | 5, 6, 7 | `4/11` |
   | 8, 9 | `62/157` |

   The minimisers are near-transitive (minimum feedback arc sets 1, 2, 3, 5, 4, 6). With S35's reduction, **HYP-9277's inequality holds for every tournament on `n <= 9` vertices** (S35 had `n <= 8`).

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

`tpo_proof_steps.py` checks every identity and inequality above on all 8,190 strong words with `m <= 14`. `tpo_u2_check.c` independently enumerates all valid `A` for every word with `m <= 20`: 524,286 words, no violation, and the extremal words are as stated. ∎

**1(d).** By S35's reduction, `δ` of a tournament is the maximum of `δ` over its strong components. Apply 1(c). ∎

**2.**
- For a pair `i < j` (in `σ`), `i` precedes `j` in exactly `(c_i − c_j) mod n` of the `n` cuts, where `c` gives the positions on `C`. So the concordance count is `n·C(n,2) − Σ_(i<j) y_ij`, with `y_ij = (c_j − c_i) mod n`.
- From a fixed `i`, the forward distances to the other vertices are exactly `1, …, n − 1`. So the `n − 1 − i` larger labels have distance sum at least `1 + … + (n − 1 − i)`.
- Summing over `i` gives `Σ_(i<j) y_ij >= (n³ − n)/6`, with equality iff each `i` is followed by `i+1, …, n−1` at distances `1, …, n−1−i`, i.e. the identity cycle. The average concordance is therefore at most `(2n − 1)/(3n)`.
- For the corollary: if no pair is balanced, every pair is concordant with `σ*` with probability `> 2/3`. Non-closable paths contribute at most 1, so `2/3 < f·(2n−1)/(3n) + (1 − f)`, where `f` is the closable fraction. ∎

**3, 4.** These are the exhaustive computations listed in `verification`.

## Reading

- **The balanced pair is special.** In the transitive-plus-one family, pairs among the `u`'s are typically unbalanced: for `w = 1^a 0^b` with `a, b` large, most pairs `i < j` have `p(u_i before u_j) ≈ 3/4`, since a wide crossing makes both free. The balance lives on the pairs `(z, ·)`. So averaging over all pairs cannot prove the conjecture. Statement 2 quantifies this for closable paths: rotations alone cannot sustain a `> 2/3` majority order, but non-closable paths can.
- **Why `u_2`?** The crossed splits starting at `β = 1` make `u_2` a fair coin, unless `u_2` is the very first `0`. The geometric weights `2^(α−β−1)` concentrate the mass on wide crossings, where `u_2` is free.
- **Towards the general strong case.** The minimisers for `n <= 9` all have small feedback arc sets. A natural induction is "transitive plus `k`", tracking how the balanced pair moves as vertices are added. The `n = 8` minimiser is transitive plus 2 (deleting two vertices leaves a transitive tournament); the `n = 9` minimisers are transitive plus at most 3.
