---
id: HYP-9277
title: "A tournament 1/3-2/3 theorem: in a uniformly random Hamiltonian path of a non-transitive tournament some pair is ordered with probability in [1/3, 2/3]; equality exactly at transitive compositions of singletons and cyclic triangles (A000930(n) - 1 classes: 128 at n = 14, the same count as the posets at balance 1/3 through 14 elements); least strong balance 4/11. The gold partition property fails once precedence answers are recorded."
status: >
  OPEN (the statement for all n).
  PROVED: the component reduction (delta(T) = max over strong components), the strong-majority lemma (a counterexample
  carries a linear order sigma* with every earlier-before-later probability > 2/3), C_3 has delta = 1/3, majority
  3-cycles force a balanced pair, and the equality-class count A000930(n) - 1.
  FINITE-EXACT: all 2^C(n,2) labelled tournaments for n <= 8 [N8]; equality-class counts 2, 16, 120, 1040, 10080 labelled (n = 3..7) as predicted;
  least balance of strong tournaments 1/3, 2/5, 4/11, 4/11, 4/11 (n = 3..7).
  REFUTED (extensions): with precedence constraints the minimum is 1/4 (tournaments n = 6, semicomplete digraphs n = 4),
  and Peczarski's t0 >= t1 + t2 fails (176 closed-family members at n = 6).
session: opus-2026-10-14-S35
source: 05-knowledge/results/one_twenty_eight_atlas_balance_gold_20261014.md
scripts:
  - 04-computation/experiments/tournament_one_third_balance_20261014.py
  - 04-computation/experiments/tournament_balance_labelled_sweep_20261014.py
  - 04-computation/experiments/tournament_balance_spectrum_gold_20261014.py
related: [HYP-2371, HYP-9273, THM-4527, THM-4526]
---

# HYP-9277: the 1/3–2/3 phenomenon for Hamiltonian paths of tournaments

**Definition.** For a tournament `T`, let `p(x, y)` be the probability that `x` precedes `y` in a uniformly random Hamiltonian path of `T`. Set `δ(T) = max_{x ≠ y} min(p, 1 − p)`. A transitive tournament has one Hamiltonian path and `δ = 0`.

## Statement

For every tournament `T` that is not transitive:
1. `δ(T) >= 1/3`;
2. `δ(T) = 1/3` iff every strong component of `T` is a single vertex or a cyclic triangle `C_3`;
3. every strong tournament on `n >= 4` vertices has `δ >= 4/11`.

## Proved parts

1. **Component reduction.** Every Hamiltonian path traverses the strong components in their (transitive) condensation order. So the uniform measure on Hamiltonian paths is the product of the uniform measures on the components' paths, and `δ(T)` is the maximum of `δ` over the strong components. Hence parts 1–2 reduce to strong tournaments: strong ones on `>= 4` vertices must have `δ > 1/3`.
2. **Strong-majority lemma** (any probability measure on linear orders, posets included). If no pair has `p ∈ [1/3, 2/3]`, then `x ≻ y ⟺ p(x, y) > 2/3` is a strict linear order.
   - It is complete, because every pair is unbalanced.
   - It is transitive, because `p(x, z) >= p(x, y) + p(y, z) − 1 > 1/3` forces `p(x, z) > 2/3`.
3. **`C_3`.** Its three Hamiltonian paths give `p = 2/3` on every arc.
4. **Majority 3-cycles.** In every linear order of three elements, exactly one or two of `x<y`, `y<z`, `z<x` hold. So `p(x,y) + p(y,z) + p(z,x) ∈ [1, 2]`, and a majority 3-cycle (all three `> 1/2`) contains a pair with `p <= 2/3`, i.e. a balanced pair.
5. **Equality classes.** By the uniqueness of the condensation, the tournaments in part 2 of the statement correspond, up to isomorphism, to compositions of `n` into parts 1 and 3 with at least one 3. There are `A000930(n) − 1` of them, Narayana's cows sequence `c(n) = c(n−1) + c(n−3)` minus one: 1, 2, 3, 5, 8, 12, 18, 27, 40, 59, 87, **128** for `n = 3 … 14`.

The same count holds for posets: `δ(P) = max` over ordinal summands, and the 3-element poset `2 + 1` plays `C_3`'s role. arXiv:2607.23926 (UNREFEREED) reports that exactly 128 classes of 14-element posets attain `1/3`, all of them ordinal sums of singletons and `2 + 1`. This is that count. For tournaments the 128 at `n = 14` is conditional on part 2 for `n <= 14`.

## Evidence (FINITE-EXACT)

**Labelled sweep** (forward/backward subset DP, numba). For every non-transitive labelled tournament with `n <= 7` [N8], `δ >= 1/3`. The labelled counts at `δ = 1/3` are 2, 16, 120, 1040, 10080, exactly `Σ n!/3^b` over the compositions. That confirms part 2 there, i.e. strong tournaments on 4–7 vertices have `δ > 1/3`.

**Strong tournaments.**
- The least `δ` is 2/5 at `n = 4` (the unique strong 4-tournament, `H = 5`), and 4/11 at `n = 5, 6, 7`. The 5-vertex extremal tournament has score sequence `(1,2,2,2,3)` and `H = 11`.
- The `δ` values at `n = 5` are 4/11, 5/13, 2/5, 4/9.
- At `n = 6` the least eight are 4/11, 17/45, 2/5, 15/37, 12/29, 18/43, 13/31, 8/19.
- Random strong tournaments drift toward `1/2` (least `δ` seen: 53/135 at `n = 7`, 118/279 at `n = 9`). Every sampled strong tournament has a balanced **arc**.

**Paley.** `p(0, 1) = 2/3, 118/189, 55553/95095, 645198853226/1172695746915` at `q = 3, 7, 11, 19`, decreasing toward `1/2`. By arc-transitivity this is the same on every arc.

## The posets' parallel, and where it stops

| | posets (1/3–2/3; arXiv:2607.23926 through 14 elements) | tournaments (here, `n <= 7`/8) |
|---|---|---|
| extremal block | `2 + 1` (3 linear extensions) | `C_3` (3 Hamiltonian paths) |
| reduction | ordinal sums | transitive composition (condensation) |
| classes at `1/3` | ordinal sums of `1` and `2+1`: 128 at order 14 | compositions of `1` and `C_3`: `A000930(n) − 1` (128 at `n = 14`) |
| least balance above `1/3` | 37/106 | 4/11 |
| least balance without a nontrivial sum | 254/725 | 4/11 |

**Sorting fails.** Peczarski's gold partition conjecture says two consecutive comparisons can be chosen with `t_0 >= t_1 + t_2` (`t_0` linear extensions, `t_1` after the first answer, `t_2` after both). That Fibonacci inequality gives sorting in `log_φ e(P)` comparisons. Its natural tournament analogue fails.
- Take the closed family `{HP(T) ∩ Ext(Q)}`: Hamiltonian paths consistent with recorded precedence answers.
- For tournaments it holds for `n <= 5` (792 members). At `n = 6` the minimum of `δ` over the family is **1/4**, and 176 members have no gold pair.
- Explicit witness: the code `(0,0,0,0,1,0,0,0,0,0,1,0,0,0,1)` tournament with the 4 surviving paths `(0,1,5,4,2,3)`, `(4,2,5,0,1,3)`, `(5,0,1,3,4,2)`, `(5,0,1,4,2,3)`, in which every pair splits 1–3 or 0–4.
- Semicomplete digraphs (2-cycles allowed, interpolating to the complete digraph = posets) fail already at `n = 4`.

So whatever proves 1/3–2/3 and the gold partition conjecture for posets must use what linear-extension sets have and Hamiltonian-path sets lack: convexity of the order polytope (Stanley, Kahn–Saks, Brightwell–Felsner–Trotter). Part 1 for tournaments, by contrast, holds without precedence constraints.

## Open

- Prove part 1. The strong-majority lemma reduces it to showing that a strong tournament admits no linear order `σ*` with all earlier-before-later probabilities `> 2/3`.
- Test `n = 9, 10` by targeted search (hill-climbing on `δ`).
- Is the least strong balance `4/11` for all `n >= 5`?
- Is there always a balanced **arc** in a strong tournament?
