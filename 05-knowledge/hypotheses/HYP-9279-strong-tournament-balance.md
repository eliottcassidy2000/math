---
id: HYP-9279
title: "Strong tournaments are strictly and vertex-wise balanced: (A) every strong tournament on n >= 4 vertices has delta(T) >= 4/11; (B) in every strong tournament EVERY vertex v has a partner w with 1/3 <= P(v before w) <= 2/3 in a uniform Hamiltonian path; (B+) the graph of balanced pairs is connected; (B_loc) every vertex lies on a forward balanced arc of a directed 3-cycle; (Z) a cut vertex z is balanced against every vertex strictly between the first and last strong components of T - z. Refuted variant (H): 'every Hamiltonian cycle carries a balanced arc' fails at n = 9"
status: >
  OPEN. (A) and (B) are FINITE-EXACT for all strong tournaments with n <= 10 (exhaustive: 1 + 6 + 35 + 353 + 6008 + 178133
  + 9355949 tournaments for n = 4..10; minima of delta 2/5, 4/11, 4/11, 4/11, 62/157, 62/157, 62/157). (B+) is FINITE-EXACT
  for n <= 9. On the transitive-plus-one family, (B) and (B+) are PROVED (THM-4655 1(g)). (B) is sharp: the vertex-wise
  constant is exactly 1/3 for every n = 5..9. (B_loc) is FINITE-EXACT for n <= 9 and exact on the transitive-plus-one family
  to n = 14. (Z) is FINITE-EXACT for n <= 9, survives random and adversarial sampling of the exact cut-vertex transfer chain
  (THM-4655 (5)) to n = 46, is sharp at 1/3, and is PROVED for singleton blocks (THM-4655 1(f)). Each of (A), (B), (B+),
  (B_loc) implies HYP-9277's inequality, and (Z) implies it for every strong tournament with a cut vertex; (A) gives it with a
  strict gap for strong components on >= 4 vertices. Variant (H) is REFUTED (explicit 9-vertex counterexample).
session: mac-mini-2026-10-14-tourn
source: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
scripts:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out)
  - 04-computation/experiments/tournament_balance_20261014/tbal2.c (+ tbal2_strong9.out, tbal2_strong10.out)
  - 04-computation/experiments/tournament_balance_20261014/tbal3.c, tbal4.c (+ tbal34_upto8.out), tbal5.c (+ tbal5_strong9.out), c_counterexample.py (+ .out)
  - 04-computation/experiments/tournament_balance_20261014/tbal6.c (+ tbal6_upto9.out), tbal7.c (+ tbal7_upto8.out, tbal7_strong9.out)
  - 04-computation/experiments/tournament_balance_20261014/cutvertex_chain.py, cutvertex_search.py (+ .out files)
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py, tpo_vertexwise.py, tpo_cycle_arcs.py, tpo_forward_balanced.py
related: [HYP-9277, THM-4655, HYP-9280 (opus S36, majority-cycle conjecture)]
---

# HYP-9279: strict and vertex-wise balance for strong tournaments

## Statement

Let `T` be strong with `n >= 4` vertices. In a uniformly random Hamiltonian path, some pair `x, y` satisfies

`4/11 <= p_T(x before y) <= 7/11`.

Equivalently, `δ(T) >= 4/11`.

- Equality would be attained exactly at the `n = 5, 6, 7` minimisers. At `n = 5` the minimiser is the transitive-plus-one tournament `w = 1100`.
- With S35's component reduction, this would imply HYP-9277's inequality for every tournament. Equality `δ = 1/3` then occurs only at strong components of size 3.

## Statement (B): vertex-wise balance

In a strong tournament on `n >= 3` vertices, every vertex `v` has a partner `w` with `1/3 <= p_T(v before w) <= 2/3`.

- **Strength.** This is much stronger than HYP-9277, which asks for one balanced pair.
  - It has no poset analogue. In the 3-element poset `{a} ∥ {b < c}`, the element `a` has balanced partners (`p = 1/3, 2/3`), but `b` and `c` are always ordered against each other.
  - Strong connectivity is what puts every vertex on cycles.
  - It fails trivially for non-strong tournaments: a singleton strong component is ordered deterministically against everything.
- **Evidence.**
  - Exhaustive for all strong tournaments with `n <= 10` (`tbal2.c`; 9,355,949 tournaments at `n = 10`).
  - PROVED on the transitive-plus-one family (THM-4655 1(g)). There `z` is balanced against every interior `u_i`.
  - Also observed for `n <= 10`: every vertex whose deletion destroys strong connectivity has a balanced partner (2,820,361 such tournaments at `n = 10`). This is implied by (B).
- **Sharpness.** The vertex-wise constant `min_v max_w min(p, 1 − p)` is exactly `1/3` for every `n = 5, …, 9` (`tbal5.c`), unlike the global `δ`. It is attained by 1, 2, 5, 17 and 82 strong tournaments respectively.
  - The tight vertices sit in a cyclic triangle hung at the end of the order by one back arc. Example at `n = 6`: `TT_3` followed by `3 → 4 → 5 → 3`, with the back arc `4 → 0`.
  - The tight vertex `5` has `p = 2/3, 1/3` against `3, 4`, exactly as in an isolated `C_3`.
  - So the equality cases of (B) are local `C_3` blocks, mirroring the global equality cases of HYP-9277.

## Statement (B⁺): the balanced graph is connected

In a strong tournament, the graph on `V` whose edges are the balanced pairs is connected.

- (B⁺) implies (B), since a connected graph on `n >= 2` vertices has no isolated vertex.
- Exhaustive for all strong tournaments with `n <= 9` (`tbal5.c`); PROVED on the transitive-plus-one family (THM-4655 1(g)).
- The balanced graph can have a vertex of degree 1 (142 of the 6,008 strong tournaments at `n = 8`; `tbal3.c`).

## Statement (B_loc): a forward balanced arc on a 3-cycle at every vertex

In a strong tournament, every vertex lies on an arc `x → y` with `1/2 < p(x before y) <= 2/3` that belongs to a directed 3-cycle.

- (B_loc) implies (B). It refines (B) in two ways: the partner shares a 3-cycle with the vertex, and the arc between them is *forward*, i.e. an arc of S36's digraph `F(T)` (HYP-9280).
- **Evidence.**
  - Exhaustive for all strong tournaments with `n <= 9` (`tbal6.c`). Each half also holds separately: a balanced partner on a common 3-cycle (R_tri), and a forward balanced arc (R_maj).
  - Exact on the transitive-plus-one family to `n = 14` (`tpo_forward_balanced.py`).
- In the refuted variant (H) below, all nine arcs of one Hamiltonian cycle are forward with `p > 2/3`. (B_loc) asks for a forward arc with `p <= 2/3` at every vertex, and that survives.

## Statement (Z): cut vertices

Let `z` be a vertex of a strong tournament `T` with `T − z` not strong, and let `D_1 ⇒ … ⇒ D_k` be the strong components of `T − z`.
- **(Z1)** If `k >= 3`, then `1/3 <= p(v before z) <= 2/3` for every `v ∈ D_2 ∪ … ∪ D_(k−1)`.
- **(Z2)** `z` has a balanced partner. For `k >= 3` this follows from (Z1); for `k = 2` it is the content.

**Why it should hold.** In the exact transfer chain of THM-4655 (5), every interior block is a fair coin against `z` inside the crossing window. Wide windows carry geometric weight, since each interior block contributes `Σ_X h(X) h(D ∖ X) >= 2h(D)`.
- For singleton blocks, `T − z` is transitive, and (Z1) is PROVED (THM-4655 1(f)).

**Evidence.**
- **Exhaustive** for `n <= 9` (`tbal7.c`): 99,695 cut-vertex instances at `n = 9`, of which 13,060 have `k >= 3`. There are no failures of (Z1) or (Z2). Moreover a balanced partner always lies in `D_1 ∪ D_2` or in `D_(k−1) ∪ D_k`.
- **Random sampling** of the exact chain (`cutvertex_chain.py`): 200 random cut-vertex tournaments with up to 9 blocks of sizes 1–5 (`n <= 46`) give no failure.
- **Adversarial hill-climbing** (`cutvertex_search.py`):
  - The minimum of (Z1) is exactly `1/3`, and never lower. It is attained by the transitive-plus-one word `1 0 1 1 0` and by block patterns such as sizes `[3, 1, 3, 1]`. Near-misses of the form `c/(3c − 1)` (for example `229376/688127`) approach `1/3` from above.
  - For `k = 2`, the best-partner balance of `z` stayed `>= 6/17`.

**Consequence.** (Z) would prove HYP-9277 for every strong tournament with a cut vertex. At `n = 10` these are 2,820,361 of the 9,355,949 strong tournaments.

**Route to a proof.** The chain reduces (Z1) to inequalities between the block parameters: the `M_j` entries normalised by `h(D_j)`.
- THM-4655 1(f) is the case where every `mid → mid` entry equals 2.
- A crude geometric bound handles blocks far from both ends.
- The delicate case is `i = 2` (and dually `i = k − 1`). There a small left weight `Σ_(X ≠ D_1) h(X) h⁺(D_1 ∖ X)/h(D_1)` must be offset, as Lemmas A and B do for singleton blocks.

## Refuted variant (H): balanced arcs on Hamiltonian cycles

**Conjecture (H):** every Hamiltonian cycle of a strong tournament carries a balanced arc, i.e. a pair `(c_i, c_(i+1))` with `1/3 <= p <= 2/3`. It would imply HYP-9277 at once, because strong tournaments are Hamiltonian (Camion).

- **True for `n <= 8`.** All 133,632 Hamiltonian cycles of the 6,008 strong 8-tournaments carry one. The minimum number of balanced arcs per cycle is exactly 1, attained for 8 tournaments.
- **True on the transitive-plus-one family**, with at least two balanced arcs (THM-4655 1(e)).
- **Shorter cycles fail earlier:** some 4-, 5- and 6-cycles have no balanced arc from `n = 7`, and some 3-cycles from `n = 8` (`tbal4.c`).
- **FALSE at `n = 9`.** Exactly one of the 178,133 strong 9-tournaments has a Hamiltonian cycle with no balanced arc. It is `110111111101001111111111111111111101` (gentourng ascii).
  - Its arcs agree with the order `0 < 1 < … < 8`, except the backward arcs `3 → 0`, `4 → 1`, `6 → 1`, `7 → 1`, `8 → 6`. Its scores are 7, 4, 6, 6, 5, 3, 2, 2, 1, and `H = 213`.
  - It has 4 Hamiltonian cycles. On the cycle `0 → 2 → 4 → 5 → 7 → 8 → 6 → 1 → 3 → 0`, the arc probabilities `p(c_i before c_(i+1))` are 143/213, 51/71, 59/71, 51/71, 143/213, 143/213, 48/71, 48/71, 143/213, all in `(0.671, 0.831)`.
  - So this cycle is a **`> 2/3` majority cycle**: on its arcs, the uniform path looks like a random rotation of the cycle.
  - The tournament itself is well balanced: `δ = 100/213`, and every vertex has a partner with balance at least `28/71`. It is no threat to HYP-9277 or to (B).
  - Confirmed independently by plain DFS enumeration of all 213 Hamiltonian paths (`c_counterexample.py`).
- **Lesson.** Hamiltonian cycles can carry a majority orientation, whereas vertices and pairs apparently cannot avoid balance. This is consistent with S36's HYP-9280, whose majority cycles live in the forward-arc digraph.

## Evidence

**Exhaustive (gentourng `-c` with an exact subset DP).** Every strong tournament with `n <= 10`. The minimum `δ` by `n`:

| `n` | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| number of strong tournaments | 1 | 6 | 35 | 353 | 6,008 | 178,133 | 9,355,949 |
| minimum `δ` | `2/5` | `4/11` | `4/11` | `4/11` | `62/157` | `62/157` | `62/157` |

The `n = 10` minimiser is `111000111101111111111111111111110001111111101` (gentourng ascii). The minimum over arcs `x → y` of `p(x before y)` falls from `17/81` at `n = 9` to `25/153` at `n = 10`: single arcs can be strongly reversed, while some pair stays balanced.

**The transitive-plus-one family.** Every strong word up to `m = 13`, i.e. `n = 14`. The family minima are `2/5, 4/11, 8/19, 2/5, 4/9, 4/9, 113/241, 8/17, 116/241, 16/33, 400/817`.

**Interpretation.** The minimisers are near-transitive, which is where balance is hardest. Away from transitivity the rotation structure of Hamiltonian cycles helps (THM-4655 (2)).

## Status of the constant

- `4/11` is the observed infimum for `n <= 10`. It is attained only at `n = 5, 6, 7`; the minimum is `62/157 = 0.3949` at `n = 8, 9, 10`.
- Whether the infimum over all strong tournaments is `4/11`, or whether larger `n` push it lower, is open.
- No strong tournament with `4 <= n <= 10` comes close to `1/3`.
