---
id: HYP-9279
title: "Strong tournaments are strictly and vertex-wise balanced: (A) [= HYP-9277(3), opus S35] every strong tournament on n >= 4 vertices has delta(T) >= 4/11, here with evidence to n = 10 and glued-block constructions to n = 15; (B) in every strong tournament EVERY vertex v has a partner w with 1/3 <= P(v before w) <= 2/3 in a uniform Hamiltonian path; (B+) the graph of balanced pairs is connected; (B_loc) every vertex lies on a forward balanced arc of a directed 3-cycle; (Z) a cut vertex z is balanced against every vertex strictly between the first and last strong components of T - z. Refuted variant (H): 'every Hamiltonian cycle carries a balanced arc' fails at n = 9"
status: >
  OPEN. Part (A) is HYP-9277 statement 3 (opus S35, evidence n <= 8), restated here with new evidence; (B), (B+), (B_loc), (Z)
  and the refutation of (H) are this file's. (A) and (B) are FINITE-EXACT for all strong tournaments with n <= 10 (exhaustive:
  1 + 6 + 35 + 353 + 6008 + 178133 + 9355949 tournaments for n = 4..10; minima of delta 2/5, 4/11, 4/11, 4/11, 62/157,
  62/157, 62/157); beyond n = 10 two extremal blocks glued at a cut vertex give delta 592/1521 (n = 11) and 1052/2721
  (n = 12, 13), so the minimum of delta is not monotone in n but stays above 4/11 in every construction tried. (B+) is FINITE-EXACT
  for n <= 9. On the transitive-plus-one family, (B) and (B+) are PROVED (THM-4655 1(g)). (B) is sharp: the vertex-wise
  constant is at most 1/3 for every n >= 5 (PROVED, THM-4655 (6), family V_n) and equals 1/3 for 5 <= n <= 10. (B_loc) is FINITE-EXACT for n <= 10 and exact on the transitive-plus-one family
  to n = 14. (Z) is FINITE-EXACT for n <= 10, survives random and adversarial sampling of the exact cut-vertex transfer chain
  (THM-4655 (5)) to n = 46, is sharp at 1/3, and (Z1) is PROVED for every interior component with at most 7 vertices, the other
  components arbitrary (THM-4655 (7), computer-assisted; singleton blocks = THM-4655 1(f)). Each of (A), (B), (B+),
  (B_loc) implies HYP-9277's inequality, and (Z) implies it for every strong tournament with a cut vertex; (A) gives it with a
  strict gap for strong components on >= 4 vertices. Variant (H) is REFUTED (explicit 9-vertex counterexample). Audited
  (audit E, MISTAKE-615: (A) re-attributed, poset example, tight-vertex structure and the k = 2 floor of (Z) corrected; audit E
  part 2, MISTAKE-616: sharpness of (B) restated as an upper bound plus the exhaustive range; glued-block values reproduced).
session: mac-mini-2026-10-14-tourn
source: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
scripts:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out)
  - 04-computation/experiments/tournament_balance_20261014/tbal2.c (+ tbal2_strong9.out, tbal2_strong10.out)
  - 04-computation/experiments/tournament_balance_20261014/tbal3.c, tbal4.c (+ tbal34_upto8.out), tbal5.c (+ tbal5_strong9.out), c_counterexample.py (+ .out)
  - 04-computation/experiments/tournament_balance_20261014/tbal6.c (+ tbal6_upto9.out, tbal6_strong10.out), tbal7.c (+ tbal7_upto8.out, tbal7_strong9.out, tbal7_strong10.out)
  - 04-computation/experiments/tournament_balance_20261014/cutvertex_chain.py, cutvertex_search.py (+ .out files)
  - 04-computation/experiments/tournament_balance_20261014/glue_blocks.py, glue_blocks_run.py (+ glue_blocks.out), glue_minimisers_check.out, check_minimiser.py, tbal_extend.c, tbal_near.c (+ tower.log), tbal_anneal.c (+ tbal_anneal.out)
  - 04-computation/experiments/audit_E_20261014/ (independent audit code and report)
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py, tpo_vertexwise.py, tpo_cycle_arcs.py, tpo_forward_balanced.py
related: [HYP-9277, THM-4655, HYP-9280 (opus S36, majority-cycle conjecture)]
---

# HYP-9279: strict and vertex-wise balance for strong tournaments

## Statement (A) = HYP-9277(3) (opus S35)

This is S35's HYP-9277 statement 3, restated with new evidence: the exhaustive `n = 9, 10` sweeps, and the glued-block constructions for `n = 11..13` below.

Let `T` be strong with `n >= 4` vertices. In a uniformly random Hamiltonian path, some pair `x, y` satisfies

`4/11 <= p_T(x before y) <= 7/11`.

Equivalently, `δ(T) >= 4/11`.

- New conjecture here: equality is attained exactly at the `n = 5, 6, 7` minimisers (one, two and one classes). At `n = 5` the minimiser is the transitive-plus-one tournament `w = 1100`.
- With S35's component reduction, this would imply HYP-9277's inequality for every tournament. Equality `δ = 1/3` then occurs only at strong components of size 3.

## Statement (B): vertex-wise balance

In a strong tournament on `n >= 3` vertices, every vertex `v` has a partner `w` with `1/3 <= p_T(v before w) <= 2/3`.

- **Strength.** This is much stronger than HYP-9277, which asks for one balanced pair.
  - It has no poset analogue (corrected by audit E).
    - In the poset `3 + 1 = {x < u_1 < u_2} + {y}`, the element `x` is incomparable only to `y`, but `P(x before y) = 3/4`. So `x`, and likewise `u_2`, has no balanced partner.
    - Four elements are needed: no poset on at most 3 elements fails, and 48 of the 219 labelled 4-element posets do (`audit_E_20261014/poset_vertexwise_check.out`).
  - Strong connectivity is what puts every vertex on cycles.
  - It fails trivially for non-strong tournaments: a singleton strong component is ordered deterministically against everything.
- **Evidence.**
  - Exhaustive for all strong tournaments with `n <= 10` (`tbal2.c`; 9,355,949 tournaments at `n = 10`).
  - PROVED on the transitive-plus-one family (THM-4655 1(g)). There `z` is balanced against every interior `u_i`.
  - Also observed for `n <= 10`: every vertex whose deletion destroys strong connectivity has a balanced partner (2,820,361 such tournaments at `n = 10`). This is implied by (B).
- **Sharpness (THM-4655 (6)).** The vertex-wise constant `min_T min_v max_w min(p, 1 − p)` is at most `1/3` for every `n >= 5` (PROVED, family `V_n`), unlike the global `δ`. It equals `1/3` wherever (B) holds, so for `5 <= n <= 10`; the wording was corrected by audit E part 2.
  - The family is `V_n = TT_(n−3) ⇒ C_3` with one back arc `b → 0`. There the order on the triangle is uniform over its three rotations.
  - Exhaustively (`tbal5.c`), `1/3` is the minimum for `n = 5, …, 9`. It is attained by 1, 2, 5, 17 and 82 strong tournaments respectively.
  - Annealing (`tbal_anneal.c`, 6 seeds per `n`) never goes below `1/3`. It reaches `1/3` exactly at `n = 11, 12, 13` (6, 3 and 2 seeds respectively); the other runs stop just above.
  - **Structure of the tight vertices** (corrected by audit E, `tight_structure.out`): each of the 32 tight vertices with `n <= 8` lies on a cyclic triangle formed with two of its exact partners (`p = 2/3, 1/3`), as in an isolated `C_3`.
    - The triangle is attached by a single back arc only in the transitive-plus-one examples: the tournament with THM-4655's equality word `1 0 1^(m−3) 0` and its reverse `V_n`. These include every tight vertex at `n = 5, 6`.
    - At `n = 7`, 5 of the 7 tight vertices need two back arcs. At `n = 8`, 18 of 20 do, and the triangle can sit in the middle of the order.

## Statement (B⁺): the balanced graph is connected

In a strong tournament, the graph on `V` whose edges are the balanced pairs is connected.

- (B⁺) implies (B), since a connected graph on `n >= 2` vertices has no isolated vertex.
- Exhaustive for all strong tournaments with `n <= 9` (`tbal5.c`); PROVED on the transitive-plus-one family (THM-4655 1(g)).
- The balanced graph can have a vertex of degree 1 (142 of the 6,008 strong tournaments at `n = 8`; `tbal3.c`).

## Statement (B_loc): a forward balanced arc on a 3-cycle at every vertex

In a strong tournament, every vertex lies on an arc `x → y` with `1/2 < p(x before y) <= 2/3` that belongs to a directed 3-cycle.

- (B_loc) implies (B). It refines (B) in two ways: the partner shares a 3-cycle with the vertex, and the arc between them is *forward*, i.e. an arc of S36's digraph `F(T)` (HYP-9280).
- **Evidence.**
  - Exhaustive for all strong tournaments with `n <= 10` (`tbal6.c`; 9,355,949 at `n = 10`). Each half also holds separately: a balanced partner on a common 3-cycle (R_tri), and a forward balanced arc (R_maj).
  - Exact on the transitive-plus-one family to `n = 14` (`tpo_forward_balanced.py`).
- In the refuted variant (H) below, all nine arcs of one Hamiltonian cycle are forward with `p > 2/3`. (B_loc) asks for a forward arc with `p <= 2/3` at every vertex, and that survives.

## Statement (Z): cut vertices

Let `z` be a vertex of a strong tournament `T` with `T − z` not strong, and let `D_1 ⇒ … ⇒ D_k` be the strong components of `T − z`.
- **(Z1)** If `k >= 3`, then `1/3 <= p(v before z) <= 2/3` for every `v ∈ D_2 ∪ … ∪ D_(k−1)`.
- **(Z2)** `z` has a balanced partner. For `k >= 3` this follows from (Z1); for `k = 2` it is the content.

**Why it holds, and where it is PROVED.** In the exact transfer chain of THM-4655 (5), every interior block is a fair coin against `z` inside the crossing window. Wide windows carry geometric weight, since each interior block contributes `Σ_X h(X) h(D ∖ X) >= 2h(D)`.
- **PROVED for every interior component with at most 7 vertices** (THM-4655 (7); computer-assisted). The other components may be arbitrary, and so may `n`.
  - The normalised chain stays in an invariant region: left `(1, g, P, Lp)` with `P >= 1` and `1 <= Lp <= P`; right `(Rp, s', Q, 1)` with `Q >= 1` and `1 <= Rp <= Q`. The bounds `>= 1` come from injections that cut a Hamiltonian path of an end block at its first out-neighbour (or last in-neighbour) of `z`.
  - So (Z1) at `v` reduces to finitely many bilinear conditions on `v`'s own block. These are checked exactly for all strong tournaments with at most 7 vertices, all `z`-patterns and all `v`: 47,642 block types.
  - It includes singleton blocks (1(f)) and `C_3` blocks.
- **Reduction in general.** All but two corner conditions (`L4`, `U4`) hold for every block. Their slack is 0 only for singletons and grows with the block (5/3, 2, 26/9, 10/3, 298/69 for sizes 3–7). Proving `L4`, `U4` for all strong tournaments would prove (Z1) in full.

**Evidence.**
- **Exhaustive** for `n <= 10` (`tbal7.c`), with no failures of (Z1) or (Z2), and a balanced partner always in `D_1 ∪ D_2` or in `D_(k−1) ∪ D_k`:
  - `n = 9`: 99,695 cut-vertex instances, 13,060 with `k >= 3`;
  - `n = 10`: 3,274,940 instances, 250,690 with `k >= 3`.
- **Random sampling** of the exact chain (`cutvertex_chain.py`): 200 random cut-vertex tournaments with up to 9 blocks of sizes 1–5 (`n <= 46`) give no failure.
- **Adversarial hill-climbing** (`cutvertex_search.py`):
  - The minimum of (Z1) is exactly `1/3`, and never lower. It is attained by the transitive-plus-one word `1 0 1 1 0` and by block patterns such as sizes `[3, 1, 3, 1]`. Near-misses of the form `c/(3c − 1)` (for example `229376/688127`) approach `1/3` from above.
  - For `k = 2` the search stayed `>= 6/17`, but the exhaustive data go lower (audit E, `z2_k2_witness_check.out`).
    - The minima of `z`'s best-partner balance at `k = 2` cut vertices are 5/13, 6/17, 10/29, 18/53, 34/101 for `n = 5..9`.
    - For `n >= 6` they are `(2^(n−4) + 2)/(3·2^(n−4) + 5)`, decreasing towards `1/3`.
    - So (Z2) at `k = 2` is asymptotically sharp, and the adversarial search under-explores two-block configurations.

**Consequence.** (Z) would prove HYP-9277 for every strong tournament with a cut vertex. At `n = 10` these are 2,820,361 of the 9,355,949 strong tournaments.

**Route to a full proof.** THM-4655 (7) reduces (Z1) to two local inequalities `L4`, `U4` on one strong tournament with a 2-colouring (`N⁺(z)`, `N⁻(z)`) and a marked vertex.
- Generic chain constraints are not enough: an abstract chain with a left weight `0.0002` violates (Z1). The end-block injections `β, γ + s >= 1` are exactly what rule this out.
- The remaining gap is a combinatorial proof of `L4`, `U4` for blocks of every size. (Z2) at `k = 2`, where there are no interior blocks, needs a different argument.

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

## Beyond `n = 10`: glued extremal blocks

Take two of the small extremal tournaments as blocks `D_1 ⇒ D_2` and attach a cut vertex `z` in every possible way (`glue_blocks.py`, `glue_blocks_run.py`; `./tbal` on all attachment patterns).

The blocks are:
- `S4`, the strong 4-tournament (`δ = 2/5`);
- `M5`, `M6a`, `M6b`, `M7`, the minimisers with `δ = 4/11`;
- `M8`, with `62/157`.

| `n` | best construction | minimum `δ` | `H` |
|---|---|---|---|
| 11 | `M5 + M5` | `592/1521 = 0.38922` | 1521 |
| 12 | `M6a + M5` | `1052/2721 = 0.38662` | 2721 |
| 13 | `M6a + M6a`, `M6b + M6b`, `M5 + M7`, `M7 + M5` | `3156/8163 = 1052/2721` | 8163 |
| 14 | `M7 + M6a` (`M8 + M5`: 0.41535; `M6a + M7`: 0.40075) | `9468/24489 = 1052/2721` | 24489 |
| 15 | best of `M7 + M7` (0.40075), `M8 + M6a` (0.41535), `M6a + M8` (0.41249) | `18261/45567 = 0.40075` | 45567 |

- Both the `n = 11` and the `n = 12` minimiser were confirmed by an independent Python DP (`glue_minimisers_check.out`).
- The extension tower (`tbal_extend.c`, `tower.log`) first found an `n = 13` tournament with `δ = 592/1521` (`H = 13689 = 117²`). It is `M6b + M6a` glued at a cut vertex. The glue search then found `1052/2721` at `n = 13`.
- Deleting suitable vertices from it keeps `δ = 592/1521` at `n = 12`.
- So the minimum of `δ` is **not** monotone: 4/11 for `n = 5, 6, 7`, then 62/157 for `n = 8, 9, 10`, then at most 0.38922 at `n = 11` and at most 0.38662 at `n = 12, 13, 14`. The value `1052/2721` repeats across `n = 12..14`, a plateau like 62/157 at `n = 8..10`.
- Every construction stays above `4/11 = 0.36364`. A random-arc annealer (`tbal_anneal.c`) stalls far above, near 0.425 at `n = 11`.

## Status of the constant

- `4/11` is the observed infimum. It is attained only at `n = 5, 6, 7` among `n <= 10` (exhaustive), and in no construction for `n = 11..13`.
- The minimum for `n = 8, 9, 10` is `62/157 = 0.3949`. For `n = 11, 12, 13, 14` it is at most 0.3892, 0.3866, 0.3866 and 0.3866.
- Whether glued constructions with larger blocks push `δ` towards or below `4/11` is open. The 62/157 block `M8` is a worse gluing ingredient (0.41535 at `n = 14`). At `n = 15` the pairs tried give at least 0.40075 (`glue_blocks2.out`), so these two-block gluings did not continue the descent.
- No strong tournament found with `4 <= n <= 14` comes close to `1/3`.
