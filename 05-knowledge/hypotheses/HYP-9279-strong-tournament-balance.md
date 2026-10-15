---
id: HYP-9279
title: "Strong tournaments are strictly and vertex-wise balanced: (A) every strong tournament on n >= 4 vertices has delta(T) >= 4/11; (B) in every strong tournament EVERY vertex v has a partner w with 1/3 <= P(v before w) <= 2/3 in a uniform Hamiltonian path"
status: >
  OPEN. (A) and (B) are FINITE-EXACT for all strong tournaments with n <= 9 (exhaustive; minima 2/5, 4/11, 4/11, 4/11,
  62/157, 62/157 for n = 4..9; every vertex of all 1 + 6 + 35 + 353 + 6008 + 178133 strong tournaments has a balanced
  partner), and for the transitive-plus-one subfamily up to n = 14 (exact). Either implies HYP-9277's inequality; (A) with a
  strict gap for strong components on >= 4 vertices.
session: mac-mini-2026-10-14-tourn
source: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
scripts:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out)
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py
related: [HYP-9277, THM-4655]
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

- This is much stronger than HYP-9277, which asks for one balanced pair. It has no poset analogue: in the 3-element poset `{a} ∥ {b < c}`, the element `a` has balanced partners, but neither does `b` against `c`. Strong connectivity is what puts every vertex on cycles.
- It fails for non-strong tournaments trivially: a singleton strong component is ordered deterministically against everything.
- Checked exhaustively for all strong tournaments with `n <= 9` (`tbal2.c`), and for every strong transitive-plus-one tournament with `n <= 14` (`tpo_vertexwise.py`).
- For the extra vertex `z` of a transitive-plus-one tournament, (B) is PROVED by THM-4655: `u_2` is always a balanced partner.
- Also observed for `n <= 9`: every vertex whose deletion destroys strong connectivity has a balanced partner. This is implied by (B).

## Evidence

**Exhaustive (gentourng `-c` with an exact subset DP).** Every strong tournament with `n <= 9`. The minimum `δ` by `n`:

| `n` | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| number of strong tournaments | 1 | 6 | 35 | 353 | 6,008 | 178,133 |
| minimum `δ` | `2/5` | `4/11` | `4/11` | `4/11` | `62/157` | `62/157` |

(`n = 10`, 9,355,949 strong tournaments, is running.)

**The transitive-plus-one family.** Every strong word up to `m = 13`, i.e. `n = 14`. The family minima are `2/5, 4/11, 8/19, 2/5, 4/9, 4/9, 113/241, 8/17, 116/241, 16/33, 400/817`.

**Interpretation.** The minimisers are near-transitive, which is where balance is hardest. Away from transitivity the rotation structure of Hamiltonian cycles helps (THM-4655 (2)).

## Status of the constant

- `4/11` is the observed infimum for `n <= 9`.
- Whether the infimum over all strong tournaments is `4/11`, or whether larger `n` push it lower, is open.
- No strong tournament with `4 <= n <= 9` comes close to `1/3`.
