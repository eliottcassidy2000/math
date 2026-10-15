---
id: HYP-9279
title: "Strong tournaments are strictly better balanced: every strong tournament on n >= 4 vertices has a pair ordered with probability in [4/11, 7/11] in a uniform Hamiltonian path (delta(T) >= 4/11), the bound attained exactly at n = 5, 6, 7"
status: >
  OPEN. FINITE-EXACT for all strong tournaments with n <= 9 (exhaustive; minima 2/5, 4/11, 4/11, 4/11, 62/157, 62/157 for
  n = 4..9); the transitive-plus-one subfamily satisfies it for n <= 14 (exact; minima rising towards 1/2). Implies
  HYP-9277's inequality with a strict gap for strong components on >= 4 vertices.
session: mac-mini-2026-10-14-tourn
source: 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md
scripts:
  - 04-computation/experiments/tournament_balance_20261014/tbal.c (+ tbal_strong9.out)
  - 04-computation/experiments/tournament_balance_20261014/trans_plus_one.py
related: [HYP-9277, THM-4655]
---

# HYP-9279: a strict balance gap for strong tournaments

## Statement

Let `T` be strong with `n >= 4` vertices. In a uniformly random Hamiltonian path, some pair `x, y` satisfies

`4/11 <= p_T(x before y) <= 7/11`.

Equivalently, `δ(T) >= 4/11`.

- Equality would be attained exactly at the `n = 5, 6, 7` minimisers. At `n = 5` the minimiser is the transitive-plus-one tournament `w = 1100`.
- With S35's component reduction, this would imply HYP-9277's inequality for every tournament. Equality `δ = 1/3` then occurs only at strong components of size 3.

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
