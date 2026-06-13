---
id: THM-322
title: Odd-Cycle Count Formulas for All-0 Staircase T_k
status: CONJECTURE (verified k=2..8)
session: opus-2026-05-23-S5
tags: [staircase, cycle-counts, combinatorics, formula]
related: [THM-318, THM-319, THM-321]
---

## Statement

In the all-0 staircase tournament T_k (n=2k vertices), the number of odd directed cycles of length 2j+1 is:

$$\#_{2j+1}(T_k) = \sum_{i=0}^{j-1} c_{j,i} \cdot \binom{k}{2j-i}$$

with the following properties and explicit formulas for j=1..4:

### j=1: Triangle (#3)
$$\#_3(T_k) = 2\binom{k}{2}$$

### j=2: Pentagon (#5)
$$\#_5(T_k) = 4\binom{k}{4} + 6\binom{k}{3} = \binom{k}{3}(k+3)$$

### j=3: Heptagon (#7)
$$\#_7(T_k) = 6\binom{k}{6} + 80\binom{k}{5} + 28\binom{k}{4}$$

### j=4: Nonagon (#9)
$$\#_9(T_k) = 8\binom{k}{8} + 672\binom{k}{7} + 1220\binom{k}{6} + 210\binom{k}{5}$$

## Coefficient Patterns

**Leading coefficient:** $c_{j,0} = 2j$ (coefficient of $\binom{k}{2j}$)

**Number of terms:** j terms for the (2j+1)-cycle formula

**Diagonal (last) coefficient:** $c_{j,j-1} = d_j$ where $d_j = \#_{2j+1}(T_{j+1})$, the number of (2j+1)-cycles in the minimum-k staircase:

| j | d_j | T_{j+1} |
|---|-----|---------|
| 1 | 2   | T_2 (4 vertices) |
| 2 | 6   | T_3 (6 vertices) |
| 3 | 28  | T_4 (8 vertices) |
| 4 | 210 | T_5 (10 vertices) |
| 5 | 2154| T_6 (12 vertices) |

## Vanishing Condition

$\#_{2j+1}(T_k) = 0$ for $k < j+1$ (a (2j+1)-cycle needs at least 2j+1 vertices, requiring $2k \geq 2j+1$, i.e., $k \geq j+1$). This is automatic since all basis elements $\binom{k}{r}=0$ for $r > k$.

## Verification

All four formulas verified at k=2..8 by exhaustive cycle enumeration (hyp1732_full_investigation.py + cycle_count_formula.py).

Data table:

| k | #3 | #5 | #7 | #9 |
|---|----|----|----|----|
| 2 | 2  | 0  | 0  | 0  |
| 3 | 6  | 6  | 0  | 0  |
| 4 | 12 | 28 | 28 | 0  |
| 5 | 20 | 80 | 220| 210|
| 6 | 30 | 180| 906| 2480|
| 7 | 42 | 350|2702|13622|
| 8 | 56 | 616|6608|51304|

## Notes

- The formula for #5 has the nice factored form $\binom{k}{3}(k+3)$
- The leading term 2j·C(k,2j) matches the structure of the α_m formula (THM-321) at leading order
- The diagonal sequence 2, 6, 28, 210, 2154, ... is an empirical integer sequence; no closed form found as of session S5
- See output file: 05-knowledge/results/cycle_count_formula.out and cycle_counts_extended.out
