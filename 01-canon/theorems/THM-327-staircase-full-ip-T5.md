---
id: THM-327
title: Full Independence Polynomial of T_5
status: VERIFIED
session: opus-2026-05-27-S6
tags: [staircase, independence-polynomial, real-rooted]
related: [THM-318, THM-319, THM-323, THM-325, THM-326]
---

## Statement

The conflict graph $\Omega(T_5)$ of the all-0 staircase $T_5$ ($n=10$) has independence polynomial:

$$I(\Omega(T_5), x) = 1 + 530x + 317x^2 + 20x^3$$

## Structure

| Cycle length | Count |
|---|---|
| 3-cycles | 20 (= 2·C(5,2) from THM-322) |
| 5-cycles | 80 (= C(5,3)·(5+3) from THM-322) |
| 7-cycles | 220 |
| 9-cycles | 210 |
| Total odd cycles | 530 |

Independence polynomial: $\alpha_0=1$, $\alpha_1=530$, $\alpha_2=317$, $\alpha_3=20$.

## Key Properties

- **Degree:** $d(\Omega(T_5)) = 3$ (consistent with THM-318: $\lfloor 2 \cdot 5/3 \rfloor = 3$)
- **$H(T_5) = 2489$**: evaluating $I$ at $x=2$ gives $1 + 1060 + 1268 + 160 = 2489$
- $H(T_5) = 2489 \ll a(10) = 15745$: staircase is not optimal for $n=10$

## Comparison across staircases

| k | n | H(T_k) | a(n) | ratio a/H |
|---|---|---------|------|-----------|
| 2 | 4  | 5    | 5     | 1.00 |
| 3 | 6  | 29   | 45    | 1.55 |
| 4 | 8  | 233  | 661   | 2.84 |
| 5 | 10 | 2489 | 15745 | 6.32 |
| 6 | 12 | 33773| ≥531205| ≥15.73 |

The gap between $H(T_k)$ and $a(2k)$ grows rapidly.
