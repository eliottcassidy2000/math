---
id: THM-328
title: Full Independence Polynomial of T_6
status: VERIFIED
session: opus-2026-05-27-S6
tags: [staircase, independence-polynomial]
related: [THM-318, THM-319, THM-325, THM-326, THM-327]
---

## Statement

The conflict graph $\Omega(T_6)$ of the all-0 staircase $T_6$ ($n=12$) has independence polynomial:

$$I(\Omega(T_6), x) = 1 + 5750x + 4244x^2 + 642x^3 + 10x^4$$

## Structure

Total odd cycles = 5750 (from 3-, 5-, 7-, 9-, 11-cycles in $T_6$).

Independence number: $\alpha(\Omega(T_6)) = 4$ (consistent with THM-318: $\lfloor 2 \cdot 6/3 \rfloor = 4$).

## Evaluation

$H(T_6) = I(\Omega(T_6), 2) = 1 + 11500 + 16976 + 5136 + 160 = 33773$

This is the number of Hamiltonian paths in the all-0 staircase on 12 vertices (verified by THM-326).

## Note

$H(T_6) = 33773 \ll a(12) \geq 531205$ (ratio ≥15.7). The staircase becomes increasingly non-optimal as k grows.
