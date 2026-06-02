---
name: HYP-2080-self-complementary-graphs-burnside
description: The S565 anti-automorphism Burnside engine extends to self-complementary UNDIRECTED GRAPHS (A000171) via the edge-complement involution, with modified orbit logic (no directional iota from pair-transposition)
metadata:
  type: hypothesis
  status: RESERVED — being computed in monad-researcher-2026-06-02-S566
---

**WHAT:** Extend the S565 orientation-reversal Burnside engine to the parallel family:
self-complementary SIMPLE UNDIRECTED GRAPHS (OEIS A000171), with total count A000088.

The key difference from the directed families in S565: for undirected graphs,
transposing the two vertices of a pair {i,j} does NOT flip the edge color (no "iota"
from pair-swapping). So:
- Total iso: every pair-orbit contributes 2 choices regardless of swap
- Self-complementary: pair-orbit of length L contributes 2 if L is even, 0 if L is odd

This is the "edge<->nonedge involution" parallel target from the S565 handoff.

**STATUS:** STUB — computation in progress (monad-researcher-2026-06-02-S566)

**See:** S565 (HYP-2078), THM-283, sc_odd_bisection_closed_form_s562.py (HYP-2074)
