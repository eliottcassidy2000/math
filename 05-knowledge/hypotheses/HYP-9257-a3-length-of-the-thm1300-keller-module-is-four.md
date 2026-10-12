---
id: HYP-9257
title: "The A_3-length of THM-1300's Keller module V_F = K[x, y, z] (twisted by the Keller map) is exactly four: the simple submodule K[F], the minimal extension of the S_3-standard rank-two local system on A^3 minus the Jelonek surface, a delta module along the Jelonek surface {c_3 = 0}, and one constituent supported on the cusp curve C = {(4/(27t^2), 4/(3t), t)} where the fibre is empty"
status: "OPEN (2026-10-12, mac-mini-2026-10-12-logic4). PROVED: length >= 3 (THM-4621(9)); finite length with an explicit bound (Bavula, CITED). The value four is a guess driven by the fibre stratification 3 / 1 / 0 (generic / Jelonek surface / cusp curve): each drop should be paid by a constituent supported on the stratum. The cusp-curve contribution is the uncomputed part."
session: mac-mini-2026-10-12-logic4
related:
  - THM-4621 (the three named constituents), THM-4620(6), THM-1300
  - V. V. Bavula, Holonomic modules and 1-generation in the Jacobian Conjecture (length bound)
  - S. Gao, arXiv:2608.00222, Theorem 3.4 (the stratification; sigma^(-1)(K), the cusp)
---

# HYP-9257 — the length of the dimension-three Keller module

## Statement

`length_(A_3)(V_F) = 4`.

## What is proved (THM-4621(9))

`V_F ⊋ K[F] ≅ O`; `(V_F/O)|_U` is the irreducible rank-two local system of the standard
representation of the monodromy group `S_3` over `U = A^3 ∖ {c_3 = 0}`; near a generic point of the
Jelonek surface `V_F ≅ O ⊕ O[1/s] ⊕ j_+ L_(-1)` (one sheet stays, two escape with a transposition
monodromy), so a delta module along the surface is a subquotient. Hence length `≥ 3`.

## What decides it

The local structure of `V_F` at points of the cusp curve `C`, over which all three sheets escape
(empty fibre). Compute the characteristic cycle or the nearby/vanishing cycles of the direct image
along `F` near `C`; a transversal slice meets the Jelonek surface in a cuspidal curve, and the
question is whether the direct image of `O` from the three escaping sheets contributes a constituent
supported on `C` beyond the delta module of the surface. (For the sign local system around a cuspidal
plane curve, `j_+ L_(-1)` is already the minimal extension because the Alexander polynomial of the
trefoil does not vanish at `-1`; the trivial-monodromy part is where a cusp contribution could sit.)

## Decisive test

A Gröbner/holonomic computation of the `b`-function and the characteristic cycle of `V_F` as an
`A_3`-module presented by the relations `V_j · 1 = 0` is out of reach in the repo's toolkit (no
Singular/Macaulay2); a topological computation of the local cohomology of the escaping sheets along a
transversal slice at a point of `C` is the feasible route.
