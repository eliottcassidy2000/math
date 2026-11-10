---
id: HYP-9656
title: "Self-feeding one-arc closures with bounded root response are finite, and complete on a growing prefix. Start from one vertex; close under one-arc joins with a fixed finite gadget library (new root inside the gadget) and sink padding; keep only rooted tournaments whose root response B has odd part at most r0. (A) For every r0 and every finite library the closure is finite. (B) With the order-<= 5 library the closure for r0 contains every odd n < N(r0) except 7 and 21, where N(r0) = 13, 35, 35, 67, 127, 431, 2303, 5153, 20921, 55343 for r0 = 3, 7, ..., 2047 (roughly r0^1.8 over r0 = 63-2047). (C) Consequently a one-arc construction reaching N needs root responses with odd part growing with N; finite banks and pure dyadic cones (r0 = 1) cannot reach every large N"
status: >
  OPEN.
  (A) is FINITE-EXACT saturation evidence, from pre-registered P1, P2, P3, P5 and P6, all CONFIRMED:
  - r0 = 1 with the 22,243 class-gadget matrices of order <= 7 gives 197 values, none above 2 * 10^6, in a run to
    10^10;
  - r0 = 63 with order <= 6 gives 12,352 values, none above 10^8, in a run to 10^9;
  - r0 = 1023 and 2047 with order <= 5 give 93,159 and 233,115 values, all below 3.41 * 10^8, in runs to 10^9 and
    4 * 10^9.
  (B) The table up to r0 = 1023 was exploratory. N(2047) = 55,343 was pre-registered (P5, window [40,000, 120,000]).
  (C) is an interpretation, not a theorem.
  For r0 = 1 and transitive gadgets the non-chaining part is PROVED (THM-5109 (E)).
session: kind-pasteur-2026-11-09-S20
source: 05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/PREREGISTRATION.md
related:
  - THM-5109
  - HYP-9453
  - HYP-9640
  - HYP-9480
  - HYP-9655
  - THM-4875
  - THM-5085
---

# HYP-9656 - bounded-response closures are finite but complete on a growing prefix

## Definitions

- **Gadget library.** The 3 x 3 one-arc matrices M(G, u, z) of THM-4875 (B), for a finite set of gadgets: all
  labelled gadgets of order <= 5, or all class gadgets of order <= 6 or <= 7. The new root z lies in the gadget.
- **Closure C(r0).** Begin at the single vertex (1, 1, 1). Repeatedly apply a matrix of the library or sink padding
  (h, e, R) -> (h, 0, 2R). Keep only the triples with odd part of R at most r0.
- **Exactness.** Every member is the triple of an actual rooted tournament. Along each step h never decreases
  (THM-4875 (C)), so the part of C(r0) below a bound is computed exactly.

## Statements

- **(A) Finiteness.** C(r0) is finite for every r0 and every finite library.
- **(B) Complete prefix.** For the order-<= 5 library, C(r0) contains every odd n < N(r0) except 7 and 21:

  | r0 | 3 | 7 | 15 | 31 | 63 | 127 | 255 | 511 | 1023 | 2047 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | N(r0) | 13 | 35 | 35 | 67 | 127 | 431 | 2,303 | 5,153 | 20,921 | 55,343 |
  | number of values | 35 | 83 | 258 | 724 | 2,103 | 5,580 | 14,392 | 36,765 | 93,159 | 233,115 |

  Conjecture: N(r0) tends to infinity. That would prove Chan-Pak for this construction class: the union over r0 is
  the full one-arc closure of a point.
- **(C) Interpretation.**
  - Finite banks have bounded root responses: compare THM-5109 (F), and HYP-9480 (C) / HYP-9655.
  - Pure dyadic cones have r0 = 1.
  - By (A), neither reaches every large N. An inductive proof must track root responses whose odd parts grow, as in
    HYP-9453's response intervals.

## Falsifiers

- **For (A):** a library and an r0 whose closure keeps producing new values beyond every bound tested.
- **For (B):** a decrease of N(r0) in r0 is impossible, because the closures are nested. A falsifier is instead
  stagnation: N(r0) constant over a long range of r0 while the number of values grows.
