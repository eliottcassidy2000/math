---
id: HYP-9656
title: "Self-feeding one-arc closures with bounded root response have finitely many values, yet are complete on a growing prefix. Start from one vertex; close under one-arc joins with a fixed finite gadget library (new root inside the gadget) and sink padding; keep only rooted tournaments whose root response B has odd part at most r0. (A) For every r0 and every finite library the VALUE set h(C(r0)) is finite (the triple set is always infinite, through sink padding). (B) With the order-<= 5 library, h(C(r0)) contains every odd n < N(r0) except 7 and 21, where N(r0) = 13, 35, 35, 67, 127, 431, 2303, 5153, 20921, 55343 for r0 = 3, 7, ..., 2047 (audit: 145,433 at r0 = 4095), growing like about r0^1.75 so far, possibly slowing. N(r0) -> infinity would make the full one-arc closure of a point complete (every odd n except 7, 21). (C) Heuristic, conditional on (A): a one-arc construction of this kind that reaches every large N needs root responses with odd part growing with N"
status: >
  OPEN.
  (A) has FINITE-EXACT saturation evidence from pre-registered P1, P2, P3, P5 and P6, all CONFIRMED and all
  reproduced by the audit:
  - r0 = 1 with the 22,243 class-gadget matrices of order <= 7 gives 197 values, max 1,995,807, unchanged up to
    10^13 (audit);
  - r0 = 63 with order <= 6 gives 12,352 values, max 86,371,467, unchanged up to 10^13;
  - r0 = 1023 and 2047 with order <= 5 give 93,159 and 233,115 values, max 340,931,297.
  Audit stress runs:
  - order <= 5, r0 = 4095: 580,505 values, max about 1.04 * 10^9, none in (1.04 * 10^9, 5 * 10^10];
  - order <= 6, r0 = 255: 104,217 values with a thin tail, 51 values in (10^8, 10^9] and 2 above 10^9 (max
    12,039,412,527). A run capped at 10^9 would have looked saturated, so saturation needs a wide margin above the
    largest value.
  (B): the table up to r0 = 1023 was exploratory; N(2047) = 55,343 was pre-registered (P5); N(4095) = 145,433 is the
  audit's.
  (C) is a heuristic, not a theorem.
  - It does not cover fixed-port tails or finite banks. Those keep the root in the core, their port responses have
    unbounded odd part, and their coverage question is HYP-9480 (C).
  - What is proved about pure cones is THM-5109's window bound and (E3).
  AUDIT: PASS WITH CORRECTIONS (MISTAKE-751);
  05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/audit/AUDIT.md.
session: kind-pasteur-2026-11-09-S20
source: 05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/PREREGISTRATION.md
related:
  - THM-5109
  - HYP-9453
  - HYP-9640
  - HYP-9480
  - THM-5070
  - THM-4875
  - THM-5085
---

# HYP-9656 - bounded-response closures have finitely many values but a growing complete prefix

## Definitions

- **Gadget library.** The 3 x 3 one-arc matrices M(G, u, z) of THM-4875 (B), for a finite set of gadgets: all
  labelled gadgets of order <= 5, or all class gadgets of order <= 6 or <= 7. The new root z lies in the gadget.
- **Closure C(r0).** Begin at the single vertex (1, 1, 1). Repeatedly apply a matrix of the library, or sink padding
  (h, e, R) -> (h, 0, 2R). Keep only the triples whose R has odd part at most r0.
  - Every member is the triple of an actual rooted tournament.
  - C(r0) is always infinite as a set of triples: it contains (1, 0, 2^k) for every k >= 1.
  - Only the value set h(C(r0)) can be finite.
- **Exactness of bounded computations.**
  - Every library matrix has E_u >= 1 and 0 <= s_u <= h_G. On realisable triples (0 <= e <= h) this gives
    h' = h_G h - s_u e + E_u R >= max(h, R).
  - So a step never lowers the value. A triple with R > HMAX produces only values above HMAX under gadget steps, and
    the same value under sink padding.
  - Hence capping both h and R at HMAX (the programs use 2 HMAX + 2 for R) leaves h(C(r0)) ∩ [1, HMAX] exact.

## Statements

- **(A) Finiteness of values.** h(C(r0)) is finite for every r0 and every finite library.
- **(B) Complete prefix.** For the order-<= 5 library, h(C(r0)) contains every odd n < N(r0) except 7 and 21:

  | r0 | 3 | 7 | 15 | 31 | 63 | 127 | 255 | 511 | 1023 | 2047 | 4095 (audit) |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | N(r0) | 13 | 35 | 35 | 67 | 127 | 431 | 2,303 | 5,153 | 20,921 | 55,343 | 145,433 |
  | number of values | 35 | 83 | 258 | 724 | 2,103 | 5,580 | 14,392 | 36,765 | 93,159 | 233,115 | 580,505 |

  - The least-squares exponent of N(r0) over r0 = 63-2047 is about 1.76. The last doubling gives about 1.4.
  - Conjecture: N(r0) tends to infinity.
    - This is equivalent to completeness of the full one-arc closure of a point: every odd n except 7 and 21. A
      closure with even one further miss keeps N(r0) bounded.
    - The r0 = infinity baseline is THM-5070 (4): the same order-<= 5 library attains every odd value up to 2 * 10^5
      except 7 and 21.
- **(C) Heuristic, conditional on (A).** A one-arc construction of this kind that reaches every large N must let the
  odd parts of its root responses grow with N.
  - Fixed-port tails are outside this scope. They keep the root in the core, and their port responses have unbounded
    odd part (for example 1, 3, 11, 87, 347, ..., 5551 from a singleton core).
  - Finite banks are outside it too. Their coverage question is HYP-9480 (C), independent of this hypothesis.
  - For pure cones the proved statements are THM-5109's window bound (no odd value beginning 11 in binary) and (E3).

## Falsifiers

- **For (A):** a library and an r0 whose value set keeps growing beyond every bound tested, with a wide margin.
- **For (B):** N(r0) cannot decrease, because the closures are nested. A falsifier is instead stagnation: N(r0)
  constant over a long range of r0 while the number of values grows. That would mean a value missed by every slice,
  that is, by the whole one-arc closure.
