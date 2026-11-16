---
id: THM-5261
title: "Gluing without mixed odd cycles. Let A, B be induced subtournaments of U covering V(U), W = A n B, and for X subset of W let H_A(X) be the odd-cycle-packing weight (THM-002) of the packings of A whose cycle-covered W-vertices are exactly X. (A) If no odd directed cycle of U meets both A - W and B - W, then H(U) = sum over disjoint X, Y subset of W of H_A(X) H_B(Y). (B) If moreover every odd cycle of A and of B meets W and no two cycles, one of A and one of B, are vertex-disjoint (e.g. all pass through one shared vertex), then H(U) = H(A) + H(B) - 1; in general H(U) - H(A) - H(B) + 1 = sum over disjoint nonempty X, Y of H_A(X) H_B(Y), which explains the rows p + q - 1, p + q + 1, p + q + 3 of THM-5312 (C). With |W| = 1 the hypothesis of (A) never holds for pieces with cycles through the shared vertex (THM-5256), so additive gluing needs |W| >= 2; the two triangles on an arc (3 + 3 - 1 = 5) are the first instance, and the census of all pairs of orders <= 4 on a shared arc finds no other additive pair."
status: >
  PROVED (klein-2026-11-16-force; (A) immediate from THM-002 and the disjointness of cycles across the overlap; (B) from
  (A); not independently audited). FINITE-EXACT: every pair of tournament classes of orders <= 4 (and the strong pairs of
  orders (5,5), (4,5), (5,6)) glued on a shared arc with every cross orientation (scripts/glue.c, receipts glue_*.out).
  CITED: THM-002, THM-5256, THM-5312, session227 multiroot (set-root grafts), session225 (single-root graft).
session: klein-2026-11-16-force
date: 2026-11-16
source: 05-knowledge/results/forcing_digits_analogy_additive_20261116/ADDITIVE.md
scripts:
  - 05-knowledge/results/forcing_digits_analogy_additive_20261116/scripts/glue.c (+ receipts/glue_*.out)
related:
  - 01-canon/theorems/THM-002-ocf.md
  - 01-canon/theorems/THM-5256-path-splice-hyp-9862-proved-one-shared-vertex-supermultiplicativity-and-the-order-statistics-bound-for-any-overlap.md
  - 01-canon/theorems/THM-5312-gluing-targets-below-125-one-shared-vertex-forces-p-plus-q-plus-7-corrected-overlap-bound-and-the-conflict-graphs-of-7-and-21.md
---
# THM-5261 — gluing without mixed odd cycles

**(A)** Under the hypothesis every odd cycle of `U` is a cycle of `A` or of `B`. A packing of `U` restricts to a
packing of `A` and a packing of `B` whose cycles are vertex-disjoint, which for cycles of different sides means
disjoint cycle-covered `W`-sets; conversely any such pair with the uncovered `W`-vertices as singletons is a packing
of `U`, and `2^c` multiplies. ∎

**(B)** With `H_A(∅) = H_B(∅) = 1` the sum in (A) is `1 + Σ_{X ≠ ∅} H_A(X) + Σ_{Y ≠ ∅} H_B(Y) + Σ_{X, Y ≠ ∅, X ∩ Y = ∅}
H_A(X) H_B(Y) = H(A) + H(B) − 1 + (cross terms)`. ∎ For two triangles `u v a`, `u v b` on the arc `u → v` the mixed
cycles `u v a b`, `u v b a` are even, every cycle passes through `u` and `v`, and `H = 5`.
