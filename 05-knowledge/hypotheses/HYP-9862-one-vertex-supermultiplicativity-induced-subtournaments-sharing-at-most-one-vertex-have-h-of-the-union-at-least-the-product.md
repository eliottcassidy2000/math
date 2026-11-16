---
id: HYP-9862
title: "One-vertex supermultiplicativity: for induced subtournaments A, B of a tournament with |A n B| <= 1, H(A u B) >= H(A) H(B)"
status: >
  OPEN (CONJECTURE), klein-2026-11-16-join.
  FINITE-EXACT: every tournament of orders 4-9 and every cover with exactly one shared vertex, no violation
  (218,925,648 instances at order 9; 139,658 equalities). |A n B| = 0 is the classical product bound.
  SHARP: false with two shared vertices (two triangles on an edge, 5 < 9; 695,729 failures among 434,403,648 at order 9).
  PROVED special case: THM-5312 (A), H >= p + q + 7 for prime counts not both 3.
session: klein-2026-11-16-join
date: 2026-11-16
source: 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/ONE_VERTEX.md
scripts:
  - 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/supmult.c, onevx.c
related:
  - THM-5312 (A), CONDITIONAL.md (conditional_determinant_vitali_20261116), fresh_unit_prime_gluing_20261116/PRIME_GLUE.md
---

# HYP-9862 - the first shared vertex buys nothing

**Conjecture.** For induced subtournaments A and B with |A ∩ B| ≤ 1, H(A ∪ B) ≥ H(A)·H(B).

**If true:**
- gluing that keeps induced witnesses needs at least two shared vertices to undercut the product;
- the minimum one-vertex excess over p + q is at least (p − 1)(q − 1) − 1. The census table `onevx.out`
  is consistent with this.

**Line of attack.** Splice two odd cycles through the shared vertex v. Each conflicting pair (C_A, C_B)
yields two mixed odd cycles: one or two long splices of length |C_A| + |C_B| − 1, or two triangles. The
open step is a weight-preserving injection, because the triangle splices forget the cycles.

**Falsifier.** One cover with |A ∩ B| = 1 and H(A ∪ B) < H(A)H(B), at any order.
