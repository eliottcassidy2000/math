---
id: HYP-9456
title: "Residual highest-prime frontier after THM-4864: in any original-module full-response collision, every highest changed node has a prime quotient Q of order >= 7 with p(Q) = 0, no singleton child, and no child with a vertex that is a source or sink of both its source and target interiors (so every child is non-transitive or a crossed TT4). Conjecture: every such child admits a one-child separator modulo 8 (named inside that child only; modulo 4 fails already for Paley-7 with unchanged C3 children), so no collision exists."
status: >
  OPEN (kind-pasteur-2026-11-04-S9; title aligned with the modulo-8 revision after audit MINOR-3, MISTAKE-703).
  - The reduction (first sentence) is PROVED: THM-4864, widened per audit NOTE-4, combined with S188 and S189.
  - Modulo 4 is FALSE as a universal one-child mechanism. The EXPLORATORY probe residual_c3 found Paley-7 with C3
    children kept unchanged blind to every query named in one or two children modulo 4.
  - Modulo 8 is SUPPORTED by a single FINITE-EXACT instance: an exact run-word computation separates that
    configuration at 4 mod 8, already on a coarse coefficient.
  - Crossed TT4 has EXPLORATORY weighted-model support: 0 blind cases among the 2,272 transitive comparisons.
session: kind-pasteur-2026-11-04-S9
source: 05-knowledge/results/reversal_cube_mod4_anchor_20261104/PROOF.md
related: [THM-4864, S188 prime_seven_frontier, S189 collision, S188 parity_profile_transfer, S188 coarse_mod8_affine]
---

# HYP-9456: residual highest-prime frontier

## What is left

After THM-4864, a still-possible original-module collision needs, at every highest changed prime node M:
- **order and parity:** quotient order >= 7 and zero endpoint parity p(Q) = 0;
- **singletons:** no child is a singleton;
- **common ends:** no child A has a vertex that is a source or a sink of both T[A] and U[A] (widened THM-4864
  part 4).

**Transitive children.** By Session189 and THM-4864, none of order 2, 3 or >= 5 survives, and a TT4 survives only with
a crossed reordering, whose target ends are {b_2, b_3}.

**Non-transitive children.** A child survives only if no vertex is an end of both interiors. Strong children, such as
C3, never have an end at all.

## Modulo 4 is not enough for C3 children (EXPLORATORY, FINITE-EXACT single instance)

The data is in `reversal_cube_mod4_anchor_20261104/residual_c3/`.
- **Blind configuration.** Take the Paley-7 quotient with C3 children: the C3 keeps its orientation, the quotient is
  conversed. Every query named in one or two C3 children has residue 0 mod 4 (weighted model).
- **Why the reviewed tools miss it.** Paley-7 has no cyclic ordered pair-module witness, so the S188 cyclic-child
  quartic does not apply.
- **Not a collision.** An exact run-word computation of T = Paley7[C3^7] against U = Paley7^op[C3^7] separates them.
  - 870 coefficient differences, with minimum 2-adic valuation 2.
  - One is the coarse coefficient [alpha_1 gamma_2], with difference -4 x odd.
  - So the configuration is separated modulo 8, not modulo 4.

## Tests

1. **Crossed TT4.**
   - Approach: compute the four-name weighted model, with the child fully named and every other child contracted.
   - Status: exploratory, 0 blind cases.
   - Goal: find an invariant functional analogous to the tau-orbit sum of THM-4864.
2. **Strong children, starting with C3.**
   - Approach: a modulo-8 one-child law. The template is the S188 five-type coarse_mod8_affine law: a C3 child is an
     odd-order child with odd Hamilton-cycle parity.
   - Data: the Paley-7 instance shows the gap already appears on coarse coefficients at 4 mod 8.
   - Next step: pre-register a census of exact one-child coefficients over children of order 3-5 (run-word DP) and
     prime quotients of order 7-9 with p = 0.
3. **Non-transitive children with an end vertex.** These survive only if U moves that end. Can hereditary profile
   equality force a common end?
   - Exploratory finite-exact census over labeled tournaments: each source child A, and each target A' with
     q_h(A[S]) = q_h(A'[S]) for every subset S of size >= 2.
   - **Order 4: yes.** All 64 such pairs share an end. So, by THM-4864, every non-transitive order-4 child with a
     source or sink is excluded, and only the 24 labeled strong order-4 children remain.
   - **Order 5: no.** 320 of 3,360 pairs have no common end, and 544 labeled strong children have no end at all.
   - The residual child list therefore starts with the strong children of order >= 3, plus specific order-5
     non-strong pairs.
