---
id: HYP-9456
title: "Residual highest-prime frontier after THM-4864: in any original-module full-response collision, every highest changed node has a prime quotient Q of order >= 7 with p(Q) = 0, and every child of it is either non-transitive (order >= 3) or a TT4 whose target order has its two middle names as ends. Conjecture: every such child also admits a one-child modulo-four separator (named inside that child only), so no collision exists."
status: >
  OPEN (kind-pasteur-2026-11-04-S9).
  The reduction itself (first sentence) is PROVED: THM-4864 combined with S188 and S189.
  The conjecture has FINITE-EXACT support in the weighted model, 2,272 comparisons with 0 blind cases:
  random Q of order <= 7, random weights, one transitive child of order <= 5, and every reordering for order <= 4
  (TT4 crossed reorderings included).
  Non-transitive children have not yet been tested by the one-child weighted model.
session: kind-pasteur-2026-11-04-S9
source: 05-knowledge/results/reversal_cube_mod4_anchor_20261104/PROOF.md
related: [THM-4864, S188 prime_seven_frontier, S189 collision, S188 parity_profile_transfer]
---

# HYP-9456: residual highest-prime frontier

## What is left

After THM-4864, a still-possible original-module collision needs, at every highest changed prime node M:
- **order and parity**: quotient order >= 7 and zero endpoint parity p(Q) = 0;
- **children**: none is a singleton;
- **transitive children**: none of order 2, 3 or >= 5;
- **TT4 children**: only with a crossed reordering, so the target order has ends {b_2, b_3}.

Every remaining child is therefore a non-transitive tournament (it contains a cyclic triangle) or a crossed TT4.

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
- **Revised conjecture for non-transitive children.** A one-child separator exists modulo 8. A modulo-4 separator
  exists for transitive children (THM-4864).

## Tests

1. **Crossed TT4.**
   - Approach: compute the four-name weighted model, with the child fully named and every other child contracted.
   - Status: exploratory, 0 blind cases.
   - Goal: find an invariant functional analogous to the tau-orbit sum of THM-4864.
2. **Non-transitive children.**
   - First case: a C3 child.
   - Kernel candidate: name the three cycle vertices separately and contract nothing inside the child.
   - Larger children: contract a module of the child when one exists.
   - Next step: pre-register a weighted-model census over all children of order 3-5, with random prime quotients of
     order 7-9.
3. **Integer layer.** If a modulo-4 one-child separator fails for some child type, the next modulus is 8. The S188
   cyclic mod-8 and mod-16 analyses are the template.
