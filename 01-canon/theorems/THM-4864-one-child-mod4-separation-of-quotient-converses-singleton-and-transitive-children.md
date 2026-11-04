---
id: THM-4864
title: "One-child modulo-four separation of quotient converses. Let W = Q[A_1..A_m] and W' = Q^op[A'_1..A'_m] on the same child sets, with |A'_i| = |A_i| (target interiors otherwise arbitrary). (A) If some child is a singleton (m >= 3), a child-uniform two-port query has F_W - F_W' = 2 mod 4. (B) If A_j is transitive of order k >= 2 and an end u of its order is also an end of the transitive target A'_j - in particular if A'_j is A_j (unchanged) or A_j^op (reversed) - a query uniform on the other children and on A_j - u, named only at u, has residue 2 mod 4 (or a coarse query does, when p_j(Q) = 1). Both lift to every fixed ambient. Consequently, in any original-module full-response collision, every child of every highest changed prime node is non-transitive or a TT4 whose target order has its two middle names as ends."
status: >
  PROVED (pre-audit), given the cited reviewed inputs:
  - the odd-cycle Hamilton expansion modulo 4 (Grinberg-Stanley Thm 1.39);
  - Redei;
  - the S188 parity_profile_transfer endpoint lemma and its ambient lift;
  - the S188 highest_quotient_reduction, prime_seven_frontier restoration and module_decomposition hereditary equality;
  - Session189 collision Section 5, which keeps a transitive child of order >= 5 unchanged or reversed.
  New ingredients: a weighted E1 identity (Proposition 4) and a swap-orbit argument at a contracted transitive child.
  FINITE-EXACT pre-registered confirmations (P1-P5) are in the packet. An independent audit is pending.
session: kind-pasteur-2026-11-04-S9
source: 05-knowledge/results/reversal_cube_mod4_anchor_20261104/PROOF.md
scripts:
  - 05-knowledge/results/reversal_cube_mod4_anchor_20261104/confirm.py (+ confirm.out, confirm.json)
  - 05-knowledge/results/reversal_cube_mod4_anchor_20261104/wmodel.py
  - 05-knowledge/results/reversal_cube_mod4_anchor_20261104/fcoef.cpp
related: [HYP-9456, S188 parity_profile_transfer, S188 prime_seven_frontier, S189 collision, two_port_parity_separation E1-E4, THM-1370-h-spectrum-omits-7-21-all-n]
---

# THM-4864: one-child modulo-four separation

## Statement

Ports satisfy x->y, a_v = [v->x] and c_v = [v->y]. F is the augmented Hamilton count. Let W = Q[A_1..A_m] and
W' = Q^op[A'_1..A'_m] share their child sets, with |A'_i| = |A_i|. Let p_j(Q) = S_j(Q) + E_j(Q) mod 2.

1. **Weighted E1.** Let Z be a tournament on n >= 3 vertices with weights w in F2. Let g_Z be the parity of
   port-containing directed cycles of odd total weight, with each port weighing 1, and let delta_Z = g_Z - g_{Z^op}.
   Then for every v

       sum_{u != v} [a_u c_v] delta_Z + sum_{v->u} [a_u a_v c_v] delta_Z = w(v) mod 2.

2. **Contraction.** On queries uniform on each child, (F_W - F_W')/2 mod 2 equals a constant plus
   g_Q - g_{Q^op}, with weights [|A_i| = 1].
3. **Singleton child.** If m >= 3 and some child is a singleton, a child-uniform query separates W from W' modulo 4.
4. **Transitive child.** Let A_j be transitive of order k >= 2, and suppose an end u of A_j's order is an end of A'_j.
   Then a query separates W from W' modulo 4. Either it is uniform on every other child and on A_j - u, with free bits
   at u, or it is child-uniform; the second happens only when p_j(Q) = 1 and no other child is a singleton.
   - The cases A'_j = A_j and A'_j = A_j^op are included, so **unchanged and reversed are both separated**.
   - Every reordering of a TT2 or a TT3 is included.
5. **Lift and corollary.** These residue-2 queries lift to every fixed ambient when the ports copy the core's exterior
   incidences. Suppose a collision F_T = F_U exists with U an original-module XOR target. Then at every highest changed
   prime node M (order >= 7 by S188):
   - no child is a singleton;
   - p(Q) = 0;
   - no child is transitive of order 2, 3 or >= 5 (order >= 5 via Session189);
   - every transitive child of order 4 is mapped to an order with its two middle names as ends.

## Proof

See the source packet: Lemmas 1-3, Proposition 4, Theorems A-B and the Corollary.
- **Unchanged case.** It is caught by the tau-invariant functional L_u + L_beta + [a_u a_beta c_beta]. Its value on the
  contrast equals its value on the whole-converse contrast, namely 1 + p_j(Q). When p_j(Q) = 1, a coarse query
  separates instead.
- **Reversed case.** It is weighted E1 at u.
- **k = 2, hard case.** The pair-flip correction has no mixed monomials, so weighted E1 at u applies directly.

## Scope

- **Not claimed.**
  - No reconstruction.
  - No integer-gap statement.
  - No separation for non-transitive children or TT4 crossed reorderings. Those are HYP-9456, with FINITE-EXACT
    evidence of 0 blind cases.
  - No collision at any order.
- **Sharpness of the coarse layer.** S188's prime-eight endpoint-obstruction quotient is coarse-blind modulo 4, but it is
  separated by the anchored queries of part 4 at every transitive child (pre-registered P4).
