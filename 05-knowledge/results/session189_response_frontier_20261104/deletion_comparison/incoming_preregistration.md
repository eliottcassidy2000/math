# Pre-registration: kind-pasteur-2026-11-04-S9

This file was written and pushed BEFORE the confirmatory runs listed below. Nothing above the "Outcomes" heading will be
edited after the runs.

## Setting

- **Objects.** These are the Session188/189 objects.
  - A named tournament has ports x->y attached, with a_v=[v->x] and c_v=[v->y].
  - F_T(a,c) is the Hamilton count of the augmented tournament.
  - A core is W = Q[A_1..A_m]. Its target is W' = Q^op[A'_1..A'_m] with |A'_i| = |A_i|.
- **Mod-4 contrast.** delta = (F_W - F_W')/2 mod 2. It is a Boolean polynomial in the query variables.
- **Coarse.** The ports are uniform on a child.
- **Kernel.** Child j is transitive, b_1 is its source, and A_j - b_1 is contracted (it is a module). Two cases:
  - k >= 3: weights (b_1, beta) = (1, 0);
  - k = 2: weights (1, 1).
- **Hard case.** No child other than A_j is a singleton.
- **Endpoint parity.** p_j(Q) = S_j(Q) + E_j(Q) mod 2 (Hamilton starts plus Hamilton ends at j).

## Exploratory findings so far (NOT pre-registered)

These come from scratch scripts and are not confirmatory.
- **Weighted parity model.**
  - The weighted odd-cycle parity model reproduces the exact named mod-4 contrasts: 6/6 unweighted pairs and 6/6 contracted-TT2 pairs.
  - The fast implementation agrees with the reference implementation in 40/40 cases.
- **Blind cases.** Universality sweep of the weighted model: 2272 comparisons, 0 blind. They cover:
  - random Q of orders 2..7, prime or not;
  - random outside weights;
  - one TT_k child with k <= 5;
  - all reorderings for k <= 4.
- **Weighted E1 identity.** Sum_{u != v} [a_u c_v] delta + Sum_{v->u} [a_u a_v c_v] delta = w(v) held in 1814/1814 random cases (n = 3..9).
- **Kernel identities**, 600/600 random cases:
  - reversed child: L_{b1} = 1;
  - unchanged child, k >= 3, hard case: Lambda' = 1 + p_j(Q);
  - swap-orbit sums of the contrast equal those of the whole-converse contrast.
- **Actual graphs**, n <= 11:
  - singleton child: coarse L_z = 1 in 8/8 cases, with arbitrary same-size target interiors;
  - TT3 kernel: functional values as predicted in 10/10 cases.

## Claims to be proved (all orders), recorded before the confirmatory runs

- **(A) Coarse reduction mod 4.** For queries uniform on child A_i (|A_i| >= 2), delta equals the weighted cycle-parity contrast with A_i contracted to one vertex of weight 0. Singletons have weight 1.
- **(B) Weighted E1.** For every weighted tournament Z on n >= 3 vertices and every v,
  L_v := Sum_{u != v} [a_u c_v](g_Z - g_{Z^op}) + Sum_{v->u} [a_u a_v c_v](g_Z - g_{Z^op}) = w(v) (mod 2).
  Here g is the weighted port-cycle parity polynomial.
- **(C) Singleton child.** If some child of the core is a singleton, a coarse query separates W from W' mod 4, for arbitrary same-size target interiors.
- **(D) Transitive child.** If A_j is transitive of order k >= 2 and A'_j is A_j or A_j^op, a query separates mod 4. That query is uniform on every other child and on A_j - b_1, and named at b_1.
- **(E) Lift.** These lift to every fixed ambient, via the reviewed S188 lift and restoration. So any original-module full-response collision has, at every highest changed prime node:
  - no singleton child;
  - no transitive child of order 2 or 3, nor of order >= 5 (via S189);
  - transitive children of order 4 only with a reordering whose ends are the two middle names.

## Predictions for the confirmatory runs

- **(P1) Weighted E1, exhaustive.**
  - Scope: every labeled tournament with n = 3, 4, 5, every weight vector w in {0,1}^n, every v.
  - Plus 300 random (Z, w) at each n = 6..9.
  - Prediction: L_v = w(v) in every case. No exception.
- **(P2) Kernel, exhaustive in the weighted model.**
  - Scope: every labeled Q with m = 2..5; every j; kernel types k=2 and k>=3; u = source or u = sink of A_j; targets "unchanged" and "reversed".
  - Outside weights: all 2^(m-1) vectors for m <= 4; for m = 5, the all-zero vector and each single-singleton vector.
  - Predictions:
    - (i) every contrast is nonzero;
    - (ii) reversed: L_u = 1;
    - (iii) hard case, unchanged, k >= 3: Lambda' = 1 + p_j(Q), where Lambda' = L_u + L_beta + the missing swap partner of the u-beta cubic;
    - (iv) hard case, unchanged, k = 2: L_u = 1;
    - (v) non-hard case: L_z = 1 at every weight-1 outside z, for both targets.
- **(P3) Random kernel, larger orders.** 400 random (Q, j) with m = 6..9. Hard case and random outside weights. Predictions (i)-(v) as in P2.
- **(P4) S188 prime-eight endpoint-obstruction quotient** Q8 = (238,4,240,118,195,146,34,74).
  - Scope: every j; one TT_k child, k = 2..5, in the weighted model; all other children nontrivial.
  - Predictions:
    - both targets have nonzero mod-4 named contrast at every j and k;
    - reversed: L_u = 1;
    - unchanged: Lambda' = 1 for k >= 3 (p(Q8) = 0), and L_u = 1 for k = 2;
    - coarse-only contrast is zero, which reproduces S188's blindness.
- **(P5) Actual graphs, no contracted model.** The exact C++ named coefficients are restricted by summing named monomials into block monomials.
  - (a) 12 random prime Q of order 5 or 6, with a singleton child and other children random tournaments of orders 2-3. Source and target interiors are independent random tournaments of equal sizes, n <= 11. Prediction: coarse L_z = 1 at every singleton z.
  - (b) 12 random cores with one TT_k child (k = 3 or 4) and all other children nontrivial random tournaments, n <= 11, both targets. Predictions: the kernel functional values of P2 (ii)-(iv), and nonzero contrast.
  - (c) 4 ambient controls. The core M is a module inside an exterior of 2-3 extra vertices, and the ports imitate M's exterior incidences. Prediction: the same functional parity as the bare core.

A failure of any prediction will be recorded below, unedited, together with its counterexample.

## Outcomes

(appended after the runs)

