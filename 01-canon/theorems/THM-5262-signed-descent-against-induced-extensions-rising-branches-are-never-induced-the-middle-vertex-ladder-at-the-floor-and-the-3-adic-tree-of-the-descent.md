---
id: THM-5262
title: "The signed descent against induced extensions. (A) In half-counts n = (H-1)/2 the increment floor reads n(T') = n(T) or n(T') >= (3n(T)+1)/2 for every T' containing T, the shortcut Collatz odd step; hence no realiser of a rising target (4d-3)/3 or (4d+1)/3 of the signed descent contains a realiser of d. (B) Middle vertex on a bridge: for the ordinal sum A -> B with the arc u -> w reversed and a new vertex x with A -> x -> B, H(T) = 2ab + Delta - E_A(u) S_B(w), H(T+x) = 3ab + 2 Delta, 2H(T+x) - 3H(T) = Delta + 3 E_A(u) S_B(w). (C) Ladder: for m = phi_a phi_b (phi_0 = 1, phi_a = 2^a+1) a tournament of count 2m-1 sits in one of count 3m; for a = 0 the half-count chain is k -> 2k -> 3k+1. (D) The backward tree of the descent: dead rays = half-counts divisible by 3; roots with j rising steps = 2^j classes mod 3^j; Haar-average level size (5/3)^k; criticality 1/2 + (2/3)(3/4) = 1; an orbit meets a dyadic window at most 3 times. (E) Finite-exact: complete parent sets of all counts up to 223; 21 is the record rising root"
status: >
  PROVED: (A) as a corollary of the increment floor (KNOWN: HYP-9469 (1), (6), proved in all orders by Session214's
  forced-shuffle theorem; sole-exit case THM-4876 (7)); (B) by direct path splitting (the one-apex case of the
  bridge is chanpak_target_ribbon_audit_20261116/TWO_APEX.md, another session); (C) from (B) and TWO_APEX (5);
  (D) all five statements (elementary).
  FINITE-EXACT: (E) the parent sets P(h) = {H(T-x) < h : H(T) = h} for the 109 counts h <= 223 (every tournament
  of order <= 10, closed under ordinal products; complete because a strong tournament of order 11 has H >= 225):
  (h+1)/2 is a parent for 38 of the 54 counts h = 1 (mod 4), not for 13, 41, 57, 69, 77, 85, 93, 105, 109, 125,
  133, 141, 153, 169, 173, 201; min |h - 2d| over parents is 1 for 79 counts h >= 5, 3 for 27, 5 for 41 and 125;
  over all tournaments of order at most 9 the nontrivial one-vertex extensions with 2H' - 3H <= 3 are
  3 -> 5 and 2m-1 -> 3m for the twelve products m = phi_a phi_b <= 65; rising-chain records n = 2, 4, 10, 244,
  2660, 4648, 5212, 34628, 58130, 502750, 737596 below 2*10^6 (10 is the count 21, nine steps).
  NUMERICAL: basin fractions 0.4811 (through 21), 0.3676 (27), 0.1136 (41) of the odd numbers below 2*10^7;
  F(t) = t rho(t) with mean 6.129 and means 3.94 / 6.6 / 7.87 over n mod 9 (the class-mean ratios 9/14, 15/14, 9/7
  are PROVED consequences of the functional equation, assuming the densities exist); the dominance-pruned
  slow-factor tree at rate 1.72 (4, 20, 24, 84, 344, 837, 2616, 4971, 11647 undominated nodes at depths 1-9; it
  does not close).
  NOT CLAIMED: Chan-Pak (OPEN); existence of the basin densities; any closing of the literal span-3 certificate;
  the persistence of near-doubling parents (HYP-9812). The signed descent, its boundary {13, 27, 41} and the
  half-count chart are from six_tutte_collatz_20261116/COUNTS.md; explicit infinite families reversing each
  branch by surgery are in inverse_fences_landing_20261116 (both other sessions).
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-descent
date: 2026-11-15
source: 05-knowledge/results/chanpak_signed_descent_lifts_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_signed_descent_lifts_20261115/scripts/lift_gap.c, lift_par.c, lift_near.c,
    gapwit.c, parents_closed.py, port_calculus_check.py, near_floor_family_check.py, ladder_parents.py,
    dbasin.py, dtree.py, domtree.py, domanalyse.py, domtree2.py
related:
  - HYP-9812, HYP-9469, THM-4876, THM-5252, THM-4506, THM-4510, HYP-9742
  - 05-knowledge/results/six_tutte_collatz_20261116/COUNTS.md
  - 05-knowledge/results/inverse_fences_landing_20261116/INVERSE.md
  - 05-knowledge/results/chanpak_target_ribbon_audit_20261116/TWO_APEX.md
  - 05-knowledge/results/session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md
---

# THM-5262 - the signed descent against induced extensions

Source: [the note](../../05-knowledge/results/chanpak_signed_descent_lifts_20261115/README.md).

`H(T)` is the number of Hamiltonian paths of a tournament, `n(T) = (H(T) - 1)/2` its half-count. `T'` contains
`T` as an induced subtournament. The signed descent is `D(h) = (h+1)/2` (`h = 1 mod 4`), `(3h+3)/4`
(`h = 3 mod 8`), `(3h-1)/4` (`h = 7 mod 8`).

**(A)** For `T` inside `T'`: `n(T') = n(T)` or `n(T') >= (3n(T) + 1)/2`. Consequently `H(T')/H(T)` is 1 or at
least `3/2`, and since `d < (4d-3)/3 < (4d+1)/3 < (3d+1)/2` for `d > 3`, a realiser of a rising target never
contains a realiser of its descent value. If `H(T') = 2H(T) - 1` then `H(T)` is the count of a one-vertex deletion
of a realiser of `H(T')`.

**(B)** Notation: `E_A(u)`, `S_A(u)` = Hamiltonian paths of `A` ending, starting at `u`; `U_A(u)` = ordered pairs
`(P1, P2)` of disjoint paths covering `A`, `P2` starting at `u`, `P1` possibly empty; `U*_A(u)` likewise with `P1`
ending at `u`, `P2` possibly empty; `D_A(u) = U_A(u) - H(A) >= H(A - u) - E_A(u) >= 0`, `D*_A(u) = U*_A(u) - H(A)`.
Let `T` be the ordinal sum `A -> B` with the arc `u -> w` reversed, `a = H(A)`, `b = H(B)`, and let `x` be a new
vertex with `A -> x -> B`. Then

    H(T) = ab - E_A(u) S_B(w) + U_A(u) U*_B(w),        H(T + x) = ab + 2 U_A(u) U*_B(w),

and with `Delta = D_A(u) b + D*_B(w) a + D_A(u) D*_B(w)`: `2H(T + x) - 3H(T) = Delta + 3 E_A(u) S_B(w)`.

*Proof.* The only arc from `B` to `A` is `w -> u`. A Hamiltonian path of `T` avoiding it is a path of `A` followed
by a path of `B` whose junction is not `(u, w)`; one using it is `P1, (Q1 w), (u P2), Q2` with all other junctions
forward. In `T + x` the vertex `x` must follow a vertex of `A` (or be first) and precede a vertex of `B` (or be
last): one place in the first kind of path, exactly two in the second.

**(C)** `Delta = 0` and `E_A(u) S_B(w) = 1` hold exactly when `A` is the fan `1 0^r` at its apex or a single
vertex, up to a transitive block placed before it, and `B` the dual fan `1^s 0` or a single vertex (then
`H(A - u) = 1` and `A(w) + 1 = h(w)`, TWO_APEX (5)).
So for `m = phi_a phi_b` the bridge has count `2m - 1` and the bridge with its middle vertex `3m`. For `a = 0`
the bridge contains the dual fan, of count `m`: half-counts `k -> 2k -> 3k + 1`.

**(D)** In half-counts the children of `n` are `2n` and, for `3` not dividing `n`, `r(n) = (4n - e)/3` with
`e = ±1`, `e = n (mod 3)`. (i) If `3 | n` the only child is `2n`; the basin of 7 is `{3·2^k + 1}`. (ii) `r` maps
the classes mod `3^j` inside a unit class mod 3 bijectively onto the classes mod `3^(j-1)`, so the `n` with at
least `j` consecutive rising steps above them form `2^j` classes mod `3^j`. (iii) The number `b_k(n)` of units at
level `k` above `n` satisfies `b_(k+1)(n) = b_k(2n) + [3 ∤ r(n)] b_k(r(n))` and its Haar average over the 3-adic
units is `(5/3)^k`. (iv) `1/2 + (2/3)(3/4) = 1`. (v) `h/2 < D(h)` and `D(h) - 3 <= (3/4)(h - 3)`, so an orbit that
starts above `2Y` has 1, 2 or 3 points in the window `[Y, 2Y)`, `Y >= 12`.

**(E)** See the status block and `out/` of the note.
