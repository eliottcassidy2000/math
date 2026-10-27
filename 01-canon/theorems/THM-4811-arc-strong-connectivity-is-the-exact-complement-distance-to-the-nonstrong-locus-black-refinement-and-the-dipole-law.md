---
id: THM-4811
title: "Arc-strong connectivity is the exact native-complement distance to the non-strong locus; the black refinement and its exceptional family; the dipole law. (A) Dipole law: a native move along a legal word Q = (a, ..., b) sends the skew matrix K to -K + 2A_Q, the score deviations d to -d + e_a - e_b, and the cut function f(X) = e(X -> V\\X) to f'(X) = f(V\\X) + [a in X] - [b in X]; hence c3' - c3 = s(a) - s(b) - 1. (B) For every tournament T, |lambda(tau_Q T) - lambda(T)| <= 1 for the arc-strong connectivity lambda, with lambda(tau_Q T) = lambda(T) - 1 iff some minimum out-cut contains a and misses b; consequently the distance from T to the non-strong locus by native moves of either colour is exactly lambda(T), at most floor((n-1)/2) with equality exactly for the (semi-)regular tournaments. (C) With black moves only: through order 10 the black distance equals lambda(T) for every strong T except C3 and R5 (no black move) and one order-10 class E10, whose black distance is 2 = lambda + 1; E10 is the first member of an infinite family E(M; Z) (doubled C3[K1, M, Z] with K1/C3-stacks M, Z, the first block of M a C3) in which every descent word is blue, so the black distance exceeds lambda by at least one."
status: >
  PROVED (elementary): (A); (B) given Camion's theorem (proof inherited in blue_quotient_20261020b/structure.md, section 2);
  (C) the family statement: E(M; Z) is strong, lambda = 1, the converse copy D* is the unique minimum out-cut, and every
  descent word is blue.
  FINITE-EXACT: black_descent.cpp over every strong class of orders 3-10 (1, 1, 6, 35, 353, 6,008, 178,133, 9,355,949): the
  descent criterion of (B) verified literally on all legal words through order 8 (2,073,672 words at order 8, no failure);
  every strong class has a descent word, and a black one except C3, R5 (pure blue) and E10; at order 10 every lambda = 2
  class adjacent to E10 (540 labelled neighbours) has a black descent avoiding E10, so the black distance is lambda for every
  other class; exc10 analysis: E10 has 657 Hamiltonian paths, 558 black, 9 descent words all blue and closed, |Aut| = 9,
  9 anti-automorphisms, 108 black moves keeping lambda = 1 to classes with black descents (black distance 2). Family
  E(M; Z): no black descent for the nine members through order 16 that were built, black distance exactly 2 (two-move
  certificates) for the six members through order 14; the variant with first block K1 is not exceptional (three minimum
  cuts, 18 black descent words). Strong classes by lambda: n=5: 5, 1; n=6: 30, 5; n=7: 266, 84, 3; n=8: 3623, 2300, 85;
  n=9: 78515, 89366, 10237, 15; n=10: 2817450, 5124698, 1400468, 13333 (lambda = 1, 2, ...); the top values are the
  (semi-)regular counts of OEIS A255599; the 2-arc-strong counts 1, 5, 87, 2385, 99618, 6538499 (n = 5..10) are not in OEIS.
  dipole_law_check.py: (A) on 2940 random path complements (n = 4..31), no failure.
  CITED / CREDITED: Camion 1959; THM-4729 (lines, blue/black), THM-4730 / blue_quotient_20261020b (pure-blue C3, R5; black
  closed paths), black_descent_20261020c eq. (1) (the c3 change in the path-reversal convention), THM-4810 Theorem L (the
  score form of the dipole law), the concurrent regular basin theorem of the aliquot_power_124 session
  (scalar_root_reparse_20261026/basins, independently audited there, uncommitted at filing: exactly m black moves from a
  regular (2m+1)-tournament to a singleton plus a strong core; its cut-charge bound is the case lambda = m of (B)), and the
  concurrent dyadic-module-spectrum session's matrix form K' = -K + 2A_P; klein's THM-4814 (D) (concurrent, pushed before
  this file: the lower bound m for regular sources and all-black geodesic witnesses at orders 7, 9 and for P_11; HYP-9422)
  and Session75 (black-core audit rounds announced).
  PROVED (D1): if lambda(T) = delta^0(T) the distance to the singleton locus is delta^0 (either colour).
  FINITE-EXACT (D2): every strong class of orders 6-9 is exactly delta^0 black moves from a singleton plus a strong core;
  (E): the black core is connected at orders 4-9 (co-basin reduction).
  OPEN: the classification of strong tournaments without a black descent (HYP-9421 (1): exactly C3, R5 and E(M; Z)); whether
  the black distance is always at most lambda + 1 (HYP-9421 (2)); black-core connectivity (HYP-9348 (b)).
session: mac-mini-2026-10-27-potentials
source: 05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md
depends_on:
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
related:
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
  - 05-knowledge/results/black_descent_20261020c/cycles.md
scripts:
  - 04-computation/experiments/black_terminal_distance_20261027/black_descent.cpp (+ black_descent_n3_9.out, black_descent_n10.out)
  - 04-computation/experiments/black_terminal_distance_20261027/black_descent_profile.cpp, score_sum_descent.cpp, all_black_cycle.cpp (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/exceptional_family.py (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/dipole_law_check.py (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/singleton_distance.cpp (+ singleton_distance_n6_8.out, singleton_distance_n9.out)
  - 04-computation/experiments/black_terminal_distance_20261027/cobasin_connectivity.cpp (+ cobasin_connectivity_n4_9.out)
  - 04-computation/experiments/black_terminal_distance_20261027/locus_distances.cpp, lambda_vs_delta0.cpp (+ .out)
---

# THM-4811 — arc connectivity is the exact distance to the non-strong locus

## 0. Conventions

`T` is a tournament on `V`, `|V| = n`, `s(v)` the out-degree, `K` the skew matrix (`K_uv = 1` if `u -> v`, `-1` if
`v -> u`). A legal word `Q = (q_1, ..., q_n)` is a directed Hamiltonian path; the native complement `tau_Q T` keeps the
arcs of `Q` and reverses every other pair (THM-4810 (C); THM-4729's lines up to converse). `Q` is blue when
`q_i -> q_(n+1-i)` is an anti-automorphism of `T`, black otherwise (THM-4729); a blue word reverses out-degrees, so its
end-point scores sum to `n - 1`. For `X` a nonempty proper subset put `f(X) = e(X -> V\X)`; the arc-strong connectivity
is `lambda(T) = min_X f(X)`, which is 0 exactly when `T` is not strong. A *minimum out-cut* is an `X` with `f(X) = lambda`.

## (A) The dipole law

Let `A_Q` be the skew matrix of the directed path `Q`. Off the path the move negates `K`; on the path it restores `K`.
So `K' = -K + 2A_Q`. Since `A_Q 1 = e_a - e_b` (`a = q_1`, `b = q_n`), the deviations `d = K1/2 = s - (n-1)/2` obey
`d' = -d + e_a - e_b`. For a cut: the arcs of `Q` keep their direction and every other pair between `X` and `V\X` is
reversed, and a path from `a` to `b` crosses from `X` to `V\X` exactly `[a in X] - [b in X]` more times than back, so

    f_(tau_Q T)(X) = f_T(V\X) + [a in X] - [b in X].                                  (A)

Finally `c3 = C(n,3) - sum C(s_v, 2)` is an affine function of `sum d_v^2`, and `sum d'^2 = sum d^2 + 2(d_b - d_a + 1)`
gives `c3(tau_Q T) - c3(T) = s(a) - s(b) - 1` (black_descent_20261020c (1), there in the path-reversal convention, which
differs from ours by the converse and has the same `c3`).

## (B) Arc connectivity is the exact distance

By (A), `lambda(tau_Q T) = min_X (f_T(V\X) + [a in X] - [b in X]) >= lambda(T) - 1`, and applying this to the inverse move
(`Q` is legal in `tau_Q T` and `tau_Q tau_Q T = T`) gives `|lambda(tau_Q T) - lambda(T)| <= 1`. Equality `lambda - 1` holds
iff `f_T(V\X) = lambda` for some `X` with `a` outside `X` and `b` inside, i.e. iff some minimum out-cut `Y = V\X` contains
the start and misses the end of `Q`. Such a *descent word* exists whenever `T` is strong: a Hamiltonian cycle (Camion)
visits both `Y` and its complement, so it has an arc `u -> w` with `u` outside `Y` and `w` in `Y`, and its cut starting at
`w` ends at `u`. Iterating, `T` reaches the non-strong locus in exactly `lambda(T)` moves, and no fewer.

Since `f({v}) = s(v)` and `f(V\{v}) = n - 1 - s(v)`, `lambda <= min(delta^+, delta^-) <= floor((n-1)/2)`, with equality in
the last step iff `T` is regular (`n` odd) or near-regular (`n` even); conversely `f(X) = k(n-k)/2 >= (n-1)/2` for every
`X` of size `k` in a regular tournament, and the analogous count holds for near-regular ones (census: the classes with
maximal `lambda` are exactly A255599's). Moreover `lambda = min(delta^+, delta^-)` whenever `n <= 4 delta^0 + 1`: a cut
`X` of size `k >= 2` has `f(X) >= k delta^+ - C(k,2)` and `f(X) >= (n-k) delta^- - C(n-k,2)`, and `f(X) < delta^0` forces
`k >= 2 delta^0 + 1` and `n - k >= 2 delta^0 + 1`. So every strong tournament of order at most 9 has `lambda = delta^0`
(census: no exception through order 9), and the first separation occurs at order 10 (two copies of `R5` joined by all
arcs one way except one: `lambda = 1`, `delta^0 = 2`). For a regular `(2m+1)`-tournament the distance is `m`: the concurrent regular
basin theorem's cut-score bound is the case `lambda = m` of (B), and (B) shows the bound is attained from every source.

## (C) Black moves

**C1 (finite census).** Through order 10 every strong class except the pure-blue `C3` and `R5` has a black descent word
except one class `E10` of order 10 (code `111001111101111111111111111111111111000111101`). Every class adjacent to `E10`
with `lambda = 2` still has a black descent into another class, so the black distance equals `lambda` for every strong
class of order at most 10 other than `C3`, `R5`, `E10`, and equals 2 for `E10`.

**C2 (the family `E(M; Z)`).** Let `M`, `Z` be `K1`/`C3`-stacks (ordered joins of `K1` and `C3` blocks), `M`'s first block
a `C3`, `Z` nonempty. Put `D = C3[x0, M, Z]` (`x0 -> M => Z => x0`), let `D*` be a copy of the converse of `D` with
`y0 = x0*`, and let `E(M; Z)` have all arcs `D -> D*` except the single arc `y0 -> x0`. Then:

1. `E` is strong: `x0` starts Hamiltonian paths of `D`, every vertex of `D*` reaches `y0`, and `y0 -> x0`.
2. `D*` is the unique minimum out-cut and `lambda = 1`. Inside `D*` a proper subset `X~*` would need exactly one arc
   entering `X~` in `D`; since `x0` dominates `M`, `M` dominates `Z` and `|M| >= 3`, that forces `X~` to be a source of
   `M`, impossible because `M` begins with `C3`. Cuts meeting `D` lose at least `|X n D| |D* \ X| - 1` arcs to `D*` plus
   at least one arc inside the strong `D`; the case analysis leaves only `X = D*`.
3. A word starting in `D*` and ending in `D` crosses `D* -> D` once, through `y0 -> x0`, and never returns, so it is
   (Hamiltonian path of `D*` ending at `y0`) `y0 -> x0` (Hamiltonian path of `D` from `x0`). The paths of `D` from `x0` are
   `x0, P_M, P_Z` with `P_M`, `P_Z` Hamiltonian paths of the stacks; `Aut(M) x Aut(Z)` acts simply transitively on them
   (THM-4741: `tc = 1` for stacks). The positional reversal of every descent word is therefore the converse identification
   `D* -> D` composed with an automorphism of `D`, together with its inverse; the cross arcs are respected because the
   reversal swaps `x0` and `y0`. So every descent word is blue and `E(M; Z)` has black distance at least 2.

The order-10 exception is `E(C3; K1)`, the only member of order 10. At order 12 the members are `E(C3; K1 => K1)` and
`E(C3 => K1; K1)`.

## (D) Singleton targets

Let `delta^0(T) = min_v min(s_v, n - 1 - s_v)`. By (A) each `|s_v - (n-1)/2|` changes by at most one per move, so
`delta^0` is a 1-Lipschitz potential for the distance to the singleton locus `SING` (some source or sink), and also for
the smaller locus `CORE` (a source or sink whose deletion is strong).

**D1 (PROVED).** If `lambda(T) = delta^0(T) = k`, then `d(T, SING) = k` with moves of either colour. Protect a vertex
`v` attaining `delta^0`, say `s(v) = k` (otherwise pass to the converse). Make `v` the start of the first word, the end
of the second, and so on: `s(v)` runs `k, n - k, k - 2, n - k + 2, ...`, so `|s(v) - (n-1)/2|` grows by one each step
and `v` is a sink or source after `k` steps. Each intermediate tournament `T_j` (`j < k`) has `lambda(T_j) >= k - j >= 1`
by (B), so it is strong and has Hamiltonian paths starting and ending at `v` (Camion). This applies to every strong
tournament of order at most 9 and to every regular or near-regular tournament (section (B)).

**D2 (FINITE-EXACT).** Every strong class of orders 6-9 (35, 353, 6,008, 178,133) reaches `CORE` in exactly `delta^0`
**black** moves (`singleton_distance.cpp`: exhaustive depth-first search over the moves that raise the maximal
deviation, which are the only moves on a route of length `delta^0`). At order 5 only `R5` fails (no black move). For a
regular `(2m+1)`-tournament this is the concurrent basin theorem's "exactly `m` black moves from a singleton plus a
strong core" (all orders there), here extended to all strong sources at these orders; HYP-9421 (6) conjectures it in
general.

## (E) Black-core connectivity through order 9

Let `B_n` be the bipartite graph between the `lambda = 1` strong classes and the non-strong classes (converse-merged),
an edge for each black descent word. Every black line at a non-strong class is such an edge (by (B) its other end has
`lambda = 1`), and by (B)-(C) every non-pure-blue strong class of order at most 10 is joined to a `lambda = 1` class by
black moves. So if `B_n` is connected, the black core of order `n` is connected. **FINITE-EXACT**
(`cobasin_connectivity.cpp`): `B_n` is connected for `n = 4, ..., 9` and covers every non-strong class with a black line
(440 of 443 merged non-strong classes at order 8, 6,741 of 6,747 at order 9; 4,216,712 black descent lines at order 9).
This reproduces THM-4731's black connectivity through order 8 by a different route and extends it to order 9.

## (F) What remains

Whether `C3`, `R5` and `E(M; Z)` are all the strong tournaments without a black descent, and whether
the black distance never exceeds `lambda + 1`, are HYP-9421 (1)-(2). Black-core connectivity (HYP-9348 (b)) is not
claimed: (B)-(C) connect each strong class to the non-strong locus, not the non-strong basins to each other.
