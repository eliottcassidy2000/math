---
id: THM-4811
title: "Arc-strong connectivity is the exact native-complement distance to the non-strong locus, and its black refinement is classified. (A) Dipole law: a native move along a legal word Q = (a, ..., b) sends the skew matrix K to -K + 2A_Q, the score deviations d to -d + e_a - e_b, and the cut function f(X) = e(X -> V\\X) to f(V\\X) + [a in X] - [b in X]; hence c3' - c3 = s(a) - s(b) - 1. (B) For every tournament, the arc-strong connectivity lambda changes by at most one per move, it drops iff the word starts in a minimum out-cut and ends outside it, and the distance to the non-strong locus by moves of either colour is exactly lambda(T); lambda = delta^0 below order 10. (C) Black moves only (Theorems A, B): a strong tournament on n >= 4 vertices has no black descent iff it is R5 or a member of the explicit family E = {E_t(S, x0)} (bottom layer S, t - 1 middle vertices, converse top layer, one upward jump chain), which has members at every order n >= 10 (1, 1, 3, 3, 8, 8, 20, 20 at n = 10..17); the black distance to the non-strong locus is lambda(T) for every strong T other than R5 (no black move) and the members of E (black distance 2 = lambda + 1). (D) Singleton targets: distance delta^0 when lambda = delta^0 (black moves suffice); every strong class of orders 6-10 is exactly delta^0 black moves from a singleton plus a strong core. (E) The black core is connected through order 10."
status: >
  PROVED (elementary): (A); (B) given Camion's theorem (proof inherited in blue_quotient_20261020b/structure.md, section 2);
  (C) Theorems A and B, all orders, given the PROVED repo inputs cited in classification/PROOF.md (THM-2477 via
  Benhocine-Wojda, LEM-003, the blue-cut coset lemma and all-blue constructions of blue_quotient_20261020b/structure.md, the
  prime-orbit block normal form and its consequences in black_descent_20261020c/cycles.md sections 3-7, the uniform
  prefix-rerouting theorem for odd circulants, Camion) and one bounded finite check (the last-block lemma in the single case
  of a 6-vertex block with a 3-fold rotation: all 35 strong 6-tournaments and all subsets, re-done by an independent checker
  for every block size <= 9); (D1) and its black version (Theorem C of the classification); (C2) my earlier subfamily proof.
  The classification proof was written in this session by a dedicated proof agent and reviewed by a second agent with its
  own code (no logical error found); not yet independently audited at the repository level.
  FINITE-EXACT: black_descent.cpp over every strong class of orders 3-10 (1, 1, 6, 35, 353, 6,008, 178,133, 9,355,949): the
  descent criterion of (B) verified on all legal words of all strong classes through order 8 (2,073,672 words at order 8,
  no failure; independent audit: every class, strong or not, through order 9, 2,143,304 words at order 8 and 270,447,284
  at order 9); every strong
  class has a descent word, and a black one except C3, R5 (pure blue) and E10 = E_1(S5); the classification checkers
  (classification/): the order-11 census of all 884,464,590 strong classes (exactly one exception, E_2(S5)), bottom-layer
  census for |S| <= 9 (1, 2, 5, 12, 28), all even shapes through order 14 and odd shapes through order 15 (exceptions
  1, 3, 8 = e(n)), exhaustive lambda = 2 < delta^0 configurations at order 14, all 44 members of E with n <= 16 (strong,
  lambda = 1, delta = 2, cut chain of length t, no black descent, two-move certificate); my own code re-checks E_1(S5),
  E_2(S5), E_3(S5) and three order-13/14 members. Strong classes by lambda: n=5: 5, 1; n=6: 30, 5; n=7: 266, 84, 3;
  n=8: 3623, 2300, 85; n=9: 78515, 89366, 10237, 15; n=10: 2817450, 5124698, 1400468, 13333 (lambda = 1, 2, ...); the top
  values are the (semi-)regular counts of OEIS A255599; the 2-arc-strong counts 1, 5, 87, 2385, 99618, 6538499 (n = 5..10)
  are not in OEIS. dipole_law_check.py: (A) on 2940 random path complements (n = 4..31), no failure.
  FINITE-EXACT (D2): every strong class of orders 6-10 (9,355,949 at order 10, among them the 16 with lambda < delta^0) is
  exactly delta^0 black moves from a singleton plus a strong core; (E): the co-basin graph is connected at orders 4-10 (order 10: 537,339,904
  black descent lines), so the black core is connected there.
  CITED / CREDITED: Camion 1959; THM-4729 (lines, blue/black), THM-4730 / blue_quotient_20261020b (pure-blue C3, R5; black
  closed paths), black_descent_20261020c (the c3 change; the tools listed above), THM-4810 Theorem L (the score form of the
  dipole law), the concurrent regular basin theorem of the aliquot_power_124 session (scalar_root_reparse_20261026/basins,
  independently audited there, uncommitted at filing: exactly m black moves from a regular (2m+1)-tournament to a singleton
  plus a strong core; its cut-charge bound is the case lambda = m of (B)), the concurrent dyadic-module-spectrum session's
  matrix form K' = -K + 2A_P (also Session75 RESULT (6.1), on main and audited before this file, with the path spectrum
  and the trace-norm tariff (6.2)-(6.3)), klein's THM-4814 (D) (pushed before this file: the bound m for regular sources and all-black
  witnesses at orders 7, 9 and for P_11; HYP-9422), Session75 (black-core audit rounds announced).
  AUDITED (independent, 2026-10-27; PASS WITH CORRECTIONS, MISTAKE-681; report
  04-computation/experiments/black_terminal_distance_20261027/audit_thm4811_4813/REPORT.md): (A), (B), (D1), (D2), (E) and the
  E(M; Z) properties, with the census through order 10 reproduced; Theorems A, B, C are under a separate independent audit.
  OPEN: the singleton-plus-strong-core distance delta^0 at every order (HYP-9421 (6)); black-core connectivity at every order
  (HYP-9348 (b)); the closed extremal-predecessor form of black_descent_20261020c/cycles.md section 7 in the case of two
  minimum-score vertices.
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
  - 04-computation/experiments/black_terminal_distance_20261027/exceptional_family.py, e10_analysis.py (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/classification/ (PROOF.md, tour.hpp, keylemma_check.cpp, census_sc.cpp, bottom_layers.cpp, lastblock_check.cpp, e_check.cpp, family*.cpp, black_distance.cpp, verify_counterexample.py, results/, audit/)
  - 04-computation/experiments/black_terminal_distance_20261027/dipole_law_check.py (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/singleton_distance.cpp (+ singleton_distance_n6_8.out, singleton_distance_n9.out, singleton_distance_n10.out)
  - 04-computation/experiments/black_terminal_distance_20261027/cobasin_connectivity.cpp (+ cobasin_connectivity_n4_9.out, cobasin_connectivity_n10.out)
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
`k >= 2 delta^0 + 1` and `n - k >= 2 delta^0 + 1`. So every strong tournament of order at most 9 has `lambda = delta^0` (for `delta^0 = 1` because strong means `lambda >= 1`;
for `delta^0 >= 2` because `n <= 9 <= 4 delta^0 + 1`)
(census: no exception through order 9), and the first separation occurs at order 10 (two copies of `R5` joined by all
arcs one way except one: `lambda = 1`, `delta^0 = 2`). For a regular `(2m+1)`-tournament the distance to the non-strong locus is
`m`: the concurrent regular basin theorem's cut-score lower bound (for its smaller target, a singleton plus a strong
core) is the case `lambda = m` of (B), and (B) shows that the non-strong-locus bound `lambda(T)` is attained from every source.

## (C) Black moves

**C1 (finite census).** Through order 10 every strong class except the pure-blue `C3` and `R5` has a black descent word,
except one class `E10` of order 10 (code `111001111101111111111111111111111111000111101`); every `lambda = 2` class
adjacent to `E10` still has a black descent into another class, so through order 10 the black distance is `lambda`
except for `C3`, `R5` and `E10`, which has black distance 2. At order 11 the census of all 884,464,590 strong classes
(classification/census_sc.cpp) finds exactly one exception.

**C2 (the family `E`).** A *bottom layer* is a pair `(S, x0)` with `S - x0 = R_1 => ... => R_q` (`q >= 2`) a stack of
`K1` and `C3` blocks whose last block `R_q` is a `C3`, `x0` uniform on every block, `x0 => R_1` and `R_q => x0` (the
relation of `x0` to the blocks between is free). For `t >= 1`, `E_t(S, x0)` has a bottom layer `S`, middle vertices
`c_1, ..., c_(t-1)` and a top layer `sigma(S)` carrying the converse of `S`; the only upward arcs are the jumps
`x0 -> c_1 -> ... -> c_(t-1) -> sigma(x0)`, and every other pair between different layers points down. Its order is
`2|S| + t - 1`; the bottom layers number `[z^a] z^5 (1 + z^2)/(1 - 2z - 2z^3)` (1, 2, 5, 12, 28, 66, 156 for
`|S| = 5..11`), so the members number 1, 1, 3, 3, 8, 8, 20, 20, 48, 48 at orders 10-19. `E10 = E_1(S5, x0)` with
`S5 = (x0 -> K1 => C3 => x0)`, and the order-11 exception is `E_2(S5, x0)`. Every member is strong with `lambda = 1`,
`delta^0 = 2`, minimum cuts exactly the `t` cuts `S`, `S + c_1`, ...; its descent words are
`(P_0, c_1, ..., c_(t-1), P_t)` with `P_0` a Hamiltonian path of `S` into `x0`, and the stabiliser of `x0` (independent
rotations of the `C3` blocks) acts simply transitively on these paths, so every descent word is the converse
identification composed with an automorphism: all are blue (classification/PROOF.md, Part VII).

The family `E(M; Z)` of the first checkpoint (doubled `C3[x0, M, Z]`, written from the top layer, i.e. in the converse
convention) is exactly the subfamily `t = 1` in which `x0` beats an initial segment of the blocks; it misses the
members with middle vertices (every odd order, `E_3(S5)` at order 12) and the patterns that are not initial segments
(first at order 14).

**C3 (Theorems A and B, PROVED).** *Theorem A.* A strong tournament on `n >= 4` vertices has no black descent word iff it
is `R5` or a member of `E`. *Theorem B.* Its black distance to the non-strong locus is `lambda(T)`, except `R5` (no black
move) and the members of `E` (black distance `2 = lambda + 1`, attained by an explicit two-move certificate).

*Proof outline* (classification/PROOF.md). Non-self-converse tournaments have only black words, and regular ones other
than `C3`, `R5` have black closed paths, so an exception is self-converse and non-regular. *Theorem C:* if
`lambda = delta^0` there is a black descent of singleton type (starting at a minimum-score vertex or ending at a
maximum-score one), by blue-cut cosets, the prime-orbit block normal form and two vertex arguments; this also proves the
open-path form of the extremal-predecessor statement left open in black_descent_20261020c/cycles.md section 7. If
`lambda < delta^0`, minimum cuts cannot cross, so they form a chain, and every anti-automorphism reverses the chain's
layers. For `2 <= lambda < delta^0` an all-blue candidate set forces a stabiliser to act simply transitively, which forces
a `C3` block with impossible degrees, so there is no exception. For `lambda = 1 < delta^0`, layers cannot be skipped,
middle layers are single vertices, and a bottom-layer lemma (THM-2477 via Benhocine-Wojda) pins down `S`: the exceptions
are exactly `E`. For Theorem B, a singleton-type descent keeps `lambda = delta^0`, a descent from
`2 <= lambda < delta^0` stays in that regime and cannot enter `E`, and for `T` in `E` a black word starting and ending in
`S` keeps `lambda = 1` and reaches `lambda = delta^0 = 1`.

## (D) Singleton targets

Let `delta^0(T) = min_v min(s_v, n - 1 - s_v)`. By (A) each `|s_v - (n-1)/2|` changes by at most one per move, so
`delta^0` is a 1-Lipschitz potential for the distance to the singleton locus `SING` (some source or sink), and also for
the smaller locus `CORE` (a source or sink whose deletion is strong).

**D1 (PROVED).** If `lambda(T) = delta^0(T) = k`, then `d(T, SING) = k` with moves of either colour. Protect a vertex
`v` attaining `delta^0`, say `s(v) = k` (otherwise pass to the converse); this is the protected-vertex route of the
concurrent regular basin theorem (B3), run here with the cut potential. Make `v` the start of the first word, the end
of the second, and so on: `s(v)` runs `k, n - k, k - 2, n - k + 2, ...`, so `|s(v) - (n-1)/2|` grows by one each step
and `v` is a sink or source after `k` steps. Each intermediate tournament `T_j` (`j < k`) has `lambda(T_j) >= k - j >= 1`
by (B), so it is strong and has Hamiltonian paths starting and ending at `v` (Camion). This applies to every strong
tournament of order at most 9 and to every regular or near-regular tournament (section (B)). *Black version:* Theorem C
gives a black singleton-type descent whenever `lambda = delta^0` (`T` not `R5`), after which `lambda = delta^0` again; at
the end `delta^0 = 0`, so the black distance to `SING` is also `delta^0`.

**D2 (FINITE-EXACT).** Every strong class of orders 6-10 (35, 353, 6,008, 178,133, 9,355,949; at order 10 including the
16 classes with `lambda < delta^0`, all with `lambda = 1`, `delta^0 = 2`, among them `E10`) reaches `CORE` in exactly
`delta^0` **black** moves (`singleton_distance.cpp`: exhaustive depth-first search over the moves that raise the maximal
deviation, which are the only moves on a route of length `delta^0`). At order 5 only `R5` fails (no black move). For a
regular `(2m+1)`-tournament this is the concurrent basin theorem's "exactly `m` black moves from a singleton plus a
strong core" (all orders there), here extended to all strong sources at these orders; HYP-9421 (6) conjectures it in
general.

## (E) Black-core connectivity through order 10

Let `B_n` be the bipartite graph between the `lambda = 1` strong classes and the non-strong classes (converse-merged),
an edge for each black descent word. Every black line at a non-strong class is such an edge (its other end is strong by THM-4729 L3, so has
`lambda = 1` by (B)), and every non-pure-blue strong class is joined by black moves to a `lambda = 1` class with a black
descent: at every order already by the inherited strong black `c3`-descent (black_descent_20261020c; the last strong
class of a black descent chain has `lambda = 1` and its last move is a black descent word), independently of Theorems A
and B, and through order 10 also by (C1). So if `B_n` is connected, the black core of order `n` is connected. **FINITE-EXACT**
(`cobasin_connectivity.cpp`): `B_n` is connected for `n = 4, ..., 10` and covers every non-strong class with a black line
(440 of 443 merged non-strong classes at order 8, 6,741 of 6,747 at order 9, 188,642 of 188,646 at order 10, the rest
pure blue; 4,216,712 and 537,339,904 black descent lines at orders 9 and 10). At order 10, 2,817,449 of the 2,817,450
`lambda = 1` strong classes have a black descent, the exception being `E10`, an independent re-count of (C1); `E10` joins
`B_10` through its black moves to `lambda = 1` classes with black descents. This reproduces THM-4731's black connectivity
through order 8 by a different route and extends it to orders 9 and 10 (HYP-9348 (b) is now FINITE-EXACT through 10).

## (F) What remains

Theorems A and B settle HYP-9421 (1)-(2) (the conjectured exceptional list `C3, R5, E(M; Z)` was incomplete: the right
family is `E`). Open: the singleton-plus-strong-core distance `delta^0` at every order (HYP-9421 (6); D1 gives the
singleton locus when `lambda = delta^0`); black-core connectivity at every order (HYP-9348 (b), FINITE-EXACT through 10 by
(E)); the closed extremal-predecessor statement when two vertices have minimum score.
