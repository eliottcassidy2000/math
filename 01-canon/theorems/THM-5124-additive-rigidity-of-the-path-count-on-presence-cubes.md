---
id: THM-5124
title: "Additive rigidity of the Hamiltonian-path count on presence cubes. (1) Pair lemma: two vertices a, b that both raise the count of the rest R (each lies on an odd cycle with R) and have zero interaction, H(R+a+b) - H(R+a) - H(R+b) + H(R) = 0, exist only in one configuration: R+a+b = Y => Q => X with Q a strong 4-tournament on {x0, y0, a, b}, x0, y0 in R its two hubs and a, b its two one-triangle vertices; the four counts are c, 3c, 3c, 5c. (2) Two modules: if disjoint vertex sets W1, W2 both raise the count of the rest B and have zero interaction, then the counts are again c, 3c, 3c, 5c, exactly one vertex of each module is active, the two active vertices are in the configuration of (1) with hubs in B, and every other vertex of W1 u W2 is a singleton strong component. (3) No additive cube of dimension three: three pairwise disjoint vertex sets cannot all be active and additive; three pairwise additive active vertices are parallel switches between one hub pair and their joint count is 9c (transitive triple) or 15c (cyclic triple), never 7c. (4) In general k parallel switches carrying a tournament A give c (H(A) + D(A)), D the number of ordered two-path covers; for A transitive this is c (2^k + 1) against the additive prediction c (2k + 1)."
status: >
  PROVED (opus-2026-11-10-S90; elementary, from the odd-cycle formula H = I(Omega, 2) and the triangle lemma).
  KNOWN and used: the Moebius positivity of the presence cube (THM-4114, eq. (7)-(9)) - the interaction of two vertex
  sets is the weight of the odd-cycle packings meeting both, hence >= 0; Claim A (THM-4876 (1), THM-070); the
  universal increment floor H(Z) = h or H(Z) >= (3h + 1)/2 for induced extensions (Session214 forced shuffle,
  HYP-9469 (1) PROVED), which gives a second, three-line proof of the dimension bound in (3); the count 2^k + 1 of
  the reversed-extreme-arc tournaments (THM-4749, THM-4876 (6)); statement 4 is the full-cut case of the inherited
  protected-closure identity (THM-4114 section 6, THM-5122) and D is Session197's count of labelled two-bin covers.
  COMPLEMENT, not a duplicate, of the orientation side: on the EAR cube of one core (single cut bits) every pair
  interaction is strictly negative (THM-4114; Session192 cut_sumsets: no additive 2-cube); along EXCHANGE directions
  (pairs of cut bits swapped) interactions are signed, THM-5122 (Session246) gives the affine criterion, the
  capacity bound k <= min(n_-(W), floor(n/2)) and one 2-cube, and THM-5117 (Session245) gives affine exchange cubes
  of every dimension k at order 7k + 1; the present theorem is about the PRESENCE cube (induced subtournaments),
  where every interaction is nonnegative. FINITE-EXACT controls
  (rigidity_check.py): statement (1) on every labelled tournament of orders 3-7 (the hypothesis and the structure
  hold on the same 24, 240, 2160, 23520 pairs at orders 4-7); statement (2) on every ordered partition of every
  labelled tournament of orders 3-6 (38,880 instances at order 6, no violation); statement (3) on every module
  triple of orders 3-6 (11,468,800 triples at order 6) and every vertex triple over every base at orders 3-7; random
  samples at orders 8-10. AUDITED (independent adversarial audit, own code,
  05-knowledge/results/s90_rigidity_redei_20261110/audit/AUDIT.md): PASS WITH CORRECTIONS - every proof step and
  number confirmed, statements 2 and 3 verified exhaustively at order 7 as well (1,095,360 module pairs; 3.57 * 10^9
  module triples); corrections applied (the description of THM-5122 / THM-5117 and of the ear cube, the empty case
  of statement 4, credits, wording; MISTAKE-752). The Chan-Pak conjecture (every odd number other than 7 and 21 is
  a Hamiltonian-path count) remains OPEN; this theorem removes one family of proof architectures; it does not
  bound the spectrum and does not reprove that 7 is not a value.
session: opus-2026-11-10-S90
source: 05-knowledge/results/s90_rigidity_redei_20261110/README.md
scripts:
  - 04-computation/experiments/s90_rigidity_redei_20261110/rigidity_check.py (+ rigidity_check.out)
related:
  - 01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md
  - 01-canon/theorems/THM-4876-* (Claim A, increments), 05-knowledge/hypotheses/HYP-9469-* and 05-knowledge/results/session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md
  - 01-canon/theorems/THM-5122-shared-core-exchange-intervals-and-negative-inertia-capacity.md (exchange cubes)
  - 01-canon/theorems/THM-4960-*, THM-5000-* (moving-core windows), THM-5117-* (sparse affine cubes)
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-4094-*, THM-4749-*
  - Grinberg-Stanley, The Redei-Berge symmetric function of a directed graph, arXiv:2307.05569 (odd-cycle formula)
---

# THM-5124 - additive rigidity of the Hamiltonian-path count on presence cubes

## Setting

`H(T)` is the number of Hamiltonian paths of a tournament `T` (`H` of the empty tournament is `1`). For `U <= V(T)`,
`T[U]` is the induced subtournament. An *odd packing* is a set of pairwise vertex-disjoint odd directed cycles;
the odd-cycle formula (Grinberg-Stanley; THM-002) is `H(T) = sum over odd packings P of 2^|P|`. Two consequences
are used throughout.

- **Claim A** (THM-4876 (1)). `H(T) - H(T - v) = 2 sum over odd cycles C through v of H(T - V(C)) >= 0`. In
  particular `H(T) > H(T - v)` iff `v` lies on an odd cycle iff `v` lies on a triangle (a shortest cycle through a
  vertex of a tournament is a triangle: if `v -> c1 -> c2 -> ... -> v` has length at least `4`, then either
  `c2 -> v`, a triangle, or `v -> c2`, a shorter cycle).
- **Interaction formula** (THM-4114, eq. (9)). For disjoint `W1, W2 <= V(T)` and `B = V(T) - W1 - W2`,

      D(W1, W2) := H(T) - H(T[B u W1]) - H(T[B u W2]) + H(T[B])
                 = sum over odd packings P of T that meet both W1 and W2 of 2^|P|  >= 0.

  So `D(W1, W2) = 0` iff no odd packing of `T` meets both sets, i.e. (a) no odd cycle meets both, and (b) every
  odd cycle meeting `W1` intersects every odd cycle meeting `W2`. If `D(W1, W2) = 0` then also `D(S1, S2) = 0`
  inside `T[B u S1 u S2]` for all `S1 <= W1`, `S2 <= W2`.

A vertex set `W` is *active over* `B` if `H(T[B u W]) > H(T[B])`. `ST_4` denotes the strong tournament of order
`4` (unique up to isomorphism, `H = 5`): it has two triangles, sharing an arc `x0 -> y0`; `x0, y0` are its *hubs* and
the other two vertices lie on one triangle each. `Y => Q => X` is the order join (all arcs from `Y` to `Q u X` and
from `Q` to `X`).

## Statements

1. **Pair lemma.** Let `Z` be a tournament, `a, b` two vertices, `R = V(Z) - {a, b}`. The following are
   equivalent.
   - (i) `H(Z - b) > H(Z[R])`, `H(Z - a) > H(Z[R])` and `H(Z) - H(Z - a) - H(Z - b) + H(Z[R]) = 0`.
   - (ii) There are `x0, y0` in `R` and a partition `R - {x0, y0} = Y u X` such that `Z = Y => Q => X`, where
     `Q = Z[{x0, y0, a, b}]` has `x0 -> y0`, `y0 -> a`, `y0 -> b`, `a -> x0`, `b -> x0` (so `Q` is `ST_4` with hubs
     `x0, y0`; the arc between `a` and `b` is arbitrary).

   In that case, with `c = H(Z[R]) = H(Z[Y]) H(Z[X])`: `H(Z - b) = H(Z - a) = 3c` and `H(Z) = 5c`.
2. **Two modules.** Let `V(T) = B u W1 u W2` (disjoint), `c = H(T[B])`, and suppose `W1` and `W2` are both active
   over `B` and `D(W1, W2) = 0`. Then `H(T[B u W1]) = H(T[B u W2]) = 3c` and `H(T) = 5c`. Moreover there are
   unique `a` in `W1`, `b` in `W2` active over `B`; `T[B u {a, b}]` has the structure (ii) with hubs in `B`; and
   every vertex of `(W1 u W2) - {a, b}` lies on no cycle of `T` (it is a singleton strong component).
3. **No additive cube of dimension three.** Let `V(T) = B u W1 u ... u Wk` (disjoint, `k >= 3`) and
   `f(S) = H(T[B u U_{i in S} W_i])`. It is impossible that `f(S) = f(empty) + sum_{i in S} d_i` for all
   `S <= [k]` with every `d_i > 0`.
   If three vertices `a1, a2, a3` outside `B` are pairwise in the situation of statement 1 over `B` (each pair
   active with zero interaction inside `T[B u {a_i, a_j}]`), then they share one hub pair `x0, y0` in `B`,
   `T[B u {a1, a2, a3}] = Y => P => X` with `P` on `{x0, y0, a1, a2, a3}`, and the eight counts are
   `c (1; 3, 3, 3; 5, 5, 5; 9)` if `{a1, a2, a3}` is transitive and `c (1; 3, 3, 3; 5, 5, 5; 15)` if it is cyclic.
   The value `7c` predicted by additivity does not occur.
4. **Parallel switches (two-hub suspension).** For a nonempty tournament `A` let `S(A)` be `A` together with two
   further vertices `x0, y0`, `x0 -> y0`, `y0 -> A`, `A -> x0`. Then

       H(S(A)) = H(A) + D(A),      D(A) = sum over U <= V(A) of H(A[U]) H(A - U)

   (`D` counts ordered pairs of disjoint paths covering `V(A)`, either possibly empty: Session197's labelled
   two-bin covers; for empty `A` the count is `1`, not `2`). This is the full-cut case of the inherited
   protected-closure identity: `S(A)` is the closure of `A` with every vertex selected, and
   `F(1) = h + 2 a.1 - sum_(i<j) W_ij = H(A) + D(A)` (THM-4114 section 6, THM-5122). For `A` transitive of order `k`
   this is `2^k + 1` (the reversed-extreme-arc tournament `R_(k+2)`), for the 3-cycle it is `15`. The presence cube
   of `k` parallel switches over `B = Y u {x0, y0} u X` is `S |-> c (H(A[S]) + D(A[S]))` for nonempty `S` and `c` at
   `S` empty; it is additive iff `k <= 2`.

## Proofs

**Statement 1, (ii) implies (i).** In `Z = Y => Q => X` the odd cycles through `a` or `b` lie in `Q`; they are the
triangles `(a, x0, y0)` and `(b, x0, y0)`, which intersect, so the interaction vanishes. By multiplicativity over
order joins, `H(Z[R]) = H(Y) H(X)` (the pair `x0 -> y0` is transitive), `H(Z - b) = H(Y) * 3 * H(X)` and
`H(Z) = 5 H(Y) H(X)`.

**Statement 1, (i) implies (ii).** By the interaction formula, (a) no odd cycle of `Z` contains both `a` and `b`, and
(b) every odd cycle through `a` meets every odd cycle through `b`. By Claim A, `a` lies on a triangle of `Z - b` and
`b` on a triangle of `Z - a`. Assume `a -> b` (otherwise exchange the names). For `r` in `R` the triple `{a, b, r}`
is not cyclic by (a), so `b -> r -> a` never happens. Hence `R = Y0 u M u X0` with

    Y0 = {r : r -> a and r -> b},   M = {r : a -> r -> b},   X0 = {r : a -> r and b -> r}.

*Step 1: `M` is empty.* Suppose `m` in `M`. If `x -> y` for some `x` in `X0`, `y` in `Y0`, then
`a -> m -> b -> x -> y -> a` is a 5-cycle through `a` and `b`, contradicting (a); so there is no arc from `X0` to `Y0`.
A triangle through `a` in `Z - b` is `a -> x -> y -> a` with `x` in `M u X0` and `y` in `Y0`, hence `x = x1` in `M`,
`y = y1` in `Y0`, `x1 -> y1`. A triangle through `b` in `Z - a` is `b -> u -> v -> b` with `u` in `X0` and `v` in
`Y0 u M`, hence `u = u1` in `X0`, `v = v1` in `M`, `u1 -> v1`. If `v1 = x1`, then `a -> b -> u1 -> x1 -> y1 -> a` is a
5-cycle through `a` and `b`, contradicting (a). If `v1 != x1`, the triangles `(a, x1, y1)` and `(b, u1, v1)` are
disjoint, contradicting (b). So `M` is empty: `a` and `b` have the same neighbours in `R`, `Y0 -> {a, b} -> X0`.

*Step 2: one back arc.* Let `F` be the set of arcs `x -> y` with `x` in `X0`, `y` in `Y0`. The triangles through
`a` in `Z - b` are `(a, x, y)` with `(x, y)` in `F`, and likewise for `b`; so `F` is nonempty.
- (C2) There is no path `x -> z -> y` with `x` in `X0`, `y` in `Y0`, `z` in `R`: otherwise
  `a -> b -> x -> z -> y -> a` is a 5-cycle through `a` and `b`.
- (C3) Any two arcs of `F` share a vertex: otherwise `(a, x, y)` and `(b, x', y')` are disjoint.
Let `(x0, y0)` in `F` and suppose `(x, y)` in `F` is a different arc; by (C3) it shares a vertex. If `x = x0`, then
`y != y0` are two vertices of `Y0`, and `x0 -> y -> y0` or `x0 -> y0 -> y` is a path forbidden by (C2). If `y = y0`,
then `x != x0` are two vertices of `X0`, and `x -> x0 -> y0` or `x0 -> x -> y0` is forbidden by (C2). So
`F = {(x0, y0)}`, and every other arc between `Y0` and `X0` points from `Y0` to `X0`. By (C2) again no vertex of
`X0` beats `x0` (else `x -> x0 -> y0`) and `y0` beats no vertex of `Y0` (else `x0 -> y0 -> y`).
With `Y = Y0 - {y0}`, `X = X0 - {x0}` and `Q = {x0, y0, a, b}` this is exactly `Z = Y => Q => X` with the stated arcs
inside `Q`. QED.

**Statement 2.** Write `f(S1, S2) = H(T[B u S1 u S2])`. Since `D(W1, W2) = 0`, no odd packing of `T` meets both
modules, so `f(S1, S2) = f(S1, empty) + f(empty, S2) - c` for all `S1 <= W1`, `S2 <= W2`. By Claim A, `f` is
nondecreasing in each argument. Call `(S, v)` an *active step* of module `i` if `v` in `S <= W_i` and the count
rises when `v` is added to `B u (S - v)`; active steps exist in both modules because both are active over `B`.

Let `(S1, a)` and `(S2, b)` be active steps of modules `1` and `2`, and `R = B u (S1 - a) u (S2 - b)`. By the
additivity just noted, `a` and `b` are active over `R` and have zero interaction in `Z = T[R u {a, b}]`. Statement 1
gives the structure (ii) with hubs `x0, y0` in `R`. A hub cannot lie in `W1 u W2`: if, say, `x0` in `W1`, the
triangle `(b, x0, y0)` would meet both modules. So `x0, y0` are in `B`, and the triangle `(a, x0, y0)` shows that `a`
is active over `B`; likewise `b`.

Statement 1 also gives `f(S1, empty) - f(S1 - a, empty) = 2 H(T[R]) = 2 (f(S1 - a, empty) + f(empty, S2 - b) - c)`.
The left side does not depend on `(S2, b)`, so `f(empty, S2 - b)` takes the same value `beta` at every active step
of module `2`. Along a chain from the empty set to `W2` the first active step starts from the value `c`; hence
`beta = c`, and by symmetry `f(S1 - a, empty) = c` at every active step of module `1`. Therefore every active step
of either module leads from `c` to `c + 2 (c + c - c) = 3c`, and along any chain there is at most one active step:
`f(S1, empty)` and `f(empty, S2)` take only the values `c` and `3c`, with `f(W1, empty) = f(empty, W2) = 3c` and
`f(W1, W2) = 5c`.

Uniqueness: if `a, a'` in `W1` were both active over `B`, apply the above to the active steps `({a}, a)`, `({a'}, a')`
and `({b}, b)`. In `T[B u {b}]` the vertex `b` lies on exactly one triangle (it is a one-triangle vertex of `Q`), so
both structures have the same hubs `x0, y0`; then `y0 -> {a, a', b} -> x0`, and the hubs together with `a, a', b`
taken along a Hamiltonian path of `T[{a, a', b}]` form a 5-cycle meeting both modules, a contradiction. Finally let
`z` in `W1 - {a}`. Since `a` lies in `W1 - z` and `f({a}, empty) = 3c`, monotonicity gives
`3c <= f(W1 - z, empty) <= f(W1, empty) = 3c`; so adding `z` to `B u (W1 - z)` does not raise the count, and Claim A shows that `z` lies on no odd cycle of
`T[B u W1]`; no odd cycle through `z` meets `W2`; so `z` lies on no odd cycle of `T`, hence on no cycle. QED.

**Statement 3.** *First proof.* Suppose the cube is additive with all `d_i > 0`, `k >= 3`. Then `W1` and
`W' = W2 u ... u Wk` are active over `B` with zero interaction, so statement 2 gives `H(T[B u W']) = 3c`. The pair
`W2, W3` is also active over `B` with zero interaction inside `T[B u W2 u W3]`, so `H(T[B u W2 u W3]) = 5c`; by
monotonicity `5c <= H(T[B u W']) = 3c`, absurd.
*Second proof (from the increment floor).* Put `D = sum d_i`. For each `i`, `T[B u U_{j != i} W_j]` is an induced
subtournament of `T` with count `c + D - d_i < c + D`, so the forced-shuffle floor (HYP-9469 (1), Session214, eq. (8))
gives `c + D >= (3 (c + D - d_i) + 1)/2`, i.e. `3 d_i >= c + D + 1`. Summing over `i`: `(3 - k) D >= k (c + 1) > 0`,
so `k <= 2`. (For `k = 2` the same inequalities give `d_1, d_2 >= (c + 1 + max(d_1, d_2))/2`, consistent with
`d_1 = d_2 = 2c`.)
*Three pairwise additive vertices.* For each pair statement 1 applies inside `T[B u {a_i, a_j}]` with `R = B`. The
hub pair of a structure containing `a_i` is determined by `B` and `a_i` (it is the unique triangle through `a_i` in
`T[B u {a_i}]`), so the three structures share `x0, y0`, and `B - {x0, y0} = Y u X` with `Y => {x0, y0, a_i} => X` for
each `i`. Hence `T[B u {a1, a2, a3}] = Y => S(A) => X` with `A = T[{a1, a2, a3}]`, and statement 4 gives `9` for the
transitive triple and `15` for the cyclic one. QED.

**Statement 4.** In `S(A)` the vertex `y0` has the single in-neighbour `x0` and `x0` has the single out-neighbour
`y0`. A Hamiltonian path either starts at `y0` - then `x0` cannot be followed by anything, so the path is `y0`, a
Hamiltonian path of `A`, `x0`: `H(A)` paths - or contains `x0 y0` consecutively: a path on some `U <= V(A)`, then
`x0 y0`, then a path on `V(A) - U`, with no condition at the two junctions: `D(A)` paths. For `A` transitive,
`H(A[U]) = 1` for all `U`, so `D = 2^k`. The cube formula follows from multiplicativity over `Y => S(A[S]) => X`, and
`H(S(A[S]))` for `|S| = 0, 1, 2` is `1, 3, 5`; for `|S| = 3` it is `9` or `15`, not `7`. QED.

## Remarks

1. **What the theorem says about 5 and 7.** `5 = 1 + 2 + 2` is the only additive identity among Hamiltonian-path
   counts of nested induced subtournaments, up to the scale factor `c`. A third parallel switch always creates the
   5-cycle `x0 -> y0 -> a1 -> a2 -> a3 -> x0` through all three (or more, if the switches are cyclic), which turns the
   additive prediction `7` into `9` or `15`. This is consistent with, and does not prove, the base-free obstruction
   "no tournament has conflict polynomial `1 + 3x`" (THM-4940 (3)); the exclusion of `7` and `21` from the spectrum
   is THM-1370 and is not reproved here.
2. **Three cubes, three signs.** For one new vertex over a fixed core the count is the quadratic
   `F(x) = h + 2 a^T x - (1/2) x^T W x` of the cut (THM-4114 section 6, THM-5122).
   - *Presence cube* (which vertices are present): interactions are packing weights, `>= 0`, and vanish between two
     active sets only in the configuration of statement 1.
   - *Ear cube* (single cut bits): the interaction of two bits is `-W_ij <= -2 H(C - i - j) < 0`, never zero
     (nonpositivity: THM-4114; strictness: THM-5122, capacity paragraph, and Session124; Session192
     `cut_sumsets/PROOF.md`: no additive cube of dimension two).
   - *Exchange directions* (a selected and an unselected vertex swapped): interactions `-d_i^T W d_j` are signed and
     can cancel. THM-5122 gives the affine criterion, the capacity bound `k <= min(n_-(W), floor(n/2))` and one
     2-cube (`345, 347, 347, 349`); THM-5117 gives affine exchange cubes of every dimension `k` at order `7k + 1`,
     with spacing `50 * 37^(k-2)`. The capacity bound is not attained in general (audit: no core of order 6 carries
     three pairwise orthogonal disjoint exchanges with nonzero increments; order 7 does).
   So additivity of `H` beyond `1, 3, 3, 5` is available only along exchange directions: it is an orientation
   phenomenon, never a presence phenomenon.
3. **Proof architectures this removes.** Any construction that asks for three or more vertex sets whose presence
   contributes independently to `H` - complete sequences in the sense of Brown's criterion, Erdos-Lewin antichain
   sums of 3-smooth numbers, sumsets of two independent rich families (Mann/Schnirelmann), affine presence cubes -
   is impossible; statement 2 shows that even two independent modules realise only `c, 3c, 5c`. Constructions must
   use interacting digits (carries), exchanges (THM-5122), or non-nested families.
4. **Width of an orthogonal exchange cube (corollary of THM-5122 and Busch's floor; scope: strong closures).** An
   interaction-orthogonal family of `k` disjoint pair exchanges in a core of order `n` has `k <= floor(n/2)`
   (THM-5122, capacity); if it is interval-complete its `2^k` counts fill an interval of at most `2^k` odd numbers
   (THM-5122, affine interval criterion; a unit-step interacting cube has at most `k + 1` values); if the closures
   are strong, of order `n + 2`, each count is at least the strong floor
   `f(n + 2) >= 5^((n+1)/3)` (Busch 2006; `f = 3, 5, 9, 15, 25, 45, 75` at orders `3-9`). Hence the interval has
   relative length at most `2^(floor(n/2)+1) / 5^((n+1)/3) < 1.17 * 0.83^n`: an interval at baseline `B` obtained this
   way has length at most about `B^0.65`. Exchange cubes with strong closures can cover a constant-ratio window
   only as a union of polynomially many cubes with different baselines.

## Non-consequences

The theorem does not bound the spectrum, does not prove that any integer is or is not a value, and says nothing
about non-induced operations (arc reversals, module replacement) beyond the cited THM-5122. The equivalence in
statement 1 needs both vertices active; a vertex on no odd cycle is additive with everything.
