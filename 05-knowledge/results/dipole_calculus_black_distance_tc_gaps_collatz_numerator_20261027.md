# Potentials and their certificates: black distance = arc connectivity, the {7, 21} presentation gaps at every order, and Collatz in the numerator coordinate

mac-mini-2026-10-27-potentials (overnight session, 2026-10-27). Theorems THM-4811, THM-4812, THM-4813; predictions
HYP-9421. **Status: pre-audit checkpoint.** Every number below is reproduced by a committed script; theorem files carry
the exact typing.

## 0. Inheritance and the concept board

The owner's dispatch merged several concurrent lanes: the regular-tournament black-core distance (`m` black moves from a
regular `(2m+1)`-tournament to a singleton plus a strong core, with a cut-score lower bound), the ROOT transport in the
coordinate `X = 3s+1`, the separation of `H(T)` from `H(T)/|Aut(T)|` in the `{7, 21}` analogy with the aliquot
untouchables `{2, 5}`, the Eisenstein reading of the directed triangle, the three-pairing repair of the dyadic ledger,
and the trace-norm obstruction for bounded module banks. Closest inherited mechanisms:

- THM-4729 / THM-4730 / black_descent_20261020c (blue and black lines; strong black descent of `c3`), and the concurrent
  regular basin theorem of the aliquot_power_124 session (`scalar_root_reparse_20261026/basins`, audited there, not yet on
  main at filing): its B1 is the cut-charge lower bound, its B3 the black route for regular sources.
- THM-4810 (score monotonicity `|L(tau_P T) - L(T)| <= 1`) and Session74 (exact Paley11 distances 6 and 15).
- THM-4790 (`tc` gaps through 11 vertices), THM-4741 (`tc = 1` iff `K1`/`C3` stack), THM-1370 (`H` gaps at every order).
- The paid-ROOT packet R1-R2 (`2^t Y = 3X + 2^a`, the price `lambda_s`), THM-4600 (anchors), THM-4633 (cycle budget).
- The concurrent dyadic-module-spectrum session (three-chart average, Eisenstein lift, trace-norm theorem T1-T6;
  internally reviewed, uncommitted at filing).

Board: **cut charge / arc connectivity / path dipole / presentation count / orbit product / numerator coordinate**.
Hostiles kept in view: the pure-blue `C3`, `R5`; the transitive-class distance of `P7` (7, above its score bound 6); the
fact that `tc` is not multiplicative under substitution.

## 1. The dipole law and the potential calculus (THM-4811 (A))

A native move along a legal word `Q = (a, ..., b)` replaces the skew matrix `K` by `-K + 2A_Q`. On score deviations
`d = K1/2` this is `d' = -d + e_a - e_b`; on cut charges `c(X) = sum_(x in X) d_x` it is
`c'(X) = -c(X) + [a in X] - [b in X]`, i.e. `f'(X) = f(V\X) + [a in X] - [b in X]` for `f(X) = e(X -> V\X)`. Every
lower bound for native-move distances in the repo is a functional of `K` that this law moves by a bounded amount:

| potential | change per move | what it bounds | attained? |
|---|---|---|---|
| arc connectivity `lambda = min_X f(X)` | at most 1 | distance to the non-strong locus | **always** (any colour); black: all but `R5` and the classified family `E` (section 2) |
| `sum d_v^2` (equivalently `c3`) | `c3' - c3 = s(a) - s(b) - 1` | triangle descent (black_descent_20261020c) | descents exist (inherited) |
| THM-4810's `L = (1/2) sum |s_(i) - i|` | at most 1 | distance to the transitive class | `P11`: 15 (Session74); `P7`: 7 > 6 |
| `(1/2) sum |d_v|` | at most 1 | distance from regular score vectors | |
| `||K||_op`, `||K||_*` | `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` | linear move bounds only (THM-4813 (C)) | far from tight |
| exceptions to a converse-closed family | at most `n - 1` | static budgets (concurrent T6) | first-order sharp for `P_q` vs `G_b` |

The cut-charge row is the new one: the "cut-score bound" of the regular basin theorem is exactly the arc-strong
connectivity, and it is attained for **every** source.

## 2. Black distance equals arc connectivity, with a classified exceptional family (THM-4811)

- **PROVED, all orders:** `|lambda(tau_Q T) - lambda(T)| <= 1`, with equality `lambda - 1` iff some minimum out-cut
  contains the start and misses the end of `Q`. Camion's Hamiltonian cycle enters every cut, so its cut at the entry
  vertex is a descent: the native-complement distance from any tournament to the non-strong locus is exactly
  `lambda(T)`; the maximum `floor((n-1)/2)` is attained exactly by the (semi-)regular tournaments (A255599 counts), and
  `lambda = delta^0` below order 10.
- **Black moves, PROVED (Theorems A, B):** a strong tournament on `n >= 4` vertices has no black descent iff it is `R5`
  or a member of the family `E = {E_t(S, x0)}`: a bottom layer `S` (a `K1`/`C3` stack ending in `C3`, with a vertex `x0`
  uniform on the blocks), `t - 1` middle vertices and a converse top layer, joined downward except for one upward jump
  chain. Members exist at every order `n >= 10` (1, 1, 3, 3, 8, 8, 20, 20 at orders 10-17). The black distance is
  `lambda` except `R5` (no black move) and `E` (black distance `2 = lambda + 1`). The proof (a dedicated proof agent, with
  an independent agent review) uses blue-cut cosets, the block normal form, a cut-chain structure theorem for
  `lambda < delta^0` and THM-2477; one bounded case (a 6-vertex block with a 3-fold rotation) is an exhaustive check.
- **How the family was found:** the order-10 census (9,355,949 strong classes) found one exception `E10`; its structure
  suggested the subfamily `E(M; Z)` (the case `t = 1` with an initial-segment pattern), the classification proof found the
  rest, and the order-11 census (884,464,590 strong classes) found exactly the one predicted member `E_2(S5)`, refuting the
  narrower pre-stated conjecture (HYP-9421 (1)).
- **Singleton targets:** when `lambda = delta^0` the (black) distance to the singleton locus is `delta^0`; every strong
  class of orders 6-9, and the 16 order-10 classes with `lambda < delta^0`, is exactly `delta^0` black moves from a
  singleton plus a strong core (HYP-9421 (6) conjectures this at every order).
- **Black-core connectivity** (HYP-9348 (b)) holds through order 10: the co-basin graph joining `lambda = 1` strong
  classes to the non-strong classes they black-descend to is connected (537,339,904 black descent lines at order 10).

## 3. The {7, 21} presentation gaps at every order (THM-4812)

The owner's observation that `H` and `tc = H/|Aut|` have different `{7, 21}` scopes is resolved: **no tournament has
`tc in {7, 21}`.** Ingredients: the order-12 extension of THM-4790 (24.9 million invariant tournaments over seven
automorphism types); the orbit product lemma `tc(T) >= prod tc(T[O])` over the `Aut(T)`-orbits, from the OCF; the
vertex-transitive lemma (only `K1, C3, R5, P7` have `tc <= 21`; complete via circulants at prime orders, the regular
census at 9, two explicit groups at 15, and a greedy triangle packing against Dixon's bound from 19 on); the counting
that leaves only `tc = 21` at orders 13-20 and 22; and an exhaustive search over the 83 (orbit configuration, index-1
or index-3 subgroup) cases, 185,642 strong tournaments, none with `tc` in `{7, 21}`.

The arithmetic side: almost all primes are touchable unconditionally (Montgomery-Vaughan plus the sieve bound for
`p = 2q + 1`); `{2, 5}` being the only prime untouchables remains conditional (THM-4805). The two gap pairs `{7, 21}` and
`{2, 5}` are still not connected by any structure-preserving map: the tournament gaps are now unconditional at every
order and in both counts, the aliquot ones are not.

## 4. Collatz in the numerator coordinate (THM-4813 (A))

`X = 3s+1` is the numerator of the next Collatz step. In it the Syracuse map is `X -> 3X/2^v + 1`, the `4x+1` lift is
multiplication by 4, and the native seed map is `X -> (3X + 2^a)/2^t`. The forcing cocycle `2 artanh(lambda) =
log(1 + 2^a/(3X))` is cohomologous to the Collatz drift cocycle `t log 2 - log 3` (transfer function `log X`). Around a
native cycle the energy identity pays exactly `(T log 2 - d log 3)/2`: the classical Steiner-Eliahou identity
(THM-4633) in forcing language. The trivial cycle pays `lambda = 1/7` per lap, `2 artanh(1/7) = log(4/3)`; the
fixed points of the `B`-free edges are THM-4600's anchors `c_a = 1/(2^a - 3)`. Native cycles and Syracuse cycles
correspond bijectively with the same `(d, T)`. Nothing here bounds the coherent ROOT readout or excludes cycles beyond
the classical results.

## 5. The Gauss-sum intertwiner and spectral budgets (THM-4813 (B), (C))

On the centred lattice with the cyclic shift acting as `zeta_q`, every circulant tournament's skew operator is
multiplication by `sum_(s in S)(zeta^s - zeta^-s)`; for Paley it is the quadratic Gauss sum with square `-q`, which is
the Gram identity `KK^T = qI - J` (klein's THM-4814 (E) gives the general doubly-regular module). At `q = 3` the concurrent session's `beta = 1 + 2 omega` appears, and the three ROOT
rows are `-R K`: the same prime ideal `(sqrt(-3))`, a different unit. The concurrent trace-norm exception bound for
bounded-arity parses is sharp to first order (Session75; `dist(P_q, G_b) <= q(q-1)/4` trivially), but spectral potentials give
only linear native-move bounds, while the score potential gives `(q^2-1)/8`.

## 6. Concurrent work credited

- aliquot_power_124 session (`scalar_root_reparse_20261026`, uncommitted at filing): the regular basin theorem B1-B4
  (exactly `m` black moves, cut-charge bound, protected-vertex route, destroyed modules), the ROOT boundary-value /
  forcing separation and the 433-versus-79 order-of-operations example.
- dyadic_module_spectrum session (uncommitted at filing): the three-chart identity `M = (E + 2S)/3`, the Eisenstein
  lattice lift, `beta = 1 + 2 omega`, the ROOT rows `(2 + omega) z`, trace-norm theorem T1-T6.
- paid_root_cycle_budget (on main): R1 corrected affine carry, R2 price bound; session73's reconciliation.
- klein-2026-10-27, THM-4814 / HYP-9422 (pushed during this session): regular sources need at least `m` native moves and
  all-black geodesic witnesses at orders 7, 9 and for `P_11` (our THM-4811 (B) is the all-source version, (C) the black
  census through order 10); nuclear-norm extremals (skew-Hadamard deletions) and the edit-count comparison with the
  score bound; the doubly-regular `Z[(1+sqrt(-q))/2]`-module (the general form of our circulant Gauss-sum statement).
- Session75 checkpoint 1 (pushed during this session, audited): the bounded-arity nuclear-norm cap and the sharp
  `(1/2 - o(1))` exception budget, which our THM-4813 (C1) only records; black-core and ROOT audit rounds continue there.
- Session74 (Paley11 distances), S63 / THM-4810 (score potential), S57 / THM-4790 (`tc` through 11).

## 7. Reproduction

- `04-computation/experiments/black_terminal_distance_20261027/`: `black_descent.cpp` (census and the Lipschitz lemma
  check), `black_descent_profile.cpp`, `score_sum_descent.cpp`, `all_black_cycle.cpp`, `dipole_law_check.py`, outputs.
- `04-computation/experiments/tc_fibres_n12_20261027/`: `tc_invariant_scan.cpp`, `cayley_tc.cpp`,
  `regular_vertex_transitive.py`, `vt15_noncayley.py`, `type_lemma_checks.py`, `window_configs.py`,
  `tc_window_scan.cpp`, outputs.
- `04-computation/experiments/collatz_numerator_carry_20261027/numerator_dictionary_check.py` (+ `.out`).
