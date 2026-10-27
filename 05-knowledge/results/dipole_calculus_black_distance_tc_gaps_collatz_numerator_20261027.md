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
| arc connectivity `lambda = min_X f(X)` | at most 1 | distance to the non-strong locus | **always** (any colour); black: all but `C3`, `R5`, family `E` (section 2) |
| `sum d_v^2` (equivalently `c3`) | `c3' - c3 = s(a) - s(b) - 1` | triangle descent (black_descent_20261020c) | descents exist (inherited) |
| THM-4810's `L = (1/2) sum |s_(i) - i|` | at most 1 | distance to the transitive class | `P11`: 15 (Session74); `P7`: 7 > 6 |
| `(1/2) sum |d_v|` | at most 1 | distance from regular score vectors | |
| `||K||_op`, `||K||_*` | `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` | linear move bounds only (THM-4813 (C)) | far from tight |
| exceptions to a converse-closed family | at most `n - 1` | static budgets (concurrent T6) | first-order sharp for `P_q` vs `G_b` |

The cut-charge row is the new one: the "cut-score bound" of the regular basin theorem is exactly the arc-strong
connectivity, and it is attained for **every** source.

## 2. Black distance equals arc connectivity (THM-4811)

See the theorem file for the complete statement; summary:

- **PROVED, all orders:** `|lambda(tau_Q T) - lambda(T)| <= 1`, with equality `lambda - 1` iff some minimum out-cut
  contains the start and misses the end of `Q`. Camion's Hamiltonian cycle enters every cut, so its cut at the entry
  vertex is a descent: the native-complement distance from any tournament to the non-strong locus is exactly
  `lambda(T)`; the maximum `floor((n-1)/2)` is attained exactly by the (semi-)regular tournaments (A255599 counts).
- **Black moves:** exhaustive through order 10 (9,355,949 strong classes at order 10): every strong class except the
  pure-blue `C3`, `R5` and **one** order-10 class has a black descent, and every descent target chain stays black, so
  `d_black = lambda` there. The exception `E10` (scores 2,2,2,2,3,6,7,7,7,7; `|Aut| = 9`; nine descent words, all blue)
  has `d_black = 2 = lambda + 1`.
- `E10` is the first member of an infinite family `E(M; Z)`: doubled copies of `D = C3[K1, M, Z]` (stacks `M`, `Z`, the
  first block of `M` a `C3`) joined by all arcs `D -> D*` except one reversed arc. **PROVED:** `lambda = 1`, the unique
  minimum cut is `D*`, every descent word is palindromic under an anti-automorphism, hence blue. Members through order 16
  have no black descent; members through order 14 have black distance exactly 2.

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
the Gram identity `KK^T = qI - J`. At `q = 3` the concurrent session's `beta = 1 + 2 omega` appears, and the three ROOT
rows are `-R K`: the same prime ideal `(sqrt(-3))`, a different unit. The concurrent trace-norm exception bound for
bounded-arity parses is sharp to first order (`dist(P_q, G_b) <= q(q-1)/4` trivially), but spectral potentials give
only linear native-move bounds, while the score potential gives `(q^2-1)/8`.

## 6. Concurrent work credited

- aliquot_power_124 session (`scalar_root_reparse_20261026`, uncommitted at filing): the regular basin theorem B1-B4
  (exactly `m` black moves, cut-charge bound, protected-vertex route, destroyed modules), the ROOT boundary-value /
  forcing separation and the 433-versus-79 order-of-operations example.
- dyadic_module_spectrum session (uncommitted at filing): the three-chart identity `M = (E + 2S)/3`, the Eisenstein
  lattice lift, `beta = 1 + 2 omega`, the ROOT rows `(2 + omega) z`, trace-norm theorem T1-T6.
- paid_root_cycle_budget (on main): R1 corrected affine carry, R2 price bound; session73's reconciliation.
- Session74 (Paley11 distances), S63 / THM-4810 (score potential), S57 / THM-4790 (`tc` through 11).

## 7. Reproduction

- `04-computation/experiments/black_terminal_distance_20261027/`: `black_descent.cpp` (census and the Lipschitz lemma
  check), `black_descent_profile.cpp`, `score_sum_descent.cpp`, `all_black_cycle.cpp`, `dipole_law_check.py`, outputs.
- `04-computation/experiments/tc_fibres_n12_20261027/`: `tc_invariant_scan.cpp`, `cayley_tc.cpp`,
  `regular_vertex_transitive.py`, `vt15_noncayley.py`, `type_lemma_checks.py`, `window_configs.py`,
  `tc_window_scan.cpp`, outputs.
- `04-computation/experiments/collatz_numerator_carry_20261027/numerator_dictionary_check.py` (+ `.out`).
