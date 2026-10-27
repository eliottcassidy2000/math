# Potentials and their certificates: black distance = arc connectivity, the {7, 21} presentation gaps at every order, and Collatz in the numerator coordinate

mac-mini-2026-10-27-potentials (overnight session, 2026-10-27). Theorems THM-4811, THM-4812, THM-4813; predictions
HYP-9421. **Status: independently audited** (three audits, all PASS WITH CORRECTIONS, corrections applied; MISTAKE-681
records them). Every number below is reproduced by a committed script; the theorem files carry the exact typing.

## 0. Inheritance and the concept board

The owner's dispatch merged several concurrent lanes: the regular-tournament black-core distance (`m` black moves from a
regular `(2m+1)`-tournament to a singleton plus a strong core, with a cut-score lower bound), the ROOT transport in the
coordinate `X = 3s+1`, the separation of `H(T)` from `H(T)/|Aut(T)|` in the `{7, 21}` analogy with the aliquot
untouchables `{2, 5}`, the Eisenstein reading of the directed triangle, the three-pairing repair of the dyadic ledger,
and the trace-norm obstruction for bounded module banks. Closest inherited mechanisms:

- THM-4729 / THM-4730 / black_descent_20261020c (blue and black lines; strong black descent of `c3`), and the concurrent
  regular basin theorem of the aliquot_power_124 session (`scalar_root_reparse_20261026/basins`, audited there,
  uncommitted at filing): its B1 is the cut-charge lower bound, its B3 the protected-vertex black route.
- THM-4810 (score monotonicity `|L(tau_P T) - L(T)| <= 1`) and Session74 (exact Paley11 distances 6 and 15).
- THM-4790 (`tc` gaps through 11 vertices), THM-4741 (`tc = 1` iff `K1`/`C3` stack), THM-1370 (`H` gaps at every order).
- The paid-ROOT packet R1-R2 (`2^t Y = 3X + 2^a`, the price `lambda_s`), THM-4600 (anchors), THM-4633 (cycle budget).
- The concurrent dyadic-module-spectrum session (three-chart average, Eisenstein lift, trace-norm theorem T1-T6;
  internally reviewed, uncommitted at filing) and Session75 (its audited spectral exception budgets and move tariffs).

Board: **cut charge / arc connectivity / path dipole / presentation count / orbit product / numerator coordinate**.
Hostiles kept in view: the pure-blue `C3`, `R5`; the transitive-class distance of `P7` (7, above its score bound 6); the
fact that `tc` is not multiplicative under substitution.

## 1. The dipole law and the potential calculus (THM-4811 (A))

A native move along a legal word `Q = (a, ..., b)` replaces the skew matrix `K` by `-K + 2A_Q` (also Session75 RESULT
(6.1)). On score deviations `d = K1/2` this is `d' = -d + e_a - e_b`; on the cut function `f(X) = e(X -> V\X)` it is
`f'(X) = f(V\X) + [a in X] - [b in X]`. Every lower bound for native-move distances in the repository is a functional of
`K` that this law moves by a bounded amount:

| potential | change per move | what it bounds | attained? |
|---|---|---|---|
| arc connectivity `lambda = min_X f(X)` | at most 1 | distance to the non-strong locus | **always** (any colour); black: all but `R5` and the classified family `E` (section 2) |
| `delta^0 = min_v min(s_v, n-1-s_v)` | at most 1 | distance to a singleton / singleton plus strong core | when `lambda = delta^0`; exhaustive through order 10 |
| `sum d_v^2` (equivalently `c3`) | `c3' - c3 = s(a) - s(b) - 1` | triangle descent (black_descent_20261020c) | descents exist (inherited) |
| THM-4810's `L = (1/2) sum |s_(i) - i|` | at most 1 | distance to the transitive class | `P11`: 15 (Session74); `P7`: 6 (corrected: THM-4810's 7 is the labelled distance to one fixed transitive tournament, forced by parity; THM-4815 (C), MISTAKE-682) |
| `||K||_op`, `||K||_*` | `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` (trace tariff: Session75 (6.2)-(6.3)) | linear move bounds only (Session75's exception tariff gives `n/4 - O(sqrt(n) log n)`) | far from tight |
| exceptions to a converse-closed family | at most `n - 1` | static budgets (Session75, dyadic T6) | first-order sharp for `P_q` vs `G_b` |

The cut-charge row is the new one: the "cut-score bound" of the regular basin theorem is exactly the arc-strong
connectivity, and it is attained for **every** source.

## 2. Black distance equals arc connectivity, with a classified exceptional family (THM-4811)

- **PROVED, all orders:** `|lambda(tau_Q T) - lambda(T)| <= 1`, with equality `lambda - 1` iff some minimum out-cut
  contains the start and misses the end of `Q`. Camion's Hamiltonian cycle enters every cut, so its cut at the entry
  vertex is a descent: the native-complement distance from any tournament to the non-strong locus is exactly
  `lambda(T)`. The maximum `floor((n-1)/2)` is attained exactly by the (semi-)regular tournaments (A255599 counts), and
  every strong tournament of order at most 9 has `lambda = delta^0` (a smaller cut needs `n >= 4 delta^0 + 2`).
- **PROVED, all orders (Theorems A, B, C):** a strong tournament on `n >= 4` vertices has no black descent iff it is `R5`
  or a member of the family `E = {E_t(S, x0)}`: a bottom layer `S` (a `K1`/`C3` stack ending in `C3`, with a vertex `x0`
  uniform on the blocks, beating the first block and beaten by the last), `t - 1` middle vertices and a converse top
  layer, joined downward except for one upward jump chain. Members exist at every order `n >= 10` (1, 1, 3, 3, 8, 8, 20,
  20, 48 at orders 10-18). The black distance is `lambda` except `R5` (no black move) and `E` (black distance
  `2 = lambda + 1`). The proof (classification/PROOF.md) uses blue-cut cosets, the prime-orbit block normal form, a
  cut-chain structure theorem for `lambda < delta^0`, and THM-2477 (via the cited Benhocine-Wojda theorem); one bounded
  case (a 6-vertex block with a 3-fold rotation) is an exhaustive check. Independently audited (own code through order
  10, all 112 members to order 18, exhaustive structured searches to order 17, 320,000 perturbed and random
  tournaments; no exception outside `E`).
- **How the family was found:** the order-10 census (9,355,949 strong classes) found one exception `E10`; its structure
  suggested the subfamily `E(M; Z)` (the case `t = 1` with an initial-segment pattern), the classification proof found
  the rest, and the order-11 census (884,464,590 strong classes, the 276,971 nonregular self-converse ones tested
  directly) found exactly the predicted member `E_2(S5)`, refuting the narrower pre-stated conjecture (HYP-9421 (1)).
- **Singleton targets:** when `lambda = delta^0` the black distance to the singleton locus is `delta^0` (Theorem C keeps
  `lambda = delta^0` along singleton-type descents); every strong class of orders 6-10 is exactly `delta^0` black moves
  from a singleton plus a strong core (all 9,355,949 at order 10, including the 16 with `lambda < delta^0`); HYP-9421 (6)
  conjectures this at every order.
- **Black-core connectivity** (HYP-9348 (b)) holds through order 10: the co-basin graph joining `lambda = 1` strong
  classes to the non-strong classes they black-descend to is connected (537,339,904 black descent lines at order 10),
  and every non-pure-blue strong class reaches it by black moves (at every order by the inherited `c3`-descent).

## 3. The {7, 21} presentation gaps at every order (THM-4812)

The owner's observation that `H` and `tc = H/|Aut|` have different `{7, 21}` scopes is resolved: **no tournament has
`tc in {7, 21}`.** Ingredients: the order-12 extension of THM-4790 (21.5 million invariant tournaments, 18.1 million
strong, over seven automorphism types); the orbit product lemma `tc(T) >= prod tc(T[O])` over the `Aut(T)`-orbits, from
the OCF; the vertex-transitive lemma (only `K1, C3, R5, P7` have `tc <= 21`; complete via circulants at prime orders, the
regular census at 9, two explicit groups at 15, and a greedy triangle packing against Dixon's bound from 19 on); the
counting that leaves only `tc = 21` at orders 13-20 and 22; and an exhaustive search over the 83 (orbit configuration,
index-1 or index-3 subgroup) cases, 185,642 strong tournaments, none with `tc` in `{7, 21}`. The audit added a
strengthening: applying the product lemma to `Aut(S)` itself forces index 1, so the window tournaments are substitutions.

Two further consequences: every **strong** presentation fibre is finite (order `O(log^2 v)`), and the fibres below 24 are
explicit: sizes 2, 1, 2, 0, 5, 4, 1, 6, 6, 4, 0, 14 for `v = 1, 3, ..., 23`, all on orders at most 8. With the
pre-registered census of section 6, the presentation spectrum below 6002 is exactly `{odd} \ {7, 21}`.

The arithmetic side: almost all primes are touchable unconditionally (a routine corollary of Montgomery-Vaughan plus
the sieve bound for `p = 2q + 1`); `{2, 5}` being the only prime untouchables remains conditional (THM-4805). The two gap
pairs `{7, 21}` and `{2, 5}` are still not connected by any structure-preserving map: the tournament gaps are now
unconditional at every order and in both counts, the aliquot ones are not.

## 4. Collatz in the numerator coordinate (THM-4813 (A))

`X = 3s+1` is the numerator of the next Collatz step. In it the Syracuse map is `X -> 3X/2^v + 1`, the `4x+1` lift is
multiplication by 4, and the native seed map is `X -> (3X + 2^a)/2^t`. The forcing cocycle
`2 artanh(lambda) = log(1 + 2^a/(3X))` is cohomologous to the Collatz drift cocycle `t log 2 - log 3` (transfer function
`log X`). Around a native cycle the forcing cocycle sums to exactly `T log 2 - d log 3` (so
`sum artanh(lambda) = (T log 2 - d log 3)/2`, and the price `sum lambda` actually paid is strictly smaller): the
classical Steiner-Eliahou identity (THM-4633) in forcing language. The trivial cycle pays `lambda = 1/7` per lap,
`2 artanh(1/7) = log(4/3)`. The fixed points of the `B`-free native edges are THM-4600's anchors `c_1 = -1` and
`c_2 = 1`; every anchor `c_a = 1/(2^a - 3)` is the fixed point of the Syracuse letter `X -> 3X/2^a + 1`. Native cycles and
Syracuse cycles correspond bijectively with the same `(d, T)`. Nothing here bounds the coherent ROOT readout or
excludes cycles beyond the classical results.

## 5. The Gauss-sum intertwiner and spectral budgets (THM-4813 (B), (C))

On the centred lattice with the cyclic shift acting as `zeta_q`, every circulant tournament's skew operator is
multiplication by `sum_(s in S)(zeta^s - zeta^-s)`; for Paley it is the quadratic Gauss sum with square `-q`, which is
the Gram identity `KK^T = qI - J` (klein's THM-4814 (E) gives the general doubly-regular module). At `q = 3` the
concurrent session's `beta = 1 + 2 omega` appears, and the three ROOT rows (root_eisenstein_first_arrival_20261018j
E1/E2) are `-R K`: the same prime ideal `(sqrt(-3))`, a different unit (equal moduli, different eigenvalues). The
bounded-arity exception budget is sharp to first order (Session75; `dist(P_q, G_b) <= q(q-1)/4` trivially), but
spectral potentials give only linear native-move bounds (Session75's `n/4 - O(sqrt(n) log n)` is the strongest), while
the score potential gives `(q^2-1)/8`.

## 6. Pre-registered tests and a dyadic remark

- **HYP-9421 (3), CONFIRMED:** every odd `v <= 2001` other than 7, 21 is `tc(T)` for some tournament on at most 10
  vertices; in fact every odd `v <= 6001` other than 7, 21 is (all 9.9 million classes of orders 1-10; no even value).
- **HYP-9421 (4), CONFIRMED:** the native price `F(s) = sum lambda` along the compressed orbit of a primitive seed has
  maximum 0.112749570730 below `10^8` (at `s = 78027147`), `2.8 10^-7` above the small-seed maximum, under the
  pre-registered 0.1127503.
- **Grassmann chart average** (generalising the dyadic session's three-chart identity `M = (E + 2S)/3`): for `2^r`
  children indexed by `F_2^r`, averaging the intermediate energy over all `k`-dimensional subspaces gives
  `S + theta_k (E - S)` with `theta_k = (2^k - 1)/(2^r - 1)` (exact; `grassmann_chart_average.py`). The averaged chain
  realises the skip-level debt `(E - S)_+` of a block as a chain of genuine averaged partitions.

## 7. Audits

Three independent audits with their own code: THM-4812 (every count of (a), (c), (d) and the (g) window reproduced; two
text majors: a non-monotone margin, an index case), THM-4811 (A, B, D, E) with THM-4813 (a factor two in the forcing
normalisation, Session75's priority for the move law and spectral tariffs, a false provenance claim, a stale exception
list), and the classification (Theorems A, B, C confirmed; a side claim about the extremal-predecessor statement
withdrawn). All corrections applied; MISTAKE-681.

## 8. Concurrent work credited

- aliquot_power_124 session (`scalar_root_reparse_20261026`, uncommitted at filing): the regular basin theorem B1-B4
  (exactly `m` black moves, cut-charge bound, protected-vertex route, destroyed modules), the ROOT boundary-value /
  forcing separation and the 433-versus-79 order-of-operations example.
- dyadic_module_spectrum session (uncommitted at filing): the three-chart identity, the Eisenstein lattice lift,
  `beta = 1 + 2 omega`, the ROOT rows `(2 + omega) z`, trace-norm theorem T1-T6.
- paid_root_cycle_budget (on main): R1 corrected affine carry, R2 price bound; session73's reconciliation.
- klein-2026-10-27, THM-4814 / HYP-9422 (pushed during this session): regular sources need at least `m` native moves to a
  singleton and all-black geodesic witnesses at orders 7, 9 and for `P_11` (our THM-4811 (B) concerns the non-strong
  locus, (D) the singleton targets); nuclear-norm extremals and the edit-count comparison with the score bound; the
  doubly-regular `Z[(1+sqrt(-q))/2]`-module.
- Session75 checkpoint 1 (pushed during this session, audited): the bounded-arity nuclear-norm cap, the sharp
  `(1/2 - o(1))` exception budget, the matrix law `-K + 2L`, the trace-norm tariff and the linear clock bound.
- Session74 (Paley11 distances), S63 / THM-4810 (score potential), S57 / THM-4790 (`tc` through 11).

## 9. Reproduction

- `04-computation/experiments/black_terminal_distance_20261027/`: `black_descent.cpp` (census and the Lipschitz lemma
  check), `black_descent_profile.cpp`, `score_sum_descent.cpp`, `all_black_cycle.cpp`, `dipole_law_check.py`,
  `exceptional_family.py`, `e10_analysis.py`, `singleton_distance.cpp`, `cobasin_connectivity.cpp`, `locus_distances.cpp`,
  `lambda_vs_delta0.cpp`, `lambda_lt_delta0_codes.cpp`, `classification/` (PROOF.md, checkers, results, the proof agent's
  internal review, and `audit_independent/`), `audit_thm4811_4813/`, outputs.
- `04-computation/experiments/tc_fibres_n12_20261027/`: `tc_invariant_scan.cpp`, `cayley_tc.cpp`,
  `regular_vertex_transitive.py`, `vt15_noncayley.py`, `type_lemma_checks.py`, `greedy_vt_bound.py`, `window_configs.py`,
  `tc_window_scan.cpp`, `fibre_configs.py`, `fibre_scan.cpp`, `tc_strong_fibres.cpp`, `strong_fibre_table.py`,
  `tc_spectrum_sweep.cpp`, `audit_thm4812/`, outputs.
- `04-computation/experiments/collatz_numerator_carry_20261027/`: `numerator_dictionary_check.py`, `forcing_tail.c`,
  `grassmann_chart_average.py`, outputs.
