---
id: THM-4813
title: "The Collatz numerator coordinate X = 3s+1 and the native ROOT seed map; the Gauss-sum intertwiner behind the shared C3 / ROOT / Paley spectra; first-order sharpness of the spectral exception budget (credited to Session75) and the native-move Lipschitz constants of spectral potentials. (A) In X = 3x+1 the Syracuse step is X -> 3X/2^v + 1, the 4x+1 lift is X -> 4X, and the native compressed seed map of the paid-ROOT packet is X -> (3X + 2^a)/2^t (t = a + 2j); its contraction b = 3X/(2^t Y) and price lambda = 2^a/(6X + 2^a) satisfy -log b = 2 artanh(lambda) = log(1 + 2^a/(3X)), so the forcing cocycle is cohomologous to the Collatz drift cocycle t log 2 - log 3 with transfer function log X; on every native cycle the forcing cocycle sums to exactly T log 2 - d log 3, i.e. sum artanh(lambda) = (T log 2 - d log 3)/2 while the price sum lambda is strictly smaller (the Steiner-Eliahou cycle identity of THM-4633, factor by factor at the lifted Syracuse elements), native cycles correspond bijectively to Syracuse cycles with the same (d, T), the trivial cycle pays lambda = 1/7 per lap (2 artanh(1/7) = log(4/3)), and the fixed points of the B-free native edges (a = 1, 2) are THM-4600's anchors c_1 = -1 and c_2 = 1 (every anchor c_a = 1/(2^a - 3) is the fixed point of the Syracuse letter X -> 3X/2^a + 1; the B-layered native edges, j >= 1, fix non-anchors such as -1/5, -3/13). (B) On the centred lattice Lambda_q = {x in Z^q : sum x = 0} with the cyclic shift R acting as zeta_q, the skew operator of the circulant tournament C(S) is multiplication by g_S = sum_{s in S} (zeta^s - zeta^-s); for the Paley tournament (q = 3 mod 4 prime) g_S is the quadratic Gauss sum G with G^2 = -q, which gives K^2 = -q on Lambda_q and K K^T = qI - J; for q = 3, G = 1 + 2 omega, and the three ROOT rows are multiplication by 2 + omega = -omega (1 + 2 omega), i.e. ROOT = -R K on Lambda_3: same prime ideal (sqrt(-3)), different unit. (C) (CITED for the exception count: Session75's audited theorem already gives (1/2 - o(1)) of all arcs with sharp leading constant) spectral norms are Lipschitz under native moves with constants 4 cos(pi/(n+1)) (operator norm) and 2 sum_k |2 cos(k pi/(n+1))| (trace norm), which bound move distances only linearly (d(P_q, transitive) >= about q/(2 pi)) where the score potential of THM-4810 gives (q^2-1)/8."
status: >
  PROVED (elementary): (A) all statements; (B) given the classical evaluation G^2 = (-1/q) q of the quadratic Gauss sum;
  (C) the operator-norm Lipschitz constant 4 cos(pi/(n+1)) from K' = -K + 2A_Q and the path spectrum 2i cos(k pi/(n+1));
  CITED: the matrix law K(tau_P T) = -K(T) + 2L, the trace-norm tariff 2||L||_* = 2 sum_k |2 cos(k pi/(n+1))| ~ 8n/pi
  and the linear clock bound t >= (n(n-1) - C_b n^(3/2) log2 n)/(4(n-1)) = n/4 - O(sqrt(n) log n) are Session75's
  audited RESULT section 6, (6.1)-(6.3) and item 5 (C), pushed (d5f10e1cf2) before this file; the
  exception-count statement (C1) is Session75's audited theorem (session75_overnight_black_root_spectral_20261027,
  pushed before this file), recorded here only as the comparison point.
  FINITE-EXACT: numerator_dictionary_check.py (D1-D7: the edge identity for every primitive seed below 2 10^6, the Syracuse
  conjugacy for every odd x below 2 10^6, the price identities, 2000 random words for the carry formula, the trivial cycle,
  the fixed points, and the cycle/orbit forcing identity at 60 digits on 3866 orbits); dipole_law_check.py (L5, the two
  spectral Lipschitz bounds, 2940 random path complements); the ROOT = -R K identity on the basis of Lambda_3.
  CITED / CREDITED: the native seed map, its corrected edge identity 2^t Y = 3X + 2^a (the three-unit correction), the
  contraction b, the price lambda and the energy identity are the paid-ROOT packet's R1-R2
  (05-knowledge/results/paid_root_cycle_budget_20261026/README.md), reconciled by session73
  (session73_paid_root_exception_budgets_20261027/session73-root/INCOMING_2F1311_RECONCILIATION.md); the cycle product
  identity 2^K/3^k = prod(1 + 1/(3 n_i)) is classical (Steiner 1977, Eliahou 1993; THM-4633); the anchors are THM-4600;
  the Eisenstein three-chart lift, the identity beta = omega - omega^2 = 1 + 2 omega for the C3 skew operator and the
  three ROOT rows (2+omega)(a + b omega) and the trace-norm theorem T1-T6 are the concurrent dyadic-module-spectrum
  session's (worktree dyadic_module_spectrum_20261026, not on main, not audited); the integral ROOT rows themselves,
  c^r = 2 Re(omega^r A) = (2a - b, -a - b, -a + 2b), are root_eisenstein_first_arrival_20261018j E1/E2 (recovered by
  Session75's inheritance report, which states that the directed-C3 map was NOT source-recovered); the audited
  bounded-arity nuclear-norm cap C_b n log2 n is Session75's own theorem
  (05-knowledge/results/session75_overnight_black_root_spectral_20261027, session75-spectral/RESULT.md); klein's THM-4814 (E) proves the
  general doubly-regular form (theta = (I + S)/2 makes the sum-zero lattice a Z[(1+sqrt(-q))/2]-module on which S acts as
  sqrt(-q)), of which (B1) for Paley circulants is the explicit Gauss-sum case; THM-4814 (C) has the nuclear-norm
  extremals and the edit-count comparison with THM-4810's score bound; the Gauss sum and the Paley Gram identity are
  classical; the transitive skew spectrum i cot((2k-1)pi/2n)
  is THM-1875; THM-4810 supplies f(P_q) and the score potential.
  AUDITED (independent, 2026-10-27; PASS WITH CORRECTIONS, MISTAKE-681; report
  04-computation/experiments/black_terminal_distance_20261027/audit_thm4811_4813/REPORT.md): A1-A4, the cycle bijection on 61
  rational cycles, the price maxima, (B1) for Paley q <= 47 and the (C) numbers reproduced; the forcing normalisation, the
  Session75 priority and the provenance of the Eisenstein material corrected.
  NOT CLAIMED: any Collatz convergence or cycle-exclusion statement beyond the cited classical identity; any bound on the
  coherent ROOT terminal readout; the second-order term of dist(P_q, G_b).
session: mac-mini-2026-10-27-potentials
source: 05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md
depends_on:
  - 01-canon/theorems/THM-4600-anchored-debts-are-run-transparent-the-borel-torus-picture-of-pair-chains.md
  - 01-canon/theorems/THM-4633-collatz-cycles-as-rational-approximations-of-log2-3-with-the-admissible-minimum-budget.md
  - 01-canon/theorems/THM-1875-transitive-skew-char-is-the-binomial-reflection-tangent-revival.md
  - 01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md
related:
  - 05-knowledge/results/paid_root_cycle_budget_20261026/README.md
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 01-canon/theorems/THM-4634-collatz-parity-number-dictionary-irrationality-exponents-of-parity-vectors-and-jarnik-strata.md
scripts:
  - 04-computation/experiments/collatz_numerator_carry_20261027/numerator_dictionary_check.py (+ .out)
  - 04-computation/experiments/black_terminal_distance_20261027/dipole_law_check.py (+ .out)
---

# THM-4813 — the numerator coordinate, the Gauss-sum intertwiner, and sharp spectral budgets

## (A) Collatz in the numerator coordinate

For odd `x` write `X = 3x + 1` (the *numerator* of the next step) and `v = v_2(X)`.

**A1 (conjugacy).** The Syracuse map `U(x) = X/2^v` becomes `3U(x) + 1 = 3X/2^v + 1`: multiply by `3/2^v`, then add the
positive carry `1`. The lift `B(x) = 4x + 1` satisfies `3B(x) + 1 = 4(3x + 1)`, so in the numerator coordinate `B` is
multiplication by `4`, and `U(B(x)) = U(x)`. The native compressed seed map of the paid-ROOT packet (R1) takes a primitive
odd seed `s` (`s != 5 mod 8`), puts `a = 2` if `s = 1 mod 8` and `a = 1` otherwise, writes `(3s+1)/2^a = B^j(y)` with `y`
primitive, and sets `P(s) = y`, `t = a + 2j`. In the numerator coordinate this is

    X -> Y = (3X/2^a + 1)/4^j,      2^t Y = 3X + 2^a.        (A1)

The identity (A1) is the packet's corrected R1 (a constant three too large in an older display); we only rederive it:
`3B^j(y) + 1 = 4^j (3y + 1)` and `3((3s+1)/2^a) + 1 = (3X + 2^a)/2^a`.

**A2 (the forcing cocycle).** With `b = 3X/(2^t Y)` and `lambda = (1 - b)/(1 + b)` (the packet's contraction and price),
(A1) gives `1/b = 1 + 2^a/(3X)` and `lambda = 2^a/(6X + 2^a)`, hence

    -log b = 2 artanh(lambda) = log(1 + 2^a/(3X)) = (t log 2 - log 3) + (log Y - log X).   (A2)

So the forcing cocycle `2 artanh(lambda)` differs from the **Collatz drift cocycle** `t log 2 - log 3` by the coboundary of
`log X = log(3s + 1)`. Over any word of `d` edges with total clock `T` and carry
`C_d = sum_i 3^(d-i) 2^(T_(i-1) + a_i)` (packet R2), telescoping gives

    sum_i 2 artanh(lambda_i) = log(1 + C_d/(3^d X_0)) = T log 2 - d log 3 + log(X_d/X_0),

and `artanh(lambda) > lambda` yields the packet's inequality (7), strictly.

**A3 (cycles).** On a closed native word (`X_d = X_0`) the forcing cocycle sums to `sum_i 2 artanh(lambda_i) = T log 2 - d log 3 > 0`, i.e.
`sum_i artanh(lambda_i)` is exactly half the logarithmic defect of `2^T/3^d`; the price actually paid by the energy
identity, `sum_i lambda_i`, is strictly smaller (`artanh(x) > x` for `x > 0`). This is the classical product
identity `2^K/3^k = prod (1 + 1/(3 n_i))` (Steiner 1977, Eliahou 1993; THM-4633) in forcing language. Native cycles and
Syracuse cycles correspond bijectively with the same `d` and `T`: the primitive-root map `rho` (strip `B`-layers) satisfies
`rho U = P rho` on odd integers, and a native cycle `(y_0, ..., y_(d-1))` lifts to the Syracuse cycle of the elements
`B^(j_(i-1))(y_i)`, whose valuations are `a_i + 2 j_(i-1)` because `v_2(3y + 1) = a` on primitive `y`. The trivial cycle
`s = 1` (`X = 4`, `a = t = 2`) has `b = 3/4`, `lambda = 1/7` and forcing `2 artanh(1/7) = log(4/3) = 2 log 2 - log 3` per lap.
Hence the per-clock bound `lambda_s <= t_s/9139` of packet R2 on the exterior of the authenticated frontier turns any
native cycle avoiding that frontier into `0 < T log 2 - d log 3 <= (2 + 10^-6) T/9139`, the convergent window of THM-4633
(a much weaker form of its budget; no new cycle bound is claimed).

**A4 (anchors).** The edge map `X -> (3X + 2^a)/2^t` has the unique fixed point `X* = 2^a/(2^t - 3)`, i.e.
`s* = (2^a - 2^t + 3)/(3(2^t - 3))`. For `j = 0` (`t = a`) this is `s* = 1/(2^a - 3) = c_a`, the anchor of THM-4600, with
`X* = 1 + 3c_a`: `c_1 = -1` (the cycle at `-1`) and `c_2 = 1` (the trivial cycle). In the numerator coordinate every anchor `c_a` (all `a >= 1`) is the
fixed point `2^a/(2^a - 3)` of the Syracuse letter `X -> 3X/2^a + 1`; the `B`-free native edges are the letters `a = 1, 2`.

**A5 (finite checks).** D1-D7 of `numerator_dictionary_check.py` (all exact or at 60 digits) confirm A1-A4. As a side
observation, the total price `F(s) = sum lambda` of the native orbit of a primitive seed `s < 2 10^6` down to `1` never
exceeds `0.11275` (its forcing `sum 2 artanh(lambda)` is about twice that, at most `0.225534` there); the maximisers share one tail, so the supremum is a property of small seeds (not claimed in general).

## (B) The Gauss-sum intertwiner

Let `q` be an odd prime, `Lambda_q = {x in Z^q : sum x = 0}` and `(Rx)_i = x_(i+1)`. On `Lambda_q`,
`1 + R + ... + R^(q-1) = 0`, so `R` acts as a primitive `q`-th root of unity `zeta` and `Lambda_q` is a `Z[zeta]`-module.
A circulant tournament `C(S)` (`S` disjoint from `-S`, `S` union `-S` = `Z_q \ 0`, `u -> v` iff `v - u in S`) has skew matrix
`K = sum_{s in S} (R^s - R^(-s))`, so on `Lambda_q`

    K = multiplication by g_S = sum_{s in S} (zeta^s - zeta^(-s)).                    (B1)

For the Paley tournament (`q = 3 mod 4`, `S` the nonzero squares) `g_S = sum_x (x/q) zeta^x` is the quadratic Gauss sum,
`g_S^2 = (-1/q) q = -q`. Therefore `K^2 = -q` on `Lambda_q`; with `K 1 = 0` this is `K K^T = qI - J` and the flat spectrum
`+-i sqrt(q)`: the Gram identity used by the concurrent trace-norm theorem is the Gauss-sum evaluation in disguise. For
`q = 3`, `g = zeta - zeta^2 = 1 + 2 omega`, the concurrent packet's `beta` with `beta^2 = -3`. In the same packet's
embedding `iota(a + b omega) = (a, -b, b - a)` the three ROOT rows are `iota((2 + omega) z)`, and

    2 + omega = -omega (1 + 2 omega),   so   ROOT = -R K   on Lambda_3                   (B2)

(checked on a basis: `-R K iota(1) = (2, -1, -1)`). The ROOT row map and the directed-triangle skew operator generate the
same prime ideal `(sqrt(-3))` and differ by the unit `-omega` of order six. This is the precise reason for the shared
singular values (`|2 + omega| = |1 + 2 omega| = sqrt(3)`; the eigenvalues differ by the unit, `(3 +- i sqrt(3))/2` against
`+- i sqrt(3)`); the marked first ROOT row is not preserved by the unit, as the concurrent packet stresses. For general circulants
(B1) gives `||K_(C(S))||_* = sum_(j=1)^(q-1) |sum_(s in S) 2 sin(2 pi j s/q)|`.

## (C) Spectral potentials: sharp for exceptions, weak for moves

Let `G_b` be the tournaments built from singletons by substitution into quotients on at most `b` vertices. Session75's
audited theorem bounds the skew nuclear norm on `G_b` by `C_b n log2 n` and deduces that a flat target needs
`(1/2 - o(1))` of all arcs as exceptions, with sharp leading constant when transitive baselines are available; the
upper bound is the trivial `dist(P_q, G_b) <= f(P_q) <= C(q,2)/2` (a ranking or its reverse has at most half the arcs
backward):

    q(q-1)/4 >= dist(P_q, G_b) >= q(q-1)/4 - O_b(q^(3/2) log q).                      (C1)

For `b = 2` (transitive targets) the spectral form
`f(P_q) >= [q(q-1) - sqrt(q) sum_k |cot((2k-1) pi/(2q))|]/4` gives `2.4, 8.9, 35.9, 57.1, 115.7` at `q = 7, 11, 19, 23, 31`
against the exact `7, 20, 64, 92, 180` of THM-4810 and the rearrangement bound `(q^2-1)/8 = 6, 15, 45, 66, 120`; the
spectral bound overtakes the rearrangement bound from `q = 43`.

By THM-4811 (A), a native move replaces `K` by `-K + 2A_Q`, and the directed path has skew spectrum `2i cos(k pi/(n+1))`.
Therefore `| ||K'||_op - ||K||_op | <= 4 cos(pi/(n+1))` and `| ||K'||_* - ||K||_* | <= 2 sum_(k=1)^n |2 cos(k pi/(n+1))|`
(about `8n/pi`). With `||K_(T_q)||_op = cot(pi/(2q))` and `||K_(P_q)||_op = sqrt(q)`,

    d(P_q, transitive class) >= (cot(pi/(2q)) - sqrt(q))/(4 cos(pi/(q+1))),  about q/(2 pi),

which is linear in `q` (precisely `q/(2 pi) - sqrt(q)/4 + O(1/q)`; `0.94` at `q = 11`) and weaker than the linear clock
bound `q/4 - O(q^(1/2) log q)` already proved by Session75 (RESULT section 6 and item 5 (C): exception tariff `n - 1` per
move; its trace-norm tariff (6.2)-(6.3) is the second Lipschitz constant above), while THM-4810's score potential gives `(q^2-1)/8` (attained at `q = 11` by Session74's fifteen
moves). The spectral potentials measure static representation complexity (exceptions, rank, trace mass); dynamic move
counts are governed by the cut and score potentials of THM-4811 (A). This is recorded as the failure anatomy of a tempting
transfer, not as a new bound.
