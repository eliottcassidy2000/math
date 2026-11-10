# The tilted tournament ensemble: a heat-flow reading of the W-polynomial, the mod-4 law of a vertex cube, and a local-limit route to the {7, 21} spectrum

**Session:** mac-mini-2026-11-10-ensemble. **Canon:** THM-5127. **Hypotheses:** HYP-9659.
**Audit:** independent, same day: PASS WITH CORRECTIONS (5 MAJOR - two missed credits, one implication, one
mis-worded hypothesis, one misattribution - and 22 minor; all applied; `audit/AUDIT.md`, MISTAKE-753).

**Owner directive (summary).** Aim to prove that 7 and 21 are the only positive odd numbers that are not the number
`H(T)` of Hamiltonian paths of a tournament; prefer reasoning and experiments that expose structure over more coverage
runs; connect to the rest of the repository and name the principles.

**Outcome.** The conjecture (Chan-Pak, arXiv:2308.10214, Conjecture 6.2) is **not proved here and remains OPEN**.
What this session adds:

1. **(PROVED, reading of old canon)** The repository's W-polynomial is the noise semigroup of the arc cube applied to
   `H`: `W_T(r) = (2r)^(n-1) (N_{1/(2r)} H)(T)`. In the coordinate `w = p(1-p)` the noised count is
   `sum_j gamma_j(T) w^j` with `gamma_0(T) = H(T)` and `sum_j gamma_j(T) 4^(-j) = n!/2^(n-1)` (Theorem A).
2. **(PROVED)** For the transitive base point this gives a one-parameter *tilted ensemble* whose mean path count is
   the gamma-polynomial of the Eulerian polynomial, `E_p[H_n] = sum_j gamma_{n,j} (p(1-p))^j`. It decreases strictly
   and continuously from `n!/2^(n-1)` to `1` as `p` runs from `1/2` to `1` (Theorem B). So every real target `>= 1` has a
   tournament order `n >= 3` and a tilt at which it is the *mean*.
3. **(PROVED)** A mod-4 law for the arcs at one vertex: for a fixed tournament `R`, the map
   `X -> (H(R + v_X) - H(R))/2 mod 2` is a quadratic form over `F_2` whose polar form is
   `beta_R(i,j) = mu_R(i,j)/2 mod 2`, where `mu_R(i,j)` is the pair term of the known cut quadratic
   (THM-4118 eq. (4)), known to be even (THM-4097 eq. (6a)) and at least `2 H(R - {i,j}) >= 2` (Session 136
   `positive_cut_geometry.md` eq. (3); THM-5122). Hence the mod-4 imbalance of `H` over that cube is exactly `0` or `2^(-rank(beta_R)/2)`, and the mod-4
   discrepancy of `H` over all tournaments of order `m+1` is at most the mean of `2^(-rank(beta_R)/2)` over order `m`
   (Theorem C). A Weyl-differencing bound covers every other frequency (Lemma C4).
4. **(FINITE-EXACT)** Censuses designed to test structure rather than coverage:
   - `beta_R` on every tournament class of orders 2-9 (198,947 classes) and every labelled tournament of orders 2-7:
     no exception to Theorem C, and the rank statistics agree with those of a uniformly random alternating matrix
     over `F_2` to within 1.2 percentage points at orders 8 and 9. A **pre-registered** out-of-sample run at order
     10 (9,733,056 classes) then agreed to within 0.08 points (P1, P2 CONFIRMED).
   - The class-count function `N_n(h)` at orders 9 and 10: in the window mean +- sd at order 10 every odd `h` has at
     least 1,504 classes, the local fluctuations are of Poisson order (variance about twice Poisson), and residue
     classes mod 3, 4, 5, 7, 8 are balanced to 1.2% or better (the only effect seen is an excess of multiples of 3).
   - A two-dimensional census of deletion pairs `(H(T), H(T + z))` over all bases `T` of order <= 10 and all `2^|T|`
     out-sets of the new vertex: for every odd `h` in `[237, 1023]` **every** odd `n` in `[3h, 5h)` occurs, and for
     every odd `h` in `[365, 1023]` every odd `n` in `[3h, 9h)` occurs.
5. **(OPEN, HYP-9659)** Three falsifiable statements, each of which is a step of a local-limit proof: (A) a tilted
   local limit law, which implies the conjecture beyond a finite range; (B) the random-alternating law for `beta_R`, in the
   quantitative form `E_R[2^(-rank(beta_R)/2)] <= C 2^(-m/2)` which gives exponential decay of the mod-4 discrepancy
   by Theorem C (4) (agreement of the rank laws in total variation alone gives only `o(1)`); (C) the window law for all `h >= 237`, which implies the
   conjecture by a one-line induction (this is the general-extension form of klein's HYP-9453 (2)).

Nothing here enlarges the range of verified values. The largest verified range is still `2^35` (THM-4885 (D); one
run, independently reproduced to `2^30`).

---

## 1. Inheritance pass and concurrency

- **Closest proved mechanism.** THM-4114 / THM-4097 (codex, 2026-08-25): a one-vertex ear response is a submodular
  directed-cut quadratic. Theorem C is a corollary of that quadratic plus Redei's parity. THM-4118 states the
  quadratic in the form (4.2) and the evenness of its coefficients; Session 136
  (`session136_full_tournament_spectrum_20261102/ear_geometry/positive_cut_geometry.md`) proves
  `mu_R(i,j) >= 2 H(R - {i,j})`. New in Theorem C: (2), (3), (4) and Lemma C4.
- **Canonical hostile example.** The one-apex family: `X ln 2` words up to `X` (Session 16; THM-4885 (A)), density
  about one half, so no single solvable family can be complete. Any argument that would also "prove" completeness of
  the one-apex spectrum is wrong.
- **Corrected near miss.** THM-4876 and Session214: the windows `[h, 2h-1]` (PROVED) and `[2h+1, 3h]` (FINITE-EXACT
  through order 8) cannot be filled from *one* realiser, and an increment `2k` from count `h` needs `h <= 4k - 1`. The two-dimensional census below is over *all* realisers and
  respects that boundary (least observed ratio 1.5058 for `h >= 237`).
- **Least-used sidecar.** THM-062 / THM-071 (opus, 2026-03): forward-edge numbers `a_k(T)`, the W-polynomial and the
  Walsh-degree decomposition of `H`. They are used here as the thermal coordinate of the problem.
- **Same-day concurrent sessions on the same directive (2026-11-10)** (read before writing; nothing below is claimed against them):
  Session246 (THM-5120..5123: lawful resets, bounded-deletion shadowing, shared-core exchange intervals), opus S90
  (THM-5124 additive rigidity on presence cubes, THM-5125, HYP-9657 mass profile of F1), kind-pasteur S21 (2-adic
  lines, in progress). S90 section 5.3 already records the Bourgain-Kontorovich analogy for the fat O/Z/J orbit and
  that positive density is unproved; opus S82 (THM-4885) records that the join families need a pointwise bound. The
  contribution here is a different analytic object: the *whole* tournament cube with a product measure, where the
  count is hyper-abundant, instead of a word semigroup.

## 2. Notation

- A tournament `T` on `n` vertices; `H(T)` its number of Hamiltonian paths. `H` is odd (Redei).
- `R + v_X`: the tournament obtained from `R` by adding a vertex `v` with out-neighbourhood `X` (and in-neighbourhood
  `V(R) - X`). This is THM-4097's ear `T + x_S`.
- For a permutation `P` of `V(T)`, `fwd_T(P)` is the number of steps `P(i) -> P(i+1)` that are arcs of `T`;
  `a_k(T) = #{P : fwd_T(P) = k}` (THM-062). So `a_{n-1}(T) = H(T)`, `sum_k a_k(T) = n!`, `a_k = a_{n-1-k}`.
- `W_T(r) = sum_k a_k(T) (r + 1/2)^k (r - 1/2)^(n-1-k)` (THM-059, THM-062).
- **Noise.** For `rho in [-1, 1]` put `p = (1 + rho)/2`, `q = (1 - rho)/2`. `(N_rho f)(T)` is the expectation of
  `f(T')`, where `T'` is `T` with each arc reversed independently with probability `q`.
- **Tilted ensemble.** Fix the order `1 < 2 < ... < n`. Under `P_p` each pair `i < j` is oriented `i -> j` with
  probability `p`, independently. `P_p` is the law of the noised transitive tournament. `P_{1/2}` is uniform.

## 3. Theorems A and B: the heat flow of H

**Theorem A (noise identity and gamma-coordinates).** For every tournament `T` on `n` vertices:

1. `(N_rho H)(T) = sum_k a_k(T) p^k q^(n-1-k)`.
2. `W_T(r) = (2r)^(n-1) (N_{1/(2r)} H)(T)` for every real `r` with `|r| >= 1/2`.
3. There are unique integers `gamma_0(T), ..., gamma_{floor((n-1)/2)}(T)` with
   `sum_k a_k(T) x^k y^(n-1-k) = sum_j gamma_j(T) (xy)^j (x + y)^(n-1-2j)`. Consequently
   `(N_rho H)(T) = sum_j gamma_j(T) w^j` with `w = pq = (1 - rho^2)/4`; `gamma_0(T) = H(T)`; and
   `sum_j gamma_j(T) 4^(-j) = n!/2^(n-1)`.

*Proof.* (1) A permutation `P` is a Hamiltonian path of `T'` iff each of its `n-1` steps is an arc of `T'`. The steps
use `n-1` distinct vertex pairs, so the events are independent. A step that is an arc of `T` survives with probability
`p`; a step whose reverse is an arc of `T` is created with probability `q`. So `P` is a path of `T'` with probability
`p^fwd q^(n-1-fwd)`; sum over `P`.
(2) Take `p = (r + 1/2)/(2r)`, `q = (r - 1/2)/(2r)`. Then `p + q = 1`, `p - q = 1/(2r)`, both lie in `[0, 1]` when
`|r| >= 1/2`, and `(2r)^(n-1) p^k q^(n-1-k) = (r + 1/2)^k (r - 1/2)^(n-1-k)`.
(3) The form `F(x, y) = sum_k a_k x^k y^(n-1-k)` is symmetric by palindromy (THM-062 (i)). The symmetric forms
`(xy)^j (x+y)^(d-2j)`, `0 <= j <= d/2`, `d = n-1`, have leading monomials `x^(d-j) y^j` with coefficient 1, so they are
a `Z`-basis of the symmetric integer forms of degree `d`; this gives existence, uniqueness and integrality. Put
`x = p, y = q` for the `w`-expansion; `(x, y) = (1, 0)` gives `gamma_0 = a_{n-1} = H`; `x = y = 1/2` gives
`sum_j gamma_j 4^(-j) = 2^(-(n-1)) sum_k a_k = n!/2^(n-1)`. QED

*Remarks.*
- In Walsh language (1)-(2) read `N_rho H = sum_j rho^(2j) H^(=2j)`, so the coefficient of `r^(n-1-2j)` in `W_T` is
  `2^(n-1-2j)` times the level-`2j` Walsh part. For every `n` this is already in canon: THM-4757 section 5
  (`sum_L (n-L)! N_L t^L = (1-t)^(n-1) F(T,(1+t)/(1-t))`, `[r^(n-1-L)] W(r) = 2^(-L)(n-L)! N_L`) with THM-4761 (1)
  (`N_L/(n)_L` is the degree-`L` Walsh chaos of `H 2^(n-1)/n!`); put `t = rho`. THM-071 (i) is the case `n = 5, 7`
  (its general case is marked CONJECTURED in that file and follows from THM-4757/4761, or from Theorem A). The
  gamma-coordinates are the triangular change of variable `rho^2 = 1 - 4w`:
  `gamma_j(T) = (-4)^j sum_{i >= j} C(i, j) H^(=2i)(T)`.
- Modulo 3 every tournament is `H` copies of the transitive point: THM-4756 (`F(T,x) = H(T) A_n(x) mod 3`) gives
  `gamma_j(T) = H(T) gamma_{n,j} (mod 3)` for every `j` (connection found by the audit, m20; checked on all labelled
  tournaments of orders 3-5).
- What is new is the reading, not the algebra: `r = 1/2` is zero temperature (the count itself, OCF lives here),
  `r -> infinity` is infinite temperature (the mean `n!/2^(n-1)`), and `W_T` interpolates along a genuine Markov
  semigroup. THM-064's imaginary evaluations are evaluations at imaginary `rho`.

**Theorem B (tilted mean).** Let `A_n(t) = sum_k A(n,k) t^k` be the Eulerian polynomial (`A(n,k)` permutations of
`[n]` with `k` descents) and `A_n(t) = sum_j gamma_{n,j} t^j (1 + t)^(n-1-2j)` its gamma-expansion. Then

    E_p[H_n] = sum_k A(n,k) p^(n-1-k) q^k = sum_j gamma_{n,j} (p(1-p))^j .

The `gamma_{n,j}` are non-negative integers with `gamma_{n,0} = 1` (Foata-Schutzenberger 1970) and, for `n >= 3`,
`gamma_{n,1} = 2^n - 2n > 0`. Therefore, for `n >= 3`, `p -> E_p[H_n]` is continuous and strictly decreasing on
`[1/2, 1]`, with `E_{1/2}[H_n] = n!/2^(n-1)` and `E_1[H_n] = 1`: every real number in `[1, n!/2^(n-1)]` is the mean
path count of exactly one tilt `p in [1/2, 1]`.

*Proof.* `P_p` is the law of `N_rho` at the transitive tournament `TT_n` with `rho = 2p - 1`. A step of `P` is an arc
of `TT_n` iff it is an ascent, so `a_k(TT_n) = A(n, n-1-k)` and Theorem A (1) gives the first equality (this is
THM-062's transitive case, i.e. THM-059's master polynomial `F_{n-1}`: `E_p[H_n] = F_{n-1}(r)/(2r)^(n-1)`). The second is
Theorem A (3) for `TT_n`, where the gamma-vector is that of `A_n`. Non-negativity is the cited theorem.
`gamma_{n,1} = A(n,1) - (n-1) = (2^n - n - 1) - (n - 1)`. A polynomial in `w` with non-negative coefficients, constant
term 1 and a positive linear coefficient is strictly increasing on `[0, 1/4]`, and `w = p(1-p)` is strictly decreasing
on `[1/2, 1]`. QED

*Checked exactly* (all labelled tournaments, `n <= 7`, `mallows.py`): the gamma-vectors `(1,2)`, `(1,8)`,
`(1,22,16)`, `(1,52,136)`, `(1,114,720,272)` reproduce `E_p[H_n]` at four rational `p` each.

*The bridging tilt (NUMERICAL, `mallows2.py`).* Let `p_n` be the tilt with `E_{p_n}[H_n] = (n-1)!/2^(n-2)`, the
uniform mean of the previous order. Then `p_n = 0.933, 0.888, 0.844, 0.804, 0.731, 0.683, 0.642` at
`n = 4, 6, 10, 16, 40, 80, 160`, and `(p_n - 1/2) sqrt(n / ln n)` is 0.7355 at `n = 4`, falls to 0.7105 at `n = 6`
and then rises slowly to 0.7954 at `n = 160`. The standard asymptotic `A_n(t)/n! ~ (1/t)((t-1)/ln t)^(n+1)` gives `E_p[H_n] ~ (n!/2^(n-1)) exp(-(4/3) n (p - 1/2)^2)` for small
`p - 1/2`, hence `p_n - 1/2 ~ sqrt(3 ln(n/2) / (4n))` (0.1433 against the computed 0.1417 at `n = 160`). This
asymptotic is not used in any proof.

*Second moments (FINITE-EXACT, observation only).* `E_p[H_n^2]` is also a polynomial in `w` (reversing the reference
order exchanges `p` and `q`). For `n = 3..7` its coefficients are non-negative integers:
`1 + 8w`; `1 + 40w + 16w^2`; `1 + 152w + 544w^2 + 56w^3 + 48w^4`; `1 + 552w + 6288w^2 + 6944w^3 + 3072w^4`;
`1 + 2040w + 55760w^2 + 189848w^3 + 199464w^4 + 72640w^5 + 4128w^6`. No interpretation is claimed. At `w = 1/4`
they give THM-589's uniform second moment.

## 4. Theorem C: the mod-4 law of a vertex cube

**Inherited cut quadratic** (THM-4097 section 2, THM-4114). For a tournament `R` on a vertex set `V` of size `m`,
let `f_R(j)` and `g_R(i)` be the numbers of Hamiltonian paths of `R` ending at `j` and starting at `i`, and `M_R[j,i]`
the number of ordered pairs `(P_1, P_2)` of disjoint non-empty paths covering `V` with `P_1` ending at `j` and `P_2`
starting at `i`. Placing `v` first, last or between two paths,

    H(R + v_X) = sum_{i in X} g_R(i) + sum_{j notin X} f_R(j) + sum_{j notin X, i in X} M_R[j,i].          (4.1)

With `a_i = H(R + v_{i}) - H(R)` and `mu_R(i,j) = M_R[i,j] + M_R[j,i]`, (4.1) is the same as

    H(R + v_X) = H(R) + sum_{i in X} a_i - sum_{{i,j} subset X} mu_R(i,j).                                 (4.2)

((4.2) is THM-4118 eq. (4) with `A_i = a_i`, `K_ij = mu_R(i,j)`; in THM-4097's notation `mu = 2w`.)
(Proof of (4.2): `X = empty` gives `H(R)`; for `X = {i}` the bracket `g(i) - f(i) + sum_{j != i} M[j,i]` is `a_i`;
in general `sum_{i in X} sum_{j notin X} M[j,i] = sum_{i in X} sum_{j != i} M[j,i] - sum_{i != j in X} M[j,i]`.)

**Theorem C.** For every tournament `R` on `m >= 2` vertices:

1. (KNOWN.) `mu_R(i,j)` is even (THM-4097 eq. (6a); THM-4118 eq. (9)) and `mu_R(i,j) >= 2 H(R - {i,j}) >= 2` for
   all `i != j` (Session 136 `positive_cut_geometry.md` eq. (3); THM-5122 "Capacity").
2. `q_R(X) = (H(R + v_X) - H(R))/2 mod 2` is a quadratic form on `F_2^V`, with polar (alternating) form
   `beta_R(i,j) = mu_R(i,j)/2 mod 2`.
3. `#{X : q_R(X) = 0} - #{X : q_R(X) = 1}` is `0` or `+-2^(m - rank(beta_R)/2)`. It is `0` exactly when `q_R` does
   not vanish identically on the radical of `beta_R`.
4. If `T` is uniform on `m+1` labelled vertices, then
   `|P(H(T) = 1 mod 4) - P(H(T) = 3 mod 4)| <= E_R[ 2^(-rank(beta_R)/2) ]`, `R` uniform on `m` labelled vertices.

*Proof.* (1) is known; the following is a second proof of `mu >= 2` that does not use the odd-cycle formula.
Parity: `H(R)`, `H(R + v_{i})`, `H(R + v_{i,j})` are odd, so `a_i` is even and by (4.2)
`mu(i,j) = a_i + a_j - (H(R + v_{i,j}) - H(R))` is even. Positivity: let `Q = q_1 ... q_r` be a Hamiltonian path of
`R - {i, j}` (`r = m - 2 >= 0`). Let `c` be the least index in `{0, ..., r}` such that `c = r` or one of `i, j` beats
`q_{c+1}`. By minimality every `q_t` with `t <= c` beats both `i` and `j`. If `c < r` let `s in {i, j}` beat `q_{c+1}`
and let `e` be the other vertex; if `c = r` take `(e, s) = (i, j)`. Then `P_1 = q_1 ... q_c e` and
`P_2 = s q_{c+1} ... q_r` are disjoint non-empty paths covering `V`, so `M[e, s] >= 1` and `mu(i,j) >= 1`. An even
integer that is at least 1 is at least 2.
(2) Divide (4.2) by 2 and reduce mod 2: `q_R(X) = sum_{i in X} (a_i/2) + sum_{{i,j} subset X} (mu(i,j)/2)`, a
polynomial of degree at most 2 in the indicator coordinates with no constant term; its polar form
`q(X + Y) - q(X) - q(Y)` has matrix `beta_R`.
(3) Classical (Dickson): a quadratic form on `F_2^m` whose polar form has rank `2r` has character sum
`sum_X (-1)^q(X)` equal to `0` or `+-2^(m-r)`; restricted to the radical the form is linear, and the sum vanishes iff
that linear form is non-zero.
(4) Write `T = R + v_X` with `v` the last vertex: `R` is uniform on `m` vertices and `X` is an independent uniform
subset. `P(H = 1) - P(H = 3) = E_R[ (-1)^((H(R)-1)/2) E_X (-1)^(q_R(X)) ]`, and by (3)
`|E_X (-1)^(q_R(X))| <= 2^(-rank(beta_R)/2)`. QED

**Lemma C4 (all frequencies).** Let `U, V'` be disjoint subsets of `V(R)` and `theta` real. With `d_j` independent,
equal to `-1, 0, 1` with probabilities `1/4, 1/2, 1/4`, and `e(x) = exp(2 pi i x)`,

    | E_X e(theta H(R + v_X)) |^2  <=  E_d  prod_{i in U} | cos( pi theta sum_{j in V'} mu_R(i,j) d_j ) | .

*Proof.* Fix the coordinates of `X` outside `U u V'`. By (4.2) the exponent is
`alpha(X_U) + gamma(X_V') - sum_{i in U, j in V'} mu(i,j) x_i x_j` in the 0/1 coordinates. Cauchy-Schwarz in `X_U`,
expand the square in two copies `X_V', X'_V'`, and take the expectation over each `x_i`, `i in U`:
`|E_{x_i} e(-theta x_i s_i)| = |cos(pi theta s_i)|` with `s_i = sum_j mu(i,j)(x_j - x'_j)`. Average over the fixed
coordinates with Jensen. QED

At `theta = 1/4` the right side is `2^(-rank of the U x V' block of beta_R)`, consistent with (3). At `theta = a/q`
with `q` an odd prime not dividing `a` it is at most `E_d cos(pi/q)^(N(d))`, `N(d)` the number of `i in U` with
`sum_j mu(i,j) d_j != 0 mod q`. So decay of the characteristic function of `H` at a rational frequency is bounded in
terms of the pair matrix `mu_R` reduced modulo the denominator (`mu_R/2` modulo 2 at `theta = 1/4`).

**What `mu_R` is.** `M_R[j,i]` counts the Hamiltonian paths of the digraph `R + (arc j -> i)` that use that arc. So
`mu_R(i,j)` is the number of Hamiltonian paths of `R` through the arc between `i` and `j`, plus the number of
Hamiltonian paths of the arc-reversed tournament `R^{ij}` through the reversed arc:
`mu_R(i,j) = U(L_i, R_j) + U(L_j, R_i)` in THM-4876 (4)'s notation (`adj(i,j) + adj'(j,i)` in THM-013's), where the
flip increment is `U(L_j, R_i) - U(L_i, R_j)`. Theorem C (1) says that for every pair at least one orientation of
the arc lies on a Hamiltonian path, and that the two counts have the same parity.

### 4.1 Exact census of `beta_R` (FINITE-EXACT)

`betarank.c` recomputes `mu_R`, checks (1), computes `rank(beta_R)` and compares (3) with the true count over all
`2^m` out-sets.

| order `m` | classes | parity / prediction failures | rank distribution (classes) | mean `abs` imbalance |
|---|---|---|---|---|
| 6 | 56 | 0 / 0 | rank 4: 30, rank 6: 26 | 0.1339 |
| 7 | 456 | 0 / 0 | rank 4: 74, rank 6: 382 | 0.0905 |
| 8 | 6,880 | 0 / 0 | rank 2: 2, rank 4: 205, rank 6: 3,795, rank 8: 2,878 | 0.0509 |
| 9 | 191,536 | 0 / 0 | rank 2: 8, rank 4: 436, rank 6: 32,535, rank 8: 158,557 | 0.0317 |
| 10 (rank only) | 9,733,056 | 0 / not run | rank 2: 7, rank 4: 668, rank 6: 210,562, rank 8: 5,444,657, rank 10: 4,077,162 | not run |

Comparison with a uniformly random alternating matrix over `F_2` (`random_alternating.py`, exact):

| `m` | rank m or m-1 (census / random) | next rank (census / random) | lower (census / random) |
|---|---|---|---|
| 8 | 0.4183 / 0.4205 | 0.5516 / 0.5585 | 0.0301 / 0.0210 |
| 9 | 0.8278 / 0.8394 | 0.1699 / 0.1593 | 0.0023 / 0.0013 |
| 10 (pre-registered) | 0.4189 / 0.4197 | 0.5594 / 0.5591 | 0.0217 / 0.0213 |

The order-10 row was predicted before it was run (`PREREGISTRATION.md`, commit 84f80f1fe2): P1 asked for the three
fractions in `[0.405, 0.430]`, `[0.550, 0.575]`, `[0.020, 0.032]`; P2 was the control "no odd `mu`, least `mu` equal
to 2". Both CONFIRMED (`out/betarank_classes_10_summary.out`; reproduced by the audit's `b_rankonly.c`). The census
columns are over isomorphism classes, the model is uniform on labelled matrices. With labelled weights `m!/|Aut|`
(audit) the fractions are `0.4205 / 0.5505 / 0.0290` at `m = 8` and `0.8286 / 0.1693 / 0.0021` at `m = 9`, and
`E_R[2^(-rank(beta_R)/2)] = 0.10239, 0.07348, 0.05077` at `m = 8, 9, 10` against `0.10135, 0.07271, 0.05072` for
random alternating matrices.

Over **all labelled** tournaments (`out/betarank_labelled_2_7.out`) the mean of the exact imbalance
`|E_X (-1)^(q_R)|` is `0.5, 0.5, 0.40625, 0.25, 0.12793, 0.08835` at `m = 2, ..., 7`. By the proof of (4) these are
rigorous bounds for the labelled mod-4 discrepancy at orders `3, ..., 8`. The true discrepancies at orders 3-7
(`out/mallows.out`, exact) are `1/2, 1/2, 13/32, 1/4, 853/8192 = 0.1041`: the bound is attained at orders 3-6 and
first has slack at order 7. For order 8 it gives `|P(H = 1 mod 4) - P(H = 3 mod 4)| <= 0.0884`.

*Exact values from the audit (class weights, two routes each).* The labelled mod-4 discrepancy is
`1361/32768 = 0.04153` at order 8, `63103/8388608 = 0.007523` at order 9 and `650471/268435456 = 0.002423` at order
10, against the bounds `0.08835, 0.05028, 0.03155` of Theorem C (4): most of the decay comes from sign cancellation
across `R`, which (4) discards. The labelled residues mod 8 at order 10 are `0.25082, 0.24901, 0.25039, 0.24978`.

### 4.2 A remark on small steps

THM-4876 (2) forbids a *vertex* from raising a count `>= 5` by exactly 2. Single *arcs* can: in the one-apex word
`1^(d-1) 0 b 1^d 0` (the word is read along the transitive chain from its source to its sink, and the apex beats
exactly the chain vertices marked 1; THM-4876 (6)) the two values of the bit `b` give `2^(2d+1) - 3*2^(d-1) - 1` and
`2^(2d+1) - 3*2^(d-1) + 1` for every `d >= 2` (PROVED: by THM-4876 (6), `H = 1 + 2N` with
`N = 2^(2d) - 3*2^(d-2) - 1` for `b = 0`, and setting `b = 1` removes `sum_{k=2..d} c(k) = 2^(d-1) - 1` and adds
`c(d+1) = 2^(d-1)`; 25 and 27, 115 and 117, 487 and 489, ...; confirmed on the actual tournaments for `d <= 8` by the
audit's `e_oneapex.c` and for `d <= 11` from the recurrence of `g` by `flip_family.py`). The reason is a balance: as a 1 the bit pairs with the final 0
(weight `2^d`), as a 0 with the `d-1` ones above it (total weight `2^d - 2`). Such balanced bits need the top half of
the word to mirror the bottom half, so heuristically they have frequency about `2^(-L/2)` among words of length `L`: unit steps
exist at every scale (compare Session246's exchange cubes, THM-5122) but do not percolate inside a solvable family.

## 5. Three regimes, and the microcanonical statement

Count structures per target.

| regime | family | structures with count `<= X` | status of coverage |
|---|---|---|---|
| critical | one-apex words | `X ln 2 + O(sqrt X)` (Session 16; THM-4885 (A)) | density about 1/2; not complete |
| supercritical | O/Z/J words (root changes) | on average at least `(3/2)^k` words per value at length `k` (S90 section 5.3) | misses exactly 29 odd values, all at most 2,647 (Session197; to `10^9` with at most seven root changes, THM-4904 / HYP-9494); a proof would be circle-method type |
| hyper-abundant | all tournaments | at least `2^C(m,2)/m!` with `m! <= X`, i.e. `exp((ln 2/2 - o(1)) (ln X / ln ln X)^2)` | every odd value to `2^35` except 7, 21 |

The last line is trivial (`H(T) <= m!` on `m` vertices) but it is the point: taken over all tournaments the problem is
not a thin-orbit problem. Every family used so far projects it to the first or second line. In the third regime the
natural statement is **microcanonical**: the energy function `H` on the arc cube takes every odd value except 7 and
21. Theorem B supplies the canonical ensemble in which a given target is the mean; what is missing is concentration
under the tilt and equivalence of ensembles at the finest scale, i.e. a local limit theorem.

**HYP-9659 (A) (tilted local limit law).** There are `n_0 >= 4` and `c > 0` such that for all `n >= n_0`, all
`p in [1/2, p_n]` and all odd `h` with `|h - E_p[H_n]| <= sd_p(H_n)`: `P_p(H_n = h) >= c / sd_p(H_n)`.

**Proposition 5.1 (PROVED implication).** HYP-9659 (A) implies that every odd `h > (n_0 - 1)!/2^(n_0 - 2)` is the
number of Hamiltonian paths of a tournament.

*Proof.* Let `mu_n = n!/2^(n-1)`; `mu_n` increases to infinity for `n >= 2`. Given odd `h > mu_{n_0 - 1}` take
`n >= n_0` with `mu_{n-1} < h <= mu_n`. By Theorem B there is `p in [1/2, p_n)` with `E_p[H_n] = h`. The hypothesis
at this `p` and this `h` gives `P_p(H_n = h) > 0`. QED

Only positivity at the mean is used; the hypothesis is stated in its natural strength because that is what a proof
by characteristic functions would give. Where the existing canon plugs in:

- small frequencies at `p = 1/2`: THM-4761 (exact second moment of every Walsh chaos,
  `Var(H)/E[H]^2 = 2/n - 14/(3 n^3) + ...`, and `sqrt(n)(H/E[H] - 1) => N(0, 2)`, CONDITIONAL on Moran's CITED central
  limit theorem for cyclic triangles) and `THM-589-metagraph-H-variance-closed-form-nud-succession`; no variance
  formula or limit law is known for `p != 1/2`;
- rational frequencies: Theorem C and Lemma C4 bound the decay through `mu_R` reduced modulo the denominator, and
  section 4.1 shows that at orders 8-10 the rank distribution of `beta_R` agrees with that of a random alternating
  matrix in three bins (HYP-9659 (B));
- the remaining range needs a decoupling of `H` over several exposed vertices; nothing here does that.

**What the data say about the hypothesis** (`analyze_hist.py`; classes, not labelled tournaments):

| order | classes | every odd `h` up to | window mean +- sd | least class count there | rms normalised second difference | spread of residue classes mod 3, 4, 5, 7, 8 |
|---|---|---|---|---|---|---|
| 9 | 191,536 | 2,881 | [757, 2066] | 108 | 1.36 | at most 4.9% (multiples of 3) |
| 10 | 9,733,056 | 14,683 | [3963, 10192] | 1,504 | 1.45 | at most 1.2% |

(1.00 would be independent Poisson counts; the measured values correspond to a variance of 1.9 and 2.1 times
Poisson.) So at order 10 the class-count function is already a smooth density with noise of Poisson order and no
visible singular series; the only residue effect seen is the excess of multiples of 3 (4.9% at order 9, 1.2% at
order 10). This is a measurement of *how* the conjecture is true in the bulk,
not of *whether*: all these values were known to occur.

Uniform labelled residues mod 8 (`out/mallows.out`), classes `1, 3, 5, 7`: order 5: `.352 .234 .352 .063`;
order 6: `.303 .176 .322 .199`; order 7: `.295 .246 .257 .202`.

## 6. The two-dimensional census

For a base `T` and a new vertex `z` with any out-set, record the pair `(h, n) = (H(T), H(T + z))`. `pairs.c` does this
for all 198,947 + 9,733,056 tournament classes of orders 2-10 and all out-sets, for `h < 1024` (`analyze_window.py`,
`out/window_le10.out`).

- Every odd `h < 1024` except 7 and 21 occurs as `H(T)` with `|T| <= 10` (510 values).
- **(W1)** For every odd `h` with `237 <= h <= 1023`, every odd `n` in `[3h, 5h)` is `H(T + z)` for some `T` with
  `H(T) = h`, `|T| <= 10`. Below 237 exactly 48 values of `h` have an incomplete window; the list is in the output
  file (largest: `h = 235`, missing `n = 839`; `h = 223`, missing `699`).
- **(W2)** For every odd `h` with `365 <= h <= 1023`, every odd `n` in `[3h, 9h)` occurs.
- For `237 <= h <= 1023` the least `n_0(h)` from which the extension set is complete up to `5h` has
  `n_0(h)/h` between 1.765 and 2.933, median 2.234. The band `[2h, 3h)` is 98.9% full and `[3h/2, 2h)` is 74% full
  (for `h >= 100` the latter is 50.5% with bases of order `<= 9` and 74.3% with order `<= 10`, audit
  `d_window_analyze.out`; so it is still filling). The least ratio is `387/257 = 1.5058`, the forced value
  `(3h+3)/(2h)` for `h = 1 mod 4`; the forced bound is `n >= (3h+1)/2` (Session214).
- Out-degree-one extensions alone fill 80% of `[3h, 5h)` and 86% of `[2h, 3h)` for `h >= 237`.

**Proposition 6.1 (PROVED implication).** Suppose that for every odd `h >= 237` that is a Hamiltonian path count and
every odd `n` with `3h <= n < 5h` there is a tournament `T` (of any order) with `H(T) = h` and a one-vertex extension
with `H(T + z) = n` (HYP-9659 (C); this is (W1) without the restriction `|T| <= 10`). Then every odd number other
than 7 and 21 is a Hamiltonian path count.

*Proof.* All odd numbers below 1024 other than 7, 21 are counts. Let `n >= 1025` be odd and let `h` be the largest
odd integer with `3h <= n`. Then `h >= 341 >= 237` and `n <= 3h + 5 < 5h`. By induction `h` is a count, and by the hypothesis
`n` is a one-vertex extension of a realiser of `h`. QED

This is klein's HYP-9453 (2) with all out-sets instead of the near-join. It is a *stronger* statement than the
conjecture and is recorded because it is two-dimensional and sharp: inside the tested range the deletion-pair set is
the full lattice cone `3h <= n < 9h`. No single base fills a window in the census: with bases
of order `<= 10` the best realiser of any `h >= 237` covers at most 44% of the odd numbers in `[3h, 5h)` (FINITE-EXACT,
`audit/d_single_realiser.out`), so the witnesses change with `n`. (PROVED is only that `(h, (3h+1)/2)` contains no
extension value, THM-4876 (2), (7) and Session214; THM-4876 (5) records FINITE-EXACT, through order 8, that no
realiser fills `[2h+1, 3h]`.) THM-5124 shows the witnesses cannot come from additive presence cubes.

## 7. Typed connections and named principles

Each connection is typed as source / target / map / preserved predicate / lost information / test.

1. **W-polynomial = heat flow (IDENTITY, Theorem A).** Source: THM-059/062, THM-4757 section 5, THM-4761 (THM-071 for `n = 5, 7`). Target: noise semigroup on the arc
   cube. Map: `r -> rho = 1/(2r)`. Preserved: every coefficient. Lost: nothing. Test: `mallows.py`.
2. **Chan-Pak = microcanonical non-emptiness (REFORMULATION).** Source: the spectrum problem. Target: local limit law
   for a polynomial of independent bits (Gilmer-Kopparty and Berkowitz for subgraph counts are the nearest proved
   templates). Preserved: the truth value beyond a finite range (Proposition 5.1). Lost: constructiveness and
   small values. Test: section 5 table; a falsifier would be a bulk value with anomalously few classes at high order.
3. **Mod-4 layer = random alternating form (MODEL, FINITE-EXACT).** Source: `beta_R`. Target: uniform alternating
   matrices over `F_2`. Preserved: rank distribution to about 1 point at orders 8-9 and 0.1 point at order 10.
   Lost: structured low-rank classes (the bin `rank <= m - 4` is over-represented by 1.43 at order 8 and 1.72 at
   order 9; rank-2 classes by 7.2, 66 and 145 at orders 8, 9, 10). Test: the pre-registered order-10 run (CONFIRMED).
4. **Zaremba / Bourgain-Kontorovich (ANALOGY, inherited from Session205, S82, S90).** It applies to the word
   semigroups (second regime), where it yields density one at best. The third regime is closer to local limit
   theorems for random-graph statistics than to thin orbits.
5. **Sun's 2-4-6-8 conjecture (CAUTION, THM-4026/4027).** Locally unobstructed and false at
   `896,315,812,331,399 = 8.96 * 10^14`. There the representation count is critical
   (`THM-4028-sun-two-four-six-eight-average-order-criticality`, cumulative count `~ V X^(25/24)`). Local solvability alone decides nothing; the regime does. For the path
   spectrum local solvability is THM-4960 (arbitrarily long runs of consecutive odd counts).
6. **Even-graph Chan-Pak (CONTRAST, THM-5040/5055).** Proved, because deficits add under joins. Here the
   join-additivity used there is not available on presence cubes (THM-5124) or on the out-set cube of one new vertex
   (`mu >= 2`); additive cubes do exist in the exchange cube of a fixed core (THM-5122). This is a typed CONTRAST,
   not an impossibility theorem.
7. **LRC good periods (SHARED PATTERN).** "Existence is a maximum or tail question, not automatically a mean
   question" (META-PATTERNS). The representation masses of THM-4885 and the class counts above are means; the
   conjecture is a statement about the minimum over targets.

**One picture.** The tournament thread began (March 2026) with three statements about one function on the arc cube:
the odd-cycle formula `H = I(Omega, 2)` (THM-002), the forward-edge deformation `a_k`, `W(r)` (THM-059, THM-062)
and the Walsh diagonalisation (THM-071). Theorem A says they are one object seen at different temperatures: `W` is
the heat flow of `H`, the Eulerian numbers are its value at the transitive point, the odd-cycle formula is its
zero-temperature fibre and `n!/2^(n-1)` its infinite-temperature fibre. What was proved later describes the
landscape of that energy: joins are independent subsystems (the count multiplies, THM-1862); an ear is a local
field with repulsive pair couplings (`+a_i`, `-mu_ij`, THM-4114; Theorem C (1) says no coupling vanishes), which is
why presence cubes are never additive (THM-5124) and why no realiser has a full window (THM-4876); the word
families are low-entropy sub-ensembles. In this language the conjecture says that the density of states has no
holes above 21, and 7 and 21 are finite-size effects below the first strong floors. The statement that is missing
in every route - join families, word semigroups, the whole cube - is the same one: pointwise control at the scale
of single integers of a sum that is only understood on average.

**Principles (recorded here, not promoted; evidence is from one problem family).**

- **Read the palindromic refinement as temperature.** Trigger: a count over orientations with a palindromic
  refinement (`a_k`, `W(r)`). Action: write it as a noise average and use `w = p(1-p)`; the zero-noise fibre is the
  exact count, the full-noise fibre is the mean. Counterindication: `|r| < 1/2` is not a probability.
- **Count the regime before choosing the proof.** Trigger: a universal coverage statement with no local obstruction.
  Action: compute structures per target for the family in hand. Bounded or logarithmic: expect a positive-density
  set of misses or a Goldbach-hard problem; polynomial: circle method, density one; super-polynomial: local limit
  law over the whole ensemble. Counterindication: a family can be fat and still miss small values (the 29 misses).
- **The arithmetic layer of a vertex cube is a quadratic form.** Trigger: a count that is a cut quadratic in the
  arcs at one vertex. Action: reduce the pair matrix modulo the frequency; rank gives the decay. Counterindication:
  the quadratic is exact only for one new vertex; two new vertices give a quartic.
- **The exceptions are finite-size effects.** 7 and 21 sit strictly between consecutive strong floors `5 < 9` and
  `15 < 25` (`THM-1370-h-spectrum-omits-7-21-all-n`) at orders where the cube has 64 and 32,768 points. No congruence sees them.

## 8. Routes tried in this session and why they stopped

- *Balanced binary doubling* (`n = 2h +- 1` from a realiser of `h`): false for single realisers at small `h`
  (FINITE-EXACT: no realiser of 13 of order `<= 10` extends to 27; audit); for all realisers it becomes the window
  law.
- *Two-dimensional induction on deletion pairs*: a pair `(h, n)` does not determine the state needed for the next
  step; closing it needs triples, then quadruples (the port state is three-dimensional and thin at computable
  heights, in agreement with THM-4875 (C)).
- *Injective sub-families for positive density*: every decoding device found (fixed low block, concatenation
  `U * 2^len + H`, slice injectivity of a two-apex readout with `gcd(x, P(x)) = 1` and `P(y) < x`) costs a square
  root. Positive density stays equivalent in difficulty to linear energy of the one-apex function (THM-4887, 4890).
- *Arithmetic progressions of counts*: in every word family the free integer `A` enters as
  `A * U - C(A) * M` with `M = H(T - u) >= 1` and `C` the pair-of-ones weight; positivity forbids cancelling `C`.
- *Unit-step flip paths (discrete intermediate value)*: balanced bits exist at every scale (section 4.2) but have
  frequency `2^(-L/2)`; no percolation in solvable families.

## 9. Reproduction

Requires `gentourng` (nauty) and a C compiler. From this directory:

    cc -O2 -o /tmp/labeled labeled.c && for n in 3 4 5 6 7; do /tmp/labeled $n > out/lab_$n.txt; done
    python3 mallows.py; python3 mallows2.py                      # Theorems A, B; tilt table; residues
    python3 verify_cut_mod4.py                                   # independent check of (4.2) and Theorem C
    cc -O2 -o /tmp/betarank betarank.c
    for m in 2 3 4 5 6 7 8 9; do gentourng -q $m | /tmp/betarank $m; done      # section 4.1 (about 10 s)
    for m in 2 3 4 5 6 7; do python3 gen_labelled.py $m | /tmp/betarank $m; done  # all labelled tournaments
    for r in 0 1 2 3 4 5 6 7; do gentourng -q 10 $r/8 | /tmp/betarank 10 fast > /tmp/b10_$r.out & done; wait   # order 10, ranks only
    python3 random_alternating.py
    cc -O2 -o /tmp/hist hist.c && gentourng -q 9 | /tmp/hist 9 > out/hist_9.txt
    for r in 0 1 2 3 4 5 6 7; do gentourng -q 10 $r/8 | /tmp/hist 10 > /tmp/h10_$r.txt & done; wait
    python3 merge_parts.py hist out/hist_10.txt /tmp/h10_*.txt   # sum the eight parts
    python3 analyze_hist.py                                       # section 5 table
    gunzip -k out/pA_le10.bin.gz out/pZ_le10.bin.gz && python3 analyze_window.py out/pA_le10.bin out/pZ_le10.bin
                                                                  # section 6 from the shipped merged bitsets, or:
    sh run_window.sh                                              # section 6: pairs.c on orders 2..10 (8 residue
                                                                  # classes at order 10), OR of the bitsets with
                                                                  # merge_parts.py, then analyze_window.py
                                                                  # (about 15 minutes on 8 cores)
    python3 flip_family.py                                        # section 4.2, d = 2..11 from the recurrence of g
    python3 ozj.py 30000                                            # non-negative O/Z/J coordinates: the known 29 misses

`ozj.py` uses the coordinates `(h, E, O - M, M)` in which the three inherited letters are non-negative matrices,
`J: (h+2E, M+O', E, h)`, `Z: (h+2E, E+O', E+O', M)`, `O: (h, E+O'+M, E+O', M)` with `O' = O - M`; it reproduces the 29
omissions of the three-letter grammar (Session197; HYP-9494) up to 30,000 (`out/ozj_30000.out`) and is included only
as a consistency check of the inherited calculus.
