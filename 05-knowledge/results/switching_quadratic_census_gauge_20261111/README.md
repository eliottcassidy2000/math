# A switching class carries a quadratic form; the census at giant orders; a reading of the Yang-Mills construction

**Session:** mac-mini-2026-11-11-gauge. **Canon:** THM-5139, THM-5140. **Hypotheses:** HYP-9669, HYP-9670.
**Audit:** AUDIT-PENDING (pre-audit checkpoint; nothing here is canonised until the independent audit is applied).

Owner directive (2026-11-11): predict A000568 at very large orders "a myriad of ways", classify the emergent
recursive tournament structures and their interactions, do the same for Chan-Pak and related sequences, metrics
and symmetry structures, and read the OpenAI Yang-Mills construction with the repository's advances in mind.

**The Chan-Pak conjecture is not proved here and remains OPEN. None of the five Yang-Mills items listed by the
owner is settled here; all five remain OPEN.**

The same directive was worked today by Codex: Session271 (`../session271_predictive_tournaments_20261111/`,
forecasts of A000568 to `n = 10^1000`, a structure atlas, fixed-core Chan-Pak growth) and Session272
(`../session272_yang_mills_bridges_20261111/`, sixteen bridges and five extension programs for the Yang-Mills
construction). Both were read before this note was written; this note does not repeat them and says where it
builds on them.

## 0. What is in this note

| | statement | status |
|---|---|---|
| A | **Endpoint graph (Theorem G).** For a tournament `T` let `e_T(v)` be the parity of the number of Hamiltonian paths that start or end at `v`. On a switching class `e` is affine: `e_{T^X} = e_T + Lambda 1_X`, where `Lambda` is the adjacency matrix of an even graph that depends only on the class. | PROVED (section 2), exhaustive control on all labelled tournaments of orders 3-7 |
| A | **Quadratic law.** The parity of the number of Hamiltonian cycles, the parity of the number of cycles of each length, and `(H - 1)/2 mod 2` are quadratic functions of the switching set; equivalently the eight switchings at any three vertices have `sum H = 0 (mod 4)`. This proves the observation (QF) of the previous session (HYP-9667, THM-5136). | PROVED |
| A | Corollaries: every switching class of order 6 has exactly 20 members with `H = 1` and 12 with `H = 3 (mod 4)`; the classes of `P_19`, `D_19`, `P_23`, `D(P_11)` with a source are constant in every cycle parity (all `2^19`, `2^23` members). | PROVED from finite data |
| A | **2-adic derivative law (HYP-9669).** `v_2(c_U(H)) >= ceil(|U|/2)` for the switching differences `c_U`; then the binary digit `k` of `H` has `F_2`-degree at most `2^k` on a class. | CONJECTURE; `|U| <= 4` proved; pre-registered run at orders 11-14 CONFIRMED |
| B | **Explicit uniform remainder** for the level expansion of A000568 (all `n`, all truncation orders), **certified enclosures** of `a(n)` up to `n = 10^1000`. | PROVED (elementary), ball arithmetic |
| B | A000568 as an expectation over **integer heights**; the moments do not determine the measure. | PROVED (elementary) |
| B | **Switching classes (A049313):** the odd branch has the tournament polynomials with the clock `n - 1`; first correction `(64/3)(n)_3 4^(-n)`. | PROVED; exact to `n = 19` |
| C | Forecasts for order 12 (least unattained value, maximum, number of values). | NUMERICAL, registered, not scorable today (HYP-9670) |
| D | Reading of the Yang-Mills construction: fourteen typed analogies, and for each of the five open items what the source contains, what blocks, and an exactly solved instance in the repository's own models; two conditional suggestions (blocked Wilson loops, a multiscale measure). | ANALOGY / SUGGESTION, typed in place |

## 1. Inheritance

Read before work started: an inheritance survey by a sub-agent (canon, hypotheses including the historical
indexes, results, backlog, reflections), the two Codex packets of today, and the files named below.

- **Switching as a gauge symmetry is old in this repository** and is not a contribution of this note: THM-474
  ("the gauge theorem": every class has one member containing a base path; Babai-Cameron 2000), THM-479
  (level permutations; Burnside for A049313), THM-1405 (star flips as gauge transformations, an `F_2` Wilson
  loop), THM-1460 / THM-1465 / THM-1470 (score parities inside a class), THM-1467 (class sum of spanning-tree
  counts, so the class mean of `H` is `n!/2^(n-1)`), THM-4524 (Walsh expansion of `H` on a class, degree at most
  `2 floor(N/3)`, Ising form for `N <= 5`), THM-4531 / THM-4721 / THM-4722 (twisted Mallows-Sloane; no
  equivariant bijection between classes and even graphs), Session34
  (`../session34_gauge_representations_20261023/gauge/SIMPLICIAL_WILSON_FANO.md`: triangle fluxes, the triangular
  Wilson action on `K_n`), Session272 (exact `Z_2` Wilson laws on a disk and a torus; "star flips become `Z_2`
  vertex gauge changes").
- **(QF)** was an observation of this machine's previous session (HYP-9667 and section 5.3 of
  `../character_of_two_records_20261111/README.md`): exhaustive for `N <= 6`, random tests to `N = 20`, "not
  proved". The vertex-cube analogue is THM-5127 (C): `X -> (H(R + v_X) - H(R))/2 mod 2` is a quadratic form with
  polar form `beta_R(i,j)` = the number of paths of `R` between `i` and `j` on an even number of vertices.
- **The level expansion of A000568** (`n! a(n)/2^C(n,2) = sum_k P_k(n) 4^(-kn)`, degree `3k`, leading term from
  `3^k`), the bit-loss lemma (a permutation moving `m` points loses at least `m(2n-m-1)/3` bits) and the
  statistical-mechanics reading are inherited: Moon 1968 and Wright 1970 (correction terms), the April
  reflections (`07-reflections/burnside-perturbation-theory.md`), THM-514, THM-516, Session98
  (`../session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md`, priority for the pole theorem,
  and `LIMITS_AND_NORMALIZATION.md` section 2 for `0 <= b_n - 1 <= q_n^3/(1-q_n)`), THM-4845 (Gaussian average
  with covariance `ln 2 gcd(p,q)`), Session249 (Mayer expansion of the cycle gas), Session271 (support envelope
  `U_M(n)`, forecasts, the factor 64 between two 3-cycles, the mixture `U = L/(1+delta) + delta D/(1+delta)`).
  Session271 lists as open a growing-order truncation theorem and interval-certified digits; sections 3.2-3.3
  answer those two points.
- **Chan-Pak data used:** THM-4889 / HYP-9481 (spectrum by order to order 11), THM-329 (lower bounds for the
  maximum at orders 12-14), MISTAKE-743 (the value "545,835" for order 12 has no source).

## 2. Part A: the endpoint graph and the quadratic law

### 2.1 Notation

`T` is a tournament on a vertex set `V`, `n = |V|`. `N^+(v)` is the out-neighbourhood of `v`, `iota_v` its indicator
vector. `s_T(v)` (`t_T(v)`) is the number of Hamiltonian paths starting (ending) at `v`; `H(T) = sum_v s_T(v)`.
`e_T(v) = s_T(v) + t_T(v) mod 2` (so `e = 0` when `n = 1`). `hc(T)` is the number of directed Hamiltonian cycles
(`n >= 3`). For `X` a set of vertices, `T^X` is `T` with the arcs between `X` and its complement reversed;
`T^X = T^(V - X)`, `(T^X)^Y = T^(X + Y)` (symmetric difference), `(T^X)[W] = T[W]^(X cap W)`, and
`(T^X) - u = (T - u)^X` when `u` is not in `X`. All congruences in this section are modulo 2 unless said otherwise.

### 2.2 Theorem G

**Lemma 1 (deletion).** For `n >= 2`: `e_T(u) = 1 + sum_{v in N^+(u)} e_{T-u}(v)`.

*Proof.* `s_T(u) = sum_{v in N^+(u)} s_{T-u}(v)` and `t_T(u) = sum_{v in N^-(u)} t_{T-u}(v) = H(T-u) - sum_{v in N^+(u)} t_{T-u}(v)`;
`H(T-u)` is odd (Redei).

**Lemma 2 (diagonal).** `e_{T^{u}}(u) = e_T(u)`.

*Proof.* `s_T(u) + s_{T^u}(u) = sum_{v != u} s_{T-u}(v) = H(T-u)`, and the same for `t`; the sum changes by `2 H(T-u)`.

**Theorem G.** For every tournament `T` there is a symmetric matrix `Lambda_T` over `F_2`, with zero diagonal and
zero row sums (the adjacency matrix of a graph in which every degree is even), such that for every `X`

    e_{T^X} = e_T + Lambda_T 1_X ,      Lambda_{T^X} = Lambda_T .

For `u != w` and `n >= 4`, with `Q = T - u - w`, `a = N^+(u) cap Q`, `b = N^+(w) cap Q`:

    (*)   Lambda_T(u,w) = e_{T-u}(w) + sum_{v in N^+(u)} Lambda_{T-u}(v,w)        (n >= 3)
    (**)  Lambda_T(u,w) = 1 + <a + b, e_Q> + a^T Lambda_Q b                      (n >= 4)

and `Lambda_T` is the complete graph for `n = 3`, zero for `n <= 2`.

*Proof.* Induction on `n`. For `n <= 2`, `e` is `(0)` or `(1,1)` for every tournament and `Lambda = 0`.
Let `n >= 3`, fix `u`, write `R = T - u`. For `X` not containing `u`, `N^+_{T^X}(u) = N^+_T(u) + X`, so by Lemma 1
and the induction hypothesis for `R`

    e_{T^X}(u) = 1 + <iota_u + 1_X , e_R + Lambda_R 1_X> = e_T(u) + <1_X , e_R + Lambda_R iota_u> + <1_X , Lambda_R 1_X> ,

and the last term is 0 because `Lambda_R` is symmetric with zero diagonal. So on the sets `X` avoiding `u`,
`e_{T^X}(u) = e_T(u) + sum_{w in X} Lambda(u,w)` with `Lambda(u,w)` given by (*). Put `Lambda(u,u) = 0`. Taking
`X = V - u`, Lemma 2 gives `sum_{w != u} Lambda(u,w) = 0`; with `T^X = T^(V-X)` the formula
`e_{T^X}(u) = e_T(u) + sum_{w in X} Lambda(u,w)` then holds for every `X`. Hence `e_{T^X} = e_T + Lambda_T 1_X`
with zero diagonal and zero row sums. `Lambda_{T^Y} = Lambda_T` because
`e_{T^(X+Y)} = e_T + Lambda_T(1_X + 1_Y) = e_{T^Y} + Lambda_T 1_X` for all `X`.

Symmetry. For `n = 3`, (*) gives `Lambda(u,w) = e_{T-u}(w) = 1`. For `n >= 4` and `u != w`: in (*),
`e_{T-u}(w) = 1 + <b, e_Q>` (Lemma 1 in `T - u`); the term `v = w` is `Lambda_{T-u}(w,w) = 0`; for `v` in `Q`,
`Lambda_{T-u}(v,w) = Lambda_{T-u}(w,v)` (symmetry at order `n - 1`) `= e_Q(v) + sum_{v' in b} Lambda_Q(v',v)` ((*) for
`T - u` at `w`). Summing over `v` in `a` gives (**), which is symmetric in `(u,a) <-> (w,b)` because `Lambda_Q`
is (symmetry at order `n - 2`). QED

**Corollary 1 (Hamiltonian cycles).** For `n >= 3` and every `X`:

    hc(T^X) = hc(T) + sum_{a in X} (1 + e_T(a)) + sum_{{a,b} in X} Lambda_T(a,b) .

In words: the first switching derivative of the parity of the number of Hamiltonian cycles is the endpoint
vector (plus one), the second is the endpoint graph, the third vanishes.

*Proof.* A Hamiltonian cycle through `a` is `a -> x ... y -> a` with a Hamiltonian path of `T - a` from `x` to `y`;
it needs `iota_a(x) = 1`, `iota_a(y) = 0` in `T` and the opposite in `T^a`. Hence
`hc(T) + hc(T^a) = sum_{x in N^+(a)} (s_{T-a}(x) + t_{T-a}(x)) = 1 + e_T(a)` by Lemma 1. Apply this to `T^X` and use
Theorem G: `hc(T^(X+a)) = hc(T^X) + 1 + e_T(a) + sum_{w in X} Lambda(a,w)`; induct on `|X|`.

**Corollary 2 (cycles of each length; `H mod 4`).** Let `c_l(T)` be the number of directed cycles of length `l`.
Since `c_l(T) = sum_{|W| = l} hc(T[W])` and `(T^X)[W] = T[W]^(X cap W)`, the parity of `c_l(T^X)` is a quadratic
function of `X`, for every `l >= 3` (odd or even), with polar form `sum_{|W| = l, W contains a, b} Lambda_{T[W]}(a,b)`.
Since `H = 1 + 2 (number of odd cycles) (mod 4)` (the odd-cycle formula, THM-002),

    (H(T^X) - H(T))/2 = sum_{a in X} lambda_a + sum_{{a,b} in X} M(a,b)        (mod 2),
    lambda_a = sum_{W contains a, |W| odd >= 3} (1 + e_{T[W]}(a)) ,     M = sum_{|W| odd} Lambda_{T[W]} ,

`M` again symmetric with zero diagonal and zero row sums, and constant on the class. Equivalently: **for any
three vertices the eight tournaments `T^S`, `S` a subset of the triple, have `sum H(T^S) = 0 (mod 4)`.**

**Corollary 3 (recursion for the Hessian).** For `n >= 4`, `a != b`, `R = T - a - b`:

    M_T(a,b) = sum_{v in R : exactly one of a -> v, b -> v} E_R(v) + iota_a^T M_R iota_b ,

where `E_R(v)` is the parity of the number of paths of `R` on an odd number (at least 3) of vertices having `v` as
an endpoint. (Sum (**) over the odd sets `W` containing `a`, `b`; there are `2^(n-3)` of them.) Also `M = 0` for
`n = 4`, and `M = J - I + Lambda_T` for `n = 5`.

**Corollary 4 (order 6).** For `n = 6`, `M(a,b) = 1 + d_a + d_b` (`d` = score): `M` is the union of the two complete
graphs on the vertices of odd and of even score. (In Corollary 3, `R` has order 4, `M_R = 0`, and
`E_R(v) = C(d,2) + C(3-d,2)` with `d` the score of `v` in `R`, which is odd for every `d`.) **Every one of the
1,024 switching classes of order 6 has exactly 20 members with `H = 1` and 12 with `H = 3 (mod 4)`**
(`out/qf_exhaustive_n3_6.txt`; this is the whole of the mod-4 discrepancy `20480 - 12288 = 1024 * 8` at order 6).

**Corollary 5 (classes that are constant).** A class on which a quadratic function takes one value at a member
and at all its single and double switchings is constant. The previous session recorded exactly that data
(`../character_of_two_records_20261111/out/cycle_parities_pairs_24.out` and the order-20 row of its table): for
the classes of `(source) + P_23` and `(source) + D(P_11)` (order 24; 300 single and double switchings each) and
of `(source) + P_19`, `(source) + D_19` (order 20; 210 each), every odd cycle parity is constant. Hence, by
Corollary 2, **all `2^23` (resp. `2^19`) members of these four classes have the cycle parities of the normalised
member**: an even number of cycles of each odd length and `H = 1 (mod 4)` at order 24; parities `C(9,k)` and
`H = 3 (mod 4)` at order 20. (At the time this was typed "constant on the members tested".) The class of
`(source) + D23_0` is not constant (MISTAKE-763); by Corollary 2 its non-constancy is visible on single and
double switchings, which is how the audit found it.

**Corollary 6 (imbalance of a class).** On a class, `#{H = 1} - #{H = 3}` (mod 4) is `0` or `+- 2^(n-1-r/2)`,
`r` the rank of `M` on `F_2^V/<1>` (the character sum of a quadratic function; Dickson). Exhaustive: order 4:
classes constant; order 5: imbalance 4 (24 classes) or 8 (40 classes), always towards `H = 1`; order 6: always
`+8`; order 7: `-16, -8, 0, +8, +16` for 1680, 720, 5040, 19280, 6048 classes (`out/qf_exhaustive_n7.txt`).

**Remarks.** (1) `Lambda` is an `S_n`-equivariant map from labelled switching classes on `n` points to even
graphs on `n` points, two sets with `2^C(n-1,2)` elements. It is not a bijection (order 3: both classes go to
the triangle; order 4: `Lambda(u,w) = d_u + d_w`, image of size 4 out of 8), in agreement with THM-4531 /
THM-4722. The number of distinct graphs is 1, 4, 22, 91, 2972 for the 2, 8, 64, 1024, 32768 classes of orders
3-7 (`lambda_image.c`, `out/lambda_image.txt`); a class and the class of its converse have the same `Lambda`, so
every fibre is even; `Lambda = 0` for 184 of the 1,024 classes of order 6. (2) Theorem G needs no parity of `n`; Corollary 2 for `H` uses odd cycles only.
(3) Relation to THM-5127: there the cube is the set of one-vertex extensions `R + v_X`; here it is the class
`{R^Y + v_(V-Y)}`. Neither statement contains the other.

**Controls** (`qf_exhaustive.c`, `out/qf_exhaustive_n3_6.txt`, `out/qf_exhaustive_n7.txt`): on every labelled
tournament of orders 3, 4, 5, 6, 7 (8; 64; 1,024; 32,768; 2,097,152 tournaments): `Lambda` symmetric with zero
diagonal and even rows, equal at `T` and at every `T^w`; `hc(T) + hc(T^a) = 1 + e_T(a)`;
`hc(T) + hc(T^a) + hc(T^b) + hc(T^ab) = Lambda(a,b)`; all triple sums of `hc` even and of `H` divisible by 4: no
violation. Random tournaments of orders 3-10 (`endpoint_graph.py`): affine law, (**), the formula for `H mod 4`
(orders at most 8), Corollary 3 (orders 5-9).

### 2.3 The 2-adic derivative law (HYP-9669)

For a set `U` of vertices put `c_U(F) = sum_{S in U} (-1)^(|U|-|S|) F(T^S)`; these are the coefficients of the
monomials when `F(T^X)` is written as a polynomial in the 0/1 indicator of `X`. For `F = H`: `c_U` is even for
`|U| = 1, 2` (Redei), divisible by 4 for `|U| = 3, 4` (Corollary 2), and zero for `|U| > 2 floor(n/3)` (THM-4524).

**Conjecture (HYP-9669).** For every tournament and every non-empty `U`:
`v_2(c_U(H)) >= ceil(|U|/2)` and `v_2(c_U(hc)) >= ceil(|U|/2) - 1`.

What is proved about it:

- **Reduction to cycles.** By the odd-cycle formula `H(T^S) = sum over families of disjoint odd cycles of 2^(number of cycles)`,
  `c_U(H) = sum_k 2^k sum over families of k disjoint odd vertex sets W_1..W_k covering U of prod_i c_{U cap W_i}(hc(T[W_i]))`
  (a family not covering `U` contributes 0). So the statement for `hc` on all odd subtournaments implies the
  statement for `H`; conversely for odd `n` the two are equivalent by induction on `n` (the term `k = 1`,
  `W_1 = V` is `2 c_U(hc(T))`).
- **Digit lemma.** If `v_2(c_U(H)) >= ceil(|U|/2)` for all `U`, then the binary digit `k` of `H` has `F_2`-degree
  at most `2^k` on the switching cube. (Write `H(T^X)` as a sum of terms `2^j x^U` with `|U| <= 2j <= 2^j`; the
  digit `k` of a sum of numbers `2^(s_i) b_i`, `b_i` Boolean, is the sum over the sets `I` with
  `sum_{i in I} 2^(s_i) = 2^k` of `prod_{i in I} b_i` - Lucas and Vandermonde - and such a product has degree at
  most `sum 2^(s_i) = 2^k`.) This is the law of a vertex cube in THM-5127 (C) (8).
- **Equivalent form.** The conjecture holds iff for every `k` and every `2k` vertices the number
  `2^(-k) c_U(H) mod 2` is the same at every member of the class ("higher Hessians are gauge invariant");
  `k = 1` is the invariance of `M`.

Evidence: every `U` in random tournaments of orders 6-10 and three structured ones (in sample); the
pre-registered run R1 at orders 11, 12, 13, 14 (ten tournaments: seven random, two circulant, one transitive;
`out/two_adic_R1.txt`): **no violation**, and digits 1, 2, 3 of `H` have degree exactly 2, 4, 8. For odd `n` the
bound is attained at every `|U| <= 6`. No proof for `|U| >= 5`.

### 2.4 Exploratory: the rank of `M`

`rank_stats.c`, 2,000-3,000 random tournaments per order (no prediction was registered; see
`PREREGISTRATION.md`): the corank of `M` is far from that of a uniformly random alternating matrix for `n <= 8`
(order 6: always 1; order 8: `M = 0` in 1.3 percent of the samples) and within 4 points of it at orders 9, 10,
within 1.2 points at order 11. The mean class imbalance reproduces the known global discrepancy (order 11:
`0.0005 +- 0.0006` against `260546929/2^39 = 0.00047`, THM-5127).

## 3. Part B: the census at giant orders

Notation: `a(n) = A000568(n)`, `b(n) = 2^C(n,2)/n!`, `f_n = a(n)/b(n) = sum_{k >= 0} L_k(n)`, `L_k(n) = P_k(n) 4^(-kn)`.
Level `k` collects the cycle types `mu` with odd parts at least 3 and defect `m - r = 2k` (`m` = number of points
moved, `r` = number of cycles); `w_mu(n) = (n)_m / z_mu * 2^(-D_mu(n))`, `D_mu(n) = (n-m)(m-r) + C(m,2) - e_mu`,
`e_mu = sum (a_i - 1)/2 + sum_{i<j} gcd(a_i, a_j)`. All of this is inherited (section 1).

### 3.1 Controls

`census_forecast.py check`: an independent Burnside sum reproduces the recorded exact values for `n <= 40`; the
level sums add up to `f_n` exactly for `3 <= n <= 60` (`out/census_check.txt`).

### 3.2 An explicit remainder, uniform in `n` and in the truncation order

**Proposition B1.** For `n >= 3` and `K >= 0` let `T_K(n) = f_n - sum_{k <= K} L_k(n)` and

    t_k(n) = C(n,2)^(2k) 2^(3k^2 + k - 2kn) / (2k)! ,
    B_K(n) = sum_{k = K+1}^{floor(n/3)} t_k(n) + n! 2^(-ceil(n(n-1)/3)) .

Then `0 <= T_K(n) <= B_K(n)`. Moreover `k -> t_k(n)` is log-convex, so
`B_K(n) <= (n/3) max(t_{K+1}(n), t_{floor(n/3)}(n)) + n! 2^(-ceil(n(n-1)/3))`; for fixed `K`,
`B_K(n) = (1 + o(1)) t_{K+1}(n)`, which has the true exponential order `4^(-(K+1)n)` of the first omitted level
and exceeds it by a factor of order `n^(K+1)`.

*Proof.* All terms are positive, so `T_K >= 0`. A permutation with odd cycles and defect `2k` moves `m` points,
`2k + 1 <= m <= min(3k, n)`, and loses `D >= 2k(n-m) + (2/3) C(m,2)` bits (each fixed point sees `m - r = 2k` fewer
orbits; pairs of moved points lie in orbits of size at least 3 - the bit-loss lemma of Session98 / Session271,
kept by level instead of by support). The right side decreases in `m` for `m <= 3k`, so `D >= 2kn - 3k^2 - k` when
`3k <= n` and `D >= n(n-1)/3` when `3k > n`. The number of permutations of `S_n` with defect `d` is the elementary
symmetric function `e_d(1, ..., n-1) <= C(n,2)^d / d!`; the number with `3k > n` is at most `n!`. Log-convexity:
`t_{k+1} t_{k-1} / t_k^2 = 64 (2k)(2k-1)/((2k+2)(2k+1)) >= 32/3`. QED

So the level expansion can be truncated at an order growing with `n` in any way, with a certified error; there
is no optimal truncation point (the terms are positive and the sum is finite: keeping all levels is the
Burnside sum). Check (`out/census_bound.txt`): `T_K(n) <= B_K(n)` in all 395 cases `6 <= n <= 60`,
`n = 80, 100, 150, 200, 250`, `K <= 6`; `B_0(n) < 1` for `n >= 8`. Examples: `n = 100`: `T_0 = 3.22e-54`,
`B_0 = 1.22e-52`; `T_3 = 1.75e-205`, `B_3 = 6.04e-201`.

### 3.3 Certified enclosures

`a(n) = b(n) (1 + sum_{k <= K} L_k(n) + theta B_K(n))` with `0 <= theta <= 1`. Evaluated with Arb ball arithmetic
(python-flint; `census_forecast.py giant`, `out/census_giant.txt`): the exponent and the 31 leading digits
below are the unique integers contained in the computed balls. Control: for `n = 50, 100, 150, 200, 250` the
enclosure with `K = 1` contains the exact value and returns its 31 leading digits (at `n = 20, 30` the width,
`10^-11.4` and `10^-21.8`, does not determine 31 digits).

| `n` | levels kept | `a(n)` | relative width below |
|---|---|---|---|
| 300 | 1 | `5.122742023739228130944393781907... x 10^12886` | `10^-337` |
| 1000 | 0 | `7.554225356196997046558502242438... x 10^147796` | `10^-587` |
| 1000 | 1 | same digits | `10^-1175` |
| `10^4` | 0 | `1.509893164724476949547293924165... x 10^15014335` | `10^-6000` |
| `10^5` | 0 | `2.340037440184861173797367161434... x 10^1504678353` | `10^-60181` |
| `10^6` | 0 | `1.190087511215613402861178947766... x 10^150509281608` | `10^-602030` |
| `10^9` | 0 | `6.012301390415698761421190338729... x 10^150514989115770076` | `10^-602059946` |
| `10^12` | 0 | `2.039399725886987711196158477614... x 10^150514997820274377090934` | `10^-602059991268` |
| `10^100` | 0 | `5.363615985283098976914061957383... x 10^E`, `E` has 200 digits | (in the output file) |
| `10^1000` | 0 | `2.078788742790582527840139178791... x 10^E`, `E` has 2000 digits | (in the output file) |

The 16-digit forecasts of Session271 at `n = 300` and `n = 1000` agree with these digits.
Reading: every level kept buys `2 n log_10 2 - O(log n)` further digits; the integer `a(n)` itself needs every
sector whose weight `2^(o(lambda))/z_lambda` is at least about 1, which includes single `n`-cycles (circulant
tournaments, weight `2^((n-1)/2)/n`); so no truncation at a fixed level returns the exact integer (the warning
of Session98 and Session270 stands), and the low digits of `a(n)` are the census of the most symmetric
tournaments. For `P_k(n) < 4^n` (roughly `k < sqrt(2n/3)`) the numbers `1, P_1(n), P_2(n), ...` are the successive
digits of `f_n` in base `4^n`: the sequence displays its own excitation polynomials.

### 3.4 Integer heights; the moments do not determine the measure

THM-4845 (A): `a(n) = E[Z_n(W)]`, `sum_n Z_n x^n = exp(sum_{p odd} W_p x^p / p)`, `W_p = 2^(-1/2) exp(X_p)`,
`X_p = sum_{d | p} Y_d` with independent centred Gaussians `Y_d` of variance `phi(d) ln 2`.

**Proposition B2.** Let `J_d`, `d` odd, be independent integer-valued random variables with
`P(J_1 = j)` proportional to `2^(-j(j+1)/2)` and, for `d >= 3`, `P(J_d = j)` proportional to `2^(-j^2/(2 phi(d)))`.
Put `H_p = sum_{d | p} J_d`. Then `E[2^(M J_1)] = 2^C(M,2)` for every integer `M >= 0` (the number of labelled
tournaments on `M` points is the `M`-th moment of `2^(J_1)`), `E[2^(M J_d)] = 2^(phi(d) M^2 / 2)` for every integer
`M`, and

    a(n) = E[ Z_n( 2^(H_p) ) ] ,

the expected cycle index of `S_n` evaluated at random powers of two. The Gaussian field of THM-4845 and this
integer-valued field give the same value to every monomial `prod_p W_p^(m_p)` with non-negative integer
exponents, and their laws are mutually singular: **the census moments do not determine the measure.**

*Proof.* Complete the square: `Mj - j^2/(2 phi) = -(j - phi M)^2/(2 phi) + phi M^2/2`, and `phi M` is an integer, so
the shifted sum over `j` is the unshifted one; for `J_1`, `Mj - j(j+1)/2 = -(j - M + 1/2)^2/2 + (M - 1/2)^2/2` and
`sum_j 2^(-(j + 1/2)^2/2)` does not depend on the integer `M`. The exponent of a cycle type is
`e = C(M_1, 2) + (1/2) sum_{d >= 3} phi(d) M_d^2` with `M_d` the number of parts divisible by `d` (THM-516). QED

This is Stieltjes' classical phenomenon (log-normal moments are indeterminate; the discrete law is the
Stieltjes-Wigert / theta one). The theta shape of the exponent is noted in HYP-2538; the probabilistic reading
and the indeterminacy were not found in the repository (searched: `discrete gaussian`, `stieltjes`,
`indetermin`, `moment problem`, `theta` near A000568 in canon, hypotheses, results). Numerical control:
`out/census_heights.txt` (moments to 40 digits; `a(n)` for `n <= 15`).

### 3.5 Classification of the excitations and of their interactions

"Emergent recursive structure" is read here as: what a rare structure is, at which exponential rate it occurs,
with which polynomial degree, and how two of them interact. Exact statements; first moments are for the
uniform labelled tournament on `n` points.

| structure | exact first moment | rate, degree | remark |
|---|---|---|---|
| twin pair (module of size 2) | `C(n,2) 2^(2-n)` | 2, 2 | the dominant reason a tournament is not prime (Session271 bounds the probability) |
| source or sink | `n 2^(2-n)` | 2, 1 | the dominant reason it is not strong: `1 - s(n)/a(n) = n 2^(2-n) (1 + o(1))` (classical; table to `n = 60` in `out/census_strong.txt`, ratio `1.000000000` at `n = 40`) |
| module of size `s` | `C(n,s) 2^(-(s-1)(n-s))` | `2^(s-1)`, `s` | |
| both a source and a sink | `8 n(n-1) 4^(-n)` | 4, 2 | |
| automorphism that is one 3-cycle = a module that is a cyclic triangle | `(16/3)(n)_3 4^(-n)` pairs `(T, sigma)` | 4, 3 | level 1 of the census; the first correction of A000568 (Moon, Wright) |
| automorphism of cycle type `mu` (odd parts) | `w_mu(n)` | `4^k`, `m = 2k + r` | level `k = (m - r)/2`; in a level, more cycles means higher degree; top: `3^k` |

- **What a sector is.** `T` is fixed by `sigma` of type `mu` iff it is a circulant tournament on each cycle of
  `sigma` and, between two cycles of lengths `a` and `b`, the arc from `x` to `y` depends only on the class of the
  difference of their positions modulo `gcd(a,b)`. The free bits are `(a-1)/2` per cycle (its connection set),
  `gcd(a,b)` per pair of cycles, one per pair (cycle, fixed point) and one per pair of fixed points; this is
  the count of pair orbits. A single cycle with fixed points is a module carrying a circulant tournament.
- **Interaction.** For two cycle types, `w_{mu + nu}(n) / (w_mu(n) w_nu(n))` tends to
  `(z_mu z_nu / z_{mu + nu}) 2^(I(mu,nu))`, `I(mu,nu) = sum_{i,j} (a_i b_j - a_i - b_j + gcd(a_i, b_j))`. Every pair
  of cycles attracts, with at least `(a-1)(b-1)` bits; two 3-cycles: `2^6 = 64` (Session271). The energy is a sum
  of one-cycle and two-cycle terms only, which is why the Burnside sum is exactly a gas with pair interactions
  (Session249) and a Gaussian average (THM-4845).
- **Selection rules.** No even cycle (a power of it reverses an arc); hence no sector of defect 1 and the
  first level costs `2n - 4` bits, against `n - 2` for graphs (transpositions): poles at `4^k` for tournaments, `2^k`
  for graphs (Session98, palette classification).
- **How the axes interact.** A 3-cycle automorphism is a module; the strong components form a transitive
  chain of modules and `Aut(T)` is the product of the automorphism groups of the components; on switching
  classes automorphisms of even order appear (level permutations, THM-479), next paragraph.

### 3.6 The census of switching classes

THM-479 (Babai-Cameron): `A049313(n) = N_odd(n) + N_lev(n)`, `N_odd(n) = sum over odd cycle types of 2^(o - l + 1)/z`
(`o` = number of pair orbits, `l` = number of cycles), `N_lev` the sum over the other level permutations.

**Proposition B3.** (i) `N_odd(n) = (2^C(n-1,2)/n!) sum_{k >= 0} P_k(n) 4^(-k(n-1))`, with the same polynomials
`P_k` as for tournaments: the census of switching classes is the tournament census with the clock `n - 1`.
(ii) `N_lev(n) = 0` for odd `n` and `0 <= N_lev(n) <= 2^(n^2/4)` for even `n`. (iii) Hence
`A049313(n) n!/2^C(n-1,2) = 1 + (64/3)(n)_3 4^(-n) + ...`, four times the first correction of the tournament
census, with remainder at most `sum_{k > K} 4^k t_k(n) + 4^((n-1)/2) n! 2^(-ceil(n(n-1)/3)) + 2^(n^2/4) n!/2^C(n-1,2)`.

*Proof.* (i) A type at level `k` has `l = n - m + r` cycles, so its weight relative to the identity is multiplied
by `2^(m - r) = 4^k`. (ii) A level permutation with an even cycle has only even cycles, so `n` is even; it has no
fixed point, so every orbit of pairs has size at least 2 except the pairs exchanged by a 2-cycle, at most
`n/2 + (C(n,2) - n/2)/2 = n^2/4` orbits; the classes fixed by `sigma` are at most `2^(number of pair orbits)` (for each
of the `2^(n-1)` sign vectors the fixed tournaments form an affine space of that dimension at most). QED

Check (`out/census_switching.txt`, values of A049313 from the OEIS): the identity (i) is exact for
`3 <= n <= 19`; `N_odd(n) = A049313(n)` for every odd `n <= 19`; `N_lev = 1, 2, 15, 280, 23464, 6819360, 6981258440,
25417204478336` for `n = 4, ..., 18`; `A049313(19) n!/2^C(18,2) = 1.000000451 = 1 + (64/3)(19)_3 4^(-19)` to the
digits shown. No asymptotic statement for A049313 was on record.

## 4. Part C: Chan-Pak side

- **What Part A adds.** The second binary digit of `H` is now known on a whole switching class from one member
  and its single and double switchings; classes are constant, balanced or of imbalance `2^(n-1-r/2)`; at order
  6 the law is the same for every class. This produces no new value of `H` and excludes none.
- **Forecasts for order 12** (registered in `PREREGISTRATION.md` R2 before any attempt to score; HYP-9670):
  least odd number other than 7, 21 that is not `H` of a 12-vertex tournament: 488,000, interval
  [476,000, 500,000]; maximum of `H`: 545,000, interval [531,205, 552,000]; number of distinct values: 262,000,
  interval [250,000, 272,000]. NUMERICAL extrapolations from the censuses of orders 8-11 (THM-4889); the method
  is in the registration.
- **Why 7 and 21 look like the last accidents (HEURISTIC, not a proof).** An odd `h` is "in the bulk" of the
  orders `n` with `n!/2^(n-1)` comparable to `h`. There the number of classes per odd value is about
  `a(n)` divided by half the width of the distribution: about 1.5 at order 5 (12 classes on the odd values up
  to 15: 7 is missed), 2.4 at order 6 (56 classes up to 45: 7, 21, 35, 39 are missed), 5 at order 7, 20 at
  order 8, and about `9 x 10^19` at order 18, the bulk order of `h` near `2^35`. Under a Poisson model the chance
  that a bulk value is missed is `exp(-density)`: 22 and 9 percent at orders 5 and 6, below `10^-8` from order 8
  on. The values 7 and 21 are bulk values of the only two orders where the density is of order 1. This is
  consistent with the verified ranges (every odd number other than 7 and 21 up to `2^32`, THM-4875, and up to
  `2^35` in one run, THM-4885) and proves nothing; the
  rigorous form of the statement is the local-limit hypothesis HYP-9659 (A).

## 5. Part D: reading the Yang-Mills construction

### 5.1 What was read, and how

Source: `yang-mills-group/OpenAI-YM-20261028` at commit `3456e5e75cf8d4cd558955a5c7638fba2d68f464` (the commit
pinned by Session272), cloned read-only into a scratch directory; nothing was built, run or verified, as the
owner asked. Read by the author in full: `README.md`, `ChallengeSpec.lean`, `PublicWilson.lean`,
`PublicTheory.lean`, `PublicView.lean`, `Construction.lean`, `Solution.lean`, `ConcreteSU2.lean`,
`WilsonBridge.lean`, `HierarchyBridge.lean`. The directories under `YM/` (9,462 files) were read through two
read-only survey agents (an architecture map and a probe of the five open items), whose reports quote
declarations with file names and say where a role is inferred from a docstring only; statements below that rest
on those reports are marked (S). File names are relative to the source root.

The statement, in one paragraph. For every compact simple gauge group there is a family of Wilson lattice
measures on periodic 4-tori of side `period(p) * denominator(N)` and mesh `1/denominator(N)`, with coupling
`beta(N,p)` such that `beta/log(1/mesh)` tends to `2 b_0`, `b_0 = 11 C_A/(48 pi^2)`, and gauge-invariant local
insertions with a fixed stencil, such that for separated test functions the mesh limit exists at each fixed
volume and, along a subsequence of any given sequence of volumes, the volume limit exists and is the
"unselected" member of a hierarchy: a compact metric space of characters with a probability measure of full
support, correlation functions continuous in the character whose mixture is the unselected limit, a Wightman
theory (W0-W3, unique vacuum) for a selected character, an operator product expansion, and the short-distance
law `r^8 log^2(1/r) <Theta Theta> -> 3 dim(g)/pi^4` for the centred stress trace. The README states what is
outside: a mass gap, the exclusion of continuous-mass generalised free fields, uniqueness over contact
prescriptions, and "the volume limit allows refinement". (S): the scheme is a Balaban-type block-spin
renormalisation group with a soft Markov blocking kernel (`YM/RG/Block.lean`), a small-field / large-field
split, polymer activities with weights `exp(kappa * cover number)`, and a coupling tuned for each `(N, p)` so that
the effective coupling at the unit scale is a fixed number (`YM/RG/History.lean`, `reaches`).

### 5.2 Fourteen analogies, typed

These are different from the sixteen connections of Session272 (which concern observers, phases, flux and
spectral visibility). Types: IDENTITY (the same mathematical mechanism), STRUCTURAL (the same algebraic or
logical shape, different objects), ANALOGY (a dictionary with a precise statement on our side), WEAK.

| | in the Yang-Mills source | in this repository | type |
|---|---|---|---|
| Y1 | large-field gas `sum over admissible families of prod K_i` (`YM/RG/Activities.lean`); Ursell coefficients and a Penrose bound by parent maps (`YM/Push/ABSOLUTE_SUMS`) (S) | `H(T) = I(Omega(T), 2)`, the hard-core gas of odd cycles (THM-002), its Mayer cumulants (THM-517) | STRUCTURAL; opposite regimes: activities `exp(-monitor^(1/5))` against fugacity 2, which is outside every cluster-expansion region (INV-028) |
| Y2 | exact step `rep(m).eval = push(rep(m+1).eval)` with a kernel of total mass 1, in a fixed format (`History.lean`, `Transfer`) (S) | the proof of Theorem G: the pair `(e, Lambda)` at order `n` from orders `n-1`, `n-2` by (*), (**); the format "affine, symmetric" is preserved | STRUCTURAL: an exact recursion that closes on a format |
| Y3 | small-field / large-field split (`flags`, `Healthy`, `forcedAllowance`) (S) | the unlabelled law is `L/(1+delta) + delta D/(1+delta)`: rigid bulk plus symmetric defects of weight `(16/3) n^3 4^(-n)` (Session271) | STRUCTURAL |
| Y4 | convergence criterion for polymer sums, `(2M+1)^4 gamma <= 1` (`YM/Ac`, `forest`) (S) | `n 2^(-(n-1)/3) < 1` (Session98) and Proposition B1: number of supports times activity | STRUCTURAL |
| Y5 | Moebius pieces of `-log Xi` on regions with exponential weights (`YM/SectorNext/VACUUM_NORMALIZATION/Sprig.lean`) (S) | the pieces `c_U` of `H` on a switching class, of 2-adic size `2^(-ceil(|U|/2))`: a cluster bound in the 2-adic absolute value (Corollary 2 for `|U| <= 4`, HYP-9669 beyond); digit `k` of `H` is an observable of range `2^k` | ANALOGY, with a theorem and a conjecture on our side |
| Y6 | characters = Gelfand spectrum of the compressed algebra on the zero-energy space; mixture identity; continuity in the character (`YM/Phase/Mixture.lean`, `YM/Ground/Spectrum.lean`) (S) | functions on a switching class = group algebra of `(Z_2)^(n-1)`; characters = Walsh characters; the mixture identity is "the class mean of `H` is `n!/2^(n-1)`" (THM-1467, THM-4524); Corollary 2: the second digit has character degree 2 | ANALOGY (finite, literal on our side) |
| Y7 | charges of the stress tensor implement the symmetry on local fields (`YM/EveryPhase/Charge.lean`, `local_action`) (S) | Corollary 1: the first gauge variation of `hc` is the endpoint operator, the second is the invariant tensor `Lambda`, the third vanishes | ANALOGY (a Ward hierarchy that stops at order 2) |
| Y8 | every gauge-invariant local field of weight below 4 is a multiple of the unit (`YM/Trace2/Transfer.lean`, `first_grade`) (S) | no sector of defect below 2: orientation forbids even cycles; poles `4^k` instead of `2^k` | ANALOGY: an invariance removes the lowest levels |
| Y9 | tree (axial) gauge inside blocks (`YM/Tree/Pines.lean`) (S) | THM-474: one member of each class contains a given spanning path; the class sum over spanning trees (THM-1467) | IDENTITY of mechanism (`Z_2` gauge on `K_n`, tree gauge); KNOWN here |
| Y10 | two limits in a fixed order, "not along a single diagonal"; coupling tuned for each volume (`YM/Extraction/Orders.lean`) (S) | the census has two parameters `(n, K)` and Proposition B1 is uniform in both, so every diagonal converges: the toy has the estimate (a "moving-volume error bound", two-clock lemma of Session272) that the construction lacks | STRUCTURAL |
| Y11 | separated correlations; contact prescriptions not unique (README) | Proposition B2: all the moments, two mutually singular measures | exact toy (PROVED) of "the correlations do not fix the measure" |
| Y12 | observables with a fixed stencil; no loop observables | `kappa(T) = H(T) 2^(n-1)/n!`: random tournaments have `kappa -> 1` in probability (variance `2/n`, THM-589, THM-4761), Paley tournaments `kappa -> e` (THM-438); both sequences are quasi-random, with the same local limit | exact toy (inherited theorems): a spanning observable is not a function of the local limit; its expansion in local counts (THM-4757) does not commute with the limit |
| Y13 | `<Theta Theta>` is not a finite sum of mass-shell Wick profiles (`YM/TraceExtension/Wick/Refute.lean`) (S) | poles of order `3k + 1`; `f_n` is not P-recursive (Session98) | WEAK: both say "not in a finite-dimensional class" |
| Y14 | reflection positivity and chessboard estimates (`YM/Slice/Reflection.lean`, `YM/Cells/Chess.lean`) (S) | the fold of a self-converse tournament (THM-5136): the crossing relation `a ~ b iff a -> psi(b)` is symmetric | OPEN question, with a negative remark: the natural kernel of the crossing arcs is `C tensor (J - C)`, which is not positive semidefinite in general (`C = I`); no chessboard inequality for `H` is claimed |

### 5.3 The five open items

For each: (a) what the source contains (S), (b) what blocks, (c) an exactly solved instance in the repository,
(d) a remark. Session272 has a program for each item; (d) is meant to be complementary to it. Nothing below is
a result about Yang-Mills.

**1. Convergence and uniqueness of the infinite-volume limit.**
(a) Full mesh limit at each fixed volume (`YM/FullVolume/Profiles.lean`, `PositionsUV` with `Sequence.full`;
`YM/FullVolume/Glide/Vac.lean`, `Winter.convergence`: effective couplings and activities at a fixed level and
period converge); bounds uniform in all periods (`YM/FullVolume/Discrete.lean`, `FixedBounds`); block steps shared
between volumes (`path_shared`); period errors `exp(-cN)` for Hessians of log-determinants
(`YM/Further/PHYSICAL_PERIOD_LIMIT_F933`). The volume limit is taken by Arzela-Ascoli
(`YM/Extraction/Lipschitz.lean`, `first_surface`).
(b) No statement compares two periods; the coupling is tuned separately for each volume; there is no analysis
above the unit scale, where the effective theory is a weakly coupled lattice gauge theory in a growing box.
(c) The census: Proposition B1 controls both parameters at once. Warning from Y12: uniqueness has to be stated
for a named class of observables - two sequences with the same local limit have different limits of `kappa`.
(d) Uniqueness for local observables would follow from a comparison of the unit-scale effective measures of two
periods on a common box, with an error decaying in the distance to the boundary; that is a clustering
statement at the unit scale, so items 1 and 3 are one problem at that scale.

**2. A probability measure for the continuum theory.**
(a) Correlations exist on separated configurations only (`YM/Euclid/Table.lean`); a contact extension as
distributions (`YM/Contacts/Solstice/Result.lean`), without positivity or uniqueness; growth in the arity
`(n/r)^(4(D+1)n)` (`YM/OS/Slabs.lean`, `AllCopy`); all fields have weight at least 4.
(b) `r^(-8) log^(-2)` is not locally integrable in dimension 4, so even second moments of one smeared field need
the contact prescription; and a growth `n^(An)` with `A > 1` is beyond Carleman's condition, so the moments
could not fix a measure even if it existed.
(c) Proposition B2 is that situation solved exactly: the measure exists, the moments grow like `2^(c k^2)`, and
there are at least two measures.
(d) SUGGESTION (conditional on the reading (S) of `Transfer` and `Winter.convergence`; not checked). The blocked
fields are honest random variables with values in compact groups. At a fixed volume, for each level `m` the
law `mu_m` of the level-`m` blocked configuration has density `rep(N,p,m).eval`, and `mu_m` is the image of `mu_(m+1)`
under the blocking kernel. The structure `Winter` (`YM/FullVolume/Glide/Vac.lean`, read by the author; its
constructor is in the cone of the audited theorems according to `verification/declaration-usage.tsv`) asserts
that, at a fixed period and level, the level-`m` coupling and the activities, in a weighted space, converge as
`N -> infinity`. If this gives weak convergence of the normalised laws `mu_m` (an additional step that was not
checked: the constant in the exponent diverges and has to cancel in the normalisation) and the kernels are
continuous, the limits satisfy the same relation, so the joint
laws of `(V_0, ..., V_M)` are consistent and Kolmogorov's theorem gives **a probability measure on the compact
space of sequences of blocked configurations of all levels**. This would be a continuum measure at finite
volume on an explicit space, with bounded coordinates (no moment problem). What it would not give: the relation
between its coordinates and the renormalised local fields (the source says its insertions are "the normalized
readout, not the actual value of a coarse field"), and the volume limit.

**3. A uniform mass gap.**
(a) Positive energy and a unique vacuum in every character; the only rates are `nu(0 < E <= u) <= M u^2` and a heat
bound `M s^(-2)` (`YM/EveryPhase/Extension.lean`, `Completed`), which follow from locality and cone support;
exponential decay exists only inside the block steps (activities, covariances), at scales up to the unit scale.
(b) Nothing above the unit scale; the limit comes from compactness, which keeps no rate. Session272 shows by a
two-level countermodel that the polynomial bounds, cyclicity and irreducibility together do not give a gap.
(c) The census has an exact spectrum: `f_n = sum_k P_k(n) e^(-k n log 4)`, gap `log 4`, uniform in `n`, because an
excitation costs 2 bits for every spectator vertex. The Yang-Mills analogue of "a cost proportional to the
size of the system" is a lower bound for the free energy of a flux sector proportional to its cross-section;
this is the known heuristic link between the gap and the area law, not a new one.
(d) No lemma is proposed. The point visible in the source: the analysis needs weak coupling at every level
including the unit scale (`delayed`, `smallField`), and a convergent expansion in infinite volume needs the
opposite.

**4. Wilson loops and the area law.**
(a) Absent as observables. Building blocks: products of links along walks (`YM/Lattice/Path.lean`), block
connectors (`YM/RG/Block.lean`); lattice reflection positivity and chessboard bounds hold for every bounded
gauge-invariant function supported in a half-space (`YM/Cells/Cut.lean`), hence for lattice Wilson loops.
(b) The insertions have a fixed stencil in lattice units and shrink to points; no perimeter renormalisation.
(c) Session34 (triangular Wilson action on `K_n`) and Session272 (disk and torus laws) are the exact toys on
record. From this note: on a switching class the count of Hamiltonian cycles is not gauge invariant, and
Corollary 1 says exactly how it fails - a quadratic law whose Hessian is invariant.
(d) SUGGESTION (same condition as in item 2). A Wilson loop of the **unit-scale blocked field** along a loop of
the level-0 lattice is a bounded gauge-invariant variable of fixed physical size (the level-0 spacing is the
unit of length), smeared by the blocking itself. Its expectation is an integral against `rep(N,p,0).eval`, so no
insertion machinery is needed, and its mesh limit at fixed volume would follow from the convergence of the
level-0 effective density. These are natural candidates for "Wilson loop observables"; their relation to the
local fields, and any area law, are untouched.

**5. Single-particle states.**
(a) A joint spectral measure with cone support and Lorentz covariance for the two-point function of every
field (`YM/TraceExtension/Damping.lean`, `Momentum.lean`), its ultraviolet law for the stress trace
(`Endpoint.lean`, `HighEnergy`), exclusions of finite Wick polynomials and of finitely many massless Fock
copies.
(b) Without a gap there is no isolated mass shell; nothing is stated about atoms in the mass variable.
Session272 exhibits spectral measures with the same ultraviolet coefficient and three different infrared
behaviours.
(c) The census has one particle: the cyclic-triangle module, of "mass" 2 bits per spectator. A cycle of
length `2h + 1` has the rate of `h` three-cycles and a lower degree (`n^(2h+1)` against `n^(3h)`): every other
excitation sits exactly at a multi-particle threshold, and the factor 64 is an interaction at threshold, not a
bound state below it.
(d) No lemma is proposed.

## 6. Pre-registration

`PREREGISTRATION.md`: R1 (2-adic law at orders 11-14) registered and pushed (9f260a24f7) before the runs,
CONFIRMED; R2 (order-12 forecasts) not scorable today; R3 exploratory. Statements typed PROVED were not
registered; their exhaustive and random checks are controls.

## 7. Reproduction

    cc -O2 -o qf_exhaustive qf_exhaustive.c && for n in 3 4 5 6 7; do ./qf_exhaustive $n; done     # order 7: about 1 minute
    python3 endpoint_graph.py                                                                      # random checks, orders 3-10
    cc -O2 -o two_adic_cube two_adic_cube.c && ./two_adic_cube 11 1 0                              # R1 (arguments in PREREGISTRATION.md)
    cc -O2 -o rank_stats rank_stats.c && ./rank_stats 9 3000 1
    cc -O2 -o lambda_image lambda_image.c && for n in 3 4 5 6 7; do ./lambda_image $n; done
    python3 census_forecast.py all                                                                 # needs python-flint, mpmath

Outputs are in `out/`.

## 8. What is not claimed

- No value of `H` is produced or excluded; Chan-Pak is OPEN.
- HYP-9669 is a conjecture beyond `|U| = 4`.
- The order-12 numbers are extrapolations.
- Section 5 contains analogies and two conditional suggestions; the source was not verified, the deeper files
  were read through survey agents, and no statement about four-dimensional Yang-Mills theory is proved here.
