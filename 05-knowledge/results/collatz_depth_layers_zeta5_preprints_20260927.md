# The depth of the rooted component, its layers and its Hankel rank; two September-2026 zeta(5) preprints read by their equations; and the themes they share with the Collatz thread

**Session:** opus, `collatz-posets-zeta5-20260927` (second note; the first is
[`collatz_posets_dags_zeta5_20260927.md`](collatz_posets_dags_zeta5_20260927.md)), 2026-09-27.
**Owner's directive:** "keep going with a long session striving for Collatz;
consider the methods in arXiv:2609.22316 and zenodo.org/records/22826419,
read equation deeply and look for creative or surprising connections to our
methods, even if they feel a bit like numerology at first; let pattern
seeking run freely; deconstruct themes and abstract strategies similar to
the pasted content" (the pasted content: the depth `d(n) = min{k : C^k(n) in
{1,2,4}}` of the rooted component, its layers `0: {1,2,4}, 1: {8}, 2: {16},
3: {5,32}, 4: {10,64}, 5: {3,20,21,128}`, the additivity `d(2^x m) = x +
d(m)`, and the target "every odd number belongs to the component").
**Inherits:** the first note (Apéry forms `A_L = 3^L n + S_(L-1) = 2^(d_L)
m_L`, the exponent `nu_L`, the series `F_n(z) = sum 2^(d_k) z^k`, Propositions
7–10), THM-4476, THM-4512, the parallel note
[`collatz_posets_dags_20260927_spine_descent_tree.md`](collatz_posets_dags_20260927_spine_descent_tree.md)
(descent tree; first-descent words as 2-adic-to-3-adic affine bijections),
the HARD-class note (Theorems S, D, Y), Krasikov–Lagarias 2003 and
Applegate–Lagarias 1995 (from memory), Zudilin 2001/2002, Hermite's formula
for the Hurwitz zeta function and the Binet integral (classical).

**Status: PROVED elementary (Propositions 1–3) + FINITE-EXACT (layers to
depth 60, bigraded layers to `d = 10`, the Hankel-rank law for odd `n <=
3000`, first-descent margins on five record orbits, thirty-digit checks of
the two identities behind the Zenodo functional) + READ (both preprints,
classified AUTHOR-CLAIMED / UNREFEREED; one refuted as written by its own
equation (13)) + THEMES (section 5) + NUMEROLOGY (section 6, tested, all
dismissed) + DIRECTION (section 7). Collatz OPEN. No finishing move.**
Script `04-computation/experiments/collatz_depth_layers_20260927.py`, output
`.out` beside it. Independent audit: OWED (the session's audit subagent was
cut off by an API rate limit; section 8).

## 0. The answer in one paragraph

The depth of the pasted content is exact and has three faces. (i) Its layers
grow like `lambda_C^k` with `lambda_C = (1 + sqrt(7/3))/2 = 1.2638` per
Collatz step (the `4/3` per Terras step of the tree-search literature, re-timed
because a `(n-1)/3` predecessor costs two Collatz steps): measured
`1.2639–1.2647` at depths 50–60, and the owner's layer table is reproduced.
(ii) Its odd-step version is a rank: the linear complexity of the sequence
`2^(d_k(n))` of 2-parts of the orbit's Apéry forms, continued along the
1-cycle, equals `d_odd(n) + 1` for every odd `n <= 3000`, so by Kronecker's
theorem the periodicity conjecture for `n` says the Hankel determinants
`det(2^(d_(i+j)(n)))` vanish from some size on, and the depth is that size
minus one. (iii) Its coverage statement is the conjecture, and the layers are
the backward (3-adic) certificates dual to THM-4512's forward (2-adic)
cylinders. The two preprints: arXiv:2609.22316 (Suman, math.GM, "at least one
of `zeta(5), zeta(7), zeta(9)` is irrational") reproduces Zudilin's 2002
well-poised construction with `q = 11`, `eta_0 = 3`, `eta_j = 1`; its decay
constant `C_0' = 5.7535` is correct (recomputed from its equations (33)–(34)),
but its prime-window product `Phi` (equation (13)) runs over `h_0 < p <= n`
with `h_0 = 3n + 2`, an empty range, so `Phi = 1`, the denominator exponent is
`9`, and the decisive inequality of its Lemma 2 reads `5.75 > 9`: false as
written; the most favourable window `(n, 3n+2]` saves `1.33–1.44`, not the
`3.25` needed. The Zenodo preprint (Fauzan, "zeta(5) is irrational", 32
pages) is a Hankel-determinant construction: a functional `mu_X` on rational
functions with poles at `t = -j^2`, whose polynomial values are Bernoulli
numbers and whose pole values are `j^4 (X - H_j^(5)) - 1/4 + 1/(2j)`; at `X =
zeta(5)` both are moments of the positive weight `(1/12) u^5 (d/du)^4 [1/(e^(2
pi u) - 1)]` in `t = u^2` (the pole values are Hermite's integral for the
Hurwitz zeta, verified to thirty digits), positivity gives nonvanishing,
Heine's multiple integral and a logarithmic-energy bound give decay
`exp(-c K^2)`, and the claim rests on a p-adic normalization (its Sections
3–5) that this note cannot check. The themes these share with the Collatz
thread are exact and all point the same way: two rates compete (Collatz:
`d_L/L` against `log_2 3`, margin `0.006–0.04` bits per odd step at the first
coefficient descent of the record orbits, against a generic drift `0.415`);
nonvanishing is the hard half of an Apéry proof and the free half of Collatz,
where `A_L > 0` trivially and smallness is everything; partial sums and
tails (`H_j` against `zeta(5)`, remainder `~ 1/(4 j^2)` after scaling) are the
shape of `r_L` against `xi` with remainder `eta_L 2^(d_L)/3^L`, except that
the Collatz remainder is the `3x-1` copy and does not shrink; determinants
amplify linear smallness to quadratic smallness in both worlds, trivially in
Collatz (entrywise divisibility) and non-trivially for zeta(5). The
numerology (`3^7 - 2^11 = 139` is the `-17` cycle's denominator and the
Zenodo decay constant is `139/5`; `37` lies on that cycle and is the Zenodo
degree slope; Suman's `q = 11` and seven maxima are the cycle's `(A, p) =
(11, 7)`) was tested and holds no content: no equation of either preprint
involves `2^11`, `3^7` or a Collatz object.

## 1. The depth and its layers (PROVED + FINITE-EXACT)

Let `C(n) = n/2` (n even), `3n+1` (n odd), and `d(n) = min{k >= 0 : C^k(n) in
{1,2,4}}` on the rooted component.

**Proposition 1.** (a) `d(2n) = d(n) + 1` for every `n` in the component
with `n not in {1, 2}`; hence `d(2^x m) = x + d(m)` for odd `m > 1`. (b) The
layer `L_k = {n : d(n) = k}` is obtained from `L_(k-1)` by `n -> 2n` always and
`n -> (n-1)/3` when `n = 4 mod 6`, `n > 4`. (c) Heuristically the layer sizes
grow like `lambda_C^k` with `1 = 1/lambda_C + (1/3)/lambda_C^2`, i.e.
`lambda_C = (1 + sqrt(7/3))/2 = 1.26376`, because a halving predecessor costs
one step and always exists while a `(n-1)/3` predecessor costs two steps and
exists for one residue class in three; this is the `4/3` per Terras step of
Applegate–Lagarias re-timed in Collatz steps.

*Proof.* (a) `2n` is not in `{1,2,4}` unless `n in {1,2}`, and `C(2n) = n`.
(b) is the definition of predecessors. (c) is the branching heuristic; it is
not a theorem (the exact growth of the pruned trees is bracketed, not known).
∎

FINITE-EXACT (script part 1): the owner's table is reproduced exactly (`L_5 =
{3, 20, 21, 128}`); `|L_k|` for `k = 0..60` is `3, 1, 1, 2, 2, 4, 4, 6, 6, 8,
10, 14, 18, 24, 29, 36, 44, 58, 72, 91, 113, 143, 179, 227, 287, 366, 460,
578, 732, 926, 1174, 1489, 1879, 2365, 2988, 3780, 4788, 6049, 7628, 9635,
12190, 15409, 19452, 24561, 31025, 39229, 49580, ...`; the ratios
`|L_k|/|L_(k-1)|` at `k = 40, 50, 60` are `1.2652, 1.2636, 1.2639` against
`lambda_C = 1.2638`; `d(2n) = d(n) + 1` holds on the whole computed component.

**The bigraded layers.** In Syracuse coding the odd-step depth alone does not
give finite layers (every odd `m` not divisible by 3 has the infinite
predecessor family `(2^v m - 1)/3`, `v` of one parity), so the finite object
is the pair (odd steps `d`, halvings `D`): odd `n` reaching 1 with word
`(v_1, ..., v_d)`, `sum v = D`, i.e. `n = (2^D - S)/3^d`. Script part 2 counts
these for `d <= 10`, `D <= 2d + 6`: `3, 4, 6, 6, 8, 11, 14, 20, 21, 26` pairs,
against the naive density heuristic `sum_D C(D-1, d-1) 3^(-d)` of `2.7, 5.0,
8.1, 12.4, 18.0, 25.5, 35.4, 48.7, 66.4, 90.0`. The count grows more slowly
than the heuristic because the congruence `2^D = S mod 3^d` is a triangular
system of parity and higher conditions on the valuations (the last valuation
has the parity of `D - d_(d-1)` forced, and so on), not an independent
`3^(-d)` event, and because words that pass through 1 before the end are
excluded. This is the backward face of the parallel note's Proposition 4: a
first-descent word maps a 2-adic source class onto a 3-adic landing class,
and the layers are the fibres over the landing point 1.

**Coverage.** The union of the layers is the tree of 1; the conjecture is
that it is everything. Krasikov–Lagarias give at least `X^0.84` elements
below `X`; the computation gives all odd `n <= 10^6` (first note, part 5).
The layers are backward certificates: every element comes with a verified
route. THM-4512's coefficient-descent cylinders are forward certificates:
every element of a certified 2-adic class above its threshold descends. The
two systems are dual (2-adic source classes, 3-adic landing classes) and
neither reaches the other's complement; "every unresolved orbit enters a
certified region" (the S11 directive) is exactly the statement that the
forward system's residual meets the backward system's union, and that is
the conjecture.

## 2. The depth is a Hankel rank (FINITE-EXACT law; Kronecker)

For odd `n` let `A_k = 3^k n + S_(k-1) = 2^(d_k) m_k` (first note, Proposition
3) and consider the sequence of 2-parts `u_k = 2^(d_k(n))`, `k >= 0`,
continued along the 1-cycle (`v = 2` for ever once `m = 1`).

**Proposition 2 (Kronecker).** `F_n(z) = sum u_k z^k` is rational iff the
Hankel determinants `det(u_(i+j))_(0 <= i,j < N)` vanish for all `N >= N_0`;
the least such `N_0` is the linear complexity `LC(u)` (the order of the minimal
linear recurrence). Hence the periodicity conjecture for `n` (first note,
Proposition 7) says `LC(u(n)) < infinity`, and a divergent orbit would have
infinite Hankel rank.

*Proof.* Kronecker's theorem on rational generating functions. ∎

**FINITE-EXACT law (script part 3).** `LC(u(n)) = d_odd(n) + 1` for every
odd `n <= 3000`, where `d_odd(n)` is the number of odd steps to 1
(Berlekamp–Massey over `F_p`, `p = 2^61 - 1`, sequences of length `3
d_odd(n) + 12`). The mechanism is generic, not arithmetic: a sequence whose
prefix of length `t` is followed by a geometric tail has linear complexity
`t` (checked on random prefixes of powers of two: `LC = t` for `t = 5, 20,
41`), and the Collatz prefix has `d_odd(n) + 1` terms `u_0, ..., u_(d_odd)`
before the tail ratio `4` sets in; the last prefix increment is `v = 4` (the
step `5 -> 1`), never `2`, so the prefix never merges with the tail. So the
depth is a rank. What this buys: a reformulation in the shape of the Zenodo
preprint's engine (Hankel determinants), with the direction inverted once
more: there, nonvanishing of the determinant is the goal and positivity
supplies it; here, vanishing is the goal, and positivity can only obstruct
it. A Collatz Hankel matrix is a moment matrix only when the sequence is
log-convex, i.e. when the valuations are non-decreasing, which happens on
no divergent word (a bounded non-decreasing integer sequence is eventually
constant, and constant valuations mean `n = -1 mod 2^infinity` or the
1-cycle).

## 3. arXiv:2609.22316 (Suman) read by its equations (READ; refuted as written)

The paper takes Zudilin's 2002 well-poised construction (Izv. Math. 66) with
`q = 11`, `r = 3`: `R(t) = n!^2 (3n + 2 + 2t) (t+1)_n^3 (t+2n+2)_n^3 / (t+n+1)_(n+1)^8`
(its equation (17)), the linear form `F(h) = (1/2) sum R''(t)` in `1,
zeta(5), zeta(7), zeta(9)` (even zetas cancel by the well-poised symmetry
`R(-t - h_0) = -R(t)`, `zeta(3)` by the residue sum), denominators `D_n^9`,
and an arithmetic saving `Phi = prod_(h_0 < p <= m_8) p^(v_p)` with `v_p` a
minimum of floor-sums (its (13)–(15)). Its decay rate is `C_0' = -Re f_0(tau_0)`
with `tau_0` the root of `(tau-3)^3 (tau-1)^11 - tau^3 (tau-2)^11` of maximal
real part below 3 (its (33)–(34)); its Lemma 2 needs `C_0' > C_2'` where `C_2'
= 9 - (saving)`.

Checked (script part 5): (a) `tau_0 = 2.86852453 + 0.11960091 i` and
`f_0(tau_0) = -5.75349395 - 8.95071225 i`, so `C_0' = 5.7535`: the paper's
numbers reproduce exactly. (b) With `eta_0 = 3`, `eta_j = 1`, the paper's own
(10) and (16) give `h_0 = 3n + 2` and `m_1 = ... = m_8 = n`, so the range `h_0 <
p <= m_8` of (13) is empty: `Phi = 1`, the denominator is `D_n^9`, `C_2' = 9`,
and Lemma 2's inequality is `5.75 > 9`, false. (c) The paper's p-adic
bookkeeping mixes ranges: its (31) uses `ord_p D_n^9 = 9` for the primes of
`Phi`, but for `p > h_0 > n` one has `ord_p D_n = 0`. (d) Even the most
favourable reading, a window `(n, 3n+2]`, gives a saving `(1/n) sum_p v_p log
p = 1.44, 1.33, 1.44` at `n = 100, 200, 400`, against the `3.25` that `C_0' >
C_2'` needs (the paper's `C_2' = 5.2769` would need `3.72`). Verdict: the
theorem is not supported by the paper's equations; Zudilin's four-number
theorem stands and the three-number statement is open. (A previous
`zeta(5)` claim by the same author is marked "found incorrect" on arXiv.)

What is worth keeping from the construction: the *shape* of the arithmetic
saving. The valuations `v_(k,p)` are sums of floor functions of linear forms
in the parameters (Legendre), the saving is a periodic sawtooth `phi(x)`
integrated against the prime density (`int phi dpsi`), and the competition
is between an analytic rate (the saddle point `tau_0`) and an arithmetic rate
(`9 - saving`). Section 5 places this against the Collatz thread.

## 4. The Zenodo preprint (Fauzan) read by its equations (READ; AUTHOR-CLAIMED)

Pages 1–4 (read in the browser's text layer). `K = 40n`, `N = 3n`, `h = 37n`,
`D_m(t) = prod_(j<=m) (t + j^2)`, `H_j^(5) = sum_(v<=j) v^(-5)`. A `Q[X]`-valued
functional on rational functions with simple poles at `t = -j^2`:

    mu_X(t^e) = (-1)^e B_(2e+2) (2e+3)(2e+4)(2e+5)/24,            (its 2.2)
    mu_X(1/(t + j^2)) = j^4 (X - H_j^(5)) - 1/4 + 1/(2j),          (its 2.3)

extended by polynomial division and partial fractions; `X` enters only through
the pole values, affinely. The Hankel matrix `G_K(X) = [mu_X(D_N(t)^6 t^(i+j) /
D_K(t))]_(0 <= i,j < h)` (poles `-j^2` for `N < j <= K`, zeros of order 5 at
`-j^2` for `j <= N`), `Delta_K = det G_K`, a Vandermonde identity for the top
coefficient (`[X^h] Delta_K != 0`, degree exactly `h`, its (2.9)), a
normalization `S_K` and a rational factor `m_(K,M)` from p-adic local
estimates giving `Q_(K,M) = m_(K,M) S_K Delta_K in Z[X]`, and the claim
(Theorem 2.1) `0 < Q_(40n,200)(zeta(5)) < exp(-139 n^2/5)`, whence
irrationality (`b^(37n) Q_n(a/b)` would be a nonzero integer below `exp(37n
log b - 139 n^2/5)`), and `|zeta(5) - a/b| > b^(-260)` for large `b`
(Corollary 1.2). Thresholds are not specified; `K >= 200 M^2` with `M = 200`
means `n >= 2·10^5` and determinants of size `7.4·10^6`, so nothing numerical
can be checked at the level of the theorem.

What can be checked is the engine (script part 7, thirty digits). The pole
values at `X = zeta(5)` are Hermite's integral for the Hurwitz zeta function:
`j^4 (zeta(5) - H_j) - 1/4 + 1/(2j) = 2 j^4 int_0^oo sin(5 arctan(u/j)) /
((j^2+u^2)^(5/2) (e^(2 pi u) - 1)) du` (values `0.2869, 0.0908, 0.0162,
0.0041` at `j = 1, 2, 5, 10`, decaying like `1/j^2`), and the polynomial values
are the Binet integral with a fourth-derivative factor: `(-1)^e B_(2e+2)
(2e+3)(2e+4)(2e+5)/24 = int_0^oo u^(2e+1) (2e+2)(2e+3)(2e+4)(2e+5) / (12
(e^(2 pi u) - 1)) du`. Since `sin(5 arctan(u/j)) / (j^2+u^2)^(5/2) = Im (j -
iu)^(-5) = (1/24) Im (d/du)^4 (j - iu)^(-1)` and `Im (j - iu)^(-1) = u/(u^2 +
j^2)`, four integrations by parts move the derivatives onto the Planck factor,
whose fourth derivative is positive (complete monotonicity of `1/(e^x - 1)`;
values `768, 0.754, 0.0026` at `x = 0.5, 2, 6`): both entry types are, up to
explicit polynomial corrections and boundary terms, moments of the positive
weight `(1/12) u^5 (d/du)^4 [1/(e^(2 pi u) - 1)] du` in the variable `t =
u^2`. That is the paper's "positive moment representation", and it is
correct in substance: positivity of a Hankel matrix of moments gives `Delta_K
(zeta(5)) > 0` (nonvanishing for free), and Heine's formula `det = (1/h!) int
prod_(i<j) (x_i - x_j)^2 prod w(x_i) dx` with a logarithmic-energy bound gives
decay on the `K^2` scale. The claim lives entirely in the p-adic normalization
`m_(K,M)` (its Sections 3–5: bases with prescribed vanishing at the pole
classes for intermediate primes, residue-valuation bases and a rank bound for
large primes, a coarse bound for small primes) and in the comparison of its
cost with the real decay. This note does not assess those sections; the
classification is AUTHOR-CLAIMED / UNREFEREED, as for the repo's earlier
p-adic zeta draft.

Three design points to record. (1) The `X`-dependence is confined to the
residues at the poles, and each pole value is a scaled tail of `zeta(5)`: the
"linear form in `1, zeta(5)`" lives in every entry, and the determinant turns
`h` such entries into a polynomial of degree `h`. (2) Nonvanishing is bought
by positivity, not by arithmetic: the well-poised symmetry of Zudilin's
forms and Nesterenko's criterion are replaced by a moment matrix. (3)
Smallness comes from two sources at once: the entries themselves are small
(the Euler–Maclaurin remainder `~ 1/(4j^2)` of the tails) and the
determinant compounds them (Vandermonde, log-energy).

## 5. Themes, deconstructed against the Collatz thread (exact statements)

| theme | in the zeta(5) constructions | in the Collatz thread | status |
|---|---|---|---|
| two rates compete | decay `C_0'` (saddle point) against denominators `C_2' = 9 - saving` (Suman); real decay `139n^2/5` against normalization cost (Fauzan) | 2-adic smallness `d_L` against height `L log_2 3 + log_2(n C_L)`: descent iff `d_L/L > log_2 3` (first note, Prop. 4) | exact; margins at the first coefficient descent of 27, 703, 6171, 77031, 837799: `0.0096, 0.0425, 0.0261, 0.0222, 0.0059` bits per odd step (script part 4) against the generic drift `2 - log_2 3 = 0.415` |
| nonvanishing | the hard half (Nesterenko; well-poised symmetry; positivity) | free: `A_L = 3^L n + S_(L-1) > 0` | inverted: Collatz's whole difficulty is smallness of a form that is trivially nonzero |
| partial sums and tails | `H_j` against `zeta(5)`, remainder `j^4(zeta(5) - H_j) - 1/4 + 1/(2j) ~ 0.41/j^2` | `r_L = n + S_(L-1)/3^L` against `xi = n C_inf`, remainder `eta_L 2^(d_L)/3^L` with `eta_L >= 1` the `3x-1` copy (S14) | same shape; the Collatz remainder does not shrink (HYP-9164: unbounded) |
| determinants amplify | entries small, `h x h` determinant `exp(-c K^2)` | `det(A_(i+j))` divisible by `2^(min_sigma sum d_(i+sigma(i)))`, quadratic in `N` | trivial in Collatz (entrywise divisibility carries no new information) |
| positivity / moments | Hankel of a positive weight (Planck factor, Hermite) | a Collatz Hankel matrix is a moment matrix iff the valuations are non-decreasing: never on a divergent word (section 2) | no positive structure to exploit |
| arithmetic saving from floor sums | `v_(k,p)` = Legendre floor sums, sawtooth `phi`, `int phi dpsi` | the only floor sums are the Beatty boundary words of `E_inf` (`d_k = ceil(k log_2 3)`, the parallel note's tight blocks); `3^L` and `S_(L-1)` are coprime, no common factor to save | no analogue; the cycle equation `n(2^A - 3^p) = S` is the one place a "saving" (a divisor of `S`) appears, e.g. `139 | 2363` for `-17` |
| rationality from a Hankel rank | Kronecker: rational iff finite rank; Fauzan aims at nonvanishing | PC(n) iff finite rank; depth = rank - 1 (section 2) | exact reformulation, direction inverted |
| depth as a grading | — | `d(2^x m) = x + d(m)`: the Euler factor at 2 of the depth Dirichlet series `Z(s,z) = sum z^(d(n)) n^(-s) = (1 - z 2^(-s))^(-1) sum_(m odd) z^(d(m)) m^(-s)`; `Z(s, 1) = zeta(s)` iff Collatz | exact, one line; the odd part has no closed functional equation (the `(n-1)/3` predecessors twist by a residue class, as in Berg–Meinardus) |

The one-paragraph reading: every mechanism in the two preprints is a way of
making a nonzero quantity small, and the Collatz thread has been the study
of a quantity (`A_L`, or the shadow error) that is nonzero for free and
refuses to be small. The shapes match exactly and the difficulty sits in the
complementary half each time.

## 6. Numerology, tested and dismissed (as the directive asked)

| coincidence | test | verdict |
|---|---|---|
| `3^7 - 2^11 = 139`; the Zenodo decay constant is `139/5`; `-17 = S/(2^11 - 3^7)`, `S = 17·139` | the Zenodo `139/5` is `C_real - C_arith` for `K = 40n`, `M = 200` (its 2.7), a difference of two energy/normalization constants; no `2^11` or `3^7` in its equations | coincidence |
| `37` lies on the `-17` cycle (`-37`) and on the `3x-1` cycle (`37`); the Zenodo degree is `37n = K - N` with `K = 40n`, `N = 3n` | `37 = 40 - 3` is a parameter choice (`alpha = 3/40`) | coincidence |
| Suman's `q = 11`, seven successive maxima `m_1..m_7` (plus `m_8`), `A = 8`; the `-17` cycle has `A = 11`, `p = 7`; `2^11 < 3^7`, `2^3 < 3^2` (the `-5` cycle) | `q = 11` is the number of hypergeometric parameters `h_1..h_11`, `A = 8` the pole order | coincidence |
| `zeta(5), zeta(7), zeta(9)` and the `3x-1` cycle `{5, 7}`, the `-5` cycle | — | coincidence |
| the negative cycles' `(A, p) = (1,1), (3,2), (11,7)` against the convergents of `log_2 3` (`1/1, 2/1, 3/2, 8/5, 19/12, 65/41, 84/53`) | `11/7` is the mediant of `3/2` and `8/5`, a semiconvergent; `1/1` and `3/2` are convergents | exact and classical (Steiner's route), not numerology; recorded in the first note as shape (B) |
| Suman's margin `C_0' - C_2' = 0.4766` against the generic Collatz drift `2 - log_2 3 = 0.415` | different units (nats per parameter `n` against bits per odd step); the Suman margin is unsupported anyway | coincidence |

## 7. Directions (DIRECTION; none proved)

* **D5. The layer growth constant.** Prove `|L_k|^(1/k) -> lambda_C = (1 +
  sqrt(7/3))/2` for the Collatz-step layers, or bracket it as
  Applegate–Lagarias do for the Terras-step trees (their `[1.302, 1.36]`
  becomes, re-timed, a bracket around `1.264`). The exact `1.2639` at depth
  60 suggests convergence from above with a `k^(-c)` correction; the
  heuristic assumes the `4 mod 6` residue is hit with density `1/3` along the
  tree, which the 3-adic structure of the tree does not guarantee.
* **D6. The Hankel rank of a word prefix as a certificate.** For a finite word
  prefix the linear complexity of `(2^(d_k))_(k < K)` is at most `K/2` and
  says nothing; the law `LC = d_odd + 1` is a property of the tail. A
  certificate would have to bound the rank of the *infinite* word, which is
  the conjecture. Not a method; recorded to close the door.
* **D7. Fauzan's design point (1) for Collatz.** Is there a natural matrix
  attached to an orbit whose entries are affine in `n` with a full-rank
  `n`-part? The Hankel matrix of the Apéry forms has a rank-one `n`-part
  (`3^i 3^j`), so its determinant is affine in `n` on each cylinder and
  carries nothing beyond the cylinder congruence. A full-rank version would
  need the orbit to enter through more than one linear form, i.e. through
  more than one place; the two-place descent tree of the parallel note
  (2-adic source, 3-adic landing) is the only such object in the thread.
* **D8. The depth Dirichlet series.** `Z(s, z)` has the Euler factor at 2
  exactly and nothing else closed; whether its odd part has a
  Berg–Meinardus-type functional equation with cube roots of unity (first
  note, D1) is the same question in Dirichlet-series clothing.

## 8. Audit status

Independent audit: OWED (the session's audit subagent was terminated by an
API rate limit before reporting; the labels above are the producer's own).
Producer self-checks, separate from the script: the layer table against the
owner's table (exact match); `d(2n) = d(n) + 1` on the computed component;
the Berlekamp–Massey implementation against random prefixes (`LC = t`); the
two Zenodo identities at thirty digits (mpmath, independent of the paper's
text); Suman's `tau_0` by two root-finders (numpy in the scratchpad,
Durand–Kerner in the script) agreeing to eight digits.
