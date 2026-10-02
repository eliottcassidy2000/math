# Edge multiset dimension of hypercubes, part 2: ln edim_m(Q_d) = Theta(d^(1/3)), sharp constants for the random method, Fourier-Hoelder forests

**Lane:** procgen_edim2, 2026-10-01 (worktree collatz-procgen-20260922, machine mac-mini). Resumption
after a reboot; started fresh.
**Source problem:** J. Allikvere, *The edge multiset dimension of hypercubes*, arXiv:2608.09983v1,
Section 10, Open Problems 2, 3, 4 (Open Problem 1 is THM-4525: edim_m(Q_6) = 15).
**Builds on:** `procgen_edim_20261001_edge_multiset_dimension.md` (lemmas L1-L5, explicit sets).
**Runner:** `04-computation/experiments/procgen_edim2_20261001_run.py` (stdout in
`05-knowledge/results/procgen_edim2_20261001.out`, a `--quick` run; see Section 7).

## Status header

| # | Claim | Label |
|---|-------|-------|
| 1 | **Theorem A (OP2, upper bound).** If `ln(q 2^d) >= (C* + eps) d^(1/3)` with `C* = (3 sqrt2 ln2)^(2/3) = 2.05262`, a random set of density `q` is edge-multiset resolving with probability `-> 1`. Hence `edim_m(Q_d) <= exp((C* + o(1)) d^(1/3))`. | PROVED |
| 2 | **Theorem B (OP2, lower bound).** `edim_m(Q_d) >= exp((c* - o(1)) d^(1/3))`, `c* = (3 ln2/(2 sqrt2))^(2/3) = 0.81458` (sharper analysis of the entropy inequality L5 of THM-4525, which was analysed as 0.6215). | PROVED |
| 3 | **Corollary (OP2).** `ln edim_m(Q_d) = Theta(d^(1/3))`; `0.8145 <= liminf ln edim / d^(1/3) <= limsup <= 2.0527`, the two constants differ by exactly `4^(2/3)`. Whether the limit exists: OPEN. | PROVED / OPEN |
| 4 | **Theorem C.** `C*` is exactly the reach of the paper's method: for every density and every choice of forests the forest-lemma union bound is `>= 1` once `ln(q 2^d) <= (C* - eps) d^(1/3)`. | PROVED |
| 5 | **Lemma A / Lemma FH.** A Fourier atom bound for `Bin(n,q) - Bin(n',q)` (Gaussian constant), and a Fourier-Hoelder multi-forest lemma that multiplies forest bounds over edge-disjoint forests. | PROVED |
| 5b | **Proposition K.** For parallel edge pairs the weighted spanning-tree count of the level graph is `2^n prod_k (C(n,k) - K_k(h))` (Krawtchouk diagonalisation; connected iff `h` odd). | PROVED + VERIFIED (`d <= 11`) |
| 6 | **Theorem D (OP3, partial).** A closed-form estimate `U_d <= d^2 2^(2d-3) e^(-Phi) + d 2^(d-2) e^(-Psi) < 1` for every `d >= 17` (no forest optimisation, no per-`d` table; three interval evaluations at `d = 17, 18, 19`). | PROVED |
| 7 | **Certificates.** Explicit resolving sets for every `6 <= d <= 16`; new: `edim_m(Q_13) <= 105`, `Q_14 <= 125`, `Q_15 <= 135`, `Q_16 <= 171`. With Theorem D: finiteness for all `d >= 6` without the union-bound computation of the paper for `11 <= d <= 50`. | VERIFIED |
| 8 | OP3 in its strict form (one analytic estimate valid from `d = 11`). | OPEN (Section 4.3: the best forest union bound has margin about 6 at `d = 11`) |
| 9 | With Lemma FH the probabilistic method works from `d = 10`: certified `U_10 <= 0.2649` at density 1/2 (the paper's bound gives `1.307`). | FINITE-EXACT (interval arithmetic) |
| 10 | **OP4.** Density 1/2 overshoots the optimum by `exp(d ln2 - O(d^(1/3)))`; sparse random sets are optimal up to the constant `4^(2/3)` in the exponent. Certified sparse bounds `edim_m(Q_d) <= M_d` for `11 <= d <= 64` (e.g. `M_11 = 361`, `M_16 = 636`, `M_32 = 3613`, `M_64 = 31808`; `ln M_d / d^(1/3) = 2.56-2.65`). | PROVED (asymptotics) / FINITE-EXACT (table) |
| 11 | EMPIRICAL: true first moments of density-1/2 sets: about 3.3, 0.8, 0.04 colliding pairs for `d = 8, 9, 10`; random `k`-sets become resolving around `k = 180` at `d = 11, 12`. | EMPIRICAL |
| 12 | **`Q_7`:** no resolving set of size `<= 10` (exhaustive, three complete normal forms, validated against THM-4525 on `Q_6` and Burnside counts), so `11 <= edim_m(Q_7) <= 19`. `k = 11`: all cases but `a = 11` done, none resolving. | FINITE-EXACT (`k <= 10`); OPEN (`k = 11`, one case left) |

## 0. Notation and cited facts

* `Q_d`, `n := d - 1`. For an edge `e` and a vertex `w`, `d(e,w) = min` over the endpoints (paper Sec. 1).
  `H_e(r) = #{s in S : d(e,s) = r}`, `r = 0..n`.
* `b_M(j) := C(M,j)/2^M`, `beta(N) := b_N(floor(N/2))` (paper's central atom).
* **Random model.** `S_q`: every vertex independently with probability `q in (0, 1/2]`.
  `m := q 2^d` (expected size), `lambda := ln m`.
* **Cells (CITED, paper Sec. 6-7).** For distinct edges `e, f`: `V_ab = {w : d(e,w) = a, d(f,w) = b}`,
  `n_ab = |V_ab|`, `N_ab = n_ab + n_ba` (`a != b`); the level graph `Gamma_{e,f}` has an edge `{a,b}` iff
  `N_ab > 0`.
* **Pair types (CITED, Lemma 13 and Prop. 14).** Every unordered pair of distinct edges is
  `Aut(Q_d)`-equivalent to exactly one representative of type `par(h)` (`1 <= h <= n`) or `crs(h)`
  (`0 <= h <= n-1`), with `#par(h) = d 2^(d-2) C(n,h)`, `#crs(h) = C(d,2) 2^d C(n-1,h)`; the total number of
  pairs is `< d^2 2^(2d-3)`. Restating the proof of Prop. 14 probabilistically: let `w` be a uniform vertex
  and `(A, B) = (d(e,w), d(f,w))`, so that `n_ab = 2^d P(A = a, B = b)`. Then
  * `par(h)`, `g := n - h`: `(A, B) = (P + R, h - P + R)`;
  * `crs(h)`, `g := n - 1 - h`: `(A, B) = (eta + P + R, eps + h - P + R)`;

  with `P ~ Bin(h, 1/2)`, `R ~ Bin(g, 1/2)`, `eps, eta ~ Bin(1, 1/2)` independent. In both cases
  `A ~ Bin(n, 1/2)` and `B ~ Bin(n, 1/2)`. (Runner: cells re-derived from these formulas agree with brute
  force for `d <= 7`.)
* **Forest lemma (CITED, Lemma 11; valid verbatim for every q).** For every forest `F` in `Gamma_{e,f}`,
  `P(H_e = H_f) <= prod_{ab in F} A_q(n_ab, n_ba)`, where `A_q(n, n') := max_t P(Bin(n,q) - Bin(n',q) = t)`.
  The proof uses only that the `X_ab = |S cap V_ab| ~ Bin(n_ab, q)` are independent; it bounds
  `P(H_e - H_f = v)` for every fixed vector `v`. (Runner: brute force over all `2^16` subsets of `Q_4`,
  `q in {1/2, 1/4, 1/10}`, every pair type.)

## 1. Two atom lemmas and binomial estimates (PROVED)

**Lemma A (Fourier atom bound).** Let `U ~ Bin(n, q)`, `V ~ Bin(n', q)` be independent, `N = n + n'`,
`x := N q (1-q) > 0`. Then

    A_q(n, n') <= (1 + 1/(4x)) / sqrt(2 pi x) + e^(-x)/2.

In particular `A_q(n, n') <= 0.69 / sqrt(x) <= x^(-1/2)` when `x >= 1`.

*Proof.* `|E e^(i theta (U-V))| = |1 - q + q e^(i theta)|^N = (1 - 4x/N sin^2(theta/2))^(N/2) <= exp(-2x sin^2(theta/2))`.
Fourier inversion gives `P(U - V = t) <= (1/2pi) int_{-pi}^{pi} exp(-2x sin^2(theta/2)) d theta
= (2/pi) int_0^1 e^(-2x s^2) (1 - s^2)^(-1/2) ds` (substituting `s = sin(theta/2)`). On `[0, 2^(-1/2)]` use
`(1 - s^2)^(-1/2) <= 1 + s^2` (convexity in `s^2`, checked at `s^2 = 1/2`) and extend to `[0, inf)`:
`int_0^inf e^(-2xs^2)(1 + s^2) ds = sqrt(pi/(8x)) (1 + 1/(4x))`. On `[2^(-1/2), 1]` use `e^(-2xs^2) <= e^(-x)`
and `int (1-s^2)^(-1/2) = pi/4`. For the last claim, `sqrt(x)` times the bound is
`(1 + 1/(4x))/sqrt(2 pi) + sqrt(x) e^(-x)/2`; both terms decrease on `[1, inf)`, and at `x = 1` the sum is
`0.4987 + 0.1840 < 0.69`. QED

At `q = 1/2` the bound is within `1 + O(1/N)` of the exact `beta(N)`; for small `q` it has the Gaussian
constant `1/sqrt(2 pi)`.

**Lemma A2 (atoms cannot be smaller).** For every integer-valued `Y`, `max_t P(Y = t) >= (12 Var Y + 1)^(-1/2)`.
Hence `A_q(n, n') >= (12 N q (1-q) + 1)^(-1/2)`.

*Proof.* With `W` uniform on `[-1/2, 1/2]` independent of `Y`, `Y + W` has a density bounded by
`M := max_t P(Y = t)` and variance `Var Y + 1/12`. A density bounded by `M` has variance at least that of
the uniform density on an interval of length `1/M` (bathtub principle), i.e. `1/(12 M^2)`. QED

**Lemma B (binomial estimates).** Let `M >= 1`, `k0 = floor(M/2)`, `0 <= t <= k0`, and let `n >= 1`.
1. `b_M(k0 +- t) >= sqrt(2 / (pi (M + 2))) * exp(-t(t+1) / (k0 + 1 - t))`.
2. For every integer `r`, with `x = r - n/2`: `b_n(r) <= sqrt(2 / (pi (n + 1/2))) * exp(-(x^2 - 1/4) / (n/2 + |x|))`.

*Proof.* The sequences `C(2k,k) 4^(-k) sqrt(k + 1/4)` and `C(2k,k) 4^(-k) sqrt(k + 1/2)` are respectively
increasing and decreasing (their squared ratios are `1 + 1/(16k^3 - 12k^2)` and `1 - 1/(4k^2)`), with limit
`pi^(-1/2)`; hence `(pi (k + 1/2))^(-1/2) <= C(2k,k)/4^k <= (pi (k + 1/4))^(-1/2)`. Since
`b_{2k+1}(k) = b_{2k+2}(k+1)`, this gives `b_M(k0) >= sqrt(2/(pi(M+2)))` and `beta(n) <= sqrt(2/(pi(n+1/2)))`.
For the decay, `b_M(k0 + t)/b_M(k0) = prod_{i=1}^t (ceil(M/2) - i + 1)/(k0 + i)`; each factor is at least
`1 - (2i-1)/(k0+i) >= exp(-(2i-1)/(k0 - i + 1))` (by `ln(1-y) >= -y/(1-y)`), and the exponents sum to at most
`t^2/(k0 + 1 - t)`. Downwards, `b_M(k0 - t)/b_M(k0) >= prod_i (1 - 2i/(k0+1+i)) >= exp(-t(t+1)/(k0+1-t))`.
For (2), the factors of `b_n(r)/beta(n)` are `<= exp(-(2j-1)/(n/2 + |x|))` (n even) or
`exp(-2j/(n/2 + |x|))` (n odd), by `1 - y <= e^(-y)`; they sum to `x^2` resp. `x^2 - 1/4`. QED

(All three statements are re-checked numerically by the runner for `M, n < 400`.)

## 2. Open Problem 2: the upper bound (Theorem A)

### 2.1 The star lemma
For a pair type `tau` put `eta_tau := h/n` for `par(h)` and `eta_tau := (h + 1/2)/n` for `crs(h)`, and
`c := n/2` (the centre of the levels). Write `x := a - c` for a level `a`.

**Lemma S.** Let `|x| >= x0(tau) := sqrt(n/2 + 1/4) / min(eta_tau, 1 - eta_tau)`. Then some level `b` with
`|b - c| < |x|` satisfies `n_ab >= (3 / (8|x|)) 2^d b_n(a)`.

*Proof.* Condition on `A = a`. The `n` fair bits behind `A` (`P`'s `h` bits and `R`'s `g` bits, plus `eta` for
`crs`) have a uniformly random set of `a` ones, so `P` is hypergeometric.
* `par(h)`: `B - A = h - 2P`, `E[P | A=a] = ah/n`, so `E[B - A | A=a] = -2x h/n`, and
  `Var(B - A | A=a) = 4 a(n-a) h(n-h) / (n^2 (n-1)) <= n^2 / (4(n-1)) <= n/2` (for `n >= 2`; for `n = 1` it is 0).
* `crs(h)`: `B - A = eps - eta + h - 2P` with `eps` independent; `E[B - A | A=a] = 1/2 - a/n + h - 2ah/n
  = -2x (h + 1/2)/n`. The term `eta + 2P` is a weighted count of a sample without replacement with weights in
  `[0, 2]`, so its variance is at most `a(n-a)/(n-1) <= n/2`; hence `Var(B - A | A=a) <= n/2 + 1/4`.

So `B - A` has conditional mean `-2 x eta_tau` and variance `<= n/2 + 1/4`. Take `x > 0` (the other sign is
symmetric). `|B - c| < x` is the event `B - A in (-2x, 0)`, whose endpoints are at distance
`2x eta_tau` and `2x (1 - eta_tau)` from the mean. By Chebyshev,
`P(|B - c| >= x | A = a) <= (n/2 + 1/4) / (4 x^2 min(eta,1-eta)^2) <= 1/4`.
So `P(A = a, |B - c| < x) >= (3/4) b_n(a)`. At most `2|x| - 1` integers `b` satisfy `|b - c| < |x|`; one of
them carries at least `3 b_n(a) / (8|x|)` of the probability, and `n_ab = 2^d P(A=a, B=b)`. QED

**Star forests.** Let `phi(a)` be a level given by Lemma S, for every `a` with `|a - c| >= x0(tau)`. The edges
`{a, phi(a)}` form a forest: orienting each from `a` to `phi(a)`, every vertex has out-degree at most one and
`|. - c|` strictly decreases along arcs, so there is no cycle (and no edge is chosen twice).

(Runner: the conditional mean and variance identities, and Lemma S itself, are checked exactly for every type,
level and `4 <= d <= 40`.)

### 2.2 The path and zigzag forests (all types except `par(n)`)
For the extreme types we use two simpler forests (as in the paper's Lemma 18, but with all indices):
* **r-path** (if `g >= h`). `par(h)`: cells `(p0, r)`, `p0 = floor((h-1)/2)`, `r = 0..g`, level pairs
  `{p0 + r, h - p0 + r}` with difference `h - 2p0 in {1, 2}`. `crs(h)`: cells `(p0, r)` with `(eps, eta) = (1, 0)`,
  `p0 = floor(h/2)`, level pairs `{p0 + r, 1 + h - p0 + r}`. Fixed nonzero difference => union of paths.
  Here `M := g`.
* **zigzag** (if `h > g`). Two matchings at consecutive `r`: `par(h)`: `r in {r0, r0+1}`, `r0 = floor((g-1)/2)`
  (needs `g >= 1`); `crs(h)`: the cells `(eps, eta) = (0,0)` and `(1,1)` at `r0 = floor(g/2)`. The two
  matchings are `a <-> s0 - a` and `a <-> s0 + 2 - a`; their composition is the shift `a -> a + 2`, so the
  union has no cycle (a union of two matchings has only paths and even cycles, and a cycle would be a finite
  orbit of a nontrivial translation). Here `M := h`.

In every case, for each `1 <= t <= floor(M/2)` the forest has two edges (indices `j = floor(M/2) +- t`) with
`N_j >= 2^d kappa_n b_M(j)`, `kappa_n := (1/4) sqrt(2/(pi(n+1)))`. (The central binomial factor of the
other index is bounded by Lemma B.1; runner: checked for every forest cell with `4 <= d <= 40`.)
Note `M >= (n-1)/2` always, and `M >= (1 - delta) n - 3/2` when `min(eta_tau, 1 - eta_tau) < delta`.

**`par(n)`** (`g = 0`): the level graph is the matching `{p, n - p}`, `N = 4 C(n,p) = 2^(d+1) b_n(p)`.

### 2.3 Exponent estimates
For a forest edge with `y := N q (1-q) >= 1`, Lemma A gives `-ln A_q >= (1/2) ln y`, and `y >= N q / 2`.

**(E1) star forests.** For `x0 <= |x| <= X := min(sqrt(n lambda'/2), 2 n^(2/3))`, Lemma S and Lemma B.1 give
`ln y >= lambda' - 2x^2/n`, where `lambda' := lambda - (7/6) ln n - K` for an absolute constant `K` and all
large `n`. (Use `y >= (3/(16|x|)) m b_n(a)`, `ln|x| <= (2/3) ln n + ln 2`, and, for `t = |a - floor(n/2)|
<= |x| + 1/2 <= 2 n^(2/3) + 1/2`,
`t(t+1)/(floor(n/2)+1-t) <= (x^2 + 2|x| + 3/4)/(n/2 - |x|) <= 2x^2/n + 8|x|^3/n^2 + 2(4|x| + 3/2)/n <= 2x^2/n + 65`
(using `1/(1-u) <= 1 + 2u` for `u = 2|x|/n <= 1/2`, and `n >= 2^12`).)
Summing the concave function `lambda' - 2x^2/n >= 0` over the levels in `x0 <= |x| <= X`:

    -ln P_tau >= (1/2) [ (4/3) lambda' X - 2 lambda' - (2 x0 + 1) lambda' ].

When `X = sqrt(n lambda'/2)` this is `(sqrt2/3) lambda'^(3/2) sqrt(n) - lambda' (x0 + 3/2)`.

**(E2) path/zigzag forests.** With `Lambda_M := lambda + ln(kappa_n/2) - (1/2) ln(pi(M+2)/2) = lambda - O(ln n)`,

    -ln P_tau >= sum_{t=1}^{floor(M/2)} ( Lambda_M - t(t+1)/(floor(M/2) + 1 - t) )_+ ,

and for `k := floor(M/2)`, `1 << Lambda << k`: taking `t <= T := floor(sqrt(k Lambda)) - 1` gives
`sum >= T Lambda - (T+2)^3 / (3(k - T)) = (2/3) Lambda^(3/2) sqrt(k) (1 - o(1)) = (sqrt2/3) Lambda^(3/2) sqrt(M) (1 - o(1))`.

**(E3) `par(n)`.** `y >= m b_n(p)`; only one side of the centre is available, so
`-ln P >= (1/2) sum_{t >= 1} (lambda - (1/2) ln(pi(n+2)/2) - t(t+1)/(floor(n/2)+1-t))_+
= (1/2)(sqrt2/3) lambda^(3/2) sqrt(n) (1 - o(1))`.

### 2.4 Theorem A
**Theorem A.** Let `eps > 0` and `C* := (3 sqrt2 ln2)^(2/3) = 2.052616...`. If `q = q(d) in (0, 1/2]`
satisfies `lambda = ln(q 2^d) >= (1 + eps)^(2/3) C* n^(1/3)`, then
`P(S_q is edge-multiset resolving and |S_q| <= 2 q 2^d) -> 1`. Hence
`edim_m(Q_d) <= exp((C* + o(1)) d^(1/3))`.

*Proof.* Fix `delta = 1/20`. Note `(sqrt2/3) (C*)^(3/2) = 2 ln 2`.
* **Bulk types** (`delta <= eta_tau <= 1 - delta`): `x0 <= sqrt(n)/delta = o(sqrt(n lambda'))`. If
  `lambda' <= 8 n^(1/3)` then `X = sqrt(n lambda'/2)` and (E1) gives
  `-ln P_tau >= (sqrt2/3) lambda^(3/2) sqrt(n) (1 - o(1)) >= 2(1 + eps)(1 - o(1)) n ln 2`.
  If `lambda' > 8 n^(1/3)`, then `X = 2n^(2/3)`, every summand is `>= lambda'/2`, and
  `-ln P_tau >= (1/2)(2X - 2x0 - 3)(lambda'/2) >= n^(2/3) lambda' (1 - o(1)) >= 8 n (1 - o(1))`.
  Either way the bulk contributes at most `d^2 2^(2d-3) 2^(-2(1 + eps/2) n) -> 0`.
* **Extreme types** (`eta_tau < delta` or `> 1 - delta`, excluding `par(n)`): by (E2) with
  `M >= (1-delta) n - 3/2`, `-ln P_tau >= 2(1+eps) sqrt(1-delta) (1 - o(1)) n ln 2 >= 1.94 n ln 2` (and more
  for larger `lambda`, since the bound increases with `lambda`). Their number is at most
  `d^2 2^(d+3) 2^(n H_2(2 delta)) <= d^2 2^(d + 3 + 0.47 n)`. Contribution `-> 0`.
* **`par(n)`**: by (E3), `-ln P >= (1 + eps)(1 - o(1)) n ln 2`, and there are `d 2^(d-2) = d 2^(n-1)` such
  pairs. Contribution `-> 0`.
* **Size**: `P(|S_q| > 2m) <= (e/4)^m -> 0`.

So with probability `1 - o(1)` all `d 2^(d-1)` histograms are distinct (which forces `S_q != empty`) and
`|S_q| <= 2m`. Letting `eps -> 0` gives the corollary. QED

At `q = 1/2` this reproves the paper's Theorem 19 for all large `d` and shows `P(resolving) -> 1` for density
`1/2`; the paper's own bounds give `U_d < B_d = C^5 d^12 2^(-3d-3) -> 0` from `d = 51` on.

### 2.5 Theorem C: C* is the exact reach of the forest-lemma union bound
**Theorem C.** Fix `eps > 0`. For `d` large, every `q in (0,1/2]` with `ln(q 2^d) <= (C* - eps) d^(1/3)` and
every choice of forests `F_tau` give `sum_tau #tau prod_{F_tau} A_q >= 1`. Consequently the union bound with
the forest lemma cannot prove anything better than `edim_m(Q_d) <= exp((C* - o(1)) d^(1/3))`.

*Proof.* Root every tree of `F_tau`; each edge `{a, parent(a)}` is charged to its child `a`, distinct edges to
distinct levels. Since `n_ab <= #{w : d(e,w) = a}` and `n_ba <= #{w : d(f,w) = a}`, Prop. 2 gives
`N_e <= 4 C(n,a) = 2^(d+1) b_n(a)`. By Lemma A2, `-ln prod_{F_tau} A_q <= (1/2) sum_a ln(1 + 24 m b_n(a))`.
With Lemma B.2 and `ln(1 + e^u) <= u_+ + e^(-|u|)`, the right side is
`(1/2)(4/3) lambda^(3/2) sqrt(n/2) (1 + o(1)) = (sqrt2/3) lambda^(3/2) sqrt(n) (1 + o(1))`
(the levels with `|u| <= 1` number `O(sqrt(n/lambda))`). This holds for every type, so the union bound is at
least `(#pairs) exp(-(sqrt2/3) lambda^(3/2) sqrt n (1 + o(1))) >= d^2 2^(2d-4) 2^(-2(1 - eps')n(1 + o(1))) -> infinity`,
where `eps' > 0` depends on `eps` (the exponent bound is increasing in `lambda`, so it suffices to take
`lambda = (C* - eps) n^(1/3)`). Any upper bound used in place of `A_q` only increases the union bound. QED

Heuristically the true collision probability of a bulk pair is also `exp(-(sqrt2/3) lambda^(3/2) sqrt(n)(1+o(1)))`,
so `C*` is presumably the true threshold of *uniformly random* landmark sets (heuristic, not proved). Indeed,
for a parallel pair the covariance of the level-divergence vector is diagonalised by the Krawtchouk transform:

**Proposition K (PROVED; runner: exact check for `3 <= d <= 11`).** For `par(h)` the weighted level graph
(weights `N_ab`) has weighted spanning-tree count `kappa = 2^n prod_{k=1}^n D_k(h)`, where
`D_k(h) = C(n,k) - K_k(h) = 2 #{xi in {0,1}^n : |xi| = k, |xi AND v| odd}` (`K_k` the Krawtchouk polynomial,
`|v| = h`). In particular the level graph is connected iff `h` is odd.
*Proof.* `sum_{a<b} N_ab (theta_a - theta_b)^2 = 2 sum_{x in Q_n} (f(x) - f(x + v))^2` with `f(x) = theta(|x|)`;
by Parseval this equals `2^(2-n) (K theta)^T diag(D_0, ..., D_n) (K theta)`, `K = (K_j(k))_{k,j}`, `D_0 = 0`. So the
Laplacian is `2^(2-n) K^T D K`; deleting level 0 and using `K^2 = 2^n I` and Jacobi's complementary-minor
formula, `det K_{00-minor} = 2^(-n) det K`, `(det K)^2 = 2^(n(n+1))`, which gives `det L_red = 2^n prod D_k`.
The matrix-tree theorem finishes. QED

The Gaussian local-limit heuristic then gives collision probability about `(2/pi)^(n/2) kappa^(-1/2)` at `q = 1/2`
(and the analogue with `q`-scaled weights in general), whose logarithm has the same leading term as the forest
bound. (Rigorous Gaussian domination on the torus fails for small cells, which is why the finite computations
use forests.)

## 3. Open Problem 2: the lower bound (Theorem B)

**Theorem B.** `edim_m(Q_d) >= exp((c* - o(1)) d^(1/3))` with `c* = (3 ln2/(2 sqrt2))^(2/3) = 0.814581...`.

*Proof.* L5 (THM-4525, CITED) states `d - 1 + log2 d <= sum_{r=0}^{n-1} g(mu_r)`, `mu_r = m b_n(r)`,
`g(mu) = (mu+1) log2(mu+1) - mu log2 mu`, for any resolving set of size `m`. Write `lambda = ln m`.
* `g(mu) <= log2 mu + 1 + log2 e` for `mu >= 1`, and `g(mu) <= mu log2(2e/mu) <= 3 sqrt(mu)` for `mu < 1`.
* By Lemma B.2, `ln mu_r <= L' - (x^2 - 1/4)/(n/2 + |x|)` with `L' := lambda - (1/2) ln(pi(n + 1/2)/2)`.
  The levels with `mu_r >= 1` have `|x| <= X2 = sqrt(n L'/2)(1 + o(1))`, and
  `sum_{mu_r >= 1} ln mu_r <= (4/3) L'^(3/2) sqrt(n/2) (1 + o(1))`.
* The remaining terms are `O(sqrt(n lambda))` (at most `2 X2 + 1` levels with `mu_r >= 1`, each with the
  additive `2.45` bits; and `sum_{mu_r < 1} sqrt(mu_r) = O(sqrt(n lambda))` by Lemma B.2).

Hence `n ln 2 <= (4/3) L'^(3/2) sqrt(n/2) (1 + o(1)) + O(sqrt(n lambda))`, i.e.
`L'^(3/2) >= (3 ln2 / (2 sqrt2)) sqrt(n) (1 - o(1))`, and `lambda = L' + O(ln n)`. QED

The previous lane's analysis of the same inequality (one `log2 m` bound per heavy level) gave `0.6215`; the
inequality itself is unchanged, and the exact values `4, 5, 5, 6, 8, 11, 15, 28, 81, 275, 1159, 50116`
(`d = 6, ..., 1024`) still stand. Convergence to `c*` is slow (`ln(50116)/1024^(1/3) = 1.07` at `d = 1024`,
because of the `ln n` terms).

## 4. Open Problem 3: a uniform estimate, and finiteness for all d >= 6 without the union-bound computation

The paper proves finiteness for `11 <= d <= 50` by evaluating, for each `d`, optimal (Kruskal) forests for all
`2d - 1` pair types and a rational certificate of the union bound `U_d` (Table 3), and for `d >= 51` by an
analytic ten-edge tail estimate. Open Problem 3 asks for one analytic estimate from `d = 11`.

**What is proved here (PARTIAL answer):**
* **Theorem D (closed-form estimate, PROVED).** For every `d >= 17`, the paper's union bound at density 1/2 is
  `U_d <= d^2 2^(2d-3) e^(-Phi(d-1)) + d 2^(d-2) e^(-Psi(d-1)) < 1`, with `Phi`, `Psi` explicit elementary
  functions (below). No forest optimisation and no per-`d` table are needed; the only numerical inputs are
  three interval evaluations of `Phi, Psi` at `d = 17, 18, 19`, exactly like the paper's evaluation of `B_51`.
* **Certificates (VERIFIED).** Explicit resolving sets for every `6 <= d <= 16`
  (`Q_6`: THM-4525; `Q_7..Q_12`: previous lane; `Q_13..Q_16`: this lane, sizes 105, 125, 135, 171, see Section 4.4),
  each checked by two independent exact verifiers.
* **Consequence.** `edim_m(Q_d) < infinity` for every `d >= 6` (and `= infinity` for `2 <= d <= 5`, CITED)
  is proved without any union-bound computation for `11 <= d <= 50`: certificates up to `d = 16`, Theorem D
  from `d = 17`.
* **OPEN:** a single analytic estimate valid from `d = 11` itself. Section 4.3 explains why this is delicate:
  `U_11 = 0.156` even with optimal forests and exact atoms (margin 6.4), and `U_10 > 1`.

### 4.1 The closed-form forests and the per-type bound
At `q = 1/2` use the forests of Section 2.2 with the following explicit cell lower bounds (all checked against
the exact cells for `4 <= d <= 40` by the runner):
* `par(h)`, path: `N_r = 4 C(h, p0) C(g, r)`; zigzag: `N = 4 C(h,p) C(g, r')`, `r' in {r0, r0 + 1}`;
* `crs(h)`, path: `N_r >= 2 C(h, floor(h/2)) C(g, r)` (the `(eps,eta) = (1,0)` cell and its mirror `(0,1)` cell);
  zigzag: `N >= 2 C(h,p) C(g, floor(g/2))` (the `(0,0)` resp. `(1,1)` cell and its mirror);
* `par(n)`: the matching, `N_p = 4 C(n,p)`.

With the paper's Lemma 11 and Lemma 17 (`beta(N) <= sqrt(2/(pi N))`, CITED) one gets, writing
`gamma(m) := ln C(m, floor(m/2))` and `ell(m) := ln prod_{r=0}^m C(m, r)`:
* path: `-ln P >= P(h, g) := (1/2)[(g+1)(ln pi + gamma(h)) + ell(g)]`,
* zigzag: `-ln P >= Z(h, g) := ceil(h/2) (ln pi + gamma(g)) + (1/2)(ell(h) - gamma(h))`,
* `par(n)`: `-ln P >= Psi(n) := (n/4) ln(2pi) + (1/4)(ell(n) - gamma(n))`.

For `par(h)` the bounds are at least `P(h, g-1)`, `Z(h, g-1)` (use `2C(h,p0) >= C(h, floor(h/2))`,
`C(g,r0) C(g,r0+1) >= C(g, floor(g/2))^2 / 2`, and monotonicity in `g`). Hence every type except `par(n)`
satisfies `-ln P_tau >= max(P(h, n-1-h), Z(h, n-1-h))` for some `0 <= h <= n-1`.

**Lemma PI.** For `m >= 1`, `ell(m) >= ell_-(m) := (m^2 - 1)/2 - ((m+1)/2) ln m`.

*Proof.* `ell(m) = sum_{k=1}^m (2k - m - 1) ln k`; pairing `k` with `m + 1 - k` gives
`ell(m) = sum f(j)` over `0 < j <= m - 1`, `j = m + 1 (mod 2)`, with `f(u) = u ln((M+u)/(M-u))`, `M = m + 1`.
`f` is convex and increasing, so the trapezoid rule (nodes spaced 2 apart) gives
`sum f(j) >= (1/2) int_{j0}^{m-1} f + (f(j0) + f(m-1))/2`, and `f(j0) >= int_0^{j0} f` (for `j0 = 2` this uses
`M >= 4`). With `F(u) = ((u^2 - M^2)/2) ln((M+u)/(M-u)) + M u` (so `F' = f`, `F(0) = 0`),
`F(m-1) = m^2 - 1 - 2m ln m` and `f(m-1) = (m-1) ln m`. QED (Runner: exact check for `m <= 2000`.)

Together with Lemma B (`gamma(m) >= gamma_-(m) := m ln2 - (1/2) ln(pi(m+2)/2)` and
`gamma(m) <= gamma_+(m) := m ln2 - (1/2) ln(pi(m + 1/2)/2)`), replace `P, Z, Psi` by the explicit lower
bounds `P_-`, `Z_-` (with `h/2` instead of `ceil(h/2)`: call it `Z~_-`) and `Psi_-`.

### 4.2 Minimising over h (tangent-line concavity) and Theorem D
Fix `n >= 16`, `h_c := floor((n-1)/2)`.
* **Path region `0 <= h <= h_c`.** By concavity of `ln`, `ln(pi(h+2)/2) <= tau(h) := ln(pi(h_c+2)/2) + (h - h_c)/(h_c + 2)`,
  so `P_-(h) >= P_tan(h)` with equality at `h_c`, where
  `P_tan(h) := (1/2)[(n - h)(ln pi + h ln2 - tau(h)/2) + ell_-(n-1-h)]`. Its second derivative is at most
  `(1/2)[-2 ln2 + 1/(h_c + 2) + 1] < 0` (use `ell_-''(g) = 1 - 1/(2g) + 1/(2g^2) <= 1`). So the minimum over the
  region is at `h = h_c` or `h = 0`.
* **Zigzag region `h_c + 1 <= h <= n - 1`.** With the tangent of `ln(pi(g+2)/2)` at `g_* := n - 2 - h_c`,
  `Z~_-(h) >= Z_tan(h)` with equality at `h_c + 1`, and `Z_tan'' <= -ln2 + 1/(2(g_* + 2)) + 1/2 < 0`. So the
  minimum is at `h = h_c + 1` or `h = n - 1`.

Hence `-ln P_tau >= Phi(n) := min(P_-(h_c), P_tan(0), Z~_-(h_c + 1), Z_tan(n-1))` for every
`tau != par(n)`, and `-ln P_{par(n)} >= Psi_-(n)`.

**Theorem D.** For every `d >= 17`, `U_d <= d^2 2^(2d-3) e^(-Phi(d-1)) + d 2^(d-2) e^(-Psi_-(d-1)) < 1`.
In particular a uniformly random subset of `V(Q_d)` is edge-multiset resolving with positive probability.

*Proof.* The number of pairs of types other than `par(n)` is `< d^2 2^(2d-3)`; `#par(n) = d 2^(d-2)`.
* `d = 17, 18, 19`: interval evaluation (runner) gives the bound `<= 0.6175`, `<= 0.01186`, `<= 0.00297`.
  (The binding term is `Z~_-(h_c + 1)`, the zigzag bound of the central cross type.)
* `d >= 20`: elementary lower bounds `E1..E4` for the four endpoint values and `E5` for `Psi_-`
  (obtained by bounding `h_c, g_c, g_*` by `(n-2)/2 <= . <= (n+1)/2` in each factor; listed in
  `thmD_crude` of the library and checked to be lower bounds for `16 <= n <= 3000`) satisfy
  `min(E1..E4) >= T(n) := 2n ln2 + 2 ln(n+1)` and `E5 >= (n-1) ln2 + ln(n+1) + ln 4`. At `n = 19` the
  margins are `5.24` (E2, binding) and larger for the others; each `E_i - T` has positive and increasing
  derivative for real `n >= 19` (the quadratic terms dominate: e.g. `E2'(19) = 4.62 > T'(19) = 1.49`).
  So the two terms are `<= 1/2` and `<= 1/4`. QED

### 4.3 Why `d = 11` is delicate; the Fourier-Hoelder lemma
* With the paper's optimal forests and exact atoms `U_10 = 1.307`, `U_11 = 0.156`, `U_12 = 0.0107`. The
  closed-form forests above give `V_14 = 2.7`, `V_15 = 0.52`: their per-type products are far from optimal at
  small `d` (the optimal forests are "potential forests": every level joins the best cell towards the centre;
  re-choosing them this way reproduces `U_11 = 0.1564`). An analytic estimate from `d = 11` would have to
  capture those cells within a factor of about 6 in total.
* **Lemma FH (Fourier-Hoelder multi-forest lemma, PROVED).** If `F_1, ..., F_k` are pairwise edge-disjoint
  forests in `Gamma_{e,f}` and `w_j > 0`, `sum w_j = 1`, then for every vector `v`,
  `P(H_e - H_f = v) <= prod_j prod_{ab in F_j} J_q(N_ab / w_j)^(w_j)`,
  `J_q(s) := (1/2pi) int_{-pi}^{pi} (1 - 4q(1-q) sin^2(t/2))^(s/2) dt <= e^(-x) I_0(x)`, `x = s q(1-q)`.
  At `q = 1/2`, `J(s) = Gamma((s+1)/2)/(sqrt(pi) Gamma(s/2 + 1))` (`= beta(s)` for even `s`).
  *Proof.* `P(D = v) <= (2pi)^(-L) int_{T^L} prod_{a<b} |phi_ab(theta_a - theta_b)| d theta` with
  `|phi_ab(t)| = (1 - 4q(1-q) sin^2(t/2))^(N_ab/2)`. Drop the factors outside `union F_j` (they are `<= 1`),
  apply Hoelder with exponents `1/w_j`, and integrate each forest leaf by leaf (translation invariance on the
  torus): `(2pi)^(-L) int prod_{ab in F} |phi_ab|^p = prod_{ab in F} J_q(p N_ab)`. QED
  (Runner: exact check on random small level graphs.) With `k` comparable forests this gains `k^(-1/2)` per
  forest slot. Certified (interval arithmetic, greedy disjoint forests, `k <= 3`): `U_10 <= 0.2649`,
  `U_11 <= 0.0448`. So the probabilistic method itself already works from `d = 10`.
* **A possible route to `d = 11` (OPEN).** The potential forests are compatible with the Pascal recursion
  `n^(d+1,h)_ab = n^(d,h)_ab + n^(d,h)_(a-1,b-1)`: shifting the upper-half edges of the `Q_d` potential forest
  by `+1` and keeping the lower-half edges gives an admissible forest for the same `h` in `Q_(d+1)` whose cells
  dominate the old ones, and similarly for `h -> h+1`. A monotonicity lemma `V_(d+1) <= V_d` for the
  potential-forest bound would reduce Open Problem 3 to the single value `V_11 = 0.1564`; the gain per step
  must come from the central levels (where the added Pascal term is comparable) and has not been made
  rigorous.
* EMPIRICAL (Monte Carlo): the true expected number of colliding pairs of a density-1/2 set is about
  `3.3`, `0.8`, `0.04` at `d = 8, 9, 10` (runner, fixed seed, 200 and 150 samples: `3.25`, `0.83`; a separate
  run with 120 samples at `d = 10`), and `P(resolving)` is about `0.27, 0.62, 0.97`. A sharp first-moment
  bound could reach `d = 9`, never `d = 8`.

### 4.4 Certificates for 13 <= d <= 16 (VERIFIED)
Found by `procgen_edim2_20261001_anneal.c` (descending annealing with on-the-fly distances), each re-verified
exactly by the annealer and independently by `is_resolving_exact` (numpy) in the runner:
`edim_m(Q_13) <= 105`, `edim_m(Q_14) <= 125`, `edim_m(Q_15) <= 135`, `edim_m(Q_16) <= 171`. The sets are
listed in the runner (section S7). They are 3.5-4.4 times smaller than the certified sparse random bounds
(`457, 490, 588, 636`) and about 6 times smaller than the previous lane's float bounds (`667, 781, 910, 1056`).

## 5. Open Problem 4: density 1/2 versus sparse random sets

**Answer.** Density 1/2 is far from optimal, and sparse random sets are optimal up to the constant in the
exponent.

* **Density 1/2 (PROVED, with CITED inputs).** A uniformly random subset `S_{1/2}` is resolving with
  probability `-> 1` (the paper's tail bound gives `U_d < C^5 d^12 2^(-3d-3) -> 0`; Theorem D gives the
  explicit bound for `d >= 17`), but `|S_{1/2}| = 2^(d-1) (1 + o(1))` with high probability. Since
  `ln edim_m(Q_d) <= (C* + o(1)) d^(1/3)` (Theorem A), the density-1/2 sets exceed the optimum by a factor
  `exp(d ln 2 - O(d^(1/3)))`: `ln |S_{1/2}| / ln edim_m(Q_d) >= (ln 2 / C* - o(1)) d^(2/3) -> infinity`.
* **Sparse random sets (PROVED).** For every density `q` with `ln(q 2^d) >= (C* + eps) d^(1/3)` (and
  `q <= 1/2`), `P(S_q resolving) -> 1` (Theorem A; the proof covers all such `q` uniformly, the exponent
  bounds being monotone in `lambda`). So sparse random sets reach `exp((C* + o(1)) d^(1/3))`, within the
  factor `C*/c* = 4^(2/3) = 2.52` (in the exponent) of the entropy lower bound, and Theorem C shows that no
  union-bound analysis of random sets (any density, any forests) can certify anything smaller.
* **Is C* the true threshold for random sets?** Heuristically yes: the covariance of the level-divergence
  vector of a pair is diagonal in the Krawtchouk basis, and a local limit theorem would give collision
  probabilities `exp(-(sqrt2/3) lambda^(3/2) sqrt(n) (1 + o(1)))` for bulk pairs, matching the forest bound.
  Not proved (OPEN).
* **EMPIRICAL: where random sets actually start to work.** Probability that a uniformly random `k`-subset is
  resolving (Monte Carlo, 30-60 samples per point): `d = 10`: 0.06 at `k = 100`, 0.22 at `k = 130`;
  `d = 11`: 0.10, 0.38, 0.62 at `k = 130, 160, 200`; `d = 12`: 0.07, 0.40, 0.63, 0.73 at `k = 130, 160, 200, 250`.
  So random sets need about `k = 180` at `d = 11, 12`: about half of the certified union bound
  (`361, 370`) and about 2.5-3 times the annealed sets (`65, 76`).
* **Better than random?** The annealed sets are 3.5-6 times smaller than the certified random bounds for
  `11 <= d <= 16` (table). Whether optimised sets beat random ones in the exponent constant is OPEN;
  the gap between `c*` and `C*` is exactly the gap between the entropy method (needs `d` bits, each heavy
  level carries at most `ln mu_r` nats) and the union bound (needs `2d` bits, each level gives only
  `(1/2) ln mu_r` nats of anti-concentration).

### 5.1 Quantitative comparison (finite d)

| d | lower bound | explicit set (VERIFIED) | certified sparse random bound `M_d` (this lane, Lemma FH) | previous float bound | density 1/2: `2^(d-1)` |
|---|---:|---:|---:|---:|---:|
| 6 | 7 | **15 = exact** (THM-4525) | - | - | 32 |
| 7 | **11** (Section 6; L4 gives 8) | 19 | - | - | 64 |
| 8 | 8 | 26 | - | - | 128 |
| 9 | 8 | 38 | - | - | 256 |
| 10 | 8 | 48 | - (but `U_10 <= 0.265` at `q = 1/2` by Lemma FH) | - | 512 |
| 11 | 8 | 65 | 361 | 511 | 1024 |
| 12 | 9 | 76 | 370 | 576 | 2048 |
| 13 | 9 | **105** | 457 | 667 | 4096 |
| 14 | 9 | **125** | 490 | 781 | 8192 |
| 15 | 10 | **135** | 588 | 910 | 16384 |
| 16 | 11 | **171** | 636 | 1056 | 32768 |
| 17..64 | 12..81 | - | see below | | |

Certified sparse bounds for `17 <= d <= 31` (interval arithmetic; `ln M_d / d^(1/3)` in brackets):
`d = 17: 772 (2.59)`, `18: 828 (2.56)`, `19: 988 (2.58)`, `20: 1058 (2.57)`, `21: 1255 (2.59)`,
`22: 1335 (2.57)`, `23: 1569 (2.59)`, `24: 1668 (2.57)`, `25: 1955 (2.59)`, `26: 2057 (2.58)`,
`27: 2374 (2.59)`, `28: 2498 (2.58)`, `29: 2889 (2.59)`, `30: 3016 (2.58)`, `31: 3470 (2.60)`,
`32: 3613 (2.58)`, `36: 5077 (2.58)`, `40: 6932 (2.59)`, `48: 12165 (2.59)`, `56: 20209 (2.59)`,
`64: 31808 (2.59)`. For comparison, density 1/2 gives about `2^63 = 9.2e18` at `d = 64`.
(The previous lane's double-precision values without Lemma FH had `ln M_d / d^(1/3) = 2.76-2.80`.)
The ratio `ln M_d / d^(1/3)` is flat near 2.58 in this range and decreases only slowly towards
`C* = 2.05` (the `ln d` corrections in Theorem A are large for moderate `d`).

## 6. Q_7 (optional item): 11 <= edim_m(Q_7) <= 19

Known before: `8 <= edim_m(Q_7) <= 19` (THM-4525). Upper bound: 40 further annealing runs (19 -> 18, three
million iterations per level) reached 19 twice and never 18 (best cost 8 colliding pairs), consistent with the
previous lane; the upper bound stays 19 (EMPIRICAL evidence that 18 is hard or impossible).

**Exhaustive search with a cheap complete normal form** (`procgen_edim2_20261001_q7search.c`).
Every `k`-set `S` of `Q_D` (`k >= 1`) is `Aut(Q_D)`-equivalent to `A x {0} u B x {1}` with
* coordinate `D-1` of maximal imbalance, `|A| = a >= b = |B|`, `a - b = beta_{D-1} >= |beta_i(S)|` (`i < D-1`);
* `0 in A` (translate by an element of `A` on coordinates `0..D-2`; only signs of `beta_i` change);
* the column sums `c_i(A) = #{x in A : x_i = 1}` non-increasing (permute coordinates `0..D-2`).

The program enumerates every such `A` (with redundancy, no canonical forms needed) and every `B` with the
imbalance constraints, and tests each leaf exactly (packed 4-bit histogram keys, `k <= 15`).

**Strengthened normal form (translate-maximal).** One may further require that for every `t in A` the sorted
column-sum vector of `A xor t` is lexicographically `<=` that of `A` (translate by a maximising `t in A`
first, then sort). This cuts the number of `A`'s by a factor 4-8.

**Canonical mode.** The automorphisms `x -> pi(x xor t)` of `Q_(D-1)` that map `A` to another set meeting all
three filters are exactly those with `t in A`, `sort(c(A xor t)) = c(A)` and `pi` sorting `c(A xor t)` (the
"residual group", usually tiny). Keeping `A` only if its mask is minimal among these images leaves exactly one
`A` per `Aut(Q_(D-1))`-orbit. Check: the numbers of accepted `A` equal the Burnside orbit counts of
`a`-subsets of `Q_6` (`1, 6, 16, 103, 497, 3253, 19735, 120843` for `a = 1..8`) and of `Q_5`
(`98804, 133576, 158658` for `a = 13, 14, 15`, as in THM-4525).

**Validation.**
* Leaf counts and `A`-counts agree with an independent Python enumeration for `(D,k) = (6,6), (6,7), (7,5)`
  (every `a`).
* Both normal forms: leaf and `A`-counts equal an independent Python enumeration for small `(D, k)`;
  300 random sets are mapped into the searched domains (runner).
* `D = 6`, `k = 15`: for `a = 14, 15` (imbalance 13, 15) no resolving set; for `a = 13` (imbalance 11)
  252 resolving leaves (plain form; 34 translate-maximal; exactly 3 in canonical mode), which canonicalise
  (`canon6` of the previous lane) to exactly the 3 orbits of THM-4525's list with maximal imbalance 11. So all
  three normal forms and the leaf test reproduce THM-4525.

**Results (FINITE-EXACT).** `Q_7` has no edge-multiset resolving set of size `k <= 10`:

| k | normal form | leaves (all `a`) | resolving |
|---|---|---:|---:|
| 1-7 | plain | 4.0e7 | 0 |
| 8 | plain | 5.4e8 | 0 |
| 9 | plain | 6.7e9 (a = 5..9: 3.3e8, 3.3e9, 2.4e9, 6.1e8, 6.0e7) | 0 |
| 9 | translate-maximal | 1.5e9 | 0 |
| 10 | translate-maximal | 1.55e10 (a = 5..10: 1.2e7, 5.5e9, 5.8e9, 3.6e9, 4.9e8, 5.7e7) | 0 |
| 10 | canonical | 9.2e8 (a = 5..10: 1.1e6, 2.0e8, 4.5e8, 2.2e8, 4.3e7, 3.6e6) | 0 |
| 11 | canonical, a = 6..10 only | 8.5e9 (2.7e8, 3.1e9, 3.6e9, 1.3e9, 2.3e8) | 0 |

Hence **`11 <= edim_m(Q_7) <= 19`** (the counting bound L4 gave 8), with `k = 10` confirmed by two different
normal forms. Running times on one shared core: `k = 9` about 11 min (plain) / 3 min (translate-maximal);
`k = 10` about 31 min (translate-maximal) / 3 min (canonical).
**`k = 11` is unfinished (OPEN):** the canonical search found no resolving 11-set with `6 <= a <= 10`
(`a = 11 - b`, imbalance `a - b`); the last case `a = 11, b = 0` (all eleven landmarks on one side of a
coordinate of maximal imbalance, 1.7e7 orbit representatives) was stopped at the session's usage limit.
Finishing it (`--q7k11` in the runner, about 30 more minutes) would give `edim_m(Q_7) >= 12`.
The `k = 10` canonical run and the partial `k = 11` run were logged outside the deposited `.out`
(same program and flags as the runner's S11).


## 6b. Audit guide (where each claim is proved and checked)

| claim | proof | runner check |
|---|---|---|
| forest lemma for general `q` | paper Lemma 11 verbatim | S2 (brute force, `Q_4`) |
| Lemma A, A2, B | Section 1 | S3 |
| Lemma S, star forests, path/zigzag forests and their cells | Sections 2.1-2.2 | S3 (exact, `d <= 40`) |
| Theorem A (C*), Theorem C | Sections 2.3-2.5 (asymptotic, no numerics) | S5 (constants) |
| Theorem B (c*) | Section 3 | S5 (L5 values) |
| Proposition K | Section 2.5 | S3, last check (exact determinants, `d <= 11`) |
| Lemma FH | Section 4.3 | S4 (exact, small graphs) |
| Lemma PI, Theorem D | Sections 4.1-4.2 | S6 (intervals at `d = 17, 18, 19`; crude bounds) |
| certificates `Q_6..Q_16` | - | S7 (exact) |
| `U_10 <= 0.2649` etc. | Lemma FH | S8 (intervals) |
| sparse bounds `M_d` | forest lemma, Lemma A / `e^-x I_0(x)`, Lemma FH | S9 (intervals, every listed `d`) |
| first moments, random thresholds | - | S10 (EMPIRICAL, fixed seeds) |
| `edim_m(Q_7) >= 11` (or more, Section 6) | normal forms (Section 6) | S11 |

## 7. Files and reproduction

* Note: this file.
* Library: `04-computation/experiments/procgen_edim2_20261001_lib.py` (cells, pair types, forests, exact
  max atoms, interval arithmetic via `mpmath.iv`, Lemma A / `e^-x I_0(x)` atom bounds, Fourier-Hoelder
  bounds, certified union bounds and sizes, Theorem D functions, exact verifier `is_resolving_exact`).
* Annealer: `04-computation/experiments/procgen_edim2_20261001_anneal.c` (memory-light; distances on the fly;
  every reported set re-verified exactly). Usage `anneal d k0 iters seed T0 T1 kmin`.
* Exhaustive search: `04-computation/experiments/procgen_edim2_20261001_q7search.c`
  (usage `q7search D k a [maxprint] [tmax]`); the runner also compiles the previous lane's
  `procgen_edim_20261001_canon6.c` for the `Q_6` control.
* Runner: `04-computation/experiments/procgen_edim2_20261001_run.py` (stdout only, ends with ALL CHECKS
  PASSED). Options: `--quick` (a few minutes), `--q7k11` (adds the unfinished `Q_7`, `k = 11` search).
  **The deposited output `05-knowledge/results/procgen_edim2_20261001.out` is a `--quick` run** (71 checks, ALL
  CHECKS PASSED, 161 s, peak RSS 251 MB; the session hit its usage limit). Quick mode re-checks everything except: Lemma FH on the larger random sample, the sparse
  bounds for `d >= 15`, the plain-form `Q_6` control, the `Q_7` searches for `k = 10` (it does `k <= 9`), and the
  random-threshold Monte Carlo. A full run (no flag, about 40 minutes) re-checks those too; an earlier full
  run of the same checks (before the canonical mode was added) passed all 83 checks up to S11 and was stopped
  during the `k = 10` search, and the `k = 10` searches were completed separately with the same program.
* Interval arithmetic: `mpmath.iv` (outward rounding) for every transcendental quantity in a certified claim;
  exact `fractions.Fraction` for exact atoms. Double precision is used only to choose `q` and forests (any
  choice is valid) and in statements labelled EMPIRICAL or "illustration".
