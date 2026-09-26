# Thin divergence: every non-periodic orbit of x -> x/2, (3x+b)/2 has at most X^(0.95+eps) elements below X; reciprocal sums converge (HYP-9160 proved); the real value of the Bernstein series is strictly below the integer

**Status: PROVED (elementary; full proof below; self-audited; independent
adversarial audit by a separate agent, verdict SOUND, recorded in the theorem
file THM-4476 with its own script and output; a first draft's extension to
non-constant sign strategies was withdrawn during self-audit, see Corollary 4) with FINITE-EXACT controls of
the two counting lemmas. No priority claim: the ingredients are Terras's
stopping-time count and an injectivity pigeonhole; the single-orbit
statement may exist in the literature under another name. Collatz, the
`3n-1` sheet and the Periodicity Conjecture remain OPEN: the theorem
constrains the shape of a hypothetical divergent orbit, it does not exclude
one. Session `collatz-squares-doubles-20260925` (opus), 2026-09-25.**

Theorem file: [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
Controls: `04-computation/experiments/collatz_thin_20260925_counts.py` and
`collatz_thin_20260925_controls2.py`, outputs `collatz_thin_20260925_counts.out`
and `collatz_thin_20260925_controls2.out`.
Origin: the owner asked to prove HYP-9160 for discrepancy `O(log l)`; the
in-house [no-bounded-strip theorem](collatz_guards_20260921_discrepancy.md)
(section 2a) supplied the density mechanism, and a pigeonhole on landing
points removed its need for a strip.

## 0. Statement

Throughout, `b` is an odd integer, `T_b(x) = x/2` for even `x` and
`(3x + b)/2` for odd `x`, acting on `Z`. An orbit `(x_i)_(i >= 0)` is
*eventually periodic* if `x_i = x_j` for some `i < j`; otherwise its terms
are pairwise distinct. Write

```text
h(p) = -p log_2 p - (1-p) log_2 (1-p),     h* = h(log_3 2) = 0.949956...,
N_T(X) = #{ i >= 0 : |x_i| <= X }.
```

**Theorem (thin divergence).** Let `b` be odd and let `(x_i)` be a
`T_b`-orbit in `Z \ {0}` that is not eventually periodic. Then for every
`eps > 0` there is a constant `C = C(eps, |b|)`, independent of the orbit,
with

```text
N_T(X) <= C X^(h* + eps)        for all X >= 1.
```

**Corollary 1 (reciprocal sums).** `sum_i 1/|x_i| < infinity` for every such
orbit. In particular every divergent positive orbit of `3n+1` or `3n-1`
has `sum_l 1/m_l < infinity` over its odd iterates: **HYP-9160 holds**.

**Corollary 2 (the two places).** Let `d = (v_1, v_2, ...)` be a halving
word (`v_l >= 1`), `d_L = v_1 + ... + v_L`, and suppose `d` is not
eventually periodic. If the 2-adic value `R_2(d) = sum_l 2^(d_l) 3^(-l-1)`
is a rational number `r` with odd denominator, then the real series
`R(d) = sum_l 2^(d_l)/3^(l+1)` converges. If moreover `r = n` is a positive
integer, then `R(d) < n` strictly, and the `3n-1` orbit of `n` diverges at
full rate: `m_L ~ c 3^L/2^(d_L)` with `c = n - R(d) = n prod_l (1 - 1/(3 m_l)) > 0`.
If `-r = n` is a positive integer, the `3n+1` orbit of `n` diverges at full
rate with `c = n + R(d) < infinity`.

**Corollary 3 (growth and discrepancy).** For every non-eventually-periodic
positive orbit of `T_(+-1)` with odd iterates `m_j` and every `a < 1/h* = 1.052681...`,
`m_j > j^a` for infinitely many `j`. Consequently:
* no such orbit has `m_j <= K j^a` for all `j` with `a < 1.05268`;
* on the plus sheet, no such orbit has `|Delta_j| <= C log_2 j` for all
  large `j` with `C < (1/h* - 1)/2 = 0.02634`, where `Delta_j = d_j - j log_2 3`;
  on the minus sheet none has `|delta_j| <= C log_2 j` with `C < 1.05268`;
* no such orbit has bounded discrepancy (the in-house theorem is recovered).

**Corollary 4 (the union of all cycles, and injective invariant sets).**
Let `A` be a subset of `Z \ {0}` with `T_b(A)` contained in `A` on which `T_b`
is injective (each element of `A` has at most one preimage in `A`). Then
`#{n in A : |n| <= X} <= C(eps,|b|) X^(h*+eps)`. In particular the set of
all periodic points of `T_b` (the union of all its cycles) is thin, and a
single cycle of `T_b` with maximum `M` has at most `C M^(h*+eps)` elements.

**Withdrawn claim (recorded as a hostile).** An earlier draft asserted the
theorem for the level-`k` sign strategies `T_sigma` of THM-4474. That is
unproved and the method fails for them: the parity map `Z/2^k -> {0,1}^k`
is a bijection only when the odd shift `b` is constant (injectivity needs
`3y + b(y) = 3y' + b(y') mod 2^K` to force `y = y'`, which fails as soon as
`b` takes both values on odd residues). Witness: for `sigma = -chi_(-4)`
every odd `n` has the all-ones word, so the no-dip count is `X/2`, not
`X^0.95`. Each single orbit of that strategy happens to be thin (it grows
like `(3/2)^i`), but no bound follows from this method for any
non-constant strategy.

**Corollary 5 (the word of a divergent plus-sheet orbit).** Let `(m_l)` be a
divergent positive `3n+1` orbit with halving word `d` and discrepancy
`Delta_l = d_l - l log_2 3`. Then `sum_l 2^(Delta_l) = 3 R(d) < infinity`; in
particular `Delta_l -> -infinity`, `2^(d_l)/3^l -> 0`, and
`m_l 2^(d_l)/3^l -> c in (n, infinity)`. On the minus sheet
`sum_l 2^(delta_l) = 3 R(d) < 3n` for every divergent positive orbit.

**Corollary 6 (uniformity).** The constant `C(eps, |b|)` does not depend on
the orbit. Hence for each odd `b` there is `K = K(b)` such that every
`T_b`-orbit with distinct terms satisfies `sum_i 1/|x_i| <= K`. For `b = +-1`
and positive orbits this gives absolute constants in the two-place
identity: on the plus sheet `n < c <= e^(K/3) n` and `R(d) <= (e^(K/3) - 1) n`;
on the minus sheet `c >= kappa n` with `kappa = exp(-K/3 - K^2/9) > 0`, so
`R(d) <= (1 - kappa) n`. So `|R(d)| <= K' |R_2(d)|` for every non-eventually-
periodic word whose 2-adic value is a nonzero integer, with `K'` absolute.

**Corollary 7 (branches of the inverse tree).** Every infinite branch
`x_0, x_1, x_2, ...` of the inverse tree (`T_b(x_(i+1)) = x_i`, terms
distinct) satisfies `#{i : |x_i| <= X} <= C(eps,|b|) X^(h*+eps)`. More
generally the theorem holds for every injective sequence, one- or
two-sided, in which consecutive terms are related by `T_b`.

**Corollary 8 (harmonic reformulation).** For every `n in Z \ {0}` and odd
`b`, the orbit of `n` under `T_b` is eventually periodic **iff**
`sum_(i>=0) 1/|T_b^i(n)| = infinity` (with the convention that an orbit
reaching `0` counts as eventually periodic). Hence the **no-divergence half**
of the Collatz conjecture (every positive orbit is eventually periodic) is
equivalent to: `sum_(i>=0) 1/T^i(n) = infinity` for every positive integer
`n`. It says nothing about the cycle half: a hypothetical nontrivial cycle
also has a divergent reciprocal sum (audit finding 13). The `3n-1`
no-divergence statement and Lagarias's Periodicity Conjecture for rationals
with odd denominators, which are no-divergence statements, have the same
harmonic form.
Moreover the Dirichlet series `sum_i |x_i|^(-s)` of a non-eventually-periodic
orbit converges for every `s > h*`: its abscissa of convergence is at most
`h(log_3 2)`.

**Corollary 8' (a uniform criterion).** With `K = K(b)` from Corollary 6: an
orbit is eventually periodic iff some partial sum `sum_(i<=N) 1/|x_i|` exceeds
`K`. So the no-divergence half of Collatz is equivalent to: for every
`n >= 1` the partial reciprocal sums of the orbit of `n` eventually exceed
the absolute constant `K(1)`.
(The proof makes `K` explicit in principle but astronomically large.)

**Hostile control for the union question.** For the level-2 strategy
`-chi_(-4)` of THM-4474 every odd `n` satisfies `v_2(3n + sigma(n)) = 1`, so
every odd orbit increases forever; each orbit is thin (exponential growth,
`N(X)` of order `log X`), yet their union is all odd integers. So the
single-orbit theorem cannot be summed over orbits, and any bound on the
density of divergent integers needs a different mechanism.

**Corollary 9 (2-adic irrationality of log-drift words).** Let `d` be a
halving word that is not eventually periodic, with discrepancy
`Delta_j = d_j - j log_2 3`. If `Delta_j >= -a log_2 j - O(1)` for all large
`j` with some `a < 1/h* = 1.052681`, then the 2-adic number
`R_2(d) = sum_l 2^(d_l) 3^(-l-1)` is irrational. This covers every word
with bounded discrepancy (recovering Proposition B of the transversality
foundry, including the Sturmian words of Theorem S, by a counting proof
instead of a Padé argument) and the words `d_j = floor(j log_2 3 - a log_2 j) + O(1)`
with `a < 1.05268`.

**The exponent ladder.** For a word `d` and a positive integer `n`, write
`m_j(n) = 2^(-Delta_j)(n + (1/3) sum_(i<j) 2^(Delta_i))` for the formal
plus-sheet orbit and let `a(d)` be its growth exponent
(`limsup log m_j / log j`). Distinctness of orbit points alone excludes
`a(d) < 1`; THM-4476 excludes `a(d) < 1/h* = 1.05268`; the Periodicity
Conjecture asserts that every non-eventually-periodic word is excluded.
Words with `Delta_j = -a log_2 j + O(1)` have `a(d) = max(1, a)`, so the new
region is exactly `1 <= a < 1.05268`; words with `Delta_j = o(log j)`
(including `Delta_j -> +infinity` slowly, and all sub-logarithmic drifts)
have `a(d) = 1` and are excluded. Using the 3-adic (backward) residue
classes instead of, or together with, the 2-adic ones does not lower the
exponent: each place costs the same `0.05` bits per bit of modulus, and the
total modulus cannot exceed `X`.

## 1. Proof of the theorem

### 1.1 Reductions

**Distinctness.** If `x_i = x_j` with `i < j` the orbit is eventually
periodic. So the `x_i` are pairwise distinct.

**One-signed segments.** If `x` is odd and `|x| > |b|/3`, then `3x + b` has
the sign of `x`; if `x` is even, `x/2` has the sign of `x`. So a sign change
`x_i x_(i+1) < 0` needs `x_i` odd with `|x_i| <= |b|/3`, and the set `S` of
such integers has at most `|b|/3 + 1` elements, each visited at most once.
Hence the orbit is a concatenation of at most `|b|/3 + 2` maximal
one-signed segments, all but the last finite. Replacing `x` by `-x`
conjugates `T_b` to `T_(-b)`, and `|b|` is unchanged, so it suffices to
prove, for a **positive** segment, the bound of section 1.5 with a constant
depending on `|b|` only; the sum over segments then gives the theorem.

So fix a segment `y_0, y_1, ..., y_M` (`M <= infinity`) of pairwise distinct
positive integers with `y_(i+1) = T_b(y_i)`, and let
`N(X) = #{ i <= M : y_i <= X }`.

### 1.2 The parity map is a bijection modulo 2^k (Terras)

For `k >= 1`, the parity word `(y mod 2, T_b(y) mod 2, ..., T_b^(k-1)(y) mod 2)`
depends only on `y mod 2^k`, and `y -> word` is a bijection
`Z/2^k -> {0,1}^k`.

*Proof.* `T_b` maps each class modulo `2^k` into a class modulo `2^(k-1)`:
for even `y`, `y/2 mod 2^(k-1)` is determined by `y mod 2^k`; for odd `y`,
`(3y+b)/2 mod 2^(k-1)` is determined by `3y + b mod 2^k`. By induction the
word of length `k` is a function of `y mod 2^k`. Injectivity: if `y, y'` have
the same word of length `k`, they have the same parity and, by induction,
`T_b(y) = T_b(y') mod 2^(k-1)`; for even `y, y'` this gives `y = y' mod 2^k`,
and for odd `y, y'` it gives `3y = 3y' mod 2^k`, hence `y = y' mod 2^k`.
There are `2^k` classes and `2^k` words. ∎

### 1.3 Affine form and the carry bound

For a word `w` of length `i` with `o` odd letters,

```text
T_b^i(y) = 3^o y / 2^i + beta_i(w),     |beta_i(w)| <= |b| ((3/2)^i - 1),
```

where `beta_i` depends only on `w`. *Proof.* `beta_0 = 0`; an even step
sends `beta -> beta/2`, an odd step `beta -> (3 beta + b)/2`, so
`|beta_(i+1)| <= (3/2)|beta_i| + |b|/2`, and the bound follows by induction. ∎

### 1.4 The counting lemma

Let `alpha = log_2 3`, `theta_0 = 1 - alpha/2 = 0.207519`, and for
`theta in (0, theta_0)` put `rho = (1 - theta)/alpha in (1/2, log_3 2)`.
For `X >= 2` let `k = floor(log_2 X)` and

```text
F_b(X, theta) = { y in [1, X] : T_b^i(y) >= y X^(-theta) for all 0 <= i <= k }.
```

**Lemma.** There are `k_0(theta)` and `A(theta) = 2 (rho/(1-rho))^2` such
that for `X >= 2^(k_0)`,

```text
#F_b(X, theta) <= 2|b| X^(log_2(3/2) + theta) + A(theta) (log_2 X + 1) X^(h(rho)).
```

*Proof.* Put `Y_0 = 2|b| X^(log_2(3/2) + theta)`. Elements `y < Y_0` of `F_b`
number at most `Y_0`. Let `y in F_b` with `y >= Y_0`, and let `o` be the
number of odd letters in its word of length `k`. Since `(3/2)^k <= X^(log_2(3/2))`,
section 1.3 gives `|beta_k| <= |b| X^(log_2(3/2)) <= y X^(-theta)/2`, so

```text
3^o / 2^k >= X^(-theta) - |beta_k|/y >= X^(-theta)/2.
```

Taking logarithms, `o alpha >= k - theta log_2 X - 1 >= (1 - theta) k - 2`
(because `log_2 X < k + 1` and `theta < 1`), i.e. `o >= rho k - 2`. The
number of words of length `k` with at least `rho k - 2` odd letters is at
most `(k+1) max_(o >= rho k - 2) C(k, o) <= (k+1) 2^(k h(rho - 2/k))` for
`k >= k_0(theta)` (so that `rho k - 2 > k/2`, using `C(k,o) <= 2^(k h(o/k))`
and the monotonicity of `h` on `[1/2, 1]`). Since `|h'(p)| = |log_2((1-p)/p)| <= log_2(rho/(1-rho))`
on `[1/2, rho]`, `2^(k h(rho - 2/k)) <= (rho/(1-rho))^2 2^(k h(rho))`. By
section 1.2 each word is one residue class modulo `2^k`, and since
`2^k > X/2` such a class meets `[1, X]` in at most two integers. Hence the
elements `y >= Y_0` of `F_b` number at most
`2 (k+1) (rho/(1-rho))^2 2^(k h(rho)) <= A(theta) (log_2 X + 1) X^(h(rho))`. ∎

The exact binomial sums are compared with the entropy bound in the control
script (`(L1)`), and the sets `F_(+-1)(X, theta)` are enumerated for
`X <= 2^20` (`(L2)`) and `2^24`: at `theta = 0.03` the value `log #F / log X`
is `0.78` at `2^20` and `0.80` at `2^24` on both sheets, but the local
doubling slopes are `0.92`-`0.96`, i.e. at the lemma's exponent `0.9635`; the
counting lemma is essentially sharp for `F` up to logarithmic factors
(audit finding 16).

### 1.5 The dichotomy and the recursion

Fix `theta in (0, theta_0)`, `X >= 2^(k_0)`, `k = floor(log_2 X)`. Every
index `i <= M` with `y_i <= X` is of one of three kinds:

* **(E)** `i > M - k` (only when `M < infinity`): at most `k` indices;
* **(D)** `i <= M - k` and `y_(i+s) < y_i X^(-theta)` for some `1 <= s <= k`;
* **(ND)** `i <= M - k` and `y_(i+s) >= y_i X^(-theta)` for all `0 <= s <= k`.

**(D) is a pigeonhole.** For a (D)-index take the least such `s`. Then
`y_(i+s) < X^(1-theta)`, and the pair `(i+s, s)` determines `i`. So
`#(D) <= k N(X^(1-theta))`: the dipping points land on the segment's own
points below `X^(1-theta)`, and each landing point serves at most `k`
dippers. This index count does not use distinctness.

**(ND) is the counting lemma.** An (ND)-index has `y_i in F_b(X, theta)`,
and distinct indices have distinct `y_i`: this, with section 1.1, is where
distinctness of the orbit is used. So `#(ND) <= #F_b(X, theta)`.

Therefore, for `X >= 2^(k_0)`,

```text
N(X) <= k + k N(X^(1-theta)) + 2|b| X^(log_2(3/2)+theta) + A(theta)(log_2 X + 1) X^(h(rho)).   (R)
```

### 1.6 Bootstrap

Fix `gamma` with `h(rho) < gamma < 1` (then also `log_2(3/2) + theta < gamma`
for the `theta` used below, since `log_2(3/2) = 0.585`). Choose `X_1 >= 2^(k_0)`
with `(log_2 X + 1) X^(-theta gamma) <= 1/4` for `X >= X_1`, and then
`C >= X_1` so large that for all `X >= 1`

```text
k + 2|b| X^(log_2(3/2)+theta) <= (C/4) X^gamma   and   A(theta)(log_2 X + 1) X^(h(rho)) <= (C/4) X^gamma,
```

which is possible because both left sides are `o(X^gamma)`. Claim:
`N(X) <= C X^gamma` for all `X >= 1`. For `X < X_1` this holds because
`N(X) <= X <= X_1 <= C X^gamma`. For `X >= X_1`, strong induction and (R)
give `N(X) <= (C/4)X^gamma + k C X^((1-theta)gamma) + (C/4)X^gamma <= (3/4) C X^gamma`,
using `k X^(-theta gamma) <= 1/4`. ∎

Since `gamma > h(rho(theta))` is arbitrary and `h(rho(theta)) -> h(log_3 2) = h*`
as `theta -> 0`, `N(X) = O_eps(X^(h*+eps))` for every `eps > 0`, with a
constant depending only on `eps` and `|b|`. Summing over the at most
`|b|/3 + 2` one-signed segments proves the theorem. ∎

### 1.7 Proofs of the corollaries

*Corollary 1.* By Abel summation, `sum_(|x_i| <= Y) 1/|x_i| = N_T(Y)/Y + int_1^Y N_T(X) X^(-2) dX`,
and `N_T(X) <= C X^gamma` with `gamma < 1` makes the right side bounded. ∎

*Corollary 2.* The map `d -> R_2(d)` inverts the parity-vector bijection of
`Z_2` onto words, and `n = -b R_2(d)` for the `T_b`-orbit of `n`
(Proposition 6(6) of the [foundry note](collatz_sqdbl_20260925_squares_doubles_foundry.md)).
So `R_2(d) = r` rational with odd denominator `D` means that `d` is the
`T_(-1)`-word of `r`, equivalently the `T_(-D)`-word of the integer `rD`,
which is not eventually periodic. By Corollary 1 its terms `m_l` satisfy
`sum 1/|m_l| < infinity`; no term is `1/3` (that orbit reaches `0` and its
word ends), and a factor `1 - 1/(3 m_l)` is negative only for `0 < m_l < 1/3`,
which happens finitely often along the distinct integer orbit `D m_l`; so
the product `prod_l (1 - 1/(3 m_l))` converges to a finite nonzero limit,
and the algebraic identity `R_L(d) = r (1 - prod_(l<L)(1 - 1/(3 m_l)))`
(Proposition 6(3), valid for every sign) shows that `R(d)` converges. If
`r = n` is a positive integer the `3n-1` orbit stays positive, every factor
lies in `(0,1)`, and `R(d) = n - c` with `c = n prod (1 - 1/(3m_l)) > 0`.
The `-r = n` case is the plus sheet with `c = n prod (1 + 1/(3 m_l)) < infinity`. ∎

*Corollary 3.* If `m_j <= j^a` for all `j >= j_1` then
`N(X) >= X^(1/a) - j_1`, contradicting `N(X) <= C X^(h*+eps)` when
`1/a > h*`. On the plus sheet `|Delta_j| <= C' log_2 j` gives, by
`m_j = 2^(-Delta_j)(n + (1/3) sum_(i<j) 2^(Delta_i))`,
`m_j <= j^(C')(n + j^(1+C')/(3(1+C'))) <= K j^(1+2C')`, excluded when
`1 + 2C' < 1/h*`. On the minus sheet `m_j = 2^(-delta_j) c_j <= n j^(C')`,
excluded when `C' < 1/h*`. A bounded strip is the case `C' = 0`. ∎

*Corollary 4.* Sections 1.4 and 1.5 use only two properties of the
sequence: distinct terms, and that the forward window of each term lies in
the set. For an invariant set `A` on which `T_b` is injective, define
`N_A(X) = #{n in A : |n| <= X}` and repeat section 1.5 with the dichotomy
applied to every `n in A`, `|n| <= X`: a dipper `n` has `z = T_b^s(n) in A`
with `|z| < X^(1-theta)` and least `s <= k`, and `(z, s)` determines `n`
because `T_b^s` is injective on `A`. Signs: discard the elements `n` with
`n in S` or `T_b^s(n) in S` for some `s <= k`, at most `(k+1)|S|` of them by
injectivity; every remaining window is one-signed, the positive non-dippers
lie in `F_b(X, theta)` and the negated negative ones in `F_(-b)(X, theta)`.
So `N_A(X) <= k N_A(X^(1-theta)) + #F_b(X,theta) + #F_(-b)(X,theta) + (k+1)(|b|/3 + 1)`,
with base case `N_A(X) <= 2X`, and the bootstrap of 1.6 applies unchanged. The periodic points form such a set: a periodic point has a
unique periodic preimage. ∎

*Corollary 5.* `R(d) = (1/3) sum_l 2^(Delta_l)` by definition, and
`R(d) < infinity` by Corollary 2 (or directly: `R(d) = n(prod(1 + 1/(3m_l)) - 1)`
with `sum 1/m_l < infinity`). Then `Delta_l = log_2 c_l - log_2 m_l` with
`c_l -> c` finite and `m_l -> infinity`. The minus-sheet statement is
Proposition 6(4) with the strictness of Corollary 2. ∎

*Corollary 6.* The constants `X_1`, `A(theta)`, `k_0(theta)` and `C` in
sections 1.4 and 1.6 depend on `theta`, `gamma` and `|b|` only, and the number
of one-signed segments is at most `|b|/3 + 2`. Abel summation gives
`sum_(|x_i| <= Y) 1/|x_i| = N_T(Y)/Y + int_1^Y N_T(X) X^(-2) dX <= C(1 + 1/(1-gamma)) =: K`.
For the products use `log(1+x) <= x` and, for `0 < x <= 1/3`,
`log(1-x) >= -x - x^2`, with `sum x_l = (1/3) sum 1/m_l <= K/3` and
`sum x_l^2 <= (sum x_l)^2`. ∎

*Corollary 7.* The proof of section 1.5 uses only that the sequence is
injective and that the `k` terms following each term (in the forward
direction of `T_b`) belong to the sequence. For a backward branch the
forward window of `x_i` is `x_(i-1), ..., x_(i-k)`, inside the branch for
`i >= k`; the fewer than `k` remaining indices are the (E) class. ∎

*Corollary 8.* An eventually periodic orbit repeats a positive term
infinitely often, so its reciprocal sum diverges; a non-eventually-periodic
orbit has distinct terms and Corollary 1 applies. For the Dirichlet series,
`sum_i |x_i|^(-s) = int_1^infinity X^(-s) dN_T(X)` converges when
`N_T(X) = O(X^gamma)` with `gamma < s`. ∎

*Corollary 9.* If `R_2(d) = r` is rational (necessarily with odd
denominator, as `r in Z_2`), then `d` is the `T_(+1)`-word of `-r` and the
`T_(-1)`-word of `r` (Proposition 6(6)), and `-r` has a non-eventually-
periodic orbit, thin by the theorem applied to `T_(+-D)` with `D` the
denominator. After discarding a finite prefix, the orbit has constant sign;
negating if necessary gives a positive orbit on one of the two sheets.
Corollary 2's convergent nonzero normalized product makes
`c_j=m_j 2^(Delta_j)` bounded on either sheet, with time and discrepancy
rebased at that tail. A finite shift changes the discrepancy condition only
by a bounded additive constant.
The assumed one-sided discrepancy bound gives `2^(-Delta_j)<=C j^a`, hence
`m_j<=C' j^a` for all large j. If a<=0 this contradicts distinctness in
the fixed rational lattice. If 0<a<1/h*, choose gamma with h*<gamma<1/a;
the first J distinct terms would lie below O(J^a), contradicting the
thinness bound `J<=O(J^(a gamma))`. This proves the claim for the entire
stated range without a separate partial-sum growth assumption. ∎

**2026-09-26 proof refinement (crossroads audit).** The previous audit in
`777e4e137` supplied the missing bounded-product justification. The rewrite
above makes it explicit and removes the unnecessary bound
`sum j^(-a-o(1))=infinity`, which is not valid at a=1 without controlling
the o(1) term. It also gives the stronger direct consequence: on either
sheet no injective positive orbit has
`Delta_j>=-a log_2 j-O(1)` eventually for any a<1/h*. The older symmetric
plus-sheet threshold 0.02634 in Corollary 3 remains true but is weaker.

## 1.8 The same proof for every contracting Collatz-like map (Matthews–Watts framework)

Let `m >= 2` and let `g(x) = (p_i x + q_i)/m` for `x = i mod m`, with integers
`p_i, q_i`, `p_i >= 1` coprime to `m`, and `q_i = -p_i i mod m` so that `g`
maps `Z` to `Z` (Conway's class; `m = 2`, `(p_0, q_0) = (1, 0)`,
`(p_1, q_1) = (3, b)` is `T_b`). Write `a_i = p_i/m` and define the Chernoff
rate

```text
I(g) = -log_m min_(lambda >= 0) (1/m) sum_i a_i^lambda,
```

which is positive exactly when the geometric mean of the `a_i` is below 1
(`prod_i p_i < m^m`, the *contracting* case of Matthews–Watts). For `T_(+-1)`,
`I = 1 - h* = 0.050044` at `lambda = 0.488`.

**Theorem (general form; PROVED, same proof).** If `prod_i p_i < m^m` and
`max_i p_i < m^2`, then every `g`-invariant subset of `Z \ {0}` on which `g`
is injective (in particular every orbit with distinct terms, every infinite
backward branch, and the union of all cycles) has at most
`C(eps, g) X^(max(1 - I(g), log_m max_i a_i) + eps)` elements of absolute
value at most `X`.

*Proof.* Replace `2` by `m` throughout sections 1.1–1.6. The residue word
`(g^j(x) mod m)_(j<k)` is a bijection `Z/m^k -> (Z/m)^k` because each
`p_i` is a unit modulo `m` (section 1.2 verbatim). The affine form is
`g^k(y) = (prod_j a_(i_j)) y + beta_k` with
`|beta_k| <= max|q_i| ((max a)^k - 1)/(max a - 1)` when `max a > 1` (and
`O(k max|q_i|)` otherwise). A no-dip word of length `k = floor(log_m X)`
from `y >= Y_0 = 2 K X^(log_m max a + theta)` has
`prod_j a_(i_j) >= X^(-theta)/2`; the number of such words is at most
`m^(k(1 - I_theta(g)))` by Chernoff, with `I_theta -> I(g)` as `theta -> 0`.
Sign changes occur only at `|x| <= max|q_i|/min p_i`, a finite set. The
dichotomy, pigeonhole and bootstrap are unchanged. ∎

So the single-orbit thinness is exactly a feature of the contracting regime
of the Matthews–Watts conjecture (their Conjecture A: all orbits of a
contracting map are eventually periodic). In the expanding regime
(`5x+1`: geometric mean `sqrt(5/4) > 1`) the rate `I` is zero and the
statement is not expected to be provable this way; there, divergent orbits
are expected to exist and to be exponentially thin for a different reason.

## 2. What the theorem says and does not say

* **Sheet.** The counting bound is the same on both sheets; the sign enters
  only in the consequence. On the minus sheet it turns Proposition 6's
  inequality `R(d) <= n` into the strict inequality `R(d) < n` for every
  divergent orbit, so the real and 2-adic values of a non-periodic word
  never coincide at a positive integer (the equality case of the two-place
  identity is exactly the eventually periodic case). On the plus sheet it
  says every divergent orbit has a finite real Bernstein value.
* **Drift.** The proof uses `log_2 3 < 2`: `rho_0 = log_3 2 > 1/2` makes the
  no-dip words rare (`h(rho_0) < 1`). For `5n+1` one has `1/log_2 5 = 0.43 < 1/2`
  and the lemma is vacuous, as it must be: `5n+1` is expected to have
  divergent orbits, and nothing here contradicts them since they are
  expected to be exponentially thin anyway.
* **Defect.** The theorem is about one orbit of one map. Density-zero
  flips of the pairing family (THM-4470 section 4) produce divergent orbits
  of *other* maps, to which Corollary 4 does not apply (they are not
  residue-determined).
* **Not a divergence exclusion.** A divergent orbit is expected to grow
  exponentially, with `N(X)` of order `log X`. The theorem excludes only the
  polynomially dense shapes, exactly the ones for which the real place
  (Proposition 6) had no answer. It does not bound the number of integers
  below `X` lying on *some* divergent orbit: the pigeonhole in section 1.5
  fails for a union of orbits, because the Collatz graph branches.
* **Exponent.** `h*` is the Hausdorff dimension of the 2-adic exceptional
  set `Bad_infinity` (choice ladder, PROVED there). The theorem says an
  integer orbit cannot be denser than that set. Nothing in the method
  reaches below `h*`.
* **Uniformity.** The constant depends only on `eps` and `|b|`; in
  particular the same bound holds for all rational orbits with a fixed odd
  denominator.

## 3. Controls (FINITE-EXACT)

From `collatz_thin_20260925_counts.out`:

| quantity | value |
|---|---|
| `h* = h(log_3 2)` | `0.949956` |
| `1/h*` | `1.052681` |
| `theta_0` | `0.207519` |
| plus-sheet log-band threshold `(1/h* - 1)/2` | `0.02634` |
| un-bootstrapped exponent (`1 - theta = h(rho(theta))`) | `0.96538` at `theta = 0.03462` |
| `(L1)` exact `sum_(o >= rho k - 2) C(k,o)` vs `2^(k h(rho))`, `k <= 320` | ratio below `k+1` in every case |
| `(L2)` `#F_(+-1)(2^20, 0.03)` | `46612` / `46611`; `log F/log X = 0.775`; local doubling slopes `0.92`-`0.96` up to `2^24`, at the bound exponent `0.9635` |
| `(L2)` `#F_(+-1)(2^20, 0.10)` | `157652` / `157648`; `log F/log X = 0.863`; local slopes `0.99`-`1.02` (bound exponent `0.9867`); at `2^20` the threshold `X^(-0.1) = 1/4` is exact and about 7,500 boundary elements with `T^i(y) = y/4` are counted |

The two sheets differ by at most four elements at every size, as the
sheet-blind class count predicts.

Further controls (`collatz_thin_20260925_controls2.py` / `.out`):

| control | result |
|---|---|
| (C1) union of all cycles of `T_b` with minimum `<= 2*10^5`, twelve odd `b` | at most 86 periodic points (b = 13, ten cycles); `count / X^(h*) <= 1.04` in every case |
| (C2) contracting `m = 3` map `x/3, (2x+1)/3, (4x+1)/3` (`prod p = 8 < 27`, `max p = 4 < 9`) | residue-word bijection mod `3^k` exact for `k <= 7`; Chernoff rate `I(g) = 0.2513`; no-dip counts at `X = 3^12` grow like `X^0.57` (`theta = 0.03`) and `X^0.70` (`theta = 0.10`) against bounds `X^0.79`, `X^0.87` |
| (C3) parity-word map mod `2^k` | a bijection for `T_(+1)` and `T_(-1)` for all `k <= 12`; for the strategy `-chi_(-4)` only `k + 1` distinct words among `2^k` residues (`k <= 10`), the witness for the withdrawn corollary |
| (C4) the log-drift word `d_j = floor(j log_2 3 - 1.03 log_2 j)` (excluded by Corollary 9) | every prefix of length `k <= 18` is realised by an odd integer, the least one being `9, 41, 169, 681, 8873, ..., 5743273` for `k = 2, 4, 6, 8, 11, ..., 18`, growing like `2^(d_k)`; no odd `n < 2^24` realises the prefix of length 19 |
| the no-dip count to `X = 2^24` | `log F/log X = 0.7754, 0.7893, 0.8001` at `2^20, 2^22, 2^24` (`theta = 0.03`) and `0.8633, 0.8753, 0.8850` (`theta = 0.10`), both sheets within four elements; these are not slopes: the local slopes sit at the lemma exponents `0.9635` and `0.9867` (`collatz_thin_20260925_counts_to_2_24.out`) |

## 4. Frontier

1. The exponent: is `N(X) = O(X^(h*+eps))` sharp for some map in the
   family, or can a second constraint (the windows of one orbit are shifts
   of one sequence) lower it?
2. The union question: bound `#{n <= X : the orbit of n diverges}`. The
   branching of the inverse tree defeats the pigeonhole; Krasikov–Lagarias
   difference inequalities are the natural tool.
3. E-SCC: the choice relaxation changes the class count; whether a
   thinness statement survives for `E`-walks is untested.
