# Thin divergence: every non-periodic orbit of x -> x/2, (3x+b)/2 has at most X^(0.95+eps) elements below X; reciprocal sums converge (HYP-9160 proved); the real value of the Bernstein series is strictly below the integer

**Status: PROVED (elementary; full proof below; self-audited, independent
audit recorded in the theorem file THM-4476) with FINITE-EXACT controls of
the two counting lemmas. No priority claim: the ingredients are Terras's
stopping-time count and an injectivity pigeonhole; the single-orbit
statement may exist in the literature under another name. Collatz, the
`3n-1` sheet and the Periodicity Conjecture remain OPEN: the theorem
constrains the shape of a hypothetical divergent orbit, it does not exclude
one. Session `collatz-squares-doubles-20260925` (opus), 2026-09-25.**

Theorem file: [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
Controls: `04-computation/experiments/collatz_thin_20260925_counts.py`,
output `collatz_thin_20260925_counts.out`.
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

**Corollary 4 (sign strategies).** The theorem holds verbatim for every
level-`k` sign strategy `T_sigma` of THM-4474 and for every map
`x -> x/2, (3x + b(x))/2` with `b(x)` odd and depending only on
`x mod 2^(k_0)`.

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
`X <= 2^20` (`(L2)`): at `theta = 0.03` their size grows like `X^0.78` on
both sheets, well below the bound `X^0.9635`.

### 1.5 The dichotomy and the recursion

Fix `theta in (0, theta_0)`, `X >= 2^(k_0)`, `k = floor(log_2 X)`. Every
index `i <= M` with `y_i <= X` is of one of three kinds:

* **(E)** `i > M - k` (only when `M < infinity`): fewer than `k` indices;
* **(D)** `i <= M - k` and `y_(i+s) < y_i X^(-theta)` for some `1 <= s <= k`;
* **(ND)** `i <= M - k` and `y_(i+s) >= y_i X^(-theta)` for all `0 <= s <= k`.

**(D) is a pigeonhole.** For a (D)-index take the least such `s`. Then
`y_(i+s) < X^(1-theta)`, and the pair `(i+s, s)` determines `i`. So
`#(D) <= k N(X^(1-theta))`: the dipping points land on the segment's own
points below `X^(1-theta)`, and each landing point serves at most `k`
dippers. This is where distinctness of the orbit is used.

**(ND) is the counting lemma.** An (ND)-index has `y_i in F_b(X, theta)`,
and distinct indices have distinct `y_i`. So `#(ND) <= #F_b(X, theta)`.

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
`sum 1/|m_l| < infinity`, and no term is `1/3` (that orbit reaches `0`), so
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

*Corollary 4.* Sections 1.2 and 1.3 only use that the odd step is
`(3x + b(x))/2` with `b(x)` odd and determined by `x mod 2^(k_0)`: the word
of length `k >= k_0` is then determined by `x mod 2^k` and the parity map is
still a bijection, and `|beta_i| <= max|b| ((3/2)^i - 1)`. ∎

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
| `(L2)` `#F_(+-1)(2^20, 0.03)` | `46612` / `46611`, growth exponent `0.775` (bound `0.9635`) |
| `(L2)` `#F_(+-1)(2^20, 0.10)` | `157652` / `157648`, growth exponent `0.863` (bound `0.9867`) |

The two sheets differ by at most four elements at every size, as the
sheet-blind class count predicts.

## 4. Frontier

1. The exponent: is `N(X) = O(X^(h*+eps))` sharp for some map in the
   family, or can a second constraint (the windows of one orbit are shifts
   of one sequence) lower it?
2. The union question: bound `#{n <= X : the orbit of n diverges}`. The
   branching of the inverse tree defeats the pigeonhole; Krasikov–Lagarias
   difference inequalities are the natural tool.
3. E-SCC: the choice relaxation changes the class count; whether a
   thinness statement survives for `E`-walks is untested.
