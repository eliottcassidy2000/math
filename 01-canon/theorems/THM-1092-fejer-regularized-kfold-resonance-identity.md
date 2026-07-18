---
id: THM-1092
title: Fejer-regularized k-fold resonance-lattice identity
status: PROVED ANALYTIC — the k-fold danger intersection is the canonical product-Fejer limit of finite resonance-lattice sums; every limit is justified in L1 and every support term has its own convergent centered-moment interpretation. No higher-cumulant bound or LRC14 closure is claimed.
source: codex-2026-07-18-S67 global-bridge audit
depends_on: [THM-965, THM-1075]
related: [THM-1070, THM-1080, THM-1085, THM-1091]
---

# THM-1092 — a rigorous meaning for the k-fold resonance sum

THM-1075 found the correct lattice but evaluated its conditionally ordered
series by symmetric truncation.  The identity does not need an unproved
ordering convention.  Fejer regularization turns it into a limit of finite
identities and gives every support term a separate `L1` limit.

## 1. Statement

Let `T=R/Z` have Haar probability measure `mu`.  Fix `0<lambda<1/2`, let

```text
D_lambda = {t in T : ||t|| < lambda},
h = 1_(D_lambda),
```

and use the Fourier convention

```text
hhat(n) = integral_T h(t) exp(-2 pi i n t) dmu(t).
```

Changing the two boundary points of `D_lambda` does not affect any integral.
For `n in Z`,

```text
hhat(0) = 2 lambda,
hhat(n) = sin(2 pi n lambda)/(pi n)       when n != 0.        (1)
```

Let `a=(a_1,...,a_k)` be nonzero integers and put

```text
Lambda(a) = {n in Z^k : sum_i a_i n_i = 0}.
w_N(m) = max(1-|m|/(N+1),0).
```

Then

```text
mu(intersection_i {t : a_i t in D_lambda})
 = lim_(N->infinity)
   sum_(n in Lambda(a)) product_i [w_N(n_i) hhat(n_i)].       (2)
```

The sum at level `N` is finite because `w_N(m)=0` for `|m|>N`.  Thus (2)
is not an assertion about an unordered conditionally convergent series.  It
is a canonical product-Fejer, or rectangular `(C,1)`, limit.

There is also a termwise support expansion.  For `S subseteq {1,...,k}`, set

```text
delta_N(S) =
  sum product_(i in S) [w_N(n_i) hhat(n_i)],                 (3)
```

where the sum is over `n_i != 0`, `|n_i|<=N`, and
`sum_(i in S) a_i n_i=0`.  Define `delta_N(empty)=1`.  Every limit

```text
delta_F(S) = lim_(N->infinity) delta_N(S)                    (4)
```

exists, `delta_F(empty)=1`, `delta_F({i})=0`, and

```text
mu(intersection_i {t : a_i t in D_lambda})
 = sum_(S subseteq [k]) (2 lambda)^(k-|S|) delta_F(S).       (5)
```

Moreover the regularized support term has the intrinsic formula

```text
delta_F(S)
 = integral_T product_(i in S) (h(a_i t)-2 lambda) dmu(t).   (6)
```

Equations (4) and (6) make the `|S|>=3` terms rigorous without claiming that
their raw lattice sums converge under every ordering.

## 2. Fejer approximation

Let `F_N` be the normalized Fejer kernel and put

```text
h_N = h * F_N.
```

The standard approximate-identity theorem (applied to the `L1` function
`h`) gives

```text
0 <= h_N <= 1,
integral h_N = 2 lambda,
||h_N-h||_1 -> 0,                                            (7)
```

and the finite Fourier expansion

```text
h_N(t) = sum_(|m|<=N) w_N(m) hhat(m) exp(2 pi i m t).         (8)
```

For every nonzero integer `a_i`, multiplication by `a_i` preserves Haar
measure on `T`.  Consequently

```text
||h_N(a_i .)-h(a_i .)||_1 = ||h_N-h||_1.                     (9)
```

All factors in the following products lie in `[0,1]`.  The elementary
telescoping inequality

```text
|product_i x_i - product_i y_i| <= sum_i |x_i-y_i|           (10)
```

therefore gives

```text
|| product_i h_N(a_i .) - product_i h(a_i .) ||_1
 <= k ||h_N-h||_1 -> 0.                                     (11)
```

The integral of the limiting product is exactly the left side of (2).

## 3. Finite resonance identity

Substitute (8) into the finite product.  Since every sum is finite, it may be
expanded and integrated without a convergence question:

```text
integral_T product_i h_N(a_i t) dmu(t)
 = sum_(|n_i|<=N) product_i [w_N(n_i) hhat(n_i)]
     * integral_T exp(2 pi i t sum_i a_i n_i) dmu(t).
```

Character orthogonality makes the last integral one when
`sum_i a_i n_i=0` and zero otherwise.  This is precisely the level-`N` sum
on the right of (2).  Taking the limit and using (11) proves (2).

For the support expansion, write

```text
g_N = h_N-2 lambda,       g = h-2 lambda.
```

The zero Fourier coefficient of `g_N` is zero, while every nonzero
coefficient is `w_N(m)hhat(m)`.  Repeating the finite character calculation
shows

```text
integral_T product_(i in S) g_N(a_i t) dmu(t) = delta_N(S).  (12)
```

Here `|g_N|,|g|<=1`, so the same telescoping argument and (9) show that the
left side of (12) converges to the integral in (6).  This proves (4) and (6).
Finally, expanding

```text
product_i (2 lambda + g(a_i t))
```

inside the integral gives the finite subset sum (5).

## 4. Consequences at the LRC threshold

At `lambda=1/14`, equation (1) becomes

```text
hhat(n) = sin(pi n/7)/(pi n),
```

so every nonzero coefficient with `7|n` vanishes.  This zero is retained at
every finite Fejer level.  The pair lattice is cyclic, recovering THM-965 and
THM-1075's exact pair-independence criterion.  For `|S|>=3`, (6) identifies
the higher term as an ordinary centered joint moment; THM-1080/1085's short-
vector and signed-cancellation questions concern its size, not its existence.

Common dilation is exact before taking a limit:

```text
sum_i (q a_i)n_i=0  iff  sum_i a_i n_i=0                     (13)
```

for every nonzero integer `q`.  Hence every finite regularized sum, every
`delta_N(S)`, and every limit above is dilation-invariant.

The finite-sheet identities in THM-1091 are the discrete analogue of the
same calculation: character orthogonality selects a zero-sum resonance
lattice, while ramification restricts the allowed Fourier support to
annihilator subgroups.

## 5. Honest boundary

This theorem repairs the analytic status of the k-fold identity; it does not
repair the quantitative certificate.  In particular it proves none of the
following:

- a useful uniform bound for `delta_F(S)` when `|S|>=3`;
- absolute convergence of the raw, unregularized lattice series in the exact
  ordering used by a separate computation;
- positivity of BONF5/BONF7 on the dense residual; or
- the n=12 rigidity theorem or LRC(14).

The gain is a clean separation: regularization and convergence are now closed,
and the remaining problem is genuinely cancellation/geometry in the explicit
resonance lattice.
