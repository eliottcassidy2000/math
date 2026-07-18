---
id: THM-1092
title: Absolute k-fold resonance-lattice identity, with a product-Fejer derivation
status: PROVED ANALYTIC — every full-support reciprocal resonance series is absolutely convergent by a uniform dyadic-shell majorant; hence the raw k-fold lattice sum is order-independent and equals both the danger intersection and its product-Fejer/L1 limit. Every support term is an ordinary centered joint moment. No useful higher-support bound or LRC14 closure is claimed.
source: codex-2026-07-18-S67 global-bridge audit
depends_on: []
related: [THM-965, THM-1075, THM-1070, THM-1080, THM-1085, THM-1091, THM-1093]
---

# THM-1092 — the raw k-fold resonance sum converges absolutely

THM-1075 found the correct lattice but did not prove convergence of its raw
series.  Two arguments now meet.  Fejer regularization turns the formula into
a limit of finite identities and gives every support term a separate `L1`
limit.  A dyadic-shell estimate then proves that every full-support resonance
series is actually absolutely convergent.  Thus the regularized limit equals
the raw, order-independent lattice sum.

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

Then the series below converges absolutely and

```text
mu(intersection_i {t : a_i t in D_lambda})
 = sum_(n in Lambda(a)) product_i hhat(n_i).                 (2)
```

It is also the product-Fejer limit

```text
mu(intersection_i {t : a_i t in D_lambda})
 = lim_(N->infinity)
   sum_(n in Lambda(a)) product_i [w_N(n_i) hhat(n_i)].       (3)
```

The sum at level `N` is finite because `w_N(m)=0` for `|m|>N`.  Formula (3)
is a canonical product-Fejer, or rectangular `(C,1)`, derivation of the raw
identity (2), not a convention needed to define it.

There is also a termwise support expansion.  For `S subseteq {1,...,k}`, set

```text
delta_N(S) =
  sum product_(i in S) [w_N(n_i) hhat(n_i)],                 (4)
```

where the sum is over `n_i != 0`, `|n_i|<=N`, and
`sum_(i in S) a_i n_i=0`.  Define `delta_N(empty)=1`, and define

```text
delta(S) =
  sum product_(i in S) hhat(n_i)                             (5)
```

where the raw sum has the same full-support resonance conditions.  It is
absolutely convergent, `delta(empty)=1`, `delta({i})=0`, and

```text
mu(intersection_i {t : a_i t in D_lambda})
 = sum_(S subseteq [k]) (2 lambda)^(k-|S|) delta(S).         (6)
```

Moreover

```text
delta(S) = lim_(N->infinity) delta_N(S)
 = integral_T product_(i in S) (h(a_i t)-2 lambda) dmu(t).   (7)
```

Thus the `|S|>=3` terms have three equivalent meanings: an absolutely
convergent lattice sum, a canonical finite-polynomial limit, and an ordinary
centered joint moment.

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
||h_N-h||_1 -> 0,                                            (8)
```

and the finite Fourier expansion

```text
h_N(t) = sum_(|m|<=N) w_N(m) hhat(m) exp(2 pi i m t).         (9)
```

For every nonzero integer `a_i`, multiplication by `a_i` preserves Haar
measure on `T`.  Consequently

```text
||h_N(a_i .)-h(a_i .)||_1 = ||h_N-h||_1.                    (10)
```

All factors in the following products lie in `[0,1]`.  The elementary
telescoping inequality

```text
|product_i x_i - product_i y_i| <= sum_i |x_i-y_i|           (11)
```

therefore gives

```text
|| product_i h_N(a_i .) - product_i h(a_i .) ||_1
 <= k ||h_N-h||_1 -> 0.                                     (12)
```

The integral of the limiting product is exactly the left side of (2).

## 3. Finite resonance identity

Substitute (9) into the finite product.  Since every sum is finite, it may be
expanded and integrated without a convergence question:

```text
integral_T product_i h_N(a_i t) dmu(t)
 = sum_(|n_i|<=N) product_i [w_N(n_i) hhat(n_i)]
     * integral_T exp(2 pi i t sum_i a_i n_i) dmu(t).
```

Character orthogonality makes the last integral one when
`sum_i a_i n_i=0` and zero otherwise.  This is precisely the level-`N` sum
on the right of (3).  Taking the limit and using (12) proves the
product-Fejer identity (3).

For the support expansion, write

```text
g_N = h_N-2 lambda,       g = h-2 lambda.
```

The zero Fourier coefficient of `g_N` is zero, while every nonzero
coefficient is `w_N(m)hhat(m)`.  Repeating the finite character calculation
shows

```text
integral_T product_(i in S) g_N(a_i t) dmu(t) = delta_N(S).  (13)
```

Here `|g_N|,|g|<=1`, so the same telescoping argument and (10) show that the
left side of (13) converges to the integral in (7).  This proves the limit
and integral equality in (7).
Finally, expanding

```text
product_i (2 lambda + g(a_i t))
```

inside the integral gives the finite subset sum (6).

## 4. Absolute convergence of the raw sums

It remains to identify the finite-polynomial limit with the raw series.  Fix
a support `S` of size `s>=2` and consider

```text
R_S = sum_{n_i != 0, sum_(i in S) a_i n_i=0}
        product_(i in S) 1/|n_i|.
```

For `r>=0`, put the vectors with
`2^r <= max_i |n_i| < 2^(r+1)` in a shell `A_r`, and cover `A_r` by the
`s` sets on which a chosen coordinate `j` realizes the maximum.  On the
`j`th set, `1/|n_j| <= 2^(-r)`, while every other coordinate is a nonzero
integer of absolute value below `2^(r+1)`.  Once those other coordinates
are fixed, the equation

```text
a_j n_j = -sum_(i != j) a_i n_i
```

determines at most one integer `n_j`, because `a_j != 0`.  Hence, writing
`H_M=sum_(m=1)^M 1/m`,

```text
sum_(n in A_r) product_i 1/|n_i|
 <= s * 2^(-r) * (2 H_(2^(r+1)-1))^(s-1).                  (14)
```

Since `H_M <= 1+log M`, the sum of the right side over `r` is bounded by a
constant times

```text
sum_(r>=0) 2^(-r) (r+1)^(s-1) < infinity.                  (15)
```

Thus `R_S` converges.  The estimate is uniform in the nonzero coefficients
`a_i`; it is deliberately crude, but it is a proof rather than the shell
heuristic previously recorded in THM-1085.

For nonzero `n`, equation (1) gives
`|hhat(n)| <= 1/(pi|n|)`.  Therefore every full-support series (5) converges
absolutely.  Splitting the full lattice by its finitely many coordinate
supports proves absolute convergence in (2).  Finally `0<=w_N<=1` and
`w_N(m)->1` for every fixed `m`, so dominated convergence for series turns
(3) into (2) and turns the limit in (7) into the raw sum (5).

## 5. Consequences at the LRC threshold

At `lambda=1/14`, equation (1) becomes

```text
hhat(n) = sin(pi n/7)/(pi n),
```

so every nonzero coefficient with `7|n` vanishes.  This zero is retained at
every finite Fejer level.  The pair lattice is cyclic, recovering THM-965 and
THM-1075's exact pair-independence criterion.  For `|S|>=3`, (7) identifies
the higher term as an ordinary centered joint moment; THM-1080/1085's short-
vector and signed-cancellation questions concern its size, not its existence.
Absolute convergence also makes THM-1093's finite partition into residue
classes modulo `14` legitimate without choosing an order of summation.

Common dilation is exact before taking a limit:

```text
sum_i (q a_i)n_i=0  iff  sum_i a_i n_i=0                     (16)
```

for every nonzero integer `q`.  Hence every finite regularized sum, every
`delta_N(S)`, and every limit above is dilation-invariant.

The finite-sheet identities in THM-1091 are the discrete analogue of the
same calculation: character orthogonality selects a zero-sum resonance
lattice, while ramification restricts the allowed Fourier support to
annihilator subgroups.

## 6. Honest boundary

This theorem repairs the analytic status of the k-fold identity; it does not
repair the quantitative certificate.  In particular it proves none of the
following:

- a useful coefficient-sensitive bound for `delta(S)` when `|S|>=3`;
- a closed-form evaluation of the higher-dimensional sums;
- positivity of BONF5/BONF7 on the dense residual; or
- the n=12 rigidity theorem or LRC(14).

The gain is a clean separation: regularization, absolute convergence, and
reordering are now closed.  The remaining problem is genuinely quantitative
cancellation/geometry in the explicit resonance lattice.
