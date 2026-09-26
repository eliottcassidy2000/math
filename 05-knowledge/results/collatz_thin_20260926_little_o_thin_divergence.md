# Thin divergence is `o(X^(h*))`: a moving-barrier ballot bound (weighted Spitzer identity) turns THM-4476's `X^(h*+eps)` into `X^(h*) (log_2 X)^a` for every `a > lambda*/h* - 3/2 = -0.9862`

**Status: PROVED (elementary, self-contained modulo THM-4476's recursion
and the counting apparatus of THM-4495) + FINITE-EXACT controls +
INDEPENDENTLY AUDITED (SOUND, 2026-09-26; cosmetic notes applied). Session
`collatz-exponent-atlas-20260926` (opus), 2026-09-26. Supersedes addendum
1.6b of the thin-divergence note (`a > 0.0138`), which used the counting
lemma with a geometric tail but no ballot factor. Collatz orbits are still
not excluded; the statement is about the size of any orbit that does
diverge (and of any injective invariant set).**

Script: `04-computation/experiments/collatz_thin_20260926_movingbarrier.py`
-> `collatz_thin_20260926_movingbarrier.out`. Canon: THM-4499.

## 0. Statement

Let `b` be odd, `T_b(x) = x/2` (`x` even), `(3x+b)/2` (`x` odd), `alpha = log_2 3`,
`rho* = log_3 2`, `h* = h(rho*) = 0.9499555`, and
`lambda* = log_2(rho*/(1-rho*))/alpha = 0.488077`.

**Theorem.** For every `a > a* := lambda*/h* - 3/2 = -0.986211` there is
`K = K(a, |b|)` such that every `T_b`-orbit `(x_i)` in `Z \ {0}` that is not
eventually periodic satisfies

```text
N(X) = #{ i : |x_i| <= X } <= K X^(h*) (log_2 X)^a     for all X >= 2.
```

In particular `N(X) = o(X^(h*))`: the Terras exponent `h*` of THM-4476 is
not attained. The same holds for every `T_b`-invariant set on which `T_b`
is injective (THM-4476, Corollary 4), with the same proof.

The two ingredients are THM-4476's recursion (R) (section 1.5 of its note)
and a new bound for the no-dip set `F_b(X, theta)` at a barrier that
moves with `X`:

**Lemma M (moving-barrier ballot bound).** Let `W_k` be the number of
positive words of length `k` (THM-4495), and for `y >= 0` let

```text
M_k(y) = #{ w in {0,1}^k : S_i(w) > -y for all 1 <= i <= k },      S_i = o_i alpha - i.
```

Then for every `s > 0` there is `D_s` such that for all `k >= 1`, `y >= 0`,

```text
M_k(y) <= D_s 2^(hk) k^(-3/2) 2^(lambda* y) e^(s y).
```

(The control suggests `M_k(y)` of order `(y+1) 2^(hk) k^(-3/2) 2^(lambda* y)` for
`y << k^(1/2)`; that linear law is an observation, not a claim, and the
lemma trades the factor `y+1` for `e^(sy)` to stay elementary.)

**Lemma 1.4c (counting lemma with the ballot factor).** For
`theta in (0, theta_1]`, `theta_1 = 0.10376`, `X >= 2^200` and every `s > 0`,

```text
#F_b(X, theta) <= 2|b| X^(log_2(3/2) + theta) + 4 D_s X^(h*) (log_2 X)^(-3/2) 2^(lambda* theta log_2 X) e^(s theta log_2 X).
```

## 1. Proof of Lemma M

### 1.1 The decomposition at the minimum

Fix `k >= 1` and `y >= 0`. Since `alpha` is irrational, no two partial sums
of a word coincide and none vanishes (THM-4495, section 1), so minima are
attained at unique indices; the barrier condition `S_i > -y` is strict and
needs no tie-breaking.

Let `w` satisfy `S_i(w) > -y` for all `i >= 1`, and let `m in [1, k]` be the
index of the minimum of `S_1, ..., S_k` (unique). Then `a = w_1..w_m` ends at
its strict minimum among its partial sums (`S_j(a) > S_m(a)` for `j < m`)
with `S_m(a) > -y`, and `b = w_(m+1)..w_k` is positive (its partial sums are
`S_(m+j)(w) - S_m(w) > 0`); conversely such a pair `(a, b)` glues to a word
counted by `M_k(y)` with minimum at `|a|`. If `S_m(a) > 0`, `a` is a positive
min-ending word (THM-4495, Step 1: these number `L_m`, the first-passage
words, and `sum_m L_m W_(k-m) = W_k`). If `S_m(a) < 0`, reverse: `a^R` has
`S_i(a^R) = S_m(a) - S_(m-i)(a) < 0` for `i < m` and `S_m(a^R) = S_m(a) in (-y, 0)`,
a *negative word* (all partial sums `< 0`) with total `> -y`, and back.
Writing `N_m(y)` for the number of negative words of length `m` with total
`> -y`,

```text
M_k(y) = W_k + sum_(m=1)^(k) N_m(y) W_(k-m).          (D)
```

(D) is checked exactly in the control for `k <= 60`, `y in {1/2, 1, 2, 4, 8}`.

### 1.2 The tilt and the exponential weight

Let `Q` be the product measure on words with `Q(1) = rho*`, `Q(0) = 1 - rho*`
(zero drift: `rho* alpha - 1 = 0`). As in THM-4495,
`Q(w) = 2^(-hm) 2^(lambda* S_m(w))` for every word `w` of length `m`. Hence,
for a negative word with total `> -y`, `Q(w) >= 2^(-hm) 2^(-lambda* y)`, and

```text
N_m(y) <= 2^(hm) 2^(lambda* y) Q( S_i < 0 for 1 <= i <= m,  S_m > -y )
       <= 2^(hm) 2^(lambda* y) e^(sy) E_Q[ e^(s S_m) ; S_i < 0 for 1 <= i <= m ]
```

for every `s > 0` (on the event, `e^(s S_m) >= e^(-sy)`).

### 1.3 The weighted Spitzer identity

For letter weights `x_0, x_1 > 0` and a word `w` put `x(w) = x_0^(#0) x_1^(#1)`.
Let `P_n^-(x)` be the total weight of the negative words of length `n` and
`B_n^-(x)` the total weight of the words of length `n` with negative total.
Then

```text
n P_n^-(x) = sum_(m=1)^(n) B_m^-(x) P_(n-m)^-(x),     P_0^- = 1,
```

i.e. `sum_n P_n^-(x) t^n = exp(sum_n B_n^-(x) t^n/n)`. *Proof.* THM-4495's
proof of Theorem A (minimum decomposition, reversal, ladder blocks,
rotation averaging) is a chain of bijections and a rotation argument; each
bijection preserves the multiset of letters, hence `x(w)`, and the rotations
of a word all carry the same weight, so every step holds with the weight
`x(w)` in place of `1`. Negating the sign convention (negative words,
negative totals, records replaced by strict minima) is the same proof with
the walk `-S`. ∎ The identity is checked exactly with rational weights
`(1,1)`, `(1/2, 2)`, `(3, 1/3)`, both signs, for `n <= 40`.

Apply it with `x_1 = rho* e^(s(alpha - 1))`, `x_0 = (1 - rho*) e^(-s)`: then
`x(w) = Q(w) e^(s S_n(w))`, so `P_n^-(x) = E_Q[e^(sS_n); S_i < 0 for i <= n] =: p_n`
and `B_n^-(x) = E_Q[e^(sS_n); S_n < 0] =: g_n`.

### 1.4 The local bound and the convolution

`g_n = sum_(o : o alpha < n) Bin(n, rho*)(o) e^(s(o alpha - n))`. The admissible
values `o alpha - n` are negative and spaced by `alpha`, so the exponential
factors are at most `1, e^(-s alpha), e^(-2 s alpha), ...`, and every
binomial probability is at most `C_0 n^(-1/2)` with `C_0 = 3.2` (the
binomial mode is within `1` of `rho* n`; Stirling for `n >= 10`, and
`Bin <= 1 <= 3.2 n^(-1/2)` for `n < 10`). Hence

```text
g_n <= C_s n^(-1/2),      C_s = 3.2/(1 - e^(-s alpha)).
```

THM-4495's convolution lemma (`(a*b)_k <= 2^(3/2)(A beta + B alpha) k^(-3/2)`
for nonnegative sequences bounded by `A n^(-3/2)`, `B n^(-3/2)` with sums
`alpha, beta`) applied to `exp(sum_n (g_n/n) t^n)` gives

```text
p_n <= C_s e^(2 sqrt2 sigma_s) n^(-3/2),      sigma_s = sum_n g_n/n <= C_s zeta(3/2) = 2.613 C_s.
```

### 1.5 Assembly

With `D'_s := C_s e^(2 sqrt2 sigma_s)`: `N_m(y) <= D'_s 2^(hm) 2^(lambda* y) e^(sy) m^(-3/2)`.
Insert in (D) with `W_(k-m) <= c_2 2^(h(k-m)) (k-m)^(-3/2)` for `k - m >= 1`
(THM-4495, `c_2 = 545`) and `W_0 = 1`:

```text
M_k(y) <= c_2 2^(hk) k^(-3/2) + D'_s 2^(lambda* y) e^(sy) 2^(hk) [ k^(-3/2) + c_2 sum_(m=1)^(k-1) m^(-3/2) (k-m)^(-3/2) ]
       <= 2^(hk) k^(-3/2) 2^(lambda* y) e^(sy) [ c_2 + D'_s (1 + 2^(3/2) * 2 * 2.613 c_2) ],
```

by the convolution lemma once more (`A = B = 1`, `alpha = beta = zeta(3/2)`).
So Lemma M holds with `D_s = c_2 + D'_s (1 + 14.8 c_2)`. ∎

## 2. Proof of Lemma 1.4c

As in THM-4476's section 1.4: elements `y_0 < Y_0 = 2|b| X^(log_2(3/2) + theta)`
of `F_b(X, theta)` number at most `Y_0`; an element `n >= Y_0` has, with
`k = floor(log_2 X)`, `3^(o_i)/2^i >= X^(-theta)/2` for every `i <= k`
(the carry bound `|beta_i| <= |b|((3/2)^i - 1)` is at most half of
`n X^(-theta)` for every `i <= k`), hence `S_i > -theta log_2 X - 1` for all
`1 <= i <= k`, strictly (the `-1` in the carry bound). The hypothesis
`theta <= theta_1` is inherited from Lemma 1.4b's form of the small-element
term and is not needed by Lemma M. So the word of `n` is counted by `M_k(y)` with
`y = theta log_2 X + 1`, each word is one class modulo `2^k` with at most
two representatives in `[1, X]`, and Lemma M gives
`2 M_k(y) <= 2 D_s 2^(hk) k^(-3/2) 2^(lambda*(theta log_2 X + 1)) e^(s(theta log_2 X + 1))`;
with `2^(hk) <= X^(h*)`, `k >= log_2 X - 1 >= 0.995 log_2 X` (`X >= 2^200`)
and `2^(lambda*) e^(s) 2 * 0.995^(-3/2) <= 4` for `s <= 0.1`, the lemma follows
(for `s > 0.1` enlarge `D_s`). ∎

## 3. Proof of the theorem

Take one positive one-signed segment of the orbit (THM-4476, 1.1; the sum
over at most `|b|/3 + 2` segments multiplies `K`). THM-4476's recursion (R)
of its section 1.5, with Lemma 1.4c in place of the counting lemma, reads:
for `theta in (0, theta_1]`, `X >= 2^200`, `L = log_2 X`, `k = floor(L)` and
`s in (0, 0.1]`,

```text
N(X) <= k + k N(X^(1-theta)) + 2|b| X^(0.585 + theta) + 4 D_s X^(h*) L^(-3/2) 2^(lambda* theta L) e^(s theta L).     (R'')
```

Fix `a > a* = lambda*/h* - 3/2` and put `eta = (a - a*)/8`, `c_1 = (1 + eta)/h*`,
`theta_X = c_1 (log_2 L)/L`, `s = min(0.1, eta ln 2/c_1)`. Then
`c_1 h* = 1 + eta`, `X^(1 - theta_X) = X L^(-c_1)`, `2^(lambda* theta_X L) = L^(lambda* c_1)`,
`e^(s theta_X L) = L^(s c_1/ln 2) <= L^(eta)`, and

```text
lambda* c_1 - 3/2 + eta = a* + eta lambda*/h* + eta <= a* + 1.52 eta = a - 6.48 eta <= a - 2 eta.
```

Choose `X_2 = X_2(eta, |b|)` such that for `X >= X_2`: `theta_X <= theta_1`,
`X >= 2^200`, `c_1 log_2 L <= L/2`, and `k + 2|b| X^(0.585 + theta_X) <= X^(h*) L^(a - 2 eta)`
(possible: the left side is `O(X^(0.7))`). For `X >= X_2`, (R'') gives

```text
N(X) <= X^(h*) L^(a - 2 eta) + L N(X L^(-c_1)) + 4 D_s X^(h*) L^(a - 2 eta).     (R''')
```

Let `X_3 >= X_2` be such that `2^(|a|) L^(-eta) <= 1/4` and `L^(-2 eta) <= 1/4`
for `X >= X_3`, and put

```text
K = max( 8 D_s + 1,  sup_(2 <= X < X_3) X^(1 - h*) (log_2 X)^(-a) ),
```

a finite number (the supremum is over a bounded range). **Claim:**
`N(X) <= K X^(h*) L^a` for all `X >= 2`. *Base.* For `2 <= X < X_3`,
`N(X) <= X = X^(1-h*) L^(-a) * X^(h*) L^a <= K X^(h*) L^a`. *Step.* Let
`X >= X_3` and assume the claim for every `Y in [2, X)`. Put `Y = X L^(-c_1)`;
since `c_1 log_2 L <= L/2`, `L^(c_1) <= X^(1/2)`, so `Y in [X^(1/2), X)` and
`log_2 Y in [L/2, L]`; hence `(log_2 Y)^a <= 2^(|a|) L^a` whatever the sign of
`a`, and `Y^(h*) = X^(h*) L^(-c_1 h*) = X^(h*) L^(-1-eta)`. Therefore

```text
L N(Y) <= L K Y^(h*) (log_2 Y)^a <= 2^(|a|) K X^(h*) L^(a - eta),
N(X)  <= X^(h*) L^a [ L^(-2 eta) + 2^(|a|) K L^(-eta) + 4 D_s L^(-2 eta) ]
      <= X^(h*) L^a [ (1 + 4 D_s)/4 + K/4 ]  <=  K X^(h*) L^a,
```

the last step because `K >= 8 D_s + 1 >= (1 + 4 D_s)/3`. The induction is
on `floor(X)`: `N` depends only on `floor(X)`, and `Y = X L^(-c_1) <= X/264`
for `X >= 2^200`, so `floor(Y) < floor(X)`; the base covers `[2, X_3)`. ∎

For an injective invariant set (Corollary 4 of THM-4476) the recursion of
its section 1.7 is `N_A(X) <= k N_A(X^(1-theta)) + #F_b + #F_(-b) + (k+1)(|b|/3+1)`
with base `N_A(X) <= 2X`; Lemma 1.4c bounds both `#F` terms, so the same
bootstrap applies with `4 D_s` replaced by `8 D_s` and the `O(k)` absorbed
in `X^(h*)`. The constant `K` is not effective (`X_3` is of the order
`2^(2^(1/eta))`); the theorem is asymptotic in `X` for each fixed `a`.

## 3b. Corollary (the discrepancy form, the owner's original question)

Let `(m_l)` be the odd iterates of a divergent positive 3n+1 orbit, `d_l`
the cumulative halving count, and `Delta_l = d_l - l log_2 3` the
discrepancy of THM-4476 / HYP-9160 (so `m_l = m_0 2^(-Delta_l) prod(1 + 1/(3m_j))`,
Proposition 6 of the squares/doubles note). Then for every `a > a*` and
all large `L`,

```text
min_(l <= L) Delta_l  <=  -(1/h*) log_2 L  +  (a/h*) log_2 log_2 L  +  O(1),
```

i.e. for every `c < -a*/h* = 1.038` (the coefficient is `-(a*+eps)/h* = 1.038 - eps/h*`, not `1.038` itself; the audit caught the display):
`min_(l <= L) Delta_l <= -1.0527 log_2 L - c log_2 log_2 L + O(1)`.

*Proof.* Put `X = max_(l <= L) m_l`. The `L + 1` odd iterates `m_0..m_L` are
distinct orbit elements `<= X`, so `L + 1 <= N(X) <= K X^(h*) (log_2 X)^a`.
Since `log_2 X = log_2 m_0 - min_(l<=L) Delta_l + O(1)` (the carry product is
bounded, THM-4476 Cor. 6), solving for `-min Delta_l` gives the claim. ∎

THM-4476's Corollary 3 excluded every log-band `m_j <= K j^a` with
`a < 1/h*`; this adds the `log log` term and, through THM-4499, replaces
"for some `eps`" by an explicit second-order coefficient. Nothing here
bears on whether divergent orbits exist.

## 4. Remarks

* **Where the exponent comes from.** `-3/2` is the ballot factor of a
  zero-drift walk pinned at a bounded height (THM-4495); `lambda*/h*` is
  the tilt price `2^(lambda* theta L)` of a barrier at depth `theta L`,
  with `theta L = c_1 log_2 L` forced by the `L` dippers per landing point
  in the recursion (`c_1 h* > 1`). The lemma's `e^(sy)` is a technical
  slack: the true dependence on the depth is linear, which would give the
  same exponent `a*` with `s = 0`.
* **What is and is not gained.** THM-4476 gave `X^(h*+eps)`; addendum 1.6b
  gave `X^(h*) (log X)^(0.014+eps)`; this note gives `X^(h*) (log X)^(-0.986+eps)`,
  so any divergent orbit has `N(X) = o(X^(h*))`. The expected truth for an
  actual divergent orbit is `N(X) = O(log X)`; the method still uses one
  free window per element (THM-4487, Theorem 4, remark (iii)), and the
  landing multiplicity `L` is the next constant to attack: a multiplicity
  `L^(1/2)` would give `a* = lambda*/(2h*) - 3/2`.
* **Constants.** `D_s` grows like `exp(c/s)` as `s -> 0`; the theorem is
  asymptotic in `X` for each fixed `a`, with `K` depending on `a` and `|b|`.
* **Controls.** The weighted identity, the decomposition (D) and the
  behaviour of `M_k(y)` are checked exactly in the script:

```text
h = 0.9499555, lambda* = 0.488077, lambda*/h - 3/2 = -0.986211
(A) weighted Spitzer identity n P_n = sum_m B_m P_(n-m), exact rationals, n <= 40:
    weights (x0, x1) = (1, 1), positive words: True; P_1..P_7 = ['1', '1', '2', '3', '4', '8', '13']
    weights (x0, x1) = (1, 1), negative words: True; P_1..P_7 = ['1', '2', '3', '6', '12', '22', '44']
    weights (x0, x1) = (1/2, 2), positive words: True; P_1..P_7 = ['2', '4', '10', '24', '56', '140', '344']
    weights (x0, x1) = (1/2, 2), negative words: True; P_1..P_7 = ['1/2', '5/4', '9/8', '45/16', '225/32', '613/64', '3065/128']
    weights (x0, x1) = (3, 1/3), positive words: True; P_1..P_7 = ['1/3', '1/9', '10/27', '19/81', '28/243', '280/729', '613/2187']
    weights (x0, x1) = (3, 1/3), negative words: True; P_1..P_7 = ['3', '10', '33', '110', '1100/3', '1222', '12220/3']
(B) moving barrier: M_k(y) / (2^(hk) k^(-3/2) 2^(lambda* y)) and the same divided by (y+1)
    y:             0        1        2        4        8       16       32
    k=  100:   7.661   15.586   28.600   48.143   59.581   32.208    0.580  | /(y+1):  7.661   7.793   9.533   9.629   6.620   1.895   0.018
    k=  200:   8.919   18.103   33.396   57.541   88.287   80.363    9.730  | /(y+1):  8.919   9.051  11.132  11.508   9.810   4.727   0.295
    k=  400:   9.840   19.566   34.837   66.560  110.165  143.160   60.588  | /(y+1):  9.840   9.783  11.612  13.312  12.241   8.421   1.836
    k=  800:  10.446   19.841   38.363   68.667  125.813  198.499  179.219  | /(y+1): 10.446   9.921  12.788  13.733  13.979  11.676   5.431
    k= 1200:  10.255   21.046   39.493   72.794  126.674  213.808  260.569  | /(y+1): 10.255  10.523  13.164  14.559  14.075  12.577   7.896
    k= 1600:  10.757   21.422   38.467   73.888  134.681  235.768  303.151  | /(y+1): 10.757  10.711  12.822  14.778  14.965  13.869   9.186
(C) decomposition M_k(y) = W_k + sum_m N_m(y) W_(k-m), exact, k <= 60:
    identity holds for all tested y and k: True
```
