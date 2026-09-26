# The no-descent count `W_k = |Bad_k|` has exact order `2^(hk) k^(-3/2)`: a self-contained Spitzer identity `k W_k = sum_n B_n W_(k-n)` and an elementary convolution bound; the strategy-cube distance, the periodic-deletion price, the expanding-necklace count and the `gamma = 1` dip count are all `Theta(2^(hk) k^(-3/2))`

**Status: PROVED (self-contained, elementary) + FINITE-EXACT (identity
checked against a ballot DP to `k = 300`, integrality of the recurrence
and all bounds to `k = 3000`). Session `collatz-exponent-atlas-20260926`
(opus), 2026-09-26. Integrates the procgen lane's THM-4479 / THM-4485 (the
exponent `1 - h` of the distance to provability, sandwiched between a
necklace count and `|Bad_k|`) with THM-4487 (the dip spectrum, whose
`gamma = 1` count is the same object): one exponent and now one polynomial
factor in four settings. Collatz itself is untouched.**

Script: `04-computation/experiments/collatz_nodescent_order_20260926.py`
-> `collatz_nodescent_order_20260926.out`.

## 0. Statements

Throughout, `T` is the Collatz shortcut (`n/2`, `(3n+1)/2`), `h = h(log_3 2) = 0.9499555`,
`rho = log_3 2`, and a residue `r mod 2^k` has the parity word
`w = (w_1, ..., w_k)`, `w_j = T^(j-1)(r) mod 2`, with `o_j = w_1 + ... + w_j`
odd letters among the first `j` and multiplier `M_j = 3^(o_j)/2^j`
(Terras: `r mod 2^k -> w` is a bijection). Put

```text
W_k = #{ w in {0,1}^k : 3^(o_j) > 2^j for every 1 <= j <= k }      (= |Bad_k| of THM-4479),
B_n = #{ w in {0,1}^n : 3^(o_n) > 2^n } = sum_(j : 3^j > 2^n) C(n, j),
N_k = # binary necklaces of length k with more than k log_3 2 ones.
```

`W_k` counts the residues mod `2^k` with no `k`-step descent (`M_j > 1` for
all `j <= k`), i.e. the residues on which every large `n` fails to fall
below itself within `k` steps. In the language of THM-4487 it is the number
of length-`k` words whose walk `S_j = o_j log_2 3 - j` stays strictly
positive; `B_n` is the number of words with positive total.

**Theorem A (Spitzer's identity, self-contained).** For every `k >= 1`,

```text
k W_k = sum_(n=1)^(k) B_n W_(k-n),      W_0 = 1,
```

equivalently `sum_k W_k t^k = exp( sum_n B_n t^n / n )`. So the
no-descent counts are determined by the binomial tails `B_n`, with no
enumeration: `W_1..W_20 = 1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328`
(the values `|Bad_k|`, `k <= 11`, of THM-4479's table, and its ballot DP).

**Theorem B (exact order).** For every `k >= 1`,

```text
0.26 * 2^(hk) k^(-3/2)  <=  B_k / k  <=  W_k  <=  545 * 2^(hk) k^(-3/2),
```

and `N_k >= C(k, floor(k rho) + 1)/k >= 0.26 * 2^(hk) k^(-3/2)`. The
normalised counts `B_k k^(1/2) 2^(-hk)` and `N_k k^(3/2) 2^(-hk)` do not
converge: they oscillate with the fractional part of `k log_3 2` inside the
asymptotic window `(2 pi rho(1-rho))^(-1/2) (rho/(2 rho - 1)) * ((1-rho)/rho)^theta`,
`theta = ceil(k rho) - k rho in (0, 1]`, whose ends differ by the factor
`rho/(1-rho) = 1.709`; `W_k k^(3/2) 2^(-hk)` is observed in the window
`[9.6648, 11.0517]` for `500 <= k <= 3000`.

**Corollaries.**

* **C1 (THM-4479 sharpened to order).** The distance to provability in the
  strategy cube satisfies `N_k <= delta_k <= |Bad_k| = W_k` (THM-4479), hence
  `0.26 * 2^(hk) k^(-3/2) <= delta_k <= 545 * 2^(hk) k^(-3/2)`: the rate
  `(1/k) log_2(delta_k/2^(k-1)) -> -(1-h)` of THM-4479 becomes
  `delta_k = Theta(2^(hk) k^(-3/2))`, with the polynomial pinned. THM-4479's
  lower bound `N_k >= 2^(hk)/(3k^2)` loses `k^(1/2)` only through the crude
  binomial bound `C(k,m) >= 2^(k h(m/k))/(k+1)`.
* **C2 (THM-4485's chain).** `N_k <= nu <= FVS_(log_3 2)(k) <= FVS^odd <= delta_k <= W_k`,
  so the cycle-packing number, the periodic-deletion price
  `FVS_(log_3 2)(k)` and its odd version are all `Theta(2^(hk) k^(-3/2))`;
  the refuted Golomb analogue of THM-4485 (ratio `FVS/N_k >= 1.35` on a
  positive-density set of `k`) lives inside this constant window.
* **C3 (THM-4487 at `gamma = 1`).** On both sheets the dip count
  `D_b(X, 1) = #{n <= X : T_b^j(n) >= n, 0 <= j <= floor(log_2 n)}` satisfies
  `D_b(X, 1) = Theta( X^h (log X)^(-3/2) )`: THM-4487's bracket
  `c X^h log^(-3/2) X <= D <= C X^h log X` closes at its lower end.
* **C4 (reading).** `1 - h = I(g)` is the Chernoff rate of the multiplier
  walk (THM-4476 section 1.8) and `k^(-3/2)` is the ballot correction for a
  zero-drift walk pinned at a bounded height; both are now proved to be the
  common order of the four counts. The only Collatz-specific inputs are
  the irrationality of `log_2 3` (no ties) and the values `rho, h`.

## 1. Setup: the walk has no ties

Give the letters the real weights `s(0) = -1`, `s(1) = alpha := log_2 3 - 1`,
so a word `w` has partial sums `S_j(w) = o_j log_2 3 - j` and `M_j = 2^(S_j)`.
Because `log_2 3` is irrational, `S_i(w) = S_j(w)` for `i < j` would force
`(o_j - o_i) log_2 3 = j - i` with `j - i >= 1`, impossible; likewise
`S_j = 0` is impossible for `j >= 1`. **No two partial sums of any word
coincide, and none vanishes**; every "minimum", "maximum" and "record"
below is attained at a unique index, and strict and weak inequalities
agree. A word is *positive* if `S_j > 0` for all `1 <= j <= k`
(`W_k` counts positive words), *first-passage* if `S_j < 0` for
`1 <= j < m` and `S_m > 0` (`L_m` counts them), and *ends at its maximum* if
`S_m > S_j` for all `0 <= j < m` (`S_0 = 0`). The reversal
`w^R = (w_m, ..., w_1)` has partial sums `S_i(w^R) = S_m(w) - S_(m-i)(w)`.

## 2. Proof of Theorem A

*Step 1 (minimum decomposition).* Let `w` be positive of length `k >= 1` and
let `m in [1, k]` be the index of its minimum partial sum. Then
`a = w_1..w_m` is positive and *min-ending* (`S_j(a) > S_m(a)` for `j < m`),
and `b = w_(m+1)..w_k` is positive (its partial sums are
`S_(m+j)(w) - S_m(w) > 0`), possibly empty. Conversely, for `a` positive
min-ending and `b` positive, `ab` is positive with minimum at `|a|`. Hence,
with `P_m` the number of positive min-ending words of length `m`,

```text
W_k = sum_(m=1)^(k) P_m W_(k-m)   (k >= 1),   i.e.   W(t) = 1/(1 - P(t)),   P(t) = sum_(m>=1) P_m t^m.
```

*Step 2 (reversal).* `a` is positive min-ending of length `m` iff `a^R` is
first-passage: `S_i(a^R) = S_m(a) - S_(m-i)(a) < 0` for `i < m` and
`S_m(a^R) = S_m(a) > 0`, and back. So `P_m = L_m` and `W(t) = 1/(1 - L(t))`.

*Step 3 (ladder blocks).* A word `y` of length `n` ends at its maximum iff it
is a concatenation `y = x_1 x_2 ... x_r` of first-passage words: cut `y` at
its strict ascending ladder epochs (the times `j` with `S_j > S_l` for all
`l < j`; the last one is `n`); between two epochs the walk stays strictly
below the previous record (no ties), so each block is first-passage, and
conversely. Writing `E_(n,r)` for the number of words of length `n` that
end at their maximum with exactly `r` ladder epochs,
`[t^n] L(t)^r = E_(n,r)`, hence

```text
[t^n] log W(t) = [t^n] log 1/(1 - L(t)) = sum_(r>=1) E_(n,r)/r.
```

*Step 4 (rotation).* Let `x` be a word of length `n` with `S_n(x) > 0`, extended
`n`-periodically to a bi-infinite sequence, so `S_(j+n) = S_j + S_n` and
`S_l -> -infinity` as `l -> -infinity`. Call `j` a *record* if `S_j > S_l` for
all `l < j`; records exist, and `j` is a record iff `j + n` is. Let
`R(x) subset Z/n` be the set of record residues and `r(x) = |R(x)| >= 1`.
For `i in Z/n` let `rot_i(x) = x_(i+1) ... x_(i+n)`, whose partial sums are
`S_(i+l) - S_i`. Then

* `rot_i(x)` ends at its maximum iff `S_(i+n) > S_(i+l)` for `0 <= l < n`
  iff `i + n` is a record (for `l' < i` compare with `S_(l'+mn) - m S_n`,
  `l' + mn in [i, i+n)`) iff `i in R(x)`;
* if `i in R(x)`, the ladder epochs of `rot_i(x)` are the `l in [1, n]` with
  `i + l` a record (for the reverse inclusion use `S_(i+l) > S_i > S_(l'')`
  for `l'' < i`), so `rot_i(x)` has exactly `r(x)` ladder epochs.

Therefore `sum_(i in Z/n) [rot_i(x) ends at its maximum] / r(rot_i(x)) = r(x)/r(x) = 1`
for every `x` with `S_n(x) > 0`. Summing over all such `x` and regrouping by
`y = rot_i(x)` (for fixed `i`, `x -> y` is a bijection of the words with
positive total):

```text
B_n = sum_(i in Z/n) sum_(y : S_n(y) > 0) [y ends at its maximum]/r(y) = n sum_(r>=1) E_(n,r)/r.
```

With Step 3, `[t^n] log W(t) = B_n/n`, i.e. `W(t) = exp(sum_n B_n t^n/n)`;
differentiating, `t W'(t) = W(t) sum_n B_n t^n`, which is the recurrence. ∎

*Remark.* This is Spitzer's 1956 combinatorial lemma specialised to two
letters; for i.i.d. steps it is the identity
`sum_n t^n P(S_1 > 0, ..., S_n > 0) = exp(sum_n (t^n/n) P(S_n > 0))`
(Feller II, XII.7). The proof above is written out so that the theorem
rests on nothing but the irrationality of `log_2 3`; the recurrence is
checked exactly against a direct ballot DP for `k <= 300` and is integral
(as it must be) for `k <= 3000`.

## 3. Proof of Theorem B

*Binomial tails.* Let `j_0 = floor(n rho) + 1`, `p = j_0/n in (rho, rho + 1/n]`.
For `j >= j_0`, `C(n, j+1)/C(n, j) = (n-j)/(j+1) < (1-rho)/rho`, so
`C(n, j_0) <= B_n <= (rho/(2 rho - 1)) C(n, j_0) = 2.4094 C(n, j_0)`. The standard
Stirling bounds (`n >= 3`, so that `p < 1`)

```text
2^(n h(p)) / sqrt(8 n p(1-p))  <=  C(n, pn)  <=  2^(n h(p)) / sqrt(2 pi n p(1-p))
```

and `h(p) <= h(rho)` (as `p > rho > 1/2` and `h` decreases on `[1/2, 1]`) give
`b_n := B_n 2^(-hn) <= 2.4094 / sqrt(2 pi n p(1-p)) <= 2.05 n^(-1/2)` for
`n >= 30` (then `p <= 0.665`, `p(1-p) >= 0.2228`, and the bound is `2.036 n^(-1/2)`),
and `b_n n^(1/2) <= 1.983` for `1 <= n < 30` by direct evaluation (the
control's maximum over all `n <= 3000`); so `b_n <= 2.05 n^(-1/2)` for all
`n >= 1`. For the lower bound, `h(j_0/k) >= h - 1.4415/k` for `k >= 10`
(`|h'| <= log_2(0.7309/0.2691) = 1.4415` on `[rho, rho + 0.1]`, as in THM-4479)
and `p(1-p) <= 1/4`, so `C(k, j_0) >= 2^(hk) 2^(-1.4415)/sqrt(2k) = 0.2603 * 2^(hk) k^(-1/2)`
for `k >= 10`; the cases `k < 10` are read off the table. Since every
necklace with `j` ones has at most `k` rotations, `N_k >= C(k, j_0)/k`.

*Lower bound for `W_k`.* The recurrence has nonnegative terms, and its
`n = k` term is `B_k W_0 = B_k`, so `W_k >= B_k/k >= C(k, j_0)/k >= 0.26 * 2^(hk) k^(-3/2)`.
(This is the cycle-lemma bound of THM-4487's lower bound and of THM-4479's
Theorem 2, now read off the identity.)

*Convolution lemma.* Let `a, b >= 0` be sequences with `a_0 = b_0 = 0`,
`a_n <= A n^(-3/2)`, `b_n <= B n^(-3/2)`, `sum a = alpha`, `sum b = beta`.
Splitting `(a*b)_k = sum_i a_i b_(k-i)` at `i = k/2`,

```text
(a*b)_k <= (k/2)^(-3/2) (B alpha + A beta) = 2^(3/2) (A beta + B alpha) k^(-3/2).
```

*Upper bound for `W_k`.* Put `w_k = W_k 2^(-hk)` and `g_n = b_n/n <= C n^(-3/2)`
with `C = 2.05`, `g_0 = 0`, `sigma = sum_n g_n = 1.8989 (+ tail <= 0.0749)` (numerically to
`n = 3000`, plus the tail `sum_(n > 3000) 2.05 n^(-3/2) <= 4.1/sqrt(3000) = 0.075`).
Theorem A reads `sum_k w_k u^k = exp(G(u))`, `G(u) = sum g_n u^n`, so
`w_k = sum_(m>=1) (g^(*m))_k / m!`. By induction with the convolution lemma,
`(g^(*m))_k <= C m (2^(3/2) sigma)^(m-1) k^(-3/2)` (the step is
`2^(3/2)(C sigma^m + C m (2^(3/2) sigma)^(m-1) sigma) <= C (m+1)(2^(3/2) sigma)^m`),
hence

```text
w_k <= C k^(-3/2) sum_(m>=1) (2^(3/2) sigma)^(m-1)/(m-1)! = C e^(2 sqrt2 sigma) k^(-3/2) = 545 k^(-3/2).
```

*Oscillation.* `B_k k^(1/2) 2^(-hk) = sum_(i>=0) C(k, j_0 + i) k^(1/2) 2^(-hk)`,
and with `theta = j_0 - k rho in (0, 1]`, `k h(j_0/k) = hk + h'(rho) theta + O(1/k)`,
`h'(rho) = log_2((1-rho)/rho) = -0.7731`, and `C(k, j_0+i)/C(k, j_0) -> ((1-rho)/rho)^i`,
one gets `B_k k^(1/2) 2^(-hk) = (2 pi rho(1-rho))^(-1/2) (rho/(2rho-1)) ((1-rho)/rho)^theta (1 + o(1))`,
a function of `theta` alone whose values fill `(1.166, 1.992]`; since
`k rho mod 1` is equidistributed, the normalised `B_k` has no limit, and
`N_k = B_k/k + O(k 2^(k/2))` inherits the same window. ∎

## 4. The corollaries

*C1, C2.* THM-4479 proves `delta_k <= |Bad_k|` (flipping `Bad_k` is a
class-(i) strategy) and `delta_k >= N_k` (every expanding necklace needs a
flip), and THM-4485 proves the chain `delta_k >= FVS^odd >= FVS >= nu >= N`.
`|Bad_k| = W_k` because both count the residues with `M_j > 1` for all
`j <= k` (THM-4479's definition; the table agrees for `k <= 11` and THM-4479's
audit DP for `k <= 200`). Theorem B bounds both ends of every sandwich by
constants times `2^(hk) k^(-3/2)`.

*C3.* Fix `t` and `n in [2^t, 2^(t+1))`, so `floor(log_2 n) = t`, and let
`w` be its length-`t` word; `T_b^j(n) = M_j n + beta_j` with `|beta_j| <= (3/2)^j - 1`
and `sign(beta_j) = b` (for `j` with an odd step). *Plus sheet, lower bound:*
if `w` is positive then `T^j(n) >= M_j n > n` for all `j`, so
`D_+(2^(t+1) - 1, 1) >= sum_(s <= t) W_s >= W_t`. *Plus sheet, upper bound:* if
`w` is not positive, some `S_j < 0`, and `T^j(n) >= n` needs
`beta_j >= (1 - M_j) n`, i.e. `1 - 2^(S_j) <= (3/2)^t/2^t`, i.e.
`S_j > -2 (3/4)^t > -1` (`t >= 3`); then `11w` (two odd letters prepended,
partial sums `0.585, 1.170, 1.170 + S_j(w) > 0.17`) is positive, and
`w -> 11w` is injective, so the non-positive non-dippers of `[2^t, 2^(t+1))`
number at most `W_(t+2)`, and `D_+(X, 1) <= sum_(t <= log_2 X) (W_t + W_(t+2)) <= 2 sum_(t <= log_2 X) W_(t+2)`.
*Minus sheet, upper bound:* `beta_j <= 0`, so a non-dipper has `M_j > 1` for
all `j`, i.e. `w` positive: `D_-(X, 1) <= sum_(t <= log_2 X) W_t`. *Minus
sheet, lower bound:* for `w` positive of length `t - 2` and `n` in the class
of `11w`, `M_1 = 3/2`, `M_2 = 9/4`, `M_j = (9/4) 2^(S_(j-2)(w)) > 9/4` for `j > 2`,
so `T_-^j(n) >= (3/2) n - (3/2)^t >= n` once `2^(t-1) >= (3/2)^t` (`t >= 3`):
`D_-(X, 1) >= sum_(3 <= t <= log_2 X - 1) W_(t-2)`. In all four cases Theorem B
and the geometric growth `W_(t+1)/W_t -> 2^h` give `D_b(X, 1) = Theta(X^h (log X)^(-3/2))`.
(FINITE-EXACT: the counts of THM-4487's control at `gamma = 1` are
`367698` on both sheets at `X = 2^24`; `sum_(t < 24) W_t = 367698`.)

*What is not claimed.* Nothing about Collatz orbits; THM-4476's `X^(h+eps)`
keeps its `eps` (it comes from the dichotomy bootstrap, not from the
no-dip count: turning it into a power of `log X` is a separate target).
For `gamma < 1` the polynomial factor of `D(X, gamma)` is not determined
here; the barrier `(gamma - 1) t` is not a constant level, and the
expected order is `X^(E(gamma)) (log X)^(-1/2)` (a zero-drift walk pinned at
the end but free at the start). The constants `0.26` and `545` are not
sharp: the observed window of `W_k k^(3/2) 2^(-hk)` is `[9.6648, 11.0517]`.

## 5. Controls (FINITE-EXACT)

```text
rho = log_3 2 = 0.6309297536, h = h(rho) = 0.9499555272, 1 - h = 0.0500444728, lambda* = log_2(rho/(1-rho))/log_2 3 = 0.488077
(S) Spitzer identity k W_k = sum_n B_n W_(k-n): ballot DP == recurrence for all k <= 300: True; recurrence integral for all k <= 3000
    W_k for k = 1..20: [1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, 1295, 2114, 4228, 7495, 14990, 27328]
    B_n for n = 1..20: [1, 1, 4, 5, 6, 22, 29, 37, 130, 176, 562, 794, 1093, 3473, 4944, 6885, 21778, 31180, 94184, 137980]
(B) for 3 <= n <= 3000 and j0 = floor(n log_3 2) + 1: min(log2 C(n,j0) - lower bound) = 0.0376 (>= 0), max(log2 C - upper bound) = 0.0000 (<= 0), max B_n / ((rho/(2rho-1)) C(n,j0)) = 0.996365 (<= 1), max h(j0/n) - h = -2.92e-08 (<= 0)
(O) k, W_k k^1.5/2^(hk), N_k k^1.5/2^(hk), B_k k^0.5/2^(hk), W_k/N_k, frac(k log_3 2)
        5    1.6622    0.8311    0.4987     2.000  0.1546
       10    2.7959    0.8300    0.7689     3.368  0.3093
       20    4.6650    1.1796    1.1777     3.955  0.6186
       30    5.5330    1.5545    1.5544     3.559  0.9279
       50    6.6436    1.3002    1.3002     5.110  0.5465
       75    7.2158    1.1853    1.1853     6.088  0.3197
      100    7.6613    1.0652    1.0652     7.192  0.0930
      150    8.5838    1.5256    1.5256     5.626  0.6395
      200    8.9188    1.1962    1.1962     7.456  0.1860
      300    9.3945    1.2891    1.2891     7.288  0.2789
      500   10.1281    1.4548    1.4548     6.962  0.4649
      700   10.2209    1.6224    1.6224     6.300  0.6508
     1000   10.6555    1.8973    1.8973     5.616  0.9298
     1500   10.6182    1.4257    1.4257     7.448  0.3946
     2000   10.9233    1.8367    1.8367     5.947  0.8595
     2500   10.5588    1.3780    1.3780     7.662  0.3244
     3000   10.8314    1.7719    1.7719     6.113  0.7893
    over 500 <= k <= 3000: W_k k^1.5/2^(hk) in [9.6648, 11.0517];  B_k k^0.5/2^(hk) in [1.1355, 1.9827]
    over 500 <= k < 560: N_k k^1.5/2^(hk) in [1.1355, 1.9395]
    b_n sqrt(n) = B_n n^0.5 2^(-hn) <= 1.9827 for n <= 3000 (the elementary bound is 2.41/sqrt(2 pi rho(1-rho)) = 1.9919); sigma = sum b_n/n = 1.898923 (+ tail <= 0.0727); explicit upper constant C e^(2 sqrt2 sigma) = 525.64
    min over 1 <= k <= 3000 of W_k k^1.5/2^(hk) = 0.5176;  min over 2 <= k < 400 of N_k k^1.5/2^(hk) = 0.6999  (both >= 0.26)
    sum_(t<T) W_t for T = 8, 12, 16, 20, 24: [32, 281, 2903, 31730, 367698]
```

Read: (S) the identity holds exactly; (B) the four inequalities behind the
constants hold at every `n <= 3000`; (O) the three normalised counts stay in
bounded windows, `W_k/N_k` stays between about `5.6` and `7.7`
(THM-4479's sandwich for `delta_k` is a constant-width window), and the
fractional part of `k log_3 2` drives the oscillation of `B_k` and `N_k`.
