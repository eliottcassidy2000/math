# The dip spectrum of 3n±1: the number of integers below X whose orbit never falls below n^γ within log₂n steps is X^{h(γ/log₂3)+o(1)}; the thin-divergence exponent h(log₃2) is sharp for its counting lemma; one entropy curve carries the repo's Collatz constants

**Status: PROVED (elementary; two-sided exponent, both sheets; proofs below;
self-audited; the `gamma = 1` boundary was repaired the same day and
independently audited by the 223 crossroads session, whose one-step buffer is
the proof of record there; the whole note was independently audited by a
separate agent, verdict SOUND, recorded in THM-4487) + FINITE-EXACT controls. Consequence for THM-4476: the
counting lemma behind the thin-divergence exponent `h* = h(log_3 2)` is
sharp up to `X^(o(1))`, so `h*` cannot be lowered by improving that lemma;
the exponent of an actual divergent orbit remains OPEN (conjecturally no
polynomially dense divergent orbit exists at all). Session
`collatz-exponent-atlas-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/collatz_dipspectrum_20260926.py`
(counts to `2^24`, both sheets); output `collatz_dipspectrum_20260926.out`.
Inheritance: the Terras bijection and carry bound (THM-4476 note, sections
1.2–1.3), the cycle-lemma rotation of THM-4478 section 4 (crossroads), Korec's
theorem as cited in the [barrier atlas](collatz_procgen_20260922_barrier_atlas.md),
and the four settings of the exponent `1 - h` in synthesis section 2k.

## 0. Statements

Throughout `b = +-1`, `T_b(x) = x/2` (`x` even), `(3x+b)/2` (`x` odd),
`alpha = log_2 3 = 1.58496`, `h(p) = -p log_2 p - (1-p) log_2 (1-p)`, and

```text
D_b(X, gamma) = #{ n in [1, X] : T_b^i(n) >= n^gamma for all 0 <= i <= floor(log_2 n) }.
```

**Theorem 1 (dip spectrum).** For every `gamma in (log_4 3, 1]` there are
constants `c, C > 0` (depending on `gamma`) with, on both sheets,

```text
c X^(h(gamma/alpha)) / log^(3/2) X  <=  D_b(X, gamma)  <=  C X^(h(gamma/alpha)) log X;
```

including `gamma = 1` on both sheets (the one-step buffer of section 1.3);
in particular `log D_b(X, gamma)/log X -> h(gamma/alpha)`;

and for `gamma in (log_2(3/2), log_4 3]` the limit is `1`. At `gamma = 1`
the exponent is `h(log_3 2) = 0.949956` (the Terras undecided count); as
`gamma` decreases to Korec's exponent `log_4 3 = 0.792481` the exponent
rises to `1`.

**Proposition 2 (sharpness of THM-4476's counting lemma).** For
`theta in (0, theta_0)`, `theta_0 = 1 - alpha/2 = 0.207519`, and
`rho = (1 - theta)/alpha`, the no-dip set
`F_b(X, theta) = {y <= X : T_b^i(y) >= y X^(-theta), 0 <= i <= floor(log_2 X)}`
satisfies

```text
c X^(h(rho)) / log^(3/2) X  <=  #F_b(X, theta)  <=  C X^(h(rho)) log^2 X.
```

The upper bound is the lemma of THM-4476 section 1.4; the lower bound is
new. Hence the recursion `N(X) <= k N(X^(1-theta)) + #F_b(X, theta)` cannot
give an exponent below `h(rho(theta))` for any `theta`, and its `theta -> 0`
limit `h* = 0.949956` is the best this dichotomy can produce.

**Corollary 3 (one curve).** The function `E(gamma) = h(max(1/2, gamma/alpha))`
carries, at explicit points, the constants that recur across the repo's
Collatz threads:

| point on the curve | value | where it appears |
|---|---|---|
| `E(1) = h(log_3 2)` | `0.949956` | Terras's undecided density `2^(-(1-h)L)`; the exceptional 2-adic dimension (choice ladder); the thin-divergence exponent (THM-4476); `1 - E(1) = 0.050044` is the sharp price exponent in the pairing family, for arbitrary edits and in the strategy cube (THM-4475/4478/4479, synthesis 2k) and the Chernoff rate `I(T)` |
| `E(log_4 3) = 1` | at `gamma = 0.792481` | the finite-window counting exponent drops below `1` exactly for `gamma > log_4 3`, matching Korec's sufficient threshold (no necessity claim for the unrestricted orbit is made); `1 - log_4 3 = theta_0 = 0.207519` is the mean drift per `T`-step in bits and the largest dip exponent the counting sees |
| `-E'(1)` | `0.48808` | exactly the Chernoff tilt `lambda*` of the price exponent: `lambda*` minimises `(2^(-lambda) + (3/2)^lambda)/2`, so `3^(lambda*) = 1/log_2(3/2)`, i.e. `lambda* = log_3(1/log_2(3/2)) = -log_2(alpha - 1)/alpha = -h'(log_3 2)/alpha = -E'(1)` (with `h'(log_3 2) = log_2(alpha - 1) = log_2(log_2(3/2)) = -0.77358`); large-deviation duality, tilt = slope of the rate |
| `1/E(1)` | `1.052681` | the growth exponent below which no divergent orbit exists (`m_j > j^a` infinitely often for `a < 1/h*`, THM-4476 Cor. 3) |
| `1 - alpha/2` | `0.207519` | `theta_0`; the drift; Korec's `1 - log_4 3` |
| `log_2(3/2) = alpha - 1` | `0.584963` | the carry exponent `(3/2)^(log_2 X) = X^(0.585)`, the lower limit of validity of the dip spectrum, and the ratio `(1 - rho_0)/rho_0` |

A numerological trap is recorded for contrast: THM-4475's lower-bound
exponent `0.7737 = (1 - h*) + log_2((3 + sqrt 13)/4)` and the slope constant
`|h'(log_3 2)| = |log_2 log_2(3/2)| = 0.7735` agree to three decimals and
are unrelated.

## 1. Proofs

### 1.1 Tools

* **Terras bijection** (THM-4476 note, 1.2): for constant odd `b` the parity
  word of length `t` is a bijection `Z/2^t -> {0,1}^t`; every interval of
  `2^t` consecutive integers contains exactly one representative of each
  word.
* **Affine form** (1.3): if the word of `n` has `o_i` odd letters among its
  first `i`, then `T_b^i(n) = 3^(o_i) n / 2^i + beta_i`, with
  `0 <= beta_i <= (3/2)^i - 1` on the plus sheet and `-((3/2)^i - 1) <= beta_i <= 0`
  on the minus sheet. Write `S_i = i - o_i alpha`, so `3^(o_i)/2^i = 2^(-S_i)`.
* **Entropy bounds.** `C(t, o) <= 2^(t h(o/t))`; for `rho > 1/2` and
  `t >= t_0(rho)`, `sum_(o >= rho t - c) C(t, o) <= (t+1) 2^(t h(rho)) (rho/(1-rho))^c`;
  and, for fixed `rho in (0,1)`, `C(t, ceil(rho t) + j) >= c_(rho,j) 2^(t h(rho))/sqrt t`
  for every fixed `j >= 0` (Stirling, and the ratio of consecutive binomials is
  bounded); the cruder form `2^(t h(rho))/(t+1)` fails for small `t`.
* **Cycle lemma.** Let `w` be a word of length `t` with partial sums `S_i`
  (any real step values) and let `p` be an index where `S_p` is maximal.
  The rotation of `w` starting after position `p` has all partial sums
  `<= max(0, S_t)`. (For `i <= t - p` the new sum is `S_(p+i) - S_p <= 0`; for
  `i > t - p` it is `S_t - S_p + S_(i-t+p) <= S_t`.) A word has at most `t`
  rotations, so among the `C(t, o)` words with `o` odd letters at least
  `C(t, o)/t` have all partial sums `<= max(0, S_t)`. This is the rotation
  step of THM-4478 section 4.

### 1.2 Theorem 1, upper bound

Let `gamma > log_4 3`, so `rho = gamma/alpha > 1/2`. Take `n in [2^t, 2^(t+1))`
with `T_b^i(n) >= n^gamma` for all `i <= t`, and let `o` be the number of odd
letters in its word of length `t`. At `i = t`, the affine form gives
`3^o n/2^t >= n^gamma - (3/2)^t` on the plus sheet and `3^o n/2^t >= n^gamma`
on the minus sheet. Since `(3/2)^t <= n^(log_2(3/2))` and `gamma > log_2(3/2)`,
for `n >= n_0(gamma)` we have `(3/2)^t <= n^gamma/2`, hence
`3^o/2^t >= n^(gamma-1)/2 >= 2^((t+1)(gamma-1) - 1)`. Taking `log_2`:
`o alpha >= t gamma + gamma - 2`, i.e. `o >= rho t - 2`. The number of words
of length `t` with at least `rho t - 2` odd letters is at most
`(t+1)(rho/(1-rho))^2 2^(t h(rho))` for `t >= t_0`, each word is one class
modulo `2^t`, and each class has exactly one representative in
`[2^t, 2^(t+1))`. So the dyadic block contributes at most
`(t+1)(rho/(1-rho))^2 2^(t h(rho)) + n_0`, and summing the
geometric progression over `t <= log_2 X` gives `D_b(X, gamma) <= C X^(h(rho)) log X`. ∎

### 1.3 Theorem 1, lower bound

Let `gamma in (log_2(3/2), 1]`, put `rho = max(1/2, gamma/alpha)` and
`o = ceil(rho t) + 1`. By the cycle lemma at least `C(t, o)/t` words of
length `t` with `o` odd letters have all partial sums `S_i <= max(0, S_t)`,
and `S_t = t - o alpha <= t(1 - rho alpha) - alpha`, so `2^(-S_t) >= 3 * 2^(-t(1 - rho alpha))`.
Take the representative `n in [2^t, 2^(t+1))` of such a word. For `i <= t`,
the affine form (both sheets) gives

```text
T_b^i(n) >= n 2^(-S_i) - (3/2)^t >= n 2^(-max(0, S_t)) - n^(log_2(3/2)).
```

*Case `gamma < 1`.* Then `rho alpha = gamma < 1` (for `rho = gamma/alpha`) and
`S_t <= t(1 - gamma) - alpha`, which is positive for `t >= t_1(gamma)`; so
`max(0, S_t) = S_t` and the first term is at least
`n 2^(-S_t) >= 3 n 2^(-t(1 - gamma)) >= 3 n^gamma` (using `2^t <= n`). Hence
`T_b^i(n) >= 3 n^gamma - n^(log_2(3/2)) >= n^gamma` once `n >= n_1(gamma)`,
because `gamma > log_2(3/2)`. (If `rho = 1/2 > gamma/alpha` the same holds with
`S_t <= t(1 - alpha/2) - alpha`.)

*Case `gamma = 1`, plus sheet.* `S_t <= -alpha < 0`, so every rotated partial sum
is `<= 0` and, the carry being nonnegative, `T_+^i(n) >= n 2^(-S_i) >= n`.

*Case `gamma = 1`, minus sheet.* Here `max(0, S_t) = 0` and a partial sum
`S_i < 0` may be arbitrarily close to `0`, while the carry is negative, so
the rotation alone does not give `T_-^i(n) >= n`. *One-step buffer* (the
223 crossroads session, same day, independently audited there with an exact
boundary check of 3,142 rotated inputs and 2,244 sheet realisations through
length 16, `crossroads223_20260926_boundary_audit.py`; MISTAKES 2026-09-26,
"223 crossroads: endpoint margin, finite lifts, and global ownership"): take a tail of length `t - 1` with
`o = ceil((t-1)/alpha) + 1` odd letters, rotate it after a maximum of its
partial sums (so its partial sums are `<= 0`, its total being `<= -alpha`),
and prepend one odd letter. Every prefix sum of the full word is then
`<= 1 - alpha`, i.e. every prefix multiplier is at least `3/2`, so on both
sheets `T_b^i(n) >= (3/2) n - (3/2)^t >= n` for `n >= 2^t`, `t >= 3`; distinct
tails give distinct words, and there are at least `C(t-1, o)/(t-1)` of them,
the same entropy exponent and polynomial prefactor. This is the proof of
record. Two further repairs are kept for the record. *Tilted rotation*
(suggested by the audit): take `o = ceil(t/alpha) + 2` and apply the cycle
lemma to the tilted steps `x_j + alpha/t`, whose partial sums are
`S'_i = S_i + i alpha/t` and whose total is `S'_t = S_t + alpha <= -alpha < 0`.
The rotation gives `S'_i <= 0`, i.e. `S_i <= -i alpha/t`, so
`2^(-S_i) >= 3^(i/t) >= 1 + (i/t) ln 3` and
`T_-^i(n) >= n + n (i/t) ln 3 - (3/2)^i > n` for `n >= 2^t`, `t >= 3`, every
`1 <= i <= t`; the count is again `C(t, o)/t`, so the polynomial-log lower
bound holds on the minus sheet as well. *Block alternative* (the first
repair, kept for the record): fix `L` and
let `W_L` be the set of words of length `L` with `o_L = ceil(L/alpha) + 1` odd
letters and all partial sums `<= 0` (at least `C(L, o_L)/L` of them, by the
cycle lemma, since their total `S_L <= -alpha < 0`). For `t = jL`, concatenate
`j` words of `W_L`. A partial sum inside the `r`-th block is at most
`-(r-1) alpha + 0`; inside the first block it is `< 0` and, since the finitely
many values `i - o alpha` with `i <= L` are nonzero, at most `-c_L` for some
`c_L > 0` depending only on `L`. Hence for `n >= n_2(L)`,
`T_-^i(n) >= n 2^(c_L) - (3/2)^t >= n` in the first block and
`T_-^i(n) >= 3n - (3/2)^t >= n` afterwards. The number of such words is at
least `(C(L, o_L)/L)^(t/L) = 2^(t h(o_L/L) - O((t/L) log L))`, and
`h(o_L/L) -> h(1/alpha) = h*` as `L -> infinity`, which gives
`D_-(X, 1) >= X^(h* - o(1))`. *Lineage.* The gap (an unjustified `3n` margin at
`gamma = 1` in the first draft) was found by the author and by the 223
crossroads session independently on 2026-09-26; the block construction was the
first repair (`X^(h* - o(1))`), the crossroads buffer and the auditor's tilted
rotation both restore the polynomial prefactor.

In every case the top complete dyadic block `[2^t, 2^(t+1))` (`2^t > X/4`)
contributes at least `C(t, o)/t >= c 2^(t h(rho))/t^(3/2)` to `D_b`, so
`D_b(X, gamma) >= c' X^(h(rho))/log^(3/2) X`. (When `S_t < 0` in the first
case the first term is `n >= 3 n^gamma` for `n >= 3^(1/(1-gamma))`.) For
`gamma <= log_4 3`, `rho = 1/2` and the exponent is `h(1/2) = 1`. ∎

### 1.4 Proposition 2

The upper bound is THM-4476's lemma. For the lower bound take `k = floor(log_2 X)`,
`o = ceil(rho (k-1)) + 2`, and the words of length `k - 1` with `o` odd letters
and all partial sums `<= max(0, S_(k-1))`, at least `C(k-1, o)/(k-1)` of them. The interval `(X/2, X]` contains at least `2^(k-1)`
consecutive integers, hence a representative `y` of every such word. For
`i <= k - 1` the affine form gives `T_b^i(y) >= y 2^(-max(0, S_(k-1))) - (3/2)^k`,
and `S_(k-1) <= (k-1)(1 - rho alpha) - 2 alpha <= (k-1) theta - 2 alpha`, so
`T_b^i(y) >= 9 y X^(-theta) - (3/2)^k` for `i <= k - 1` (using `2^(k-1) <= X`), and the
last step `i = k` is at worst a halving, so `T_b^k(y) >= (9/2) y X^(-theta) - (3/2)^k`.
For `X >= X_1(theta)` both are at least `y X^(-theta)` for every `y > X/2` (as
`theta < theta_0 < 1 - log_2(3/2)`). So
`#F_b(X, theta) >= C(k-1, o)/(k-1) >= c X^(h(rho))/log^(3/2) X`. ∎

**Why the exponent of THM-4476 stops at `h*`.** THM-4476 bounds the orbit
points below `X` by `k N(X^(1-theta)) + #F_b(X, theta)`. The first term is
negligible under any bootstrap, so the exponent produced is
`inf_theta h(rho(theta)) = h*`, and Proposition 2 shows `#F_b` really is that
large. A better exponent for a divergent *orbit* needs a constraint on
no-dip points beyond their `log_2 X`-window, which means information about
small integers in residue classes modulo more than `X`: the transversality
barrier of the synthesis. The auditor's numerics for THM-4476 (local
slopes at the lemma exponent) are the finite-size face of Proposition 2.

### 1.5 The general form: the dip spectrum is a constrained maximum entropy

For a contracting Conway/Matthews–Watts map `g(x) = (p_i x + q_i)/m` on
`x = i mod m` (`p_i` coprime to `m`, `a_i = p_i/m`, `prod p_i < m^m`,
`max p_i < m^2`), the same two proofs give, for the window
`floor(log_m n)`,

```text
#{ n <= X : g^j(n) >= n^gamma for all j <= floor(log_m n) } = X^(E_g(gamma) + o(1)),
E_g(gamma) = max { H_m(pi) : pi a probability law on Z/m with sum_i pi_i log_m a_i >= gamma - 1 },
```

with `H_m` the entropy in base `m`: the upper bound counts residue words by
their letter frequencies (multinomial entropy), the lower bound rotates a
word of the optimal frequency after its maximal partial sum. For `T_(+-1)`
(`m = 2`, `a = (1/2, 3/2)`) the constraint is `pi(odd) >= gamma/log_2 3` and
`E_g = h(max(1/2, gamma/alpha))`, Theorem 1. The proof is the one above with
`2` replaced by `m` and the binomial by the multinomial; it is not written
out separately here.

**Status of the audit:** SOUND (theorem file, audit field).

## 2. Controls (FINITE-EXACT)

`collatz_dipspectrum_20260926.py 24` iterates every `n <= 2^24` for
`floor(log_2 n)` steps on both sheets and buckets by the minimum ratio
`log T^i(n)/log n`. Single-doubling slopes carry a sawtooth (the window
length jumps at powers of two), so four-doubling slopes are reported.

| `gamma` | 0.792 | 0.820 | 0.850 | 0.880 | 0.910 | 0.940 | 0.970 | 1.000 |
|---|---|---|---|---|---|---|---|---|
| predicted exponent `h(rho)` | 1.0000 | 0.9991 | 0.9962 | 0.9912 | 0.9841 | 0.9749 | 0.9635 | 0.9500 |
| `log D/log X` at `2^20` (plus) | 0.9361 | 0.9216 | 0.9025 | 0.8859 | 0.8597 | 0.8273 | 0.7890 | 0.7477 |
| `log D/log X` at `2^24` (plus) | 0.9474 | 0.9351 | 0.9182 | 0.9002 | 0.8770 | 0.8497 | 0.8069 | 0.7703 |
| `log D/log X` at `2^24` (minus) | 0.9474 | 0.9351 | 0.9182 | 0.9002 | 0.8770 | 0.8497 | 0.8069 | 0.7703 |
| slope `2^16 -> 2^20` | 1.0043 | 1.0065 | 0.9794 | 1.0031 | 0.9639 | 0.9709 | 0.9196 | 0.8626 |
| slope `2^20 -> 2^24` | 1.0038 | 1.0027 | 0.9965 | 0.9717 | 0.9632 | 0.9618 | 0.8963 | 0.8837 |
| `h(rho) - 1/(k ln 2)`, `k = 22` | 0.9344 | 0.9336 | 0.9306 | 0.9256 | 0.9185 | 0.9093 | 0.8979 | 0.8844 |

The counts carry a polynomial prefactor (the prefix condition costs a factor
of order `1/t` at word length `t`, see the cycle lemma), so the local slope at
`k` steps sits below `h(rho)` by a term of order `1/k`; the last two rows
compare one such heuristic with the data. An earlier version of the control script counted `n < 2^(t+1)` into
the row for `2^t`; the corrected counts are the ones shown.

The two sheets agree to within 196 integers at every size to `2^24` (and
exactly at `gamma = 1`), as the sheet-blind class count predicts. Conventions
(audit finding 11): the control script does not count `n = 1`, so its
values are the exact `D_b` minus one, and the column headed `0.792` is
`gamma = 0.7925`, slightly above `log_4 3 = 0.792481`. The `1/(k ln 2)`
correction matches the observed slopes for `gamma >= 0.97` and undershoots
by `0.04`-`0.07` for smaller `gamma`, where the polynomial prefactor is
different; it is a heuristic for the top row only.

## 3. Reading

* **Korec's theorem is the endpoint of Terras's count.** Korec proved that
  almost every `n` dips below `n^c` for `c > log_4 3`; Theorem 1 says the
  exceptional set has exponent `h(c/log_2 3) < 1`, decreasing to
  `h(log_3 2) = 0.95` at `c = 1`, and that below `log_4 3` the exponent is `1`.
  The whole family of constants in Corollary 3 is one function of the
  critical density.
* **What the curve does not know.** It is a statement about residue classes
  modulo `2^(log_2 n)`; it is sheet-blind (identical exponents for `3n+1` and
  `3n-1`), it is blind to defects, and its `q`-analogue for `qn+1` has
  `rho = gamma/log_2 q`, which is below `1/2` for all `gamma <= 1` once
  `q >= 5`: no Korec-type dip theorem exists for `5n+1`, in agreement with
  the expected divergent orbits there.
* **For the atlas.** The recurrence of `0.95`, `0.05`, `0.7925`, `0.2075`,
  `1.0527`, `0.488` across the procgen, crossroads and opus lanes is
  structural: they are values, slopes and reciprocals of `E(gamma)` at two
  points. Closed forms in `alpha = log_2 3`: `h* = log_2 alpha - (1 - 1/alpha) log_2(alpha - 1)`,
  `1 - h*` the Chernoff rate, `lambda* = log_3(1/(alpha - 1))`, `theta_0 = 1 - alpha/2`. The near-coincidence `0.7737 ~ 0.7735` is not.
