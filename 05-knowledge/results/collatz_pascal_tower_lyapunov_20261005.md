# The Lyapunov exponent of the step-1 Pascal tower (2026-10-05): rate `0.56988 +- 0.00004` against the universal `0.5727 +- 0.0001`, on an exact twisted-Pascal reformulation (the tower is the Fourier transform of the random backward Syracuse law at a Haar frequency); the excess is a bulk effect of the per-level increments while the second moment is carried by ever rarer frequencies; `q = -3` has a larger excess and `q = 1` the fully coherent value; the "closed families with a common decreasing rank" question is the left tail of the same renewal process (no-descent rate `a^a/(3 (a-1)^(a-1)) = 0.946505`, `a = log_2 3`, PROVED), and `{2,3,11}` does not enter

**Session:** opus, `pascal-lyapunov-20261005` (worktree `codex/session-pascal-lyapunov-20261005`), 2026-10-05.
Owner's directive: "compute the Lyapunov exponent of the step-1 Pascal tower; assemble families whose dependency
rules close their entire union, with a common decreasing rank across family changes; consider deeply recursive
`p`-ary trees; the excess as a heavy-tail effect from arithmetic deep digits favoring small primes; how `{2,3,11}`
and `eta(tau)^2 eta(11 tau)^2` factor in."
**Inherits (cited):** the Fourier note
[`collatz_fourier_experiments_20261004.md`](collatz_fourier_experiments_20261004.md) sections 4j (Pascal tower
identity, PROVED), 4k (uniform-start second moment exactly `3^-n`, PROVED), 4m (uniform-start rates `0.5698 +-
0.0008` / `0.5726 +- 0.0005`, 24 seeds, OBSERVED); [HYP-9176](../hypotheses/HYP-9176-h1-decomposition-digit-incoherence-plus-ridge-control-universal-cold-rate.md);
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md) (coefficient-descent
cylinders); the [partition-cover obstruction](collatz_partition_cover_20261004.md); the finite-seed packages
([kernel](finite_seed_kernel_20261005.md), [refuel tree](finite_seed_refuel_tree_20261005.md),
[bootstrap boundary](finite_seed_bootstrap_boundary_20261005.md)); the [level-eleven note](eta11_lossless_coordinates_20261005.md).

**Status: PROVED (the twisted-Pascal reformulation and its column recursion; the no-descent large-deviation
rate in closed form); FINITE-EXACT / VERIFIED (every number below is from the two scripts named in section 7,
validated against brute force, against the Fourier note's E1 values to four digits, and by two independent
kernels agreeing to ten digits); OBSERVED (the rates; `200` seeds at level `3000`, `100` at `6000`);
CONJECTURED (the mechanism, section 4, with its decisive test recorded); OPEN (an exact value of the exponent;
H1; Collatz). NOT independently audited; audit OWED. Nothing here is a Collatz step.**

---

## 0. What is new, in one screen

1. **The object is a twisted Pascal triangle (PROVED, section 1).** Writing the uniform-start `q`-tower of the
   Fourier note as a path sum over the valuation words `c = (c_1, ..., c_N)` gives `f_N(0) = E_c e(R Phi_N)`,
   `Phi_N = sum_i q^(i-1) 2^(-S_i)`, `S_i = c_1 + ... + c_i`: the Fourier transform, at the Haar-random integer
   frequency `R`, of the law of the `N`-fold random backward Syracuse iterate `y -> (1 + q y)/2^c` of `0`. In the
   coin picture (`S_i` = position of the `i`-th one of a fair coin sequence) this is Pascal's rule with a unimodular
   twist on every up-step, `H_m(j) = (H_(m-1)(j) + omega_(m,j) H_(m-1)(j-1))/2`, `omega_(m,j) = e((R q^(j-1) mod
   2^m)/2^m)`, and `f_j(0)` is the total flux into column `j`. Computing column by column (a one-dimensional
   recursion in the depth with a per-column scale) costs `5 N^2` twists per sample against `N^2 A^2 / 2` for the
   window recursion (`A = 40`): one seed at `N = 1500` takes `1.0 s` against `5 s` in the October 4 log, and one
   sample gives every level `j <= N` at the same `R`.
2. **The Lyapunov exponent of the step-1 tower (OBSERVED, section 2).** Over levels `1000..6000` (100 seeds)
   and `1000..3000` (200 seeds): `lambda_3 = -0.56234 +- 0.00006`, i.e. **rate `0.56988 +- 0.00004`**. The step-2
   tower (`q = 5`), the i.i.d.-digit model and the random-odd-multiplier tower all give **`0.5727 +- 0.0001`**
   (`-0.5573 +- 0.0001`). The excess is `0.0049` per level in the log, more than twenty standard errors, and the
   Jensen gaps against the exact mean-square rate `1/sqrt3` are `0.0130` and `0.0081` per level. The step-1
   tower has a longer transient (local rate `0.5684` on levels `0..500`, `0.5692` on `500..1000`, flat at
   `0.5699 +- 0.0002` from `1000` on); the October 4 value `0.5698 +- 0.0008` on `300..1500` was pulled down by it.
3. **The `q`-landscape (OBSERVED, section 2c).** At level `1000` (60 seeds each): `q = 1` gives `0.5068` (the
   fully coherent tower, where `Phi` is the coin string itself), `q = -3` gives **`0.5647`** (a larger excess
   than `q = 3`), `q = 3` gives `0.5693`, and `q = 5, 7, 9, ..., 257` all sit in `0.5720..0.5729`. So the step-1
   tower is not the extreme case: the sign of the multiplier matters more than its binary length.
   [[PENDING-Q3000]]
4. **Where the excess lives (section 3).** (a) In the per-level increments `log|f_j/f_(j-1)|` the step-1 tower is
   not heavier-tailed: its left tail (increments below `-1`) contributes `+0.004` relative to the i.i.d. model,
   the deficit `-0.005` comes from the band `[-1, -0.5)` (`-0.0033`) and from fewer large positive increments
   (`-0.0045` from increments above `0.5`). (b) Over frequencies, the second moment (exactly `3^-N` for every
   model) is carried by ever rarer `R`: at level `100` the top `0.1%` of frequencies carry `77%` of it for
   `q = 3`, `45%` for `q = 7`, `30%` for `q = 5`, `27%` for `q = 9`, `34%` for i.i.d. digits; at level `200`,
   `85%` for `q = 3`. The pressure `P(s) = (1/N) log E|f_N|^(2s)` of the `q = 3` tower lies below the others for
   `s < 1`, meets them at `s = 1` (exactly) and lies above them for `s > 1`: same mean square, a more
   multifractal distribution. (c) [[PENDING-A]] (d) [[PENDING-P]]
5. **Mechanism (CONJECTURED, section 4).** The real backward map `y -> (1 + q y)/2^c` with `c ~ geometric(p)`
   has Lyapunov exponent `gamma(q, p) = log|q| - (log 2)/p`: contracting iff `|q| < 2^(1/p)`, i.e. for `p = 1/2`
   iff `|q| <= 3`. The exact pair-coherence census (section 3e) shows the law of `Phi mod 1` concentrated above
   the uniform law at coarse scales exactly for `q = +-3` (and not for `q = +-1`, whose law is uniform). The
   candidate mechanism is that the typical Fourier decay at a Haar frequency is faster when the real law of the
   path values is stationary (contracting multiplier) than when it spreads (expanding), with the sign of `q`
   setting the shape of the stationary law. Decisive test: the biased-coin scan moves the threshold
   (`q = 3` expanding for `p >= 0.64`, `q = 5` contracting for `p <= 0.43`). [[PENDING-THRESHOLD]]
6. **Families whose dependency rules close their union, with a common decreasing rank (section 5; typed).**
   With the integer itself as the rank, a union of smaller-child families is closed iff its leak is empty; for
   every bounded-depth portfolio the leak after `k` odd steps has density exactly `P(Bin(floor(k log_2 3), 1/2)
   >= k)`, which decays at the rate `e^(-I) = a^a / (3 (a-1)^(a-1)) = 0.946505...`, `a = log_2 3` (PROVED,
   Cramer for the geometric renewal; this is the repo's "no-descent rate `3^(h*-1)`" with `h* = 1 - I/log 3 =
   0.94996`). It is the left tail of the SAME coin process whose Fourier transform is the Pascal tower. The
   partition-cover theorem (CITED) says one arithmetic cell escapes every bounded-depth portfolio, and the
   depth-growing valuation-fuel families (the `H` kernel) close only until their fuel is spent. No assembly with
   a bounded dependency depth closes; this note adds no closure theorem and claims none.
7. **`{2, 3, 11}` and `eta(tau)^2 eta(11 tau)^2` do not enter (section 6).** The Pascal tower uses the coin
   (binary renewal) and the multiplier `q` only; the level-eleven bridges of the repo are local matrix identities
   (Frobenius-at-2 on 3-torsion, the index-2 lattice cube root with local factor `1 + T + 3T^2`), typed there as
   coordinate maps, not conjugacies. The only place `11` appears here is `q = 11`, whose rate is universal.

---

## 1. The object and its exact reformulation (PROVED)

**The uniform-start `q`-tower** (Fourier note 4k/4m). Level-`n` phases `theta_(n,d) = (x_n mod 2^d)/2^d`,
`x_n = q^(N-n) x_N`, `x_N = R` a Haar-random 2-adic unit (uniform odd residue to the window depth); the window
recursion `f_n(k) = sum_(c>=1) 2^-c e(theta_(n,c-k)) f_(n-1)(k-c)`, `f_0 = 1`; the tower is "step-1" for
`q = 3 = 2 + 1` (consecutive levels are cube roots, 4c; every level is a step-1 binomial transform of the deepest,
4j) and "step-2" for `q = 5 = 4 + 1`.

**Path sum.** Expanding the recursion from level `N` down to level `1`, with `D_n = c_n + ... + c_N` the depth
reached at level `n`,
```
f_N(0) = sum_(c_1..c_N >= 1) 2^(-sum c_n) prod_(n=1)^N e(theta_(n, D_n)) = E_c e(R Phi_N),
Phi_N = sum_(i=1)^N q^(i-1) 2^(-S_i),   S_i = c_1 + ... + c_i   (reindexed from the deepest level),
```
because `theta_(n,d) = frac(q^(N-n) R / 2^d)`. So `f_N(0) = hat nu_N(R)` is the Fourier transform at the integer
frequency `R` of the law `nu_N` of `Phi_N`, and `Phi_N = T_(c_1) o ... o T_(c_N)(0)` with `T_c(y) = (1 + q y)/2^c`:
the `N`-fold random backward Syracuse iterate of `0` on the real line with i.i.d. geometric(1/2) exponents.
Only `Phi_N mod 1` matters (`R` is an integer), but the circle projection is not Markov, so the real-valued
process is the natural state.

**Twisted Pascal triangle.** The partial sums `S_i` of i.i.d. geometric(1/2) variables are the positions of the
ones in a fair coin sequence. Let
```
H_m(j) := 2^-m sum_(coin words of length m with j ones) prod_(i<=j) e( R q^(i-1) / 2^(S_i) ).
```
Then Pascal's rule holds with a twist on the up-step,
```
H_m(j) = (1/2) [ H_(m-1)(j) + omega_(m,j) H_(m-1)(j-1) ],   omega_(m,j) = e( (R q^(j-1) mod 2^m) / 2^m ),
```
and `f_j(0) = sum_(m>=j) (1/2) omega_(m,j) H_(m-1)(j-1)` (the total flux into column `j`). Without the twist
`H_m(j) = 2^-m C(m,j)` -- this is why the Fourier note's identity (iii) is a Pascal identity. The twist field has
the two local relations `omega_(m-1,j) = omega_(m,j)^2` (odometer) and, for `q = 2^s + 1`,
`omega_(m,j+1) = omega_(m,j) omega_(m-s,j)` (`s = 1` for `q = 3`, `s = 2` for `q = 5`); for `q = 2^s - 1`,
`omega_(m,j+1) = omega_(m-s,j) conj(omega_(m,j))`; for `q = -3`, `omega_(m,j+1) = conj(omega_(m,j) omega_(m-1,j))`.

**Column recursion (what the script computes).** Column `j` is a function of the depth `m` only; it is
computed from column `j-1` by the one-dimensional recursion above, to depth `M_max = 5N`, with a per-column
scale (so there is no over- or underflow and no window truncation `A`). The neglected flux is at most
`P(S_j > 5j) ~ e^(-0.964 j)` against `|f_j(0)| ~ e^(-0.5625 j)`. Cost `5 N^2` twists per sample; every level
`j <= N` is obtained at the same `R` (Haar invariance makes each level's law the uniform-start law).

**Twist models.** `tower q`: `x_j = R q^j`; `iid`: independent Haar columns (the i.i.d.-digit model of E5);
`randmult`: `x_j = R Q^j` with a random odd 40-bit `Q`; `collatz`: the exact coefficient `mu_hat_N(u)` of the
3-adic Syracuse law, whose level-`n` phase `(u 2^-d mod 3^n)/3^n` equals `frac((-u 3^-n mod 2^d)/2^d) + (u mod
3^n)/(2^d 3^n)` (identity checked exactly for `u <= 13`, `n <= 8`, `d < 30`; the correction is the term the
uniform-start model drops, negligible except at the lowest levels).

**Validation.** (i) Brute-force path sums at `N = 6` (valuations `<= 14`) agree with the triangle to the
truncation error (`0.0210230` vs `0.0210207`, `q = 3`; `0.0284779` vs `0.0284789`, `q = 5`; cap error
`4 10^-4`). (ii) The exact-Collatz mode reproduces the Fourier note's E1 values `3^(n/2)|mu_hat_n(1)| = 1.341,
2.043, 3.882, 7.221, 10.377, 2.624` at `n = 5, 17, 127, 129, 131, 135` as `1.3409, 2.0434, 3.8819, 7.2211,
10.3773, 2.6235`, and `|mu_hat_7(7)| = 0.0129696632` against `0.0129696631` -- the window recursion of October 4
(with `A = 40`) is an independent implementation. (iii) A second kernel implementing the level/window recursion
directly (with a truncation `A`) agrees with the column kernel to ten digits on the same digit strings
(`log|f_60| = -33.9863248230` both).

## 2. The Lyapunov exponent

### 2a. Values (OBSERVED; `collatz_pascal_tower_lyapunov_20261005.py`)

`lambda := lim (1/N) log|f_N(0)|`; estimator `(L_N - L_(N0))/(N - N0)` per seed with `L_j = log|f_j(0)|`,
mean and standard error over seeds; the level `N0` cuts the transient (2b).

| model | `N`, seeds | levels | `lambda` | rate `e^lambda` | Jensen gap `lambda + (log 3)/2` |
|---|---|---|---|---|---|
| step-1 tower, `q = 3` | 6000, 100 | 1000..6000 | `-0.56229 +- 0.00008` | `0.56990 +- 0.00005` | `-0.01298` |
| step-1 tower, `q = 3` | 3000, 200 | 1000..3000 | `-0.56241 +- 0.00010` | `0.56983 +- 0.00006` | `-0.01310` |
| **step-1, combined** | | | **`-0.56234 +- 0.00006`** | **`0.56988 +- 0.00004`** | `-0.0130` |
| step-2 tower, `q = 5` | 3000, 200 | 500..3000 | `-0.55724 +- 0.00008` | `0.57279 +- 0.00005` | `-0.00794` |
| i.i.d. digits | 3000, 200 | 500..3000 | `-0.55749 +- 0.00008` | `0.57265 +- 0.00004` | `-0.00818` |
| random odd multiplier | 3000, 100 | 500..3000 | `-0.55753 +- 0.00013` | `0.57262 +- 0.00007` | `-0.00823` |

The three non-step-1 models agree within one standard error (`0.5727 +- 0.0001`); the step-1 tower is
`0.0049 +- 0.0001` lower in the log, i.e. `0.0028` in the rate. The per-seed standard deviation of the
estimator over `2400` levels is `0.0013` (`q = 3`) and `0.0011` (others); the increment variance per level is
`0.73` (`D = 1`), falling to `0.005` per level for blocks of `1000` levels, i.e. the increments are strongly
anticorrelated over short ranges (a cancellation at one level is followed by partial recovery) and the
asymptotic variance per level is about `0.004-0.005` for every model.

### 2b. The transient (OBSERVED)

Local rates `exp(mean increment)` over windows of `500` levels, standard errors `0.0002`:

| levels | 0-500 | 500-1000 | 1000-1500 | 1500-2000 | 2000-2500 | 2500-3000 | 3000-3500 | ... | 5500-6000 |
|---|---|---|---|---|---|---|---|---|---|
| `q = 3` (N = 6000) | 0.5684 | 0.5692 | 0.5700 | 0.5698 | 0.5697 | 0.5697 | 0.5699 | 0.5700, 0.5702, 0.5704, 0.5696 | 0.5698 |
| `q = 3` (N = 3000) | 0.5682 | 0.5693 | 0.5696 | 0.5700 | 0.5698 | 0.5699 | | | |
| `q = 5` | 0.5719 | 0.5726 | 0.5729 | 0.5728 | 0.5729 | 0.5728 | | | |
| i.i.d. | 0.5714 | 0.5724 | 0.5727 | 0.5726 | 0.5730 | 0.5727 | | | |

The step-1 tower needs about `1000` levels to reach its asymptotic rate; the others about `500`. The October 4
estimates on `300..1500` (`0.5698`, `0.5726`) mixed the transient in; the asymptotic values are `0.5699` and
`0.5727`. The windows `4000..5000` of the `N = 6000` run sit `2` standard errors high and the next window
`2` low: there is no trend past level `1000` at the `0.0002` level.

### 2c. The `q`-landscape (OBSERVED; level 1000, 60 seeds, levels 200..1000, standard errors `0.0002`)

| `q` | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 | 25 | 27 | 29 | 31 | 33 | 63 | 65 | 127 | 129 | 257 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rate | 0.5068 | 0.5693 | 0.5727 | 0.5720 | 0.5724 | 0.5722 | 0.5724 | 0.5726 | 0.5729 | 0.5728 | 0.5728 | 0.5721 | 0.5726 | 0.5722 | 0.5722 | 0.5721 | 0.5725 | 0.5726 | 0.5723 | 0.5722 | 0.5722 | 0.5725 |

and `q = -3`: **`0.5647 +- 0.0002`** (the alternating-sign step-1 tower, `x_j = R (-3)^j`, i.e. the `3x-1`
side). [[PENDING-NEGQ]] At this precision every positive `q >= 5` is within `0.0007` of the universal value
(the `q = 7, 23, 27, 29, 31, 127, 129` entries at `0.5720-0.5722` are two to three standard errors low on
levels `200..1000`, where the transient is not yet over). [[PENDING-Q3000-TABLE]]

`q = 1` is the degenerate control: `Phi_N = 0.epsilon_1 epsilon_2 ... epsilon_(S_N)` is the coin string read as
a binary fraction, the law of `Phi mod 1` is uniform at every coarse scale, and the transform at a Haar
frequency is smaller by `0.063` per level in the log than for any `q >= 5`, with the SAME exact second moment
`3^-N` (injectivity holds for `q = 1`): a Jensen gap of `0.13` per level, the extreme of the multifractal
picture of 3b.

## 3. Anatomy of the excess

### 3a. Per-level increments (OBSERVED; levels `>= 1000`, 200 seeds)

Increments `Delta_j = L_j - L_(j-1)` at the same Haar `R` (in the window formulation of October 4 consecutive
levels use `R_n = -u 3^-n`, so this decomposition is formulation-specific; the mean is not).

| model | mean | sd | skew | quantiles 0.1%, 1%, 5%, 25%, 50%, 75%, 95%, 99%, 99.9% |
|---|---|---|---|---|
| `q = 3` | `-0.5624` | `0.856` | `-0.006` | `-3.88, -2.78, -1.96, -1.065, -0.560, -0.060, 0.83, 1.64, 2.80` |
| `q = 5` | `-0.5572` | `0.869` | `-0.009` | `-3.97, -2.79, -1.97, -1.070, -0.556, -0.043, 0.85, 1.67, 2.84` |
| i.i.d. | `-0.5574` | `0.869` | `-0.007` | `-3.96, -2.79, -1.97, -1.069, -0.555, -0.045, 0.85, 1.68, 2.82` |

Contribution of each band to the mean, `q = 3` minus i.i.d.: `(-inf, -3)`: `+0.0017`; `[-3, -2)`: `+0.0008`;
`[-2, -1)`: `+0.0016`; `[-1, -0.5)`: **`-0.0033`**; `[-0.5, 0)`: `-0.0008`; `[0, 0.5)`: `-0.0005`;
`[0.5, 1)`: `-0.0016`; `[1, 2)`: `-0.0016`; `[2, inf)`: `-0.0013`; total `-0.0051`. **The step-1 tower's
increment distribution is slightly narrower, its left tail is lighter, and its deficit comes from the
moderate-negative bulk and from fewer large positive increments (fewer recoveries).** "Heavy-tail effect" is
false at the level of per-level increments.

### 3b. The pressure and the frequency tail (FINITE-EXACT samples; `collatz_pascal_tower_pressure_20261005.py`)

`2 10^5` samples per model at level `100` (`10^5` at `200`, `2 10^5` at `50`); `P_N(s) = (1/N) log E_R
|f_N|^(2s)`; exact control `P(1) = -log 3 = -1.0986` (the empirical values `-1.095 .. -1.100` are within the
sampling noise of a mean carried by `0.1%` of the samples).

| level 100 | `P(0.25)` | `P(0.5)` | `P(0.75)` | `P(1)` | `P(1.25)` | `P(1.5)` | share of `E|f|^2` in top 0.1% / 1% / 10% |
|---|---|---|---|---|---|---|---|
| `q = 3` | `-0.2832` | `-0.5619` | `-0.8343` | `-1.0954` | `-1.3466` | `-1.5947` | **`0.771 / 0.885 / 0.974`** |
| `q = 5` | `-0.2806` | `-0.5575` | `-0.8307` | `-1.0997` | `-1.3647` | `-1.6264` | `0.304 / 0.585 / 0.885` |
| `q = 7` | | | | | | | `0.453 / 0.690 / 0.919` |
| `q = 9` | | | | | | | `0.273 / 0.563 / 0.877` |
| i.i.d. | | | | | | | `0.341 / 0.632 / 0.903` |

[[PENDING-PRESSURE-ROWS]] At level `50` the shares are `0.295 / 0.545 / 0.842` (`q = 3`), `0.180 / 0.400 /
0.763` (`q = 5`), `0.154 / 0.391 / 0.769` (i.i.d.); at level `200`, `0.851 / 0.951 / 0.994` for `q = 3`
[[PENDING-200]]. **Reading.** All models share `E|f_N|^2 = 3^-N` exactly; the step-1 tower distributes it over
far fewer frequencies, and the concentration grows with the level. Its pressure is below the others for
`s < 1` (lower typical value), equal at `s = 1`, above for `s > 1` (heavier high moments). This is the
"heavy-tail effect" that is true: the Jensen gap `lambda - (-(log 3)/2)` is the price of a second moment
carried by resonant frequencies, and the step-1 tower pays `0.013` per level against `0.008`.

### 3c. Truncating the valuations (OBSERVED; window kernel, `N = 1500`, 100 seeds)

[[PENDING-A]]

### 3d. Changing the coin (OBSERVED; biased valuations `c ~ geometric(p)`, weights `p (1-p)^(c-1)`)

[[PENDING-P]]

### 3e. Exact pair-coherence census (FINITE-EXACT; `collatz_pascal_tower_pair_coherence_20261005.py`)

For `N = 6`, valuations `<= 8` (`97.7%` of the weight), `C_v = sum_(classes) mass^2` over the classes of the
first `v` binary digits of `Phi mod 1` = the mean square of `f_6` over the low frequencies `R < 2^v`; for the
uniform law `C_v = 2^-v C_0`. Values of `C_v / 3^-6` for `v = 1..9`:

| `q` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| uniform (`q = 1`) | 347.8 | 173.9 | 86.96 | 43.49 | 21.75 | 10.88 | 5.54 | 3.00 | 1.84 |
| `q = -3` | **356.3** | **195.3** | **100.2** | **51.3** | **26.0** | **13.3** | **7.1** | 3.97 | 2.41 |
| `q = 3` | **352.6** | **176.5** | **89.8** | **45.1** | **22.9** | **12.0** | **6.5** | 3.71 | 2.30 |
| `q = 5` | 347.8 | 173.9 | 87.3 | 44.0 | 22.5 | 11.7 | 6.2 | 3.50 | 2.19 |
| `q = 7` | 348.0 | 174.0 | 87.5 | 44.1 | 22.4 | 11.5 | 6.2 | 3.50 | 2.19 |
| `q = 9` | 348.2 | 174.3 | 87.3 | 43.9 | 22.1 | 11.6 | 6.1 | 3.57 | 2.20 |
| `q = -5` | 347.8 | 177.2 | 89.1 | 44.9 | 22.8 | 11.7 | 6.3 | 3.58 | 2.25 |
| `q = 11, 13, 15, 17` | 347.8-349.6 | 174.2-174.8 | 87.2-87.7 | 44.1-44.2 | 22.5-22.6 | 11.7-11.8 | 6.1-6.4 | 3.54-3.68 | 2.20-2.25 |

The law of `Phi mod 1` is uniform at the coarse scales for `q = +-1` and for every `|q| >= 5` to within `0.3%`
at `v <= 3`, and visibly concentrated for `q = +-3` (`+1.4%`, `+1.5%`, `+3.3%` at `v = 1, 2, 3` for `q = 3`;
`+2.4%`, `+12%`, `+15%` for `q = -3`). The ordering `-3 > 3 > rest` is the ordering of the excess. The dual
profile `D_v` (classes of the digits beyond `v`, i.e. the mean square over the frequencies `2^v R'`, the
weighted mass of path pairs with `2^v (Phi_a - Phi_b)` integral) is **the same for every `q`** (`1.000, 1.000,
1.375, 2.125, 3.461, 5.804, 9.92, 17.7, 34.3, ...` for `v = 0, 1, 2, ...`, identical to three decimals for all
thirteen multipliers): two paths of the same level `N` agree in the digits beyond `v` iff they have the same set
of ones beyond position `v` (same deep set implies the same number of shallow ones, so the multiplier powers
cancel; PROVED by the deepest-differing-depth argument of 4k). So the dyadic shells of the frequency carry no
`q`-information at all; what differs between multipliers is only the coarse-scale shape of the circle law.

## 4. Mechanism (CONJECTURED; decisive test recorded)

The real backward map `T_c(y) = (1 + q y)/2^c`, `c ~ geometric(p)`, has Lyapunov exponent
`gamma(q,p) = E log|q 2^-c| = log|q| - (log 2)/p`: for `p = 1/2`, `gamma = log|q| - 2 log 2`, negative
(contracting, the law of `Phi_N` converges to a stationary law with Kesten tail exponent `kappa = 1`, since
`E (3 2^-c)^kappa = 3^kappa/(2^(1+kappa) - 1) = 1` at `kappa = 1`) for `|q| <= 3` and positive (the real values
spread like `(|q|/4)^N`) for `|q| >= 5`. The census of 3e shows the circle law concentrated exactly in the
contracting cases with `|q| = 3` (for `|q| = 1` the stationary law happens to project to the uniform law).

**Conjecture (mechanism).** The typical Fourier decay at a Haar frequency is governed by the real-line
dynamics of the backward map: a contracting multiplier gives a stationary, non-uniform circle law whose
transform at random integer frequencies is typically smaller (and at rare frequencies larger) than that of a
spreading law, which behaves like the i.i.d. model; the sign of `q` selects the stationary law (`q = -3`:
fixed points `1/(2^c + 3)`; `q = 3`: `1/(2^c - 3)`, the rational cycles of the forward map), and `-3` is the
more concentrated one. The exact second moment `3^-N` is blind to this (injectivity), the typical value is not.

**Decisive test.** Change the coin: the threshold is `|q| < 2^(1/p)`. Predictions: `q = 3` loses its excess
for `p >= 0.64` (`gamma(3, 0.64) = +0.016`) and keeps it for `p <= 0.62`; `q = 5` acquires an excess for
`p <= 0.43` (`gamma(5, 0.43) = -0.003`) and `q = 7, 9` for `p = 0.3` (`gamma = -0.37, -0.11`); the i.i.d. model
at the same `p` is the reference at every `p`. Outcome: [[PENDING-THRESHOLD]]

**What this does and does not touch.** Nothing above uses the digits of `3^-n`; the integer-start Collatz
family (unit `u` fixed, `x_n = -u 3^-n`) [[PENDING-UNITS]]. The ridges of HYP-9166 (coincidences `u 2^Q = -+1
mod 3^k`) are a property of integer starts and are not the `0.0049`: the uniform-start tower has no such
coincidences and the full excess.

## 5. Families whose dependency rules close their union, with a common decreasing rank (typed)

**The object.** A family is a set `F` of positive odd integers with a dependency rule `n -> d(n)` certified by
a common-future receipt `U^r(n) = U^s(d(n))` (finite-seed kernel, section 2; partition cover, (1)); a rank is a
map `rho` to a well-ordered set with `rho(d(n)) < rho(n)`. A union `D = union F_i` is **closed** when
`d_i(n) in D` for every `n in F_i`, and then `Root(seeds) => Root(D)` by induction on `rho`, across family
changes included, provided ONE rank works for every rule (the finite-seed kernel's "decreasing rank" and the
bootstrap note's "destination obligation"). The hostile `D = {n = 1 mod 4}` (`U(n) < n` but `9 -> 7` leaves
`D`) is exactly a failure of closure, not of rank.

**With the integer as the rank, the leak of every bounded-depth portfolio is the left tail of the Pascal
tower's coin (PROVED, elementary).** The valuation word `(c_1, ..., c_k)` of the first `k` odd steps has density
`2^(-S_k)`, `S_k = c_1 + ... + c_k` (Terras), so the first `k` odd steps give an affine map with coefficient
`3^k / 2^(S_k)` and the class pays (all members but at most one bounded head, THM-4512) iff `S_k > k log_2 3`.
The density of odd `n` not paid within `k` odd steps is therefore exactly
```
P(S_k <= floor(k log_2 3)) = P( Bin(floor(k log_2 3), 1/2) >= k ),
```
(`7.48 10^-2, 1.19 10^-2, 5.21 10^-4, 1.33 10^-6, 1.61 10^-11` at `k = 20, 50, 100, 200, 400`), and by Cramer's
theorem for the geometric renewal (`Lambda(t) = log(e^t/(2 - e^t))`, `I(a) = sup_t (a t - Lambda(t))`,
`a = log_2 3`, optimum at `e^(t*) = 2(a-1)/a`):
```
I(log_2 3) = log 3 + (a-1) log(a-1) - a log a = 0.054979...,
e^(-I) = a^a / (3 (a-1)^(a-1)) = 0.946505...   (= 3^(h*-1) with h* = 1 - I/log 3 = 0.949956).
```
This is the "no-descent rate `3^(h*-1) = 0.9465`" of the renewal notes, now in closed form, and it is a
statement about the same coin sequence `epsilon` whose Fourier transform `E_epsilon e(R Phi(epsilon))` is the
tower: descent coverage is the lower large deviation of `S_k`, the cold-frequency rate is its Fourier dual.
The two are not the same problem (one is a tail, one a transform), and the `0.0049` excess is not visible in
the tail (it has no `q`-dependence at all).

**What closes and what does not (CITED + typed).**
* The [partition-cover theorem](collatz_partition_cover_20261004.md) (mixed congruence obstruction) gives, for
  every pair of depth bounds, a whole arithmetic cell none of whose members has a smaller join within those
  depths: no bounded-depth union is closed, whatever the guards.
* The depth-growing families close within themselves while their fuel lasts: the `H` kernel on `n = 155 mod
  2048` repeats exactly `floor((v_2(295 n - 669) - 1)/10)` times, each step strictly decreasing the integer and
  staying on the cylinder; the exit `c = H^m(n) >= 111` is a destination obligation outside the family
  ([recursive dependency kernel](collatz_recursive_dependency_kernel_20261004.md), section 2). The
  [refuel tree](finite_seed_refuel_tree_20261005.md) closes an injective countably branching family from the
  single seed `3` with the rank `6d + 2` -- a closed union with a common rank, but not a union of residue
  classes and with no density.
* Hence the honest status of the directive: a union closed under its own rules with one decreasing rank exists
  (the refuel tree) and has density zero; every union of residue classes with bounded-depth rules leaks at the
  rate `0.9465` per odd step; and unbounded-depth families leak at fuel exhaustion. Assembling more families
  changes the constant, not the shape. No closure theorem is added here.

**`p`-ary trees.** The only `p`-ary recursion in these objects with `p` an odd prime is the ternary fuel of the
refuel constructor (`4` has order `3^6` in `(1 + 3Z)/3^7`, so the branch index `kappa(u)` lives on a 3-adic
tree). The Pascal tower is binary (the coin) with a `q`-adic multiplier; its "step" is the binary length of
`q - 1`, and the measured excess does not follow that length (`q = 7 = 2^3 - 1` and `q = 9 = 2^3 + 1` are
universal, `q = -3` exceeds `q = 3`). The structure that the excess follows is the real contraction of the
multiplier (section 4), not a `p`-ary depth.

## 6. `{2, 3, 11}`, the eta product, and the pasted items (typed)

`eta(tau)^2 eta(11 tau)^2 = q - 2q^2 - q^3 + 2q^4 + q^5 + 2q^6 - 2q^7 - 2q^9 - 2q^10 + q^11 + ...` is the
weight-two newform of `X_0(11)` (CITED in the level-eleven note). What the repo has established about it:
Frobenius-at-2 on the 3-torsion realizes the golden `F_9` clock (`phi` of order `8`); a marked modular dessin
gives the Paley tournament on `F_11`; the denominator-9 fan return of the centroid note has an index-2-lattice
cube root with the local factor `1 + T + 3T^2` of `a_3 = -1` -- all typed there as local coordinate maps with
loss boundaries, none as a conjugacy of integer Collatz maps, and the level `11` is NOT derived from any Collatz
object (the `LB = 4 + 7` cost and the odd-square-bracket set `{2, 3, 11}` are separate coincidences, said so
there). **None of this enters the Pascal tower.** The tower is determined by the coin, the multiplier `q` and
the Haar frequency; its exact second moment is the collision count `3^-N`, its typical rate is the Lyapunov
exponent above, and the only `11` in this note is the multiplier `q = 11`, whose rate is the universal one. A
connection would have to map the stationary law of `y -> (1 + 3y)/2^c` (section 4) to a modular object; no
such map is known here, and the Kesten tail exponent `kappa = 1` of that law is a renewal fact, not a modular one.

The pasted items are repo results and are cited, not re-derived: translation modulo `9` identifying the last
letter of the controller and the single-translation recovery of a word ([translation decoder](translation_phase_decoder_20261005.md));
canonical source/target representatives `c, d` with the Fourier summand recording `d mod 3^m` and the payment
comparing the actual source with its endpoint ([carry interfaces](collatz_carry_interfaces_20261004.md),
[two clocks](source_address_two_clocks_20261004.md)); the fixed centroid encoding a finite refinement word of
exact denominator `3^(len+1)` and the period-three decoding cycle of denominator `9` ([centroid
membership](centroid_membership_20261005.md)); the golden field `F_9 = Z[phi]/(3)` with `phi` of order `8`
([five-eight-nine transfer](five_eight_nine_transfer_20261005.md), section 2); the eight fixed-path
presentations of four-vertex tournaments collapsing to four classes with the strong class of five
presentations and five Hamiltonian paths, and path count nine for the joined diamonds ([tournament
codec](tournament_recursive_four_state_20261004.md)); and `G(n) = (9n + 5)/8` with `G(n) + 5 = (9/8)(n + 5)`,
the actual word `(1, 2)` with anchor `-5` and the credit account `E(x, k) = (x + 5)(9/8)^k` ([paid guard
budget](paid_guard_budget_20261005.md)). The Pascal tower's analogue of the anchor is the fixed point
`1/(2^c - 3)` of the backward map at valuation `c` (the rational cycles), which is where the sign of `q`
enters the stationary law of section 4.

## 7. Reproduction

```bash
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --selftest                     # 35 s
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model tower --q 3 --N 3000 --seeds 200 --workers 10   # 34 s
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model tower --q 3 --N 6000 --seeds 100 --workers 10   # 90 s
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model iid --N 3000 --seeds 200 --workers 10
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model tower --q -3 --N 1500 --seeds 100 --workers 10
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model tower --q 3 --N 1500 --seeds 100 --A 4            # truncation
python 04-computation/experiments/collatz_pascal_tower_lyapunov_20261005.py --model tower --q 3 --N 1500 --seeds 100 --p 0.66        # biased coin
python 04-computation/experiments/collatz_pascal_tower_pressure_20261005.py 100 200000 10                                            # 6 min
python 04-computation/experiments/collatz_pascal_tower_pair_coherence_20261005.py 6 8                                                # 90 s
```
The batch logs are kept as `collatz_pascal_tower_lyapunov_20261005_batch*.out` beside the scripts (seeds
`0, 1, ...` per run; `numpy.random.default_rng`).

## 8. Verdicts

| claim | status |
|---|---|
| the uniform-start `q`-tower is a twisted Pascal triangle; `f_N(0)` is the transform at a Haar frequency of the random backward Syracuse law | PROVED (identity; brute-force, E1 and two-kernel checks) |
| step-1 tower rate `0.56988 +- 0.00004`; step-2, i.i.d., random multiplier `0.5727 +- 0.0001` | OBSERVED (levels to `6000`, `200` seeds) |
| the excess is `0.0049 +- 0.0001` per level in the log; not a transient (flat from level `1000` to `6000`) | OBSERVED |
| `q = -3`: `0.5647`; `q = 1`: `0.5068`; positive `q >= 5`: universal within `0.0007` at level `1000` | OBSERVED |
| the step-1 deficit is a bulk effect of the increments (lighter left tail, fewer recoveries) | OBSERVED |
| the second moment of the step-1 tower is carried by the top `0.1%` of frequencies (`77%` at level `100`, `85%` at `200`) against `30%` for step 2 | FINITE-EXACT samples |
| the circle law of `Phi` is concentrated at coarse scales exactly for `q = +-3` (`N = 6`, exact) | FINITE-EXACT |
| mechanism: real contraction `|q| < 2^(1/p)` of the backward map | CONJECTURED; test [[PENDING-THRESHOLD-VERDICT]] |
| no-descent rate `a^a/(3 (a-1)^(a-1)) = 0.946505`, `a = log_2 3`; bounded-depth unions never close | PROVED (Cramer) + CITED (partition cover) |
| `{2, 3, 11}` / `eta(tau)^2 eta(11 tau)^2` enter the tower | NOT FOUND (no object in common beyond the integer `11` as a multiplier) |
| an exact value of `lambda_3`; H1; Collatz | OPEN |
