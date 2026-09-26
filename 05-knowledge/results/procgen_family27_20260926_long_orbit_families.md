# Numbers beyond 27: five families, the rates at which they occur, and the one Moran function of the owner's backward recursion that governs them

**Status.**
* **PROVED** (full hand proofs in §§1, 3, 4; every finite ingredient is re-checked by the runner):
  1. **Proposition M (one function).** The Moran (Malthusian) function of the owner's backward recursion `A <- 2A`, `A <- (2A-1)/3` is `g(s) = 2^-s + (1/3)(3/2)^s`, and `g(s) = phi(s-1)` with `phi(t) = ((3/2)^t + (1/2)^t)/2` the pressure of the forward parity walk (this identity is Lagarias-Weiss's duality `M_BP(t) = M_RRW(t+1)`, CITED). `g` is convex and `g = 1` exactly at `s = 1, 2`; `-g'(1) = ln(2/sqrt 3)`, `g'(2) = (3/4)ln 3 - ln 2`; `min g = 2^-(1-h)` at `s = 1 + lambda*`; and `beta := max_s (-ln g(s))/s` equals `ln 2 - p* ln 3` where `p*` is the unique root of `H(p) = p ln 3` on `(1/2, 1)`. Numerically `1/beta = 41.6776476556544`, `p* = 0.6090898`.
  2. **Theorem R (rise law, words).** For the fair parity walk and `M_j = 3^(o_j)/2^j`, `tau = min{j : M_j >= W}`, every `k` and `W > 1`: `P(tau <= k) = (1 - eps_k(W)) / E[M_tau | tau <= k]` exactly, with `W <= E[M_tau | tau <= k] < 3W/2` and `eps_k(W) = E[M_k; tau > k] <= min_(0<=t<=1) W^(1-t) phi(t)^k`. Hence `2/(3W) < P(W) := P(sup_j M_j >= W) <= 1/W`. The exponent `-1` is the root `t = 1` of `phi`, i.e. the second root `s2 = 2` of `g`: AM-fairness `3 + 1 = 4`.
  3. **Theorem R' (rise law, integers).** For every block `[2^k, 2^(k+1))` and `W > 1`: `N_k(W) <= #{n : max_(j<=k) T^j n >= W n} <= N_k(W - (3/4)^k)`, `N_k` the exact word count; the block density tends to `P(W)`. **Corollary:** for every `W > 1` the set `{n : t(n) >= W n}` (`t(n)` = trajectory maximum) has **lower density `>= P(W) > 2/(3W)`**.
  4. **Theorem S (window rise spectrum).** `#{n in [2^k, 2^(k+1)) : max_(j<=k) T^j n >= n^beta} = 2^(k E_win(beta) + o(k))` with `E_win(beta) = 2 - beta` on `[1, 1 + zeta]` and `h(beta/log_2 3)` on `[1 + zeta, log_2 3]`, `zeta = (3/4)log_2 3 - 1 = 1 - h(3/4) = 0.188722`. The line `2 - beta` is the **slope -1 tangent** to the dip-spectrum curve `F(x) = h(x/log_2 3)` of THM-4487, touching at `x = (3/4)log_2 3` and hitting `0` at `beta = 2` (the Lagarias-Weiss path-record exponent).
  5. **Theorem G (glide, exact).** For `L - 1 <= k log_3 2`: `#{n in [2^k, 2^(k+1)) : glide_T(n) >= L} = 2^(k-L+1) W_(L-1)` exactly, `W_m = |Bad_m|`.
  6. **Theorem D (rate = covering number = dimension).** The number of residue classes mod `2^m` meeting the closed exceptional set `Bad = {x in Z_2 : 3^(o_j) > 2^j for all j}` is exactly `W_m`; with Theorem G, the block count of the long-glide family is `2^(k-m)` times the `2^(-m)`-covering number of `Bad`, so its occurrence exponent is `dim_B Bad = dim_H Bad = h(log_3 2) = 0.949956 = 1 + log_2 min g`. The W-riser set `{sup_j M_j >= W}` has Haar measure `P(W)`, and the infinite-riser set `{sup_j M_j = infinity}` has measure 0 and Hausdorff dimension `h(log_3 2)`.
  7. **Proposition B (exponent-branching identity).** If a backward-invariant set of integers (e.g. any backward tree) has a regularly varying counting function of index `s > 0` and a limiting fraction `kappa` of members `= 2 mod 3`, then `2^-s + kappa (3/2)^s = 1`; so `kappa = 1/3` forces `s in {1, 2}` and `s = 1` forces `kappa = 1/3`.
  8. **Proposition F (record count in the Frechet-scale model).** If `t(n)` are independent with `P(t(n) <= y) = exp(-C n^theta y^-theta)`, then `P(n is a record) = n^theta / sum_(m<=n) m^theta`; for `theta = 1` the expected number of records up to `X` is `2(H_(X+1) - 1)`, **independent of `C`**.
* **FINITE-EXACT** (all `n <= 2^32` by an exhaustive C scan validated against an independent brute-force Python reference; all OEIS b-file terms by exact big-integer orbits):
  * all 148 A006877 / A006878, 97 A006885, 35 A060413 and A217934 values re-derived from their starting values; the scan reproduces every term `<= 2^32` of A006877, A006884, A060412; the records of `gamma(n) = sigma_T(n)/ln n` are K-L Table 1 minus its row `5649499` (a ones-ratio record, `gamma = 24.699176 < 24.714906`);
  * `rho(n) = ln t(n)/ln n` is maximal at `n = 27` over all `3 <= n <= 2^32` (`rho(27) = 2.559982`); the `n <= 2^32` with `t(n) > n^2` are exactly 21 numbers in 3 clusters, the first being `{27, 31, 41, 47, 54, 55, 62, 63}` = 27's branch below `sqrt 4616`;
  * Theorem R' checked in all 597 (block `k <= 31`, `W = 2^(i/2)`) cases: the block count **equals** the exact word count every time; Theorem G checked in all 327 admissible (block, `L`) cases up to `k = 31`;
  * **27's branch** `B27 = {n : trajectory of n meets 27's trajectory before its peak} = Pred*(3077)` has block densities in `(0.3925, 0.3929)` for every block from `2^24` to `2^31` (overall `0.392660` at `2^32`; an independent plain-Python Monte Carlo of 20000 random `n` in `[2^31, 2^32)` gives `0.3941 +- 0.0035` against the scan's `0.392634` for that block), 87 times the Haar-model value `2R/3077 = 0.0045`; its members are equidistributed mod 3 to `3e-4`; densities of the nested trees along 27's orbit, of the owner's `4x` ladder rungs and of the trunk `(4^j-1)/3` (e.g. `Pred*(5)` has density `0.9379`, `Pred*(16)` `0.0621`);
  * of the OEIS records, 104/148 delay records, 36/98 path records and 13/35 glide records lie in `B27`.
* **EMPIRICAL:** the full-orbit W-riser density equals `P(W)` within 0.3% at `2^31` for `W <= 2^10` (the window alone gives only 20% there); numbers with delays like 27's (`sigma_T >= 21.24 ln n`) grow with exponent `0.656` (Lagarias-Weiss model `0.6688`); path records number 98 up to `2.36e21` against `2(H - 1) = 97.6`; glide records satisfy `glide/L* in [0.80, 1.16]` against the exact `W_k` model (27 is the outlier, 3.7 times its model value); the conditioned glide excursion is Brownian (`E u = sqrt(pi/2) - 0.885/sqrt L`) and 27's is typical; 27's branch holds 66-73% of the numbers with `sigma_T(n) >= 8..24 ln n` although it has density 39%.
* **VERIFIED** (numerical reproduction of cited constants from `g` alone): K-L's `41.677647`, `beta_BP = 0.02399`, ones-ratio `0.609091`, typical slope `6.95212`, time-to-peak slope `7.645`; THM-4487's tilt `0.488077` and `1 - h = 0.050044`; the Brownian excursion maximum's mean `sqrt(pi/2)` and second moment `pi^2/6`.
* **CITED:** Lagarias-Weiss 1992 and Kontorovich-Lagarias 2009/10 (arXiv 0910.1944, read in full: the RRW and branching-random-walk models, `gamma = 41.677647`, `rho = 2`, the spectra `x^(1 - a g(1/a))` and `x^(2-beta)`, the duality); OEIS A006877, A006878, A006884, A006885, A060412, A060413, A217934 (b-files, sha256 in the .out); Terras 1976 (parity bijection); Chung 1976 / Kennedy 1976 (Brownian excursion maximum); Nevzorov's `F^alpha`-scheme (independent record indicators; from memory); the renewal theorem for the existence of `lim W P(W)`; in-repo THM-4476, THM-4480, THM-4487, THM-4495 and the choice-ladder and inverse-tree notes.
* **OPEN:** an upper bound `#{n <= X : t(n) >= Wn} <= (P(W) + o(1)) X` (it would give density zero to the divergent integers); the existence of the density of any single backward tree, including `B27`; the limsup constants `rho = 2`, `gamma = 41.68`, glide `19.98 log_2 n` (model predictions, not theorems).
* **CORRECTED / REFUTED:** the Haar tree law `c R/a` fails for `B27` by a factor 87 and along 27's orbit by orders of magnitude (§4.3); the naive "pre-peak of any earlier record" genealogy makes 145/148 delay records "descendants" and is useless (§6).
* **NUMEROLOGY:** `delta(B27) = 0.39266` is within `4e-5` of `pi/8 = 0.392699`; the last blocks drift from `0.39274` to `0.39263` and no mechanism is known. Recorded only so that nobody builds on it.
* **Not claimed:** anything about Collatz itself; novelty of the heuristic dictionary (it is Lagarias-Weiss's); that any single tree has a density.

Session `collatz-procgen-20260922`, lane "family27", 2026-09-26. No HYP or THM file was created.
* **Scripts:** `04-computation/experiments/procgen_family27_20260926_{run,theory,bfiles,reference}.py` and `procgen_family27_20260926_scan.c`.
* **Output:** [`procgen_family27_20260926.out`](procgen_family27_20260926.out) (one runner; every claim is a `check` that raises).

## 0. The question, and the answer in brief

The owner: *"think about numbers beyond 27 in the same family, probably with very long collatz orbits or some similar property, and see how the rate at which that family occurs in the naturals governs fractal recursion".*

Throughout `T(x) = x/2` (even), `(3x+1)/2` (odd), trajectories stop at 1; `sigma_T(n)` is the total stopping time in `T`-steps, `t(n)` the trajectory maximum, `glide_T(n)` the first `j` with `T^j n < n`, `o_j` the number of odd steps among the first `j`, `M_j = 3^(o_j)/2^j`, `h` the binary entropy, `alpha = log_2 3`. For 27: `sigma_T = 70` (41 odd steps; 111 standard steps), `t = 4616` at step 45 (9232 standard), glide 59 (96 standard).

| family | 27's place | how often it occurs | status | feature of `g` |
|---|---|---|---|---|
| (a) delay records `A006877`; "delays like 27": `sigma_T(n) >= 21.24 ln n` | 9th record; the big jump (23 -> 111 standard steps) | records `~3.4` per unit of `ln X`; the like-27 family `X^(0.669+o(1))` | model CITED; exponent `0.656` measured at `2^32` | `beta = max_s (-ln g(s))/s = 1/41.6776` |
| (b) path records `A006884`; "rises like 27": `t(n) >= 171 n`; "super-quadratic": `t(n) > n^2` | 6th record; `rho(27) = 2.560` is the maximum of `rho` over `3 <= n <= 2^32` | like-27 risers: density `P(171) = 0.00486`; records `2 ln X`; `t > n^2`: 21 numbers up to `2^32`, in 3 tree clusters | density: lower bound PROVED, equality EMPIRICAL; record count model PROVED + EMPIRICAL | the second root `s2 = 2` |
| (c) glide records `A060412`; "glides like 27": `glide_T >= 12.41 log_2 n` | 4th record; 3.7 times the model's expected record at 27 | `2^(k-m) W_m` exactly in the window (PROVED); records follow `W_k` | PROVED + EMPIRICAL | `min g = 2^-(1-h)`: dimension `h` |
| (d) 27's branch `B27 = Pred*(3077)` | the family itself | density `0.3927`, exponent 1, one third in each class mod 3 | FINITE-EXACT; identity PROVED | the first root `s1 = 1` |
| (e) near-critical words (long glides) | a typical long glide: normalised height at the 72nd percentile | as (c) | EMPIRICAL (shape), PROVED (counts) | the Chernoff tilt `lambda*` |

**The precise sense in which the rate governs the fractal recursion** (§4.4): each family is a large-deviation family of the parity walk, its occurrence exponent is a Legendre-type transform of `ln g`, and `g` is the Moran function of the owner's recursion. Four statements are proved: the long-glide family's block counts are `2^(k-m)` times the covering numbers of the exceptional fractal, so its exponent is that fractal's dimension `h = 1 + log_2 min g` (Theorems G, D); the W-riser density `P(W)` lies in `(2/(3W), 1/W]` because `g(2) = 1` (Theorem R); the window rise spectrum is the slope -1 tangent (slope = `1 - s2`) to the dip-spectrum curve (Theorem S); a backward tree of exponent `s` has a fraction `kappa = (1 - 2^-s)(2/3)^s` of members `= 2 mod 3` (Proposition B). The owner's `4x` recursion is exact as a decomposition of sets, but the densities of individual trees are **not** governed by the averaged recursion: they are dominated by the smallest numbers the tree contains (§4.3).

## 1. One function: the Moran function of the backward recursion

The owner's inverse tree (inverse-tree note §0): every `x` has the predecessor `D(x) = 2x`, and `E(x) = (2x-1)/3` when `x = 2 mod 3`. Weighting a `D`-child by its size ratio `2` and an `E`-child by `2/3` with the 3-adic legality probability `1/3` gives the Moran (Malthusian) function

`g(s) = 2^-s + (1/3)(3/2)^s`,

the Lagarias-Weiss branching-process moment function `M_BP(-s)`.

**Proposition M.** (i) `g(s) = phi(s-1)` for all real `s`, `phi(t) = ((3/2)^t + (1/2)^t)/2`. (ii) `g` is strictly convex, `g(1) = g(2) = 1`, so `g < 1` exactly on `(1,2)`. (iii) `-g'(1) = ln(2/sqrt 3) = 0.143841` and `g'(2) = (3/4)ln 3 - ln 2 = 0.130812`. (iv) `g` is minimal at `s = 1 + lambda*`, `lambda* = log_3(ln 2/ln(3/2)) = 0.488077`, with `min g = 2^-(1-h(log_3 2))`. (v) `beta := max_(s>0) (-ln g(s))/s` is attained at a unique `s_opt` (`= 1.40368`), and `beta = ln 2 - q ln 3`, where `q = (1/3)(3/2)^s_opt / g(s_opt)` is the unique root in `(1/2, 1)` of `H(q) = q ln 3` (`H` the natural entropy). (vi) `zeta := (3/4)alpha - 1 = 1 - h(3/4) = phi'(1)/ln 2`.

*Proof.* (i) `(1/2)(1/2)^(s-1) = 2^-s` and `(1/2)(3/2)^(s-1) = (1/3)(3/2)^s`. (ii) A positive combination of exponentials is strictly convex; `g(1) = 1/2 + 1/2`, `g(2) = 1/4 + 3/4`. (iii) Differentiate. (iv) `g'(s) = 0` iff `2^-s ln 2 = (1/3)(3/2)^s ln(3/2)` iff `3^s = 3 ln 2/ln(3/2)`. At that point the tilted law `q_s = (1/3)(3/2)^s/g(s)` of the odd child satisfies `q/(1-q) = 3^(s-1) = ln 2/ln(3/2)`, i.e. `q = log_3 2`, the odd frequency of zero drift; by Cramer's theorem `inf_t phi(t) = exp(-D(log_3 2 || 1/2)) = 2^-(1-h(log_3 2))`. (v) `ln g` is strictly convex (a log-sum of exponentials with distinct rates), so `u(s) = s (ln g)'(s) - ln g(s)` has `u'(s) = s (ln g)''(s) > 0` for `s > 0`: `-ln g(s)/s` has at most one critical point on `(0, infinity)`. It is positive exactly on `(1,2)` and vanishes at 1 and 2, so its maximum is attained at a unique `s_opt in (1,2)`, characterised by `s (ln g)'(s) = ln g(s)`. With `q = q_s`, `(ln g)'(s) = q ln 3 - ln 2`, hence `beta = -ln g(s)/s = ln 2 - q ln 3`. For the entropy: `ln q = s ln(3/2) - ln 3 - ln g` and `ln(1-q) = -s ln 2 - ln g`, so `H(q) = -s(q ln 3 - ln 2) + q ln 3 + ln g = q ln 3` at `s_opt`. `H(p) - p ln 3` is strictly concave, positive at `p = 1/2` and negative near `p = 1`, so its root in `(1/2,1)` is unique, and `q(s_opt) > 1/2` because `s_opt > 1`. (vi) `h(3/4) = 2 - (3/4)alpha` and `phi'(1) = (3/4)ln(3/2) - (1/4)ln 2`. ∎

**The dictionary** (values checked in `.out` section A; interpretations as indicated):

| feature of `g` | value | what it governs | status of the interpretation |
|---|---|---|---|
| root `s1 = 1` | `g(1) = 1` | a backward tree has counting exponent 1 | model theorem (K-L Thm 6.5), K-L's `x^0.84` CITED; Prop. B |
| residue at `s1` | `R1 = 1/ln(2/sqrt3) = 6.952119` | typical `sigma_T(n) = 6.952 ln n`; Haar tree density `cR1/a` | CITED (K-L), inverse-tree note Prop. 10 |
| root `s2 = 2` (`phi(1) = 1`) | `3 + 1 = 4` | W-riser density `P(W) ~ C/W`; `rho = 2`; path-record count `2 ln X` | Theorem R, R' (PROVED); `rho = 2` CITED model; Prop. F |
| residue at `s2` | `R2 = 1/g'(2) = 7.644557` | time to reach the peak per unit of `ln(t/n)` | CITED (K-L §4.3: 7.645); EMPIRICAL median 7.07, mean 7.28 on 84 path records |
| minimum | `2^-(1-h)`, `1-h = 0.050044`, at `1 + lambda*` | dimension `h` of the exceptional fractal; glide records `log_2 n /(1-h) = 19.98 log_2 n` | Theorem D (PROVED); record constant heuristic |
| tangent from 0 | `beta = 0.0239937`, `1/beta = 41.677648` | delay records `sigma_T <= 41.68 ln n`; ones-ratio `p* = 0.609090` | CITED (Lagarias-Weiss / K-L conjecture 4.1) |
| `phi'(1)/ln 2` | `zeta = 0.188722` | the window rise spectrum switches from the tangent to the curve at `beta = 1 + zeta` | Theorem S (PROVED) |

The three computations of the delay constant (the tangent, the entropy root, and the Lagarias-Weiss fixed point `gamma g_LW(1/gamma) = 1`) agree to `1e-20`. They reproduce K-L's `41.677647`, `beta_BP = 0.02399` and ones-ratio `0.609091` (ours: `0.6090898`), which fixes K-L's convention: `gamma` counts `T`-steps.

## 2. The families, defined exactly, with the data

All statistics below are for `T`; the standard map counts an odd step twice.

**(a) Delay records** (`A006877`, values `A006878`, standard steps; b-file: 148 terms to `1.47e19`, from Roosendaal as of 2024-08-06): `1, 2, 3, 6, 7, 9, 18, 25, 27, 54, 73, 97, 129, 171, 231, 313, 327, 649, 703, ...`. The records of the ratio `gamma(n) = sigma_T(n)/ln n` over odd `n <= 2^32` are exactly `3, 7, 9, 27, 230631, 626331, 837799, 1723519, 3732423, 6649279, 8400511, 63728127`, with the values of K-L's Table 1 (`gamma(27) = 21.238915`, `gamma(63728127) = 32.943545`); K-L's row `5649499` is a ones-ratio record, not a `gamma` record. The model constant is `41.68`; the largest value in K-L's Table 1 is `36.72` (`n = 7.2e21`, not certified as a record).

**(b) Path records** (`A006884`, maxima `A006885` in the standard map; b-file: 98 terms to `2.36e21`): `1, 2, 3, 7, 15, 27, 255, 447, 639, 703, ...`. Among all 98, exactly 8 have `t(n) > n^2`: `27, 319804831, 1410123943, 3716509988199, 9016346070511, 1254251874774375, 10709980568908647, 1980976057694848447`. (K-L's Table 3 lists the other 7; the term `1.07e16`, with `t/n^2 = 1.528`, is missing from it.) `rho(27) = ln 4616/ln 27 = 2.559982` exceeds `rho(n)` for every other `3 <= n <= 2^32` (FINITE-EXACT; even `n` never exceed their odd part, since `t(2^a m) = max(2^a m, t(m))`).

**(c) Glide records** (`A060412`; glides `A060413` in `T`-steps and `A217934` standard; b-file 35 terms to `2.6e18`): `2, 3, 7, 27, 703, 10087, 35655, 270271, ...`. Up to `2^32` the `T`-glide and standard-glide records have the same 23 starting values.

**(d) 27's branch.** The integers whose trajectory meets 27's trajectory before 27's peak `4616` (standard `9232`) are those meeting one of `27, 41, 62, ..., 2051, 3077`; since every earlier orbit point flows into `3077`, this is the backward tree `B27 = Pred*(3077)` (the odd point `3077` is the same in both maps). More generally the merge index `mu(n)` = first point of 27's orbit on the orbit of `n` classifies all integers; `Pred*(T^j 27) = {n : mu(n) <= j}`.

**(e) Near-critical words.** 27's word stays above the critical line for 59 steps: its walk `S_j = ln(T^j 27/27)` climbs to `5.14` and returns. The family is `{n : glide_T(n) >= L}`; the critical band of THM-4480 (walks confined to a tube of height `~L^(1/3)`) is its thin core. §3.5 shows that 27 is a typical member.

## 3. Rates (T2)

### 3.1 Risers: the exact rise law

**Theorem R.** Let `u` be uniform on `{0,1}^k`, `M_j = 3^(o_j)/2^j`, `W > 1`, `tau = min{j <= k : M_j >= W}` (`infinity` if none), `eps_k(W) = E[M_k ; tau > k]`. Then
`P(tau <= k) = (1 - eps_k(W)) / E[M_tau | tau <= k]`,  `W <= E[M_tau | tau <= k] < (3/2) W`,  `eps_k(W) <= min_(0<=t<=1) W^(1-t) phi(t)^k`,
so `(2/3)(1 - eps_k)/W < P(tau <= k) <= (1 - eps_k)/W`, and `2/(3W) < P(W) := P(sup_j M_j >= W) <= 1/W`.

*Proof.* `M_(j+1) = M_j (3/2)` or `M_j/2` with probability `1/2` each, so `(M_j)` is a martingale: `E[M_(j+1) | F_j] = M_j (3/2 + 1/2)/2 = M_j` (AM-fairness, `3 + 1 = 4`, i.e. `phi(1) = 1`). Optional stopping at the bounded time `tau ∧ k` gives `1 = E[M_(tau∧k)] = E[M_tau ; tau <= k] + eps_k`, i.e. `E[M_tau; tau <= k] = 1 - eps_k`, which is the identity. At `tau` we have `M_tau >= W > M_(tau-1)` (`tau >= 1` since `M_0 = 1 < W`) and `M_tau <= (3/2)M_(tau-1) < (3/2)W`. On `{tau > k}`, `M_k < W`, so `M_k = M_k^t M_k^(1-t) <= M_k^t W^(1-t)` and `E[M_k^t] = phi(t)^k`. For `t in (0,1)`, `phi(t) < 1` by strict convexity, so `eps_k -> 0`. Letting `k -> infinity` (monotone convergence) gives `E[M_tau ; tau < infinity] = 1`, i.e. `P(W) = P(tau < infinity) = 1/E[M_tau | tau < infinity]`, and `W <= E[M_tau | tau < infinity] < 3W/2` gives `2/(3W) < P(W) <= 1/W`. ∎

**Theorem R'.** For `k >= 1` and `W > 1`, with `N_k(V) = #{u in {0,1}^k : max_(0<=j<=k) M_j(u) >= V}`:
`N_k(W) <= #{n in [2^k, 2^(k+1)) : max_(0<=j<=k) T^j(n) >= W n} <= N_k(W - (3/4)^k)`,
and `2^-k` times the middle term tends to `P(W)`.

*Proof.* The block is a complete residue system mod `2^k`, so `n -> (parity word of length k)` is a bijection onto `{0,1}^k` (Terras). By induction `T^j(n) = M_j n + c_j` with `0 <= c_j <= (3/2)^j - 1` (odd step: `c -> (3/2)c + 1/2`; even step: `c -> c/2`). For `n >= 2^k` and `j <= k`, `M_j <= T^j(n)/n < M_j + (3/4)^k`, which gives the two inclusions. `N_k(W)/2^k = P(tau <= k)` increases to `P(W)`; `N_k(W - (3/4)^k)/2^k <= P(W - (3/4)^k)`, which decreases to `P(W)` as `k -> infinity` by continuity of `P(sup M >= V)` from the left in `V`. ∎

**Corollary (PROVED).** For every `W > 1`, `liminf_X X^-1 #{n <= X : t(n) >= W n} >= P(W) > 2/(3W)`: at least a `2/(3W)` fraction of all integers rise by a factor `W`. *Proof.* Fix `k`. Any `2^k` consecutive positive integers form a complete residue system mod `2^k` and so realise every word of length `k` once; since `c_j >= 0`, each `n` whose word has `max_(j<=k) M_j >= W` satisfies `t(n) >= max_j T^j n >= W n`. `[1, X]` contains `floor(X/2^k)` such runs, so the count is at least `floor(X/2^k) N_k(W)`, and the liminf is at least `N_k(W)/2^k = P(tau <= k)` for every `k`; let `k -> infinity`. ∎ The matching upper bound is OPEN: it would imply that the integers with divergent trajectories have density zero, which THM-4476 explicitly does not give.

**Data** (`.out` §B, §E):

| `W` | `P(W)` | `W P(W)` | window density, block 31 | full-orbit density, block 31 | full/`P(W)` |
|---|---|---|---|---|---|
| `2` | `0.4273692` | `0.8547` | `0.41554` | `0.42738` | `1.0000` |
| `4` | `0.2034377` | `0.8138` | `0.19028` | `0.20341` | `0.9999` |
| `16` | `0.0532728` | `0.8524` | `0.04213` | `0.05329` | `1.0003` |
| `128` | `0.0065110` | `0.8334` | `0.00337` | `0.00651` | `0.9995` |
| `1024` | `8.1222e-4` | `0.8317` | `1.64e-4` | `8.1455e-4` | `1.0029` |
| `2^20` | `7.9784e-7` | `0.8366` | `0` | `8.52e-7` | `1.07` |

* `P(W)` is bracketed to relative width `1e-9` (pruned DP plus a Ville tail bound) for `W = 2^(i/2)`, `i <= 120`; `2/3 < W P(W) <= 1` in all cases, and `W P(W) in [0.8274, 0.8366]` for `W in [2^20, 2^60]`. The Cramer-Lundberg constant `C = lim W P(W)` (it exists by the renewal theorem, since `ln 3/ln 2` is irrational) is therefore about `0.83`; the slow oscillation comes from the convergents of `log_2 3`.
* Theorem R' was checked in all 597 cases (`k <= 31`, `W = 2^(i/2)`, `i <= 1.17k + 1`). The lower bound is attained **every time**: at these `W` the block counts are the exact word counts.
* The window captures only 20% of the 1024-risers, yet the full-orbit density equals `P(W)` to 0.3% (EMPIRICAL). The rise usually takes `~ R2 ln W = 53` steps, far beyond `log_2 n = 31`; the actual trajectory keeps behaving like the fair walk there.
* **27's rise:** `t(n) >= (4616/27) n` holds for a fraction `0.004857` of the top block, against `P(170.96) = 0.004858` (ratio `0.9999`; window alone `0.00205`). The "numbers rising like 27" are a positive-density family.

### 3.2 The window rise spectrum and its relation to the dip spectrum

**Theorem S.** For `1 < beta < alpha`,
`lim_k k^-1 log_2 #{n in [2^k, 2^(k+1)) : max_(j<=k) T^j(n) >= n^beta} = E_win(beta)`,
with `E_win(beta) = 2 - beta` for `beta <= 1 + zeta` and `h(beta/alpha)` for `beta >= 1 + zeta`.

*Proof.* Put `b = beta - 1`. For `n` in the block, `n^b in [2^(kb), 2^((k+1)b))` and `T^j(n)/n - M_j in [0, (3/4)^k)`, so up to a change `b -> b + O(1/k)`, which does not affect the exponent, we count words `u in {0,1}^k` with `max_(j<=k) 3^(o_j) 2^-j >= 2^(kb)`.
*Upper bound.* The number of such words is at most `sum_(j<=k) 2^(k-j) B(j, (j+kb)/alpha)` with `B(j,m) = sum_(o>=m) C(j,o) <= 2^(j h(m/j))` for `m/j >= 1/2` (here `m/j > 1/alpha > 1/2`), and `B = 0` when `m > j`. With `t = j/k` and `q(t) = (t+b)/(t alpha)`, the `j`-th term is at most `2^(k f(t))`, `f(t) = 1 - t + t h(q(t))`, so the count is at most `k 2^(k max f)`. Since `t q'(t) = -(q - 1/alpha)`, `f'(t) = psi(q(t))` with `psi(q) = -1 + h(q) - h'(q)(q - 1/alpha)`. Then `psi'(q) = -h''(q)(q - 1/alpha) > 0` for `q > 1/alpha`, and `psi(3/4) = -1 + (2 - 3alpha/4) + alpha(3/4 - 1/alpha) = 0` (using `h'(3/4) = -alpha`). Since `q(t)` decreases in `t`, `f` increases while `q > 3/4` and decreases afterwards. If `b >= zeta` then `q(1) = (1+b)/alpha >= 3/4`, the maximum is at `t = 1`, and `f(1) = h(beta/alpha)`. If `b < zeta` the maximum is at `t* = b/zeta` (where `q = 3/4`), and `f(t*) = 1 - t*(1 - h(3/4)) = 1 - b`.
*Lower bound.* If `b < zeta`: take a prefix of length `j = ceil(t* k)` with `o = ceil((j + (k+1)b)/alpha)` odd letters (`o/j -> 3/4`), in any order, followed by `k - j` free letters. This gives `C(j,o) 2^(k-j) = 2^(k(1-b) - O(log k))` words with `M_j >= 2^((k+1)b)`, and their integers satisfy `T^j n >= M_j n >= n^beta`. If `b >= zeta`: take all `k` letters with `o = ceil((k + (k+1)b)/alpha)` odd letters, giving `2^(k h(beta/alpha) - O(log k))` words (for `beta < alpha`). ∎

**Reading.** Let `F(x) = h(x/alpha)` be THM-4487's dip-spectrum curve (`F(gamma)` = exponent of the `n` staying above `n^gamma` through the window, `gamma <= 1`). Its tangent of slope `-1` touches at `x = 3alpha/4 = 1 + zeta` and is the line `2 - x`, which passes through `(1,1)` and `(2,0)`. So:
* the **window** rise spectrum is this tangent on `[1, 1+zeta]` followed by the curve `F` itself on `[1+zeta, alpha]`;
* the **Lagarias-Weiss full-orbit** spectrum `x^(2-beta)` (K-L Thm 4.4, model) is the same tangent continued to `beta = 2`, where it hits 0. That is the path-record exponent `rho = 2`.

The slope `-1` is the Cramer root `theta = 1 = s2 - 1`. Every rise rate in this note is a point of the dip-spectrum curve or of its slope `-1` tangent. The glide spectrum `1 - (1 - F(1)) c` of §3.3 uses the same curve through its value `F(1) = h`.

Corollary (window inside full orbit, PROVED): `#{n <= X : t(n) >= n^beta} >= X^(E_win(beta) - o(1))`. So the Lagarias-Weiss value `2 - beta` is proved **as a lower bound** for full trajectories when `beta <= 1 + zeta = 1.1887`; the upper bound is OPEN.

### 3.3 Glides

**Theorem G.** If `k >= 1` and `1 <= L <= 1 + k log_3 2`, then `#{n in [2^k, 2^(k+1)) : glide_T(n) >= L} = 2^(k-L+1) W_(L-1)`, where `W_m = #{u in {0,1}^m : 3^(o_j) > 2^j, 1 <= j <= m}`.

*Proof.* `glide_T(n) >= L` iff `T^j n > n` for `1 <= j <= L-1`. If `M_j > 1` then `T^j n >= M_j n > n`. If `M_j < 1` (never `= 1` for `j >= 1`), then `1 - M_j >= 2^-j`, so `(1 - M_j) n >= 2^(k-j) >= (3/2)^j > c_j` whenever `3^j <= 2^k`, i.e. `j <= k log_3 2`; then `T^j n = M_j n + c_j < n`. So for `L - 1 <= k log_3 2` the event depends only on the first `L - 1` parity letters, whose words run over `{0,1}^(L-1)`, each exactly `2^(k-L+1)` times in the block (Terras). ∎

Checked exactly in all 327 admissible (block, `L`) pairs up to `k = 31` (`.out` §E). `W_m` (exact big integers to `m = 1100`) reproduces THM-4479's table and THM-4495's sums `281, 2903, 31730, 367698`, and `W_m m^(3/2) 2^(-hm) in [9.66, 10.81]` for `500 <= m <= 1100` (THM-4495's window is `[9.66, 11.05]` to 3000). Hence `#{n in block k : glide_T(n) > c k} = 2^(k(1 - (1-h)c) + O(log k))` for `c < log_3 2`, PROVED; THM-4487/4495 give the same exponent at `c = 1`.

**Beyond the window (heuristic, EMPIRICAL test).** If parity letters stay fair beyond the window, `#{n <= X : glide >= L} ~ X W_(L-1) 2^-(L-1)`, and the record glide up to `X` is about `L*(X) = max{L : X W_(L-1) 2^-(L-1) >= 1}`, asymptotically `log_2 X/(1-h) = 19.98 log_2 X`. The actual records (`.out` §D) satisfy `glide/L* in [0.80, 1.16]` for all 30 records `n >= 1e4`, and `1.012` at the largest (`2.6e18`: glide 1005, `L* = 993`). The asymptotic `19.98 log_2 n = 1222` is far off at that size because of the `m^(-3/2)` ballot factor that `W_m` carries. **27 is the outlier: glide 59 against `L*(27) = 16`, a factor 3.69.** Glide records occur at about `0.83` per unit of `ln X` (35 up to `ln X = 42.4`), somewhat below the rate `~1` of an i.i.d. sequence.

### 3.4 Delays

Lagarias-Weiss's model (K-L Thm 4.2) gives `E #{n <= x : sigma_T(n) >= a ln n} = x^(1 - a gLW(1/a) + o(1))` with `gLW(y) = sup_t (t y - ln M_RRW(t))`. We checked that this equals `a (H(p) - p ln 3)`, `p = (ln 2 - 1/a)/ln 3`, to `1e-25`. It is 1 at `a = R1 = 6.952` and 0 at `a = 41.68`, so the records sit at its zero.
* **Delays like 27's** (`a = gamma(27) = 21.2389`): model exponent `0.6688`; the top four blocks give the least-squares exponent `0.656` (`.out` §E). EMPIRICAL agreement at the 2% level. Unlike §§3.1-3.3, these delay statements have no proved window version: a delay is decided far beyond the window.
* **Delay records** occur at about `3.4` per unit of `ln X` (148 up to `ln X = 44`), well above the i.i.d. rate. The reason is in §5: most delay champions are relatives of 27 in the backward tree.

### 3.5 The near-critical family: 27 is a typical long glide

Conditioned on `glide_T = L`, the walk `S_j = ln(T^j n/n)` is an excursion of the walk tilted to zero drift (odd frequency `log_3 2`, the Chernoff tilt `lambda*`), with step variance `sigma^2 = log_3 2 (1 - log_3 2)(ln 3)^2`, `sigma = 0.530138`. Its normalised height `u = max_j S_j/(sigma sqrt L)` should tend to the Brownian excursion maximum: `P(M <= x) = 1 + 2 sum_k (1 - 4k^2x^2) e^(-2k^2x^2)` (Chung, Kennedy, CITED), whose mean `sqrt(pi/2)` and second moment `pi^2/6` we verified numerically. Over all odd `n <= 2^32` with `glide >= 20` (EMPIRICAL, `.out` §E7): `E u = sqrt(pi/2) - (0.885 +- 0.015)/sqrt L` uniformly for `20 <= L < 200` (offsets `0.899, 0.890, 0.886, 0.882, 0.872, 0.884, 0.876, 0.876` in the eight bins, each with at least 94000 samples), i.e. Brownian plus a constant discrete offset. **27** has `u = ln(4616/27)/(sigma sqrt 59) = 1.2626`: the 55.6th percentile of the Brownian law and the 72nd percentile of the actual `u` of all `n <= 2^32` with glide in `[40,60)`. 27's famous peak is simply the typical height of a 59-step critical excursion. THM-4480's peak discount `1/w*` weights exactly these excursions down, and its critical band (height `~L^(1/3)`) is their thin core.

### 3.6 Path records: the count law `2 ln X`, and the clusters

**Proposition F.** If `t(1), t(2), ...` are independent with `P(t(n) <= y) = exp(-C n^theta y^-theta)`, then `P(t(n) > max_(m<n) t(m)) = n^theta / sum_(m<=n) m^theta`. *Proof.* Substitute `v = C y^-theta` in `int P(max_(m<n) t(m) < y) dP(t(n) <= y) = int_0^infinity n^theta e^(-v sum_(m<=n) m^theta) dv`. ∎

The Cramer law `P(t(n)/n >= w) ~ C/w` makes the extreme tail of `t(n)` Frechet with `theta = 1` and scale `C n`. The record probability `2/(n+1)` does not involve `C`, so the expected number of path records up to `X` is `2(H_(X+1) - 1) ~ 2 ln X` (HEURISTIC for Collatz; in the model the record indicators are also independent (Nevzorov's `F^alpha`-scheme, CITED from memory), so the standard deviation is about `sqrt(2 ln X)`, which is the scale used for `z` below).

| `X` | `10^3` | `10^6` | `10^9` | `10^12` | `10^15` | `10^18` | `10^21` | last term `2.36e21` |
|---|---|---|---|---|---|---|---|---|
| path records `<= X` | 10 | 25 | 44 | 61 | 73 | 87 | 94 | 98 |
| `2(H_(X+1) - 1)` | 13.0 | 26.8 | 40.6 | 54.4 | 68.2 | 82.0 | 95.9 | 97.6 |

Over all 98 records the largest deviation is `|z| = 1.34` standard deviations, and the least-squares slope of the count against `ln n` is `2.10` (`.out` §D). **The occurrence rate of path records (`2/n`) is the second Moran root `s2 = 2`.**

The model gets the count right and the values wrong. Records have median `t/n^2 = 0.07`, while independent Frechet records would sit at `t/n^2 ~ C/(2E)`, about `0.6`. The reason is clustering. Every `n` whose trajectory joins a champion's trajectory before the champion's peak shares that peak. So the super-quadratic family `{n : t(n) > n^2}` is the union, over peak values `P`, of the finite pieces `{n < sqrt P : t(n) = P}` of the backward trees `Pred*(P)`: up to `2^32` it consists of 21 numbers in 3 clusters (`.out` §E2):
`{27, 31, 41, 47, 54, 55, 62, 63}` (peak 4616; this is exactly `B27 cap [1, sqrt 4616]`),
`{319804831, 379027947, 426406441, 479707247, 568541921, 598957743, 639609662, 639609663, 719560871, 758055894, 758055895}` (peak `7.07e17`),
`{1410123943, 1880165257}` (peak `3.56e18`).
The expected count `C ln X ~ 18` at `2^32` matches (21); the number of independent champions (3) is far smaller.

## 4. Fractal recursion (T3)

### 4.1 Occurrence rate = covering number = dimension

**Theorem D.** Let `Bad = {x in Z_2 : 3^(o_j(x)) > 2^j for all j >= 1}` (closed), `R_W = {x : sup_j M_j(x) >= W}` (open), `R_inf = {x : sup_j M_j(x) = infinity}`.
1. For every `m`, the number of classes mod `2^m` meeting `Bad` is exactly `W_m`. Hence `dim_B Bad = lim m^-1 log_2 W_m = h(log_3 2)` (THM-4495), and `dim_H Bad = h(log_3 2)` (choice-ladder note, Besicovitch-Eggleston).
2. `mu_Haar(R_W) = P(W) in (2/(3W), 1/W]`; `mu(R_inf) = 0`; `dim_H R_inf = h(log_3 2)`.

*Proof.* The parity-vector map `Phi` is a 2-adic isometry: `x = y mod 2^m` iff their first `m` parity letters agree (Terras, Bernstein-Lagarias). It carries Haar measure to the fair coin measure, so everything can be done on parity sequences. (1) A cylinder `[u]`, `u in {0,1}^m`, meets `Phi(Bad)` iff `u in Bad_m`: if `u in Bad_m` then `u 1 1 1 ...` stays above the line, because odd letters only raise `M_j`. (2) The measure is Theorem R's `P(W)`. `R_inf` is null by Ville. For the upper bound on the dimension: `R_inf` lies in `{x : 3^(o_j) >= 2^j for infinitely many j}`, which for every `J` is covered by the cylinders of length `j >= J` with at least `j log_3 2` ones. There are at most `2^(j h(log_3 2))` of them (binomial tail, `log_3 2 > 1/2`), so the `s`-dimensional Hausdorff sum is at most `sum_(j>=J) 2^(j(h-s)) -> 0` for `s > h`. For the lower bound: for `p' in (log_3 2, 1)`, Bernoulli(`p'`) sequences have `S_j/j -> p' ln 3 - ln 2 > 0`, so `R_inf` has full `mu_(p')`-measure, and `mu_(p')` has local dimension `h(p')` almost everywhere. Hence `dim_H R_inf >= h(p') -> h(log_3 2)` as `p'` decreases to `log_3 2`. ∎

**The exact identity asked for in T3.** For `m = floor(k log_3 2)` and `L = m + 1`, Theorem G gives
`#{n in [2^k, 2^(k+1)) : glide_T(n) > m} = 2^(k-m) N_(2^-m)(Bad)`.
The occurrence count of the long-glide family is, exactly, a covering number of the exceptional fractal. So its counting exponent is `dim Bad = h = 1 + log_2 min g`. The rise family is the dual case: its density `P(W) ~ 0.83/W` decays with the exponent of the second root `s2 = 2`, while its limit fractal `R_inf` has the same dimension `h`. Both are read off `g`: one at a root, the other at the minimum.

### 4.2 The exponent-branching identity, tested on 27's branch

**Proposition B.** Let `F` be a set of positive integers with `T^-1(F) ⊂ F` and `F \ T^-1(F)` finite; every backward tree `Pred*(a)` of a non-periodic `a` qualifies (`F \ T^-1 F = {a}`). Let `N(X) = #F cap [1,X]` be regularly varying of index `s > 0`, and let `kappa = lim N_2(X)/N(X)`, `N_2` counting the members `= 2 mod 3`. Then `2^-s + kappa (3/2)^s = 1`.

*Proof.* `T^-1(F) = D(F) ⊔ E(F cap G)` (images of `D` even, of `E` odd; both injective), so `N(X) = O(1) + N(X/2) + N_2((3X+1)/2)`. Divide by `N(X)`; regular variation gives `N(X/2)/N(X) -> 2^-s` and `N_2((3X+1)/2)/N(X) -> kappa (3/2)^s`. ∎

With `kappa = 1/3` the identity reads `g(s) = 1`, so `s in {1, 2}`, and `s = 1` since `N(X) <= X`. **Test on 27's branch** (`.out` §E5): local exponents `log_2(N(2^(b+1))/N(2^b))` lie in `[0.9999, 1.0002]` for `b = 24..30`, the class fractions are `1/3` within `3e-4`, and `|2^-s + kappa(3/2)^s - 1| < 3e-3`. So `B27` is a family of exponent `s1 = 1` in branching balance. A tree growing like `x^0.84` (Krasikov-Lagarias's proved lower bound) would have to carry `kappa = 0.314`: the identity turns an exponent into a checkable class statistic.

### 4.3 The owner's 4x recursion: exact for sets, not for densities

For `x = 2 mod 3` (inverse-tree note, Prop. 11), `Pred*(x) = {x 2^i} ⊔ ⋃_(j>=0) Pred*(S^j E(x))` with `S(p) = 4p + 1`. The rungs `E(x), 4E(x)+1, ...` are the owner's `(4A-1)/3 = A + (A-1)/3` ladder. As counting identities these hold exactly, and the runner verifies them for `3077` (rungs `2051, 8205, 32821, ...`) and for the trunk `(4^j-1)/3`. In the Haar-averaged model (inverse-tree note, Prop. 10), a class-2 node's tree splits `1/4 : 3/4` between its doubling child and its odd child, and rung `j` carries `(3/4)(1/4)^j`. The ratio `1/4` is `4^-s` with `s = s1 = 1`: in the averaged model the rate exponent is literally the `4x` scaling exponent.

**Actual trees are not averaged** (FINITE-EXACT at `2^32`):

| tree | density | Haar `c R1/a` |
|---|---|---|
| `Pred*(41)` | `0.0566` | `0.339` |
| `Pred*(47)` | `0.1864` | `0.296` |
| `Pred*(182)` | `0.2829` | `0.076` |
| `Pred*(911)` | `0.3919` | `0.0153` |
| `Pred*(3077) = B27` | `0.39266` | `0.0045` |
| `Pred*(4616)` (the peak) | `0.39273` | `0.0030` |
| `Pred*(20)` | `0.9379` | `0.695` |
| `Pred*(5)` | `0.9379` | `2.78` |

* Along 27's orbit the densities increase (the trees are nested) from 0 at 27 (a multiple of 3: its tree is a doubling chain) to 1 at the root.
* Rung 0 (`2051`) carries `0.9998` of `3077`'s tree, against `3/4` in the model.
* At the root the owner's split `Pred*(8) = {8} ⊔ Pred*(16) ⊔ Pred*(5)` is `0.062 : 0.938`, against `1/4 : 3/4`.
* The rung densities `4^j delta` of the trunk scatter over `0.3 .. 9.7`.

**The reason is the second root again.** The tree of `3077` contains `27`, which is 114 times smaller. The Haar tree of a root `a` has few nodes below `a/K`: the pole of `F(s) = 1/(1 - g(s))` at `s2 = 2` (residue `R2`) suggests about `(R2/2) K^-2` of them (a Tauberian heuristic on the small side, not proved here), i.e. about `3e-4` for `K = 114`. A tree containing such a node is the backward image of a large rise: 27's orbit climbs by a factor 114 from 27 to 3077 and by 171 to its peak, and forward a rise by 171 has probability `0.83/171`. Such a tree inherits the fat trees of the small numbers on the rising orbit (`31, 41, 47, 62, 71, ...`). **The density of a single tree is governed by the smallest numbers it contains, not by the averaged recursion.** That is why 27's branch has density `0.39` instead of `0.0045`.

### 4.4 The precise answer to "the rate governs fractal recursion"

1. **(PROVED)** The family "no descent for `m` steps" occurs in every dyadic block with count exactly `2^(k-m)` times the `2^-m`-covering number `W_m` of the exceptional fractal `Bad`, for `m <= k log_3 2` (Theorems G, D). Its rate exponent is the fractal's dimension `h(log_3 2) = 1 + log_2 min g`.
2. **(PROVED)** The family "rises by a factor `W`" occurs with density between `2/(3W)` and `1/W` in the window, and with lower density `>= P(W)` for full trajectories. The exponent `-1` is `1 - s2`, where `s2 = 2` is the second root of the owner's Moran function (`g(2) = 1/4 + 3/4`, AM-fairness). The window spectrum of rises beyond `n^beta` is the slope `1 - s2 = -1` tangent to the dip-spectrum curve, then the curve itself (Theorem S).
3. **(PROVED identity + FINITE-EXACT test)** A backward-closed family of exponent `s` has branch fraction `kappa = (1 - 2^-s)(2/3)^s`; 27's branch has `s = 1`, `kappa = 1/3`.
4. **(Model PROVED, EMPIRICAL match)** Records of the rise family occur at rate `s2/n = 2/n`: 98 path records up to `2.36e21` against `97.6`.
5. **(CORRECTED)** The averaged recursion does **not** fix the density of an individual family such as 27's branch. The owner's `4x` self-similarity holds for sets and on average over 3-adic classes, not for single trees.

## 5. What "numbers beyond 27 in the same family" are (T4)

* **In the delay sense** they are, overwhelmingly, **27's relatives in the backward tree**. Of the 148 delay records, 104 lie in `B27`, as do 73% of all `n in [2^20, 2^32]` with `sigma_T(n) >= 16 ln n`, against 39% of all integers (`.out` §E8). Joining 27's trajectory early buys 27's long tail: from `31` or `47`, 66 or 67 more steps remain, against 6 from `20`. Delay records are decided by margins of a few steps. The numbers as long-lived as 27 relative to their size (`sigma_T >= 21.24 ln n`) occur at rate `X^0.67`.
* **In the height sense** they are the champions of a family occurring at the critical rate: `t(n) > n^2` holds for 21 integers up to `2^32`, in 3 tree clusters, and the path records number `2 ln X`, the second Moran root. 27 is the unique maximiser of `ln t(n)/ln n` up to `2^32` and among all 98 known path records. Rising by 27's factor 171 is common: the measured density is `0.486%`, and the lower density is PROVED to be at least `P(170.96) = 0.486%` (at least `2/(3 x 170.96) = 0.39%` in closed form).
* **In the glide sense** 27 is the largest outlier relative to the exact `W_k` model (3.7 times its expected record). Later glide records sit within 16% of the model, with an exact window theory (Theorem G) and a fractal of dimension `h` behind it.
* **In the parity-word sense** 27 is a typical critical excursion.
* **As a set**, 27's branch is a positive-density family: 39.27% of all integers join 27's trajectory before its peak.

## 6. Failures and corrections

* **Genealogy.** My first genealogy ("record `r` descends from an earlier record whose orbit it meets before that record's peak") declared 145 of 148 delay records descendants: small early records (`25`, `18`) sit on everyone's final approach. It was replaced by membership in `B27`, a fixed set.
* **Local exponents.** A first local-exponent formula for block counts had a spurious `+1`, since block counts already scale like `2^(bE)`. It was fixed before any claim was made.
* **Excursion law.** The first numerical check of the excursion law (mpmath `nsum` inside a quadrature) gave a wrong mean `0.666`, because the series was evaluated near 0 where it converges slowly. A float series with `F = 0` below `0.3` (true value `< 1e-20`) gives `sqrt(pi/2)` to 10 digits.
* **Haar tree law.** Its use for 27's branch is off by a factor 87 (§4.3); it is a statement about 3-adic averages.
* **Tolerances.** The EMPIRICAL checks are calibrated for `X = 2^32`. At `X = 2^28` the "delays like 27" exponent is still `0.77` (finite-size drift), and the runner prints such checks as data below `X = 2^31`.
* **K-L tables.** K-L Table 1's row `5649499` is not a `gamma` record, and K-L Table 3 omits the path record `10709980568908647`, which has `t > n^2`.

## 7. Reproduction

```bash
cd <worktree>
python3 -u 04-computation/experiments/procgen_family27_20260926_run.py 4294967296 \
    > 05-knowledge/results/procgen_family27_20260926.out
```

* The runner compiles `procgen_family27_20260926_scan.c` three times: production `N0 = 2^24`, and validation builds with `N0 = 2^12` and `N0 = 2^20`. It checks the `N0 = 2^12` build against the brute-force `procgen_family27_20260926_reference.py` at `X = 2^17` and the `N0 = 2^20` build against the production build at `X = 2^25` (identical reports), then runs the `2^32` scan. A plain-Python Monte Carlo cross-checks the scan's top block.
* It fetches the OEIS b-files once into `scratch/procgen_family27/oeis_cache` with `curl` and a generic User-Agent, and prints their sha256.
* Apple M2, 2026-09-26: 596 s wall time (the scan alone 562 s). Peak RSS 33.6 MB for the runner and 241.9 MB for the scan child (limit 700 MB). One heavy process at a time. Result: `ALL 84 CHECKS PASSED`.
* The EMPIRICAL tolerances are calibrated for `X = 2^32`; at `X < 2^31` those checks are printed as data (a run at `X = 2^25` passes its 77 checks in 35 s).

| file | sha256 |
|---|---|
| `04-computation/experiments/procgen_family27_20260926_run.py` | `83cb891b4b501b926f050d2d6ece82310bee70ae76ea05ba1ba45f04652b049a` |
| `04-computation/experiments/procgen_family27_20260926_theory.py` | `9180a2d746811d077bbaed3eb306d5588bf57bfa66c9aa79ad995db83deced06` |
| `04-computation/experiments/procgen_family27_20260926_bfiles.py` | `49c50504d574e2d116dea01ef43ea5dbfde6ec6373cb4ade4f9c13dcda357bc3` |
| `04-computation/experiments/procgen_family27_20260926_reference.py` | `5e6e5c9a8b6fa384d9458e7e838d0f8821ae5d9939cdad87b8bc951980d0dbd5` |
| `04-computation/experiments/procgen_family27_20260926_scan.c` | `74ef55ba059ffe12ea166c2259d060c0160e6d9cbad7cf3c429f043d09948466` |
| `05-knowledge/results/procgen_family27_20260926.out` | `19c8ea0f09fd33fe9e693fbbdf547145236ae132b62b7f0d2c50e58439783365` (raw bytes; contains timing lines) |

## 8. Parents and links

* [THM-4487](../../01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md) (dip spectrum `F`; §3.2 adds the slope -1 tangent). [THM-4495](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md) (`W_k`; Theorem D uses its order). [THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md) (the critical band and the peak discount; §3.5). [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md) (why the full-orbit rise upper bound is out of reach). [THM-4477](../../01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md) (AM-fairness `3 + 1 = 4`).
* [inverse tree mod 192](collatz_procgen_20260924_inverse_tree_mod192.md) (the owner's recursion, Prop. 10 Haar densities, Prop. 11 ladders), [choice ladder](collatz_procgen_20260922_choice_ladder.md) (`dim_H Bad`), [barrier atlas](collatz_procgen_20260922_barrier_atlas.md) (K-L and Krasikov-Lagarias entries).
* Lagarias, Weiss, *The 3x+1 problem: two stochastic models*, Ann. Appl. Probab. 2 (1992) 229-261 (via K-L). Kontorovich, Lagarias, *Stochastic models for the 3x+1 and 5x+1 problems*, arXiv:0910.1944 (read: Thm 4.1-4.4, 6.1-6.5, Tables 1-3).
