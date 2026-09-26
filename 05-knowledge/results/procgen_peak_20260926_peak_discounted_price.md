# The peak-discounted provability price: catch every undecided orbit at its highest point

Lane `peak`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_peak_20260926_{lib,integers,pairing,run}.py`.
Output: [procgen_peak_20260926.out](procgen_peak_20260926.out) (one runner; every `ok:` line is a `check(...)` that aborts on failure).

## Status

- **PROVED (full elementary hand proofs below).**
  - **Theorem 1 (two-sided peak theorem; every odd `q >= 3`, every `L >= 1`).**
    `rho^peak_L / M*_L <= eps_L(q) <= rho^peak_L`, where `rho^peak_L = 2^-L sum_{u in Bad_L} 1/w*(u)`,
    `w*(u) = max_{0<=j<=L-1} q^(e_j)/2^j`, and `M*_L <= L(L+1)(2L-2+3q)/(6q)`, `M*_L ~ (1-c^2)L^3/(6q)`.
    The upper bound is attained by *sending every undecided source to 1 at the first maximiser of its first `L` points*.
    The lower bound holds for the lower density of every `L`-step descent modification and, verbatim, for the
    flipped pair indices of the pairing family. The stratified form requested in T1 holds with factor
    `2((L-1)log_2(q/2)+1) * L * (L + L(L-1)/(2q)) = O(L^4)`. Definitions are fixed consistently: peak over times
    `0..L-1`, badness over `1..L`; the upper slope condition at `j = L` in THM-4478's band is never used (§2.4).
  - **Theorem 2 (exponent; every odd `q >= 3`).** For every `L >= 1`,
    `eps_L(q) <= rho^peak_L <= (q/2) sum_{e > cL} C(L,e) q^-e <= (q/2) 2^(-(1-H(c))L)`, `c = log_q 2`,
    and `eps_L(q) >= 2^(-(1-H(c))L - O(L^(1/3)))`. Hence `-(1/L) log_2 eps_L(q) -> 1 - H(log_q 2)`:
    `0.050044` (q=3), `0.013911` (q=5), `0.060510` (q=7), `0.100619` (q=9).
    For `q = 5` this closes the orchestrator's OPEN interval `[0.0119, 0.0139]` at its top, `0.013911`, although
    `rho_L(5) -> 0.1760` does not tend to 0. The arbitrary-edit price exponent is the entropy deficit of the
    critical density for every multiplier, independent of the drift sign.
  - **Theorem 3 (second order, explicit elementary constants).** With `sigma^2 = c(1-c)`, `lambda = (1-c)/c`,
    `z = q/max(1, lambda)`:
    `exp(-(C_q + o(1)) L^(1/3)) <= 2^((1-H)L) rho^peak_L <= C'_q(L) exp(-c_q (L-1)^(1/3))`, `C'_q(L) = O(L^(1/3))`,
    `c_3 = 0.3605, C_3 = 31.75`; `c_5 = 0.4167, C_5 = 33.16`; `c_7 = 0.4123, C_7 = 36.00`.
    Consequently `eps_L(q) = 2^(-(1-H(c))L) exp(-Theta(L^(1/3)))` for every odd `q`, and for `q = 3`
    `rho^peak_L / rho_L = exp(-Theta(L^(1/3)))`, `eps_L(3) <= rho_L * poly(L) * exp(-0.36 L^(1/3))`.
    **So the board's open question P1 ("is the price within a polynomial factor of `rho_L`?") has a NEGATIVE
    answer for arbitrary edits.**
  - **Lemma F (post-flip freshness).** In the pairing family, the value reached by a flip at time `k` has its next
    `l` Collatz parities uniformly distributed, independently of the first `k+1` letters (as `n` runs over residue
    classes).
- **PROVED modulo one CITED theorem** (Mogul'skii 1974, *Small deviations in the space of trajectories*, Theory
  Probab. Appl. 19, 726–736; statement used as verified in Gantert–Hu–Shi 2011, Ann. IHP 47, arXiv:0811.0262,
  Lemma 2.1, eqs. (2.12)–(2.13)).
  - **Theorem 4 (sharp second-order constant).**
    `ln rho^peak_L(q) = -(1-H(c)) L ln 2 - kappa_q L^(1/3) (1+o(1))`, and the same for `ln eps_L(q)`, with
    `kappa_q = (3/2) (pi^2 c(1-c))^(1/3) (ln z_q)^(2/3)`, `z_q = q min(1, c/(1-c))`:
    `kappa_3 = 2.107580`, `kappa_5 = 2.435966`, `kappa_7 = 2.410440`.
    For `q = 3`: `ln(rho^peak_L / rho_L) = -kappa_3 L^(1/3) (1 + o(1))`.
- **FINITE-EXACT.**
  - Enumeration of all `2^L` words equals the exact layer-cake DP (exact rationals) for `q = 3,5,7`, `L <= 16`
    (and `q = 3, L = 18`); the float DP matches the exact rationals to `1.3e-15` (q=3 up to `L = 64`, q=5,7 up to
    `L = 48`).
  - The digest's numbers `E[1/w* | bad] = 0.0825, 0.0102, 3.03e-3, 4.03e-4, 2.79e-5` (`L = 16, 60, 100, 200, 400`)
    are reproduced to 3 significant figures.
  - The peak-catch modification on actual integers: for `q = 3,5,7`, `L = 8,10,12,14,16`, every `2 <= n <= 10^6`
    descends within `L`; exceptional sets are empty (`q=3,7`) and `{13, 17}` (`q=5`, the 5n+1 cycles).
  - Pairing family: the partner-isolated peak catch gives valid sections (every `3 <= n <= 10^6` descends within
    `L`) for `L = 8..32`.
- **EMPIRICAL.**
  - The fitted second-order constants: q=3 `2.095` (prefactor `L^(2/3)` fixed) / `2.158` (free),
    q=5 `2.418` / `2.517`, q=7 `2.421` / `2.405`, against `kappa_q` above. The naive `z = q` constants
    (`2.766`, `3.070` for q=5,7) are excluded.
  - The measured density of the peak edit set is `0.59–0.78 rho^peak` (several bad sources share one peak point).
  - **T5 (pairing family):** the partner-isolated greedy peak catch has flip density `1.19–1.51 x 2rho^peak` for
    `L = 8..32`, i.e. `0.21 -> 0.047 x 2rho_L`: on the tested range the pairing price behaves like `rho^peak`,
    not like `rho_L`.
- **OPEN.**
  - P1 for the pairing family (the board's literal P1, `delta_L`): is `delta_L <= rho_L exp(-c L^(1/3))`?
    Proved: `rho^peak_L/M*_L <= delta_L <= 2 rho_L`. Numerics say yes-like (§7); no proof.
  - Whether the greedy pairing peak catch ever gets stuck for larger `n`; the asymptotics of the shared-peak
    factor of `E`; the polynomial prefactors (heuristics: `L^(2/3)` for the q=3 ratio, `L^(-5/6)` for the
    normalised q>=5 quantity).
- Nothing here bears on the truth of the Collatz conjecture: these are fixed-horizon modification prices. No novelty
  or priority claim; the probabilistic input is classical (Chernoff, Kolmogorov, Paley–Zygmund, Mogul'skii).

## 1. Setting and definitions

Fix an odd `q >= 3` and put `T(n) = n/2` (n even), `(qn+1)/2` (n odd), `c = log_q 2` (irrational, as
`q^a != 2^b` for `a, b >= 1`), `H = H_2`. An *`L`-step descent modification* is any `G: N_(>0) -> N_(>0)` with: for
every `n >= 2` some `1 <= j <= L` has `G^j(n) < n`. `E(G) = {v : G(v) != T(v)}`; `eps_L(q)` is the infimum of the
upper density of `E(G)`.

For a word `u in {0,1}^L`: `e_j(u)` = number of ones among `u_0..u_(j-1)`; `w_j(u) = q^(e_j)/2^j`;
`S_j = e_j - jc`, so that `w_j = q^(S_j)`. Then

- `Bad_L = {u : w_j(u) > 1 for 1 <= j <= L}` (`w_j = 1` is impossible for `j >= 1`), `rho_L = |Bad_L|/2^L`;
- `w*(u) = max_{0<=j<=L-1} w_j(u) = q^M`, `M = max_{j<L} S_j`; `rho^peak_L = 2^-L sum_{u in Bad_L} 1/w*(u)`;
- strata `S_j = {u in Bad_L : 2^j <= w*(u) < 2^(j+1)}`.

**Terras and the affine form (standard; any odd q).** The parity word of `(n, Tn, ..., T^(L-1) n)` depends only on
`n mod 2^L`, and this is a bijection `Z/2^L -> {0,1}^L`. On the class of a word `u`,
`T^j(n) = w_j (n + h_j)`, `h_j = sum_{i<j, u_i=1} 1/(q w_i) >= 0` (an odd step is `(q/2)(x + 1/q)`). If all
`w_i >= 1` then `h_j <= e_j/q <= j/q`.

**Why the index ranges are what they are.** An edit at the time-`j` point acts at time `j+1`, so only the points at
times `0..L-1` can be edit points for descent by time `L`: the peak is taken over `0..L-1`. A source needs help iff
`T^j(n) >= n` for all `1 <= j <= L`: badness is over `1..L` (the lower slope condition at `j = L` is essential).

**Exceptional set.** `Bad^act_L = {n >= 2 : T^j(n) >= n, 1 <= j <= L}` contains every lift of a bad word
(`T^j n >= w_j n > n`). Conversely `X_L = Bad^act_L \ (bad classes)` is finite: if `u` is not bad, some `w_j < 1`,
and `T^j(n) >= n` forces `n <= w_j h_j/(1-w_j)`. FINITE-EXACT: `X_L cap [2,10^6]` is empty for `q = 3, 7` and equals
`{13, 17}` for `q = 5` (`L = 8..16`), the minima of the 5n+1 cycles `13 -> 33 -> 83 -> 208 -> 104 -> 52 -> 26` and
`17 -> 43 -> 108 -> 54 -> 27 -> 68 -> 34`.

## 2. Theorem 1: the price is `rho^peak_L` up to a cubic factor

**Theorem 1.** For every odd `q >= 3` and `L >= 1`:

- (upper) the modification `G_peak` below satisfies `upper-density(E(G_peak)) <= rho^peak_L`;
- (lower) every `L`-step descent modification satisfies `lower-density(E(G)) >= rho^peak_L / M*_L`, where
  `M*_L = sum_{k=0}^{L-1} sum_{e=e_min(k)}^{k} (floor(e/q) + 1)`, `e_min(0) = 0`, `e_min(k) = floor(kc) + 1`.

Hence `rho^peak_L / M*_L <= eps_L(q) <= rho^peak_L`.

### 2.1 Upper bound: catch at the peak

For `n in Bad^act_L` let `j*(n)` be the least maximiser of `j -> T^j(n)` on `0 <= j <= L-1`, and `P(n) = T^(j*(n))(n)`.
Put `E = {P(n) : n in Bad^act_L}`, `G = 1` on `E`, `G = T` elsewhere. (`E` avoids `{1,2}`: `P(n) >= n >= 2`, and
`2` is not bad.)

*`G` is an `L`-step descent modification.* Let `n >= 2` and `tau = min{k >= 0 : T^k(n) in E}` (possibly infinite).
By induction `G^i(n) = T^i(n)` for `0 <= i <= tau` (G = T off E), and `G^(tau+1)(n) = 1 < n` if `tau` is finite.
If `n in Bad^act_L`, then `P(n) in E`, so `tau <= j*(n) <= L-1` and `G^(tau+1)(n) = 1 < n` with `tau + 1 <= L`.
Otherwise let `d = min{j in [1,L] : T^j(n) < n}`. If `tau >= d` then `G^d(n) = T^d(n) < n`; if `tau < d` then
`G^(tau+1)(n) = 1 < n` with `tau + 1 <= d <= L`. ∎

*Density.* `E(G) ⊆ E` and `|E cap [1,Y]| <= #{n in Bad^act_L : P(n) <= Y}`. For `n` in the class `r_u` of a bad word
`u`, `P(n) = max_{j<L} w_j(n + h_j) >= w*(u) n`, so `P(n) <= Y` forces `n <= Y/w*(u)`: at most
`Y/(2^L w*(u)) + 1` such `n`. Summing over `Bad_L` and adding `X_L`:

    |E cap [1, Y]| <= Y rho^peak_L + |Bad_L| + |X_L|.

Hence `upper-density(E) <= rho^peak_L`. ∎ (The runner checks this finite inequality with `Y = N = 10^6`.)

### 2.2 Lower bound: single-scale integer capacity

Let `G` be any `L`-step descent modification and `Y >= 1`. Let `Src(Y)` be the set of `n >= 2` lying in the class
of a bad word `u` with `n <= Y/w*(u) - (L-1)/q`. Each class contributes at least
`(Y/w* - (L-1)/q)/2^L - 2 >= Y/(2^L w*) - 3` points, so `|Src(Y)| >= Y rho^peak_L - 3|Bad_L|`.

Take `n in Src(Y)`. All `T^j(n) = w_j(n+h_j) > n` (`1 <= j <= L`). Let `k` be the first index with
`T^k(n) in E(G)`; as above `G^i(n) = T^i(n)` for `i <= k`, so if `k >= L` (or no such `k`) then
`G^j(n) = T^j(n) > n` for all `j <= L`, a contradiction. Thus `k <= L-1`, and `v = T^k(n) <= w*(n + (L-1)/q) <= Y`.
Assign `n` to `v`.

*Capacity.* Fix `v`, a depth `k <= L-1`, and `e = e_k(n)`. Then `n = v 2^k/q^e - h_k` with `0 <= h_k <= e/q`
(all `w_i >= 1`), so `n` lies in an interval of length `e/q`: at most `floor(e/q) + 1` integers. A source has one
trajectory, hence one `(k, e)`; and `e_min(k) <= e <= k` (`q^e > 2^k` for `k >= 1`). So `v` receives at most `M*_L`
sources, and

    Y rho^peak_L - 3|Bad_L| <= |Src(Y)| <= M*_L |E(G) cap [1, Y]|.

Dividing by `Y` and letting `Y -> infinity` along any sequence gives `lower-density(E(G)) >= rho^peak_L / M*_L`. ∎

*Explicit factor.* `e` takes at most `k` values (one for `k = 0`) and `floor(e/q) + 1 <= k/q + 1`, so
`M*_L <= sum_{k<L} (k+1)(k/q+1) = L(L+1)(2L-2+3q)/(6q)`; asymptotically `M*_L ~ (1-c^2)L^3/(6q)` (ratio 1.009,
1.013, 1.018 at `L = 512`, q = 3,5,7).

*Pairing family.* For a member `F` of the pairing family (q = 3) the edit set is exactly the union of the flipped
pairs `{2i-1, 2i}`, so `|E(F) cap [1,Y]| <= 2 #{flipped i <= (Y+1)/2}`; bad sources are odd and `>= 3`, so the
descent condition of `P_L` (`n >= 3`) suffices. Hence the lower density of flipped pair indices is also
`>= rho^peak_L / M*_L`. ∎

### 2.3 The stratified form (as requested in T1)

**THM-4478 Theorem A with the band `B'_L(K)`.** Put `B'_L(K) = {u : w_j >= 1 (0<=j<=L), w_j <= K (0<=j<=L-1)}`.
THM-4478's proof (§§2–3) of `lower-density(E(G)) >= |B_L(K)|/(2^L K M_L(K))`,
`M_L(K) = (floor(log_q K)+1) R_L`, `R_L = sum_{k<L}(floor(k/q)+1)`, uses the upper slope bound only at the first-hit
depth `k <= L-1` (to bound `T^k(n) <= K(X + L/q)` and to count the admissible `e`); it never uses `w_L <= K`. So it holds
verbatim with `B'_L(K) ⊇ B_L(K)`.

Every `u in S_j` lies in `B'_L(2^(j+1))`, and also in THM-4478's own band `B_L(q 2^j)`, because
`w_L <= (q/2) w_(L-1) < q 2^j` (checked by enumeration: `max_{bad} w_L/w* <= q/2` at `L = 12`). Hence

    lower-density(E(G)) >= |S_j| / (2^L 2^(j+1) M_L(2^(j+1)))    for every j.

Since `2^-L sum_{S_j} 1/w* <= |S_j| 2^(-L-j)`, summing over the `J+1` strata (`J = floor((L-1) log_2(q/2))`) gives
`lower-density >= rho^peak_L / (2(J+1) M_L(2^(J+1)))`, and `floor((J+1)c) + 1 <= L - (L-2)c <= L`,
`R_L <= L + L(L-1)/(2q)`, so the factor is at most `2((L-1)log_2(q/2)+1) L (L + L(L-1)/(2q)) = O(L^4)`. Using THM-4478's
band verbatim costs `K = q 2^j` instead of `2^(j+1)` (a factor `q/2` and one more admissible `e`).

**How good are the three lower bounds?** (runner §2c; `best lower / rho^peak`)

| q | L | rho^peak | single `rho^peak/M*_L` | stratified | THM-4478 Theorem A (best `K = 2^t`) | best / rho^peak |
|---|---|---|---|---|---|---|
| 3 | 16 | 2.661e-3 | 1.470e-5 | 8.134e-6 | 6.775e-6 (K=2^4) | 1/181 |
| 3 | 32 | 3.012e-4 | 2.383e-7 | 1.266e-7 | 1.518e-7 (K=2^5) | 1/1264 |
| 3 | 64 | 1.302e-5 | 1.382e-9 | 1.041e-9 | 1.395e-9 (K=2^6) | 1/9339 |
| 3 | 128 | 1.125e-7 | 1.547e-12 | 1.263e-12 | 1.937e-12 (K=2^9) | 1/58083 |
| 5 | 16 | 3.188e-3 | 2.018e-5 | 8.527e-6 | 9.983e-6 (K=2^6) | 1/158 |
| 5 | 128 | 1.171e-6 | 1.955e-11 | 1.681e-11 | 2.359e-11 (K=2^12) | 1/49653 |

All lower bounds are rigorous; so is `rho^peak`; the true `eps_L` lies between. The peak-weighted strata show where
`rho^peak` lives: at `q = 3, L = 128` the shares of the strata `j = 5..9` (`2^j <= w* < 2^(j+1)`) are
`0.068, 0.198, 0.261, 0.214, 0.134`.

### 2.4 Answer to "does the stratification need the slope condition at `j = L`?"

The *lower* condition `w_L > 1` is needed and is part of badness (otherwise the source may descend at time `L`
without any edit). The *upper* condition at `j = L` is not needed: Theorem A is valid for `B'_L(K)`. With THM-4478's
band verbatim use `K = q 2^j`. The peak must be taken over `0..L-1` (not `0..L`); the two versions differ by at most
the factor `q/2`, but only the `0..L-1` version is realised by a construction.

## 3. Theorem 2: the exponent is `1 - H(log_q 2)` for every odd `q`

**Theorem 2.** For every odd `q >= 3` and `L >= 1`,

    eps_L(q) <= rho^peak_L <= (q/2) sum_{e > cL} C(L,e) q^-e <= (q/2) 2^(-(1-H(c))L),

and `eps_L(q) >= rho^peak_L / M*_L >= 2^(-(1-H(c))L - O(L^(1/3)))` (Theorem 3; the block construction gives the weaker
`O(sqrt(L log L))`). So `lim -(1/L) log_2 eps_L(q) = 1 - H(log_q 2)`.

*Majorant.* For `u in Bad_L`, `w_L <= (q/2) w_(L-1) <= (q/2) w*`, so `1/w* <= (q/2)/w_L`, and `Bad_L ⊆ {w_L > 1}`:

    rho^peak_L <= (q/2) 2^-L sum_{e: q^e > 2^L} C(L,e) 2^L q^-e = (q/2) sum_{e > cL} C(L,e) q^-e.

*Chernoff.* For `z >= 1`, `sum_{e >= cL} C(L,e) q^-e <= z^(-cL) (1 + z/q)^L`. Take `z = qc/(1-c)`; `z >= 1` iff
`c(q+1) >= 1` iff `(q+1) ln 2 >= ln q`, true for every `q >= 1`. Then `1 + z/q = 1/(1-c)` and the bound is
`[(1-c)^(-(1-c)) c^(-c) q^(-c)]^L = [2^(H(c)) / 2]^L` since `q^c = 2`. ∎

*Remark (the monotonicity asked for in T2).* With `f(sigma) = sigma + 1 - H(c(1+sigma))` (final slope
`<= 2^(sigma L)`), `f'(sigma) = 1 - c log_2((1-p)/p)`, `p = c(1+sigma)`, is increasing, and `f'(0) >= 0` iff
`((1-c)/c)^c <= 2 = q^c` iff `c(q+1) >= 1`: exactly the condition `z >= 1` above. So `sigma + rate(sigma)` is
increasing on `sigma >= 0` for every odd `q`, and the one-line Chernoff bound with `z = qc/(1-c)` packages it.

*Minorant (block construction for `c < 1/2`).* THM-4478 §4 transfers: take `k_b = ceil(cb)`; a block with `k_b`
ones has total log-slope `(k_b - cb) ln q in (0, ln q)`; rotating after the last minimum of the cumulative
log-slopes makes every prefix slope `>= 1` (cycle lemma), at most `b` rotations per output, so at least
`C(b,k_b)/b` blocks; internal slopes `<= (q/2)^(k_b) <= q^b`. Concatenating `t` blocks (`L = tb + r`) and `r` ones
gives words with all slopes in `[1, q^(t+b+r)]`, hence `rho^peak_L >= 2^-L (C(b,k_b)/b)^t q^-(t+b+r)`. Nothing uses
`c > 1/2`: for `c < 1/2` and `b >= 1/(1-2c)`, `c <= k_b/b < c + 1/b <= 1-c`, so `H(k_b/b) >= H(c)` (rounding up moves
towards 1/2, no entropy loss), and the error is `O((L/b) log b + b)`, i.e. `O(sqrt(L log L))`. The runner checks the
block minorant against the exact `rho^peak` for `q = 3,5,7`, `L = 100, 400, 1000` (it is weak at these `L`:
`log_2 = -335` against `-83.7` at `q = 3, L = 1000`).

*Error terms (T2).* Elementary (Theorems 1–3):

    -(1-H)L - (C_q/ln 2 + o(1)) L^(1/3) - log_2 M*_L  <=  log_2 eps_L(q)  <=  -(1-H)L - (c_q/ln 2) L^(1/3) + O(log L),

together with the clean `log_2 eps_L(q) <= -(1-H)L + log_2(q/2)` for every `L`, and `log_2 M*_L = 3 log_2 L + O(1)`.
Sharp, with the cited lemma (Theorem 4): `log_2 eps_L(q) = -(1-H)L - (kappa_q/ln 2) L^(1/3)(1+o(1))`, where
`kappa_q/ln 2 = 3.0406` (q=3), `3.5144` (q=5), `3.4775` (q=7) bits.

*The q = 5 conclusion.* `eps_L(5) = 2^(-0.013911 L + o(L))` exactly, with `rho_L(5) -> 0.1760`. The peak
construction beats the catch-high construction (orchestrator's Proposition 4, balance exponent `theta_5 = 0.011921`)
by factors `3.0e-7`, `4.4e-9`, `1.7e-11` at `L = 200, 400, 800` (e.g. `rho^peak_800(5) = 2.527e-14` against the
catch-high bound `1.497e-3`).

## 4. Theorem 3: the second order is `exp(-Theta(L^(1/3)))` (elementary)

**Tilting identity.** Let `Q` make the letters iid Bernoulli(`c`); then `S_j = e_j - jc` is a centred walk with
steps `1-c` (prob. `c`) and `-c`, variance `sigma^2 = c(1-c)` per step. For every word, `2^-L = Q(u) 2^-L c^-e (1-c)^-(L-e)`
and `e = cL + S_L` give

    2^-L = Q(u) 2^(-(1-H(c))L) lambda^(S_L),        lambda = (1-c)/c,

hence `rho^peak_L = 2^(-(1-H)L) E_Q[lambda^(S_L) q^-M ; bad]` and `rho_L = 2^(-(1-H)L) E_Q[lambda^(S_L) ; bad]`
(checked by enumeration at `L = 12`). For `q = 3`, `lambda = 0.585 < 1`; for `q >= 5`, `lambda > 1`.

**Lemma C (confinement of the Q-walk; elementary).**

1. `E S_n^4 = 3n^2 sigma^4 + n sigma^2 (1 - 6 sigma^2)`. (`E xi^4 = sigma^2(1 - 3 sigma^2)`; expand.)
2. If `n sigma^2 >= 1` then `Q(S_n > 0) >= 1/16` and `Q(S_n < 0) >= 1/16`. *Proof:* for centred `X`, Hölder gives
   `E X^2 <= (E|X|)^(2/3) (E X^4)^(1/3)`; `E X^+ = E|X|/2`; Cauchy–Schwarz `E X^+ <= (E X^2 Q(X>0))^(1/2)`; so
   `Q(X > 0) >= (E X^2)^2/(4 E X^4) >= 1/16` by (1) (`E S_n^4 <= 4 n^2 sigma^4`).
3. Kolmogorov: `Q(max_{i<=n} |S_i| >= t) <= n sigma^2 / t^2` (CITED, classical).
4. If `a >= 1` and `n sigma^2 >= 2a^2` then `Q(|S_n| > a) >= 1/14`. *Proof:* Paley–Zygmund for `Z = S_n^2` with
   `theta = a^2/(n sigma^2) <= 1/2`: `Q(Z > theta EZ) >= (1/4)(EZ)^2/EZ^2 >= (1/4)/3.5`.
5. (upper) For `a >= 1`, `N >= 1`: `Q(0 < S_j < a, 1 <= j <= N) <= (13/14)^floor(N/n_a)`, `n_a = ceil(2a^2/sigma^2)`.
   *Proof:* cut `[0, N]` into `floor(N/n_a)` blocks; confinement forces every block increment into `(-a, a)`, which by
   (4) and independence has probability `<= 13/14` per block.
6. (lower) For `a >= 38` with `n = floor(a^2/(1152 sigma^2))` and `n sigma^2 >= 1`, and `N >= 1`:
   `Q(0 < S_j <= a, 1 <= j <= N) >= c^ceil(a/(2(1-c))) 32^-ceil(N/n)`.
   *Proof:* first `r = ceil(a/(2(1-c)))` letters are ones: `S` climbs to `S_r in [a/2, a/2+1) ⊂ [a/3, 2a/3]`.
   Then blocks of length `n`: from `x in [a/3, 2a/3]` require (i) `max_{i<=n} |S'_i| <= a/6` and (ii) the block
   increment `>= 0` if `x < a/2`, `<= 0` if `x >= a/2`. On this event the block stays in `[a/6, 5a/6]` and ends in
   `[a/3, 2a/3]`. By (2) and (3), its probability is `>= 1/16 - 36 n sigma^2/a^2 >= 1/16 - 1/32`. (A final partial
   block only needs (i).) ∎

Runner checks: (1) to 40 digits (`n <= 60`, `q in {3,5,7,9,27,101}`); the minima of the probabilities in (2) and (4)
over the tested ranges are `0.3512 >= 1/16` and `0.4795 >= 1/14`.

**Theorem 3.** Let `z = q/max(1, lambda)`, `C_lambda = max(1, lambda^(1-c))`. For every `L >= 2`,

    2^((1-H)L) rho^peak_L <= C_lambda sum_{m>=0} z^-m (13/14)^floor((L-1)/ceil(2(m+1)^2/sigma^2))
                          <= C'_q(L) exp(-c_q (L-1)^(1/3)),
    C'_q(L) = C_lambda (14/13) z (c_q (L-1)^(1/3)/ln z + 1 + z/(z-1)) = O(L^(1/3)),

    2^((1-H)L) rho^peak_L >= max_a (q/min(1,lambda))^-a c^ceil(a/(2(1-c))) 32^-ceil(L/floor(a^2/(1152 sigma^2))),

the max over integers `a >= 38` with `floor(a^2/(1152 sigma^2)) sigma^2 >= 1`; this is `exp(-(C_q + o(1)) L^(1/3))`, where
`c_q = (3/2^(2/3)) (ln(14/13) sigma^2/3)^(1/3) (ln z)^(2/3)` and
`C_q = (3/2^(2/3)) (1152 sigma^2 ln 32)^(1/3) (ln(q/min(1,lambda)) + ln(1/c)/(2(1-c)))^(2/3)`.

*Proof.* Upper: on `bad`, `lambda^(S_L) <= 1` if `lambda <= 1`; if `lambda > 1`, `S_L <= S_(L-1) + 1 - c <= M + 1 - c`,
so `lambda^(S_L) q^-M <= lambda^(1-c) z^-M`. In both cases the integrand is `<= C_lambda z^-M`. Split on
`M in [m, m+1)`: then `0 < S_j < m+1` for `1 <= j <= L-1`, and Lemma C(5) applies. The closed form uses
`n_a <= 3a^2/sigma^2`, `min_a (a ln z + B/a^2) = (3/2^(2/3)) B^(1/3) (ln z)^(2/3)`, and
`sum_a exp(-f(a)) <= sum_a min(e^(-f_min), z^-a)`. Lower: on the event of Lemma C(6) with `N = L`, the word is bad,
`M <= a`, `S_L <= a`, so the integrand is `>= min(1,lambda)^a q^-a`. The asymptotic form takes `a ~ theta L^(1/3)`. ∎

**Consequences.** (a) `eps_L(q) = 2^(-(1-H(c))L) exp(-Theta(L^(1/3)))` for every odd `q` (Theorems 1, 3).
(b) For `q = 3`: `rho_L <= 2^(-eta L)` (as `lambda < 1`) and `rho_L >= 2^(-eta L)/(2.72 L(L+1))` for `L >= 10`
(cycle lemma: at least `C(L, ceil(cL))/L` bad words; `H(ceil(cL)/L) >= H(c) - 1.4415/L`; as in THM-4479 §3), so

    exp(-(C_3 + o(1)) L^(1/3)) <= rho^peak_L / rho_L <= 2.72 L(L+1) C'_3(L) exp(-0.3605 (L-1)^(1/3)),

and `eps_L(3) <= rho^peak_L = o(rho_L / L^A)` for every `A`: **P1 is answered negatively for arbitrary edits.**
The constants are far from sharp (`0.36` and `31.8` against the true `2.108`; the runner prints the elementary
bracket, e.g. `-319.9 <= -23.33 <= -2.09` for `ln[2^(eta L) rho^peak_1000]`); the sharp constant is Theorem 4.
(c) THM-4478's Theorem A with `K = q^(theta L^(1/3))` already gives `eps_L >= 2^(-eta L - O(L^(1/3)))`, improving its
stated `O(sqrt(L log L))` error to the right order: `rho_L(K) >= 2^(-(1-H)L) min(1,lambda)^a Q(0 < S_j <= a, j <= L)`.

## 5. Theorem 4: the sharp constant (Mogul'skii)

**Cited lemma** (Mogul'skii 1974; in the form of Gantert–Hu–Shi 2011, Lemma 2.1, (2.12)–(2.13), read from
arXiv:0811.0262). Let `X_i` be iid, centred, with a finite `(2+eta)`-moment and variance `sigma^2 > 0`; let
`g_1 < g_2` be continuous on `[0,1]` with `g_1(0) < 0 < g_2(0)`; let `a_n -> infinity`, `a_n^2/n -> 0`. Then

    lim (a_n^2/n) log P(g_1(i/n) <= S_i/a_n <= g_2(i/n), 1 <= i <= n) = -(pi^2 sigma^2/2) int_0^1 dt/(g_2 - g_1)^2,

and the same limit holds with the extra condition `S_n/a_n >= g_2(1) - b`, any `b > 0` (by reflection, also with
`S_n/a_n <= g_1(1) + b`).

**Theorem 4.** For every odd `q >= 3`, with `theta* = (pi^2 sigma^2/ln z)^(1/3)` and
`kappa_q = (3/2)(pi^2 sigma^2)^(1/3)(ln z)^(2/3) = min_theta [theta ln z + pi^2 sigma^2/(2 theta^2)]`,

    ln E_Q[lambda^(S_L) q^-M ; bad] = -kappa_q L^(1/3) (1 + o(1)).

Hence `ln rho^peak_L = -(1-H(c)) L ln 2 - kappa_q L^(1/3)(1+o(1))`; the same holds for `ln eps_L(q)` (Theorem 1, the
factor `M*_L` is polynomial); and for `q = 3`, `ln(rho^peak_L/rho_L) = -kappa_3 L^(1/3)(1+o(1))` (Consequence (b)).

*Proof, upper bound.* Integrand `<= C_lambda z^-M` (Theorem 3). Fix small `delta` and `eps` with
`ln z/eps > kappa_q` and `ln(14/13) sigma^2/(3 eps^2) > kappa_q`. Terms with `m+1 < eps L^(1/3)` are
`<= exp(-(kappa_q + const) L^(1/3))` by Lemma C(5); terms with `m+1 > L^(1/3)/eps` sum to
`O(exp(-(ln z/eps) L^(1/3)))`. For the middle range use the finite grid `a_i = eps(1+delta)^i L^(1/3)`,
`0 <= i <= ceil(ln(eps^-2)/ln(1+delta))`. For `m+1 in [a_i, a_(i+1)]`,
`Q(0 < S_j < m+1, j <= L-1) <= Q(-delta a_(i+1) <= S_j <= a_(i+1), j <= L-1)`, and the lemma with `n = L-1`,
`a_n = eps(1+delta)^(i+1)(n+1)^(1/3)`, `g_1 = -delta`, `g_2 = 1` gives `-(pi^2 sigma^2/2)(1+delta)^-2` as the limit of
`(a_n^2/n) log`. With finitely many `i`, for `L >= L_0(eps, delta)` every middle term is
`<= z exp(-(m+1) ln z - (1-delta)(1+delta)^-4 (pi^2 sigma^2/2)(L-1)/(m+1)^2)`, whose maximum over `m` is
`exp(-kappa_q L^(1/3)(1 - O(delta)))`; there are `O(L^(1/3))` terms. Let `delta -> 0`.

*Proof, lower bound.* Let `a = theta* L^(1/3)`, `r = ceil(delta a/(1-c))`. Event: the first `r` letters are ones
(`S_r in [delta a, delta a + 1)`, cost `c^r = exp(-O(delta a))`); then the next `n = L - r` increments
`S'_i = S_(r+i) - S_r` satisfy `-delta a'_n <= S'_i <= (1-2delta) a'_n` with `a'_n = theta* n^(1/3)` and the endpoint
condition `S'_n <= -delta a'_n/2` if `lambda <= 1` (end low), `S'_n >= (1-3delta) a'_n` if `lambda > 1` (end high).
Since `a'_n <= a` and `a >= 1/delta` for large `L`, the whole walk stays in `(0, a]` (`S_j != 0` for `j >= 1`), so the
word is bad, `M <= a`, and the integrand is `>= exp(-a ln q - O(delta a))` (`lambda <= 1`: `S_L = O(delta a)`) or
`>= z^-a lambda^(-O(delta a))` (`lambda > 1`: `S_L >= (1 - 3 delta - o(1)) a`). The lemma
along the sequence `n` (width `1-delta`) gives probability `exp(-(pi^2 sigma^2/2)(1+o(1)) n/((1-delta)^2 a'^2_n))`.
Altogether `>= exp(-L^(1/3)[theta* ln z + pi^2 sigma^2/(2 theta*^2)](1 + O(delta)))`. Let `delta -> 0`. ∎

*Why `z = qc/(1-c)` for `q >= 5`:* the tilt factor `lambda^(S_L)` rewards ending at the top of the strip, so each unit
of strip height costs `ln q - ln lambda`, not `ln q`. The data exclude the naive `z = q` constants
(`2.766`, `3.070`).

*Heuristic prefactors (not proved).* A walk confined to `(0, a)` for time `L`, entering and leaving at distance `O(1)`
from the boundary, has probability `~ a^-3 exp(-pi^2 sigma^2 L/(2a^2))`; the Laplace sum over `a ~ L^(1/3)` adds
`L^(1/6)`. So the numerator is `~ L^(-5/6) e^(-kappa L^(1/3))` and, as `E_Q[lambda^(S_L); bad] ~ L^(-3/2)`, the q=3 ratio
is `~ L^(2/3) e^(-kappa_3 L^(1/3))`.

## 6. Numerics (T4)

All tables are copied from [procgen_peak_20260926.out](procgen_peak_20260926.out).

**Exact values, `q = 3`** (exact rationals, printed in floating point):

| L | rho_L | rho^peak_L | E[1/w* \| bad] |
|---|---|---|---|
| 8 | 7.421875e-02 | 1.310299e-02 | 0.176546 |
| 16 | 3.225708e-02 | 2.660874e-03 | 0.082490 |
| 24 | 1.708156e-02 | 8.333989e-04 | 0.048789 |
| 32 | 9.626961e-03 | 3.012321e-04 | 0.031290 |
| 48 | 3.537286e-03 | 5.431076e-05 | 0.015354 |
| 64 | 1.493065e-03 | 1.302336e-05 | 0.008723 |

**Float layer cake, `q = 3`** (bracket relative width `<= 1.2e-15`):

| L | rho_L | rho^peak_L | ratio | ln ratio | ln ratio / L^(1/3) |
|---|---|---|---|---|---|
| 100 | 2.3868e-04 | 7.2427e-07 | 3.0345e-03 | -5.7977 | -1.2491 |
| 200 | 3.0604e-06 | 1.2340e-09 | 4.0322e-04 | -7.8160 | -1.3365 |
| 400 | 1.1587e-09 | 3.2360e-14 | 2.7929e-05 | -10.4858 | -1.4231 |
| 600 | 6.4455e-13 | 2.6028e-18 | 4.0381e-06 | -12.4197 | -1.4725 |
| 1000 | 2.9019e-19 | 6.3523e-26 | 2.1890e-07 | -15.3347 | -1.5335 |
| 1500 | 4.6194e-27 | 6.4561e-35 | 1.3976e-08 | -18.0859 | -1.5800 |
| 2000 | 9.0581e-35 | 1.3817e-43 | 1.5254e-09 | -20.3010 | -1.6113 |

Fits of `ln ratio` on `L >= 100` (model `-kappa L^(1/3) + beta ln L + gamma + delta L^(-1/3)`):
`beta = 2/3` fixed gives `kappa = 2.0954` (within 0.6% of `kappa_3 = 2.1076`); free gives `kappa = 2.158`,
`beta = 1.00`; `kappa = kappa_3` fixed gives `beta = 0.733`. Local slopes `d ln ratio / d L^(1/3)` between `L = 400` and
`2000`: `-1.805, -1.823, -1.853, -1.872, -1.901, -1.901, -1.923`, against `-kappa_3 + 2/L^(1/3)`:
`-1.846 ... -1.941` (all within 0.05).

**Normalised `ln[2^((1-H)L) rho^peak_L]`, `q = 5, 7`:**

| L | q=5: rho_L | q=5: norm | /L^(1/3) | q=7: rho_L | q=7: norm | /L^(1/3) |
|---|---|---|---|---|---|---|
| 100 | 0.1806 | -11.3711 | -2.4498 | 0.3008 | -11.8345 | -2.5497 |
| 200 | 0.1769 | -14.5731 | -2.4920 | 0.3008 | -15.1276 | -2.5868 |
| 400 | 0.1761 | -18.5808 | -2.5218 | 0.3008 | -19.2146 | -2.6078 |
| 600 | 0.1760 | -21.3777 | -2.5346 | 0.3008 | -22.0513 | -2.6145 |
| 1000 | 0.1760 | -25.4621 | -2.5462 | 0.3008 | -26.1781 | -2.6178 |

Fits on `L >= 200`: q=5 free `kappa = 2.517` (`beta = -0.32`), `beta = -5/6` fixed `2.418`, against `2.436`;
q=7 free `2.405` (`beta = -0.92`), `beta = -5/6` fixed `2.421`, against `2.410`.

**The peak construction on actual integers** (`N = 10^6`; all `2 <= n <= N` descend within `L`):

| q | L | rho_L | rho^peak_L | #{bad n: P(n) <= N}/N | \|E cap [1,N]\|/N | / rho^peak |
|---|---|---|---|---|---|---|
| 3 | 8 | 0.07422 | 0.013103 | 0.013102 | 0.009957 | 0.760 |
| 3 | 12 | 0.05518 | 0.006681 | 0.006681 | 0.004250 | 0.636 |
| 3 | 16 | 0.03226 | 0.002661 | 0.002659 | 0.001586 | 0.596 |
| 5 | 8 | 0.27344 | 0.014295 | 0.014297 | 0.011060 | 0.774 |
| 5 | 16 | 0.23376 | 0.003188 | 0.003171 | 0.002063 | 0.647 |
| 7 | 8 | 0.34375 | 0.008858 | 0.008859 | 0.006921 | 0.781 |
| 7 | 16 | 0.31839 | 0.001519 | 0.001524 | 0.001055 | 0.694 |

The count of bad sources with peak `<= N` equals `N rho^peak` to within 0.6% in every run (the counting in the
upper bound is tight); the edit set itself is 22–40% smaller because bad sources on one orbit share its peak point.
(At `L = 16` there are only about 15 lifts per residue class below `10^6`; larger `L` would need a larger `N`.)

## 7. T5: the pairing family can catch high too (EMPIRICAL)

Pairing family (THM-4470/4475): one bit per pair `{2i-1, 2i}`; a flip sends `2i-1 -> i-1` and `2i -> 3i`. Unlike an
arbitrary edit, a flip cannot jump to 1: it replaces the step factor `3/2` by `1/2` at an odd point, and pushes the
partner up.

**Constructions** (processing `n = 3, 4, ..., N` as in THM-4475; a successful path is frozen):

- *THM-4475 rule* (A = own pair, F, B; re-implemented from the note): reproduces the audited densities
  `0.06806` (L=8) and `0.02850` (L=16) on pairs `<= N/4` — a check of the simulator.
- *Partner-isolated peak catch:* try the odd path points `v` with `v = 2 (mod 3)` in *decreasing* order of value; accept
  the first flip after which `n` descends within `L` and the partner `v+1` (now sent to `3(v+1)/2`) descends below
  itself; freeze both paths. `v = 2 (mod 3)` makes the partner a multiple of 3, which has no odd Collatz preimage
  (THM-4475's isolation device). Fallbacks:
  - a THM-4475-style source-level flip: the own pair, the first odd path point `w = 3 (mod 4)` below `2n`, or the
    pair of `T(n)`; used 13 times in all ten runs together;
  - a depth-4 search, never used.

| L | THM-4475 rule | peak catch | 2 rho_L | 2 rho^peak | peak / 2rho^peak | peak / 2rho_L | mean rank of accepted flip |
|---|---|---|---|---|---|---|---|
| 8 | 0.068052 | 0.031274 | 0.148438 | 0.026206 | 1.193 | 0.211 | 1.29 |
| 12 | 0.049190 | 0.017886 | 0.110352 | 0.013363 | 1.339 | 0.162 | 2.01 |
| 16 | 0.028426 | 0.006902 | 0.064514 | 0.005322 | 1.297 | 0.107 | 2.43 |
| 20 | 0.022694 | 0.004658 | 0.052124 | 0.003414 | 1.364 | 0.089 | 3.02 |
| 24 | 0.014986 | 0.002284 | 0.034163 | 0.001667 | 1.370 | 0.067 | 3.39 |
| 28 | 0.011532 | 0.001542 | 0.026260 | 0.001032 | 1.495 | 0.059 | 3.84 |
| 32 | 0.008398 | 0.000908 | 0.019254 | 0.000602 | 1.507 | 0.047 | 4.18 |

(densities on pair indices `<= N/2`, `N = 10^6`; `L = 10, 14, 18` are in the output.) Every `3 <= n <= 10^6` descends
within `L` in every run (FINITE-EXACT: a valid section of a member of `P_L`).

**Reading.** On `L = 8..32` the pairing price `delta_L` is at most the peak-catch density, which stays within
`1.19–1.51` of `2 rho^peak` while its ratio to `2 rho_L` falls from `0.21` to `0.047`, tracking `rho^peak/rho_L`
(`0.177 -> 0.031`). So delta_L behaves like `rho^peak`, not like `rho_L`, on the tested range. The constant drifts
upward and the accepted flip moves down the candidate list (mean rank `1.3 -> 4.2`): the highest odd point is often
too late in the window for the halved continuation to come down in time. The data cannot exclude a slowly growing
factor.

**Lemma F (post-flip freshness; PROVED).** Fix a prefix `u_0..u_k` with `u_k = 1`, class `r mod 2^(k+1)`. For
`n = r + 2^(k+1) t`, `T^k(n) = T^k(r) + 2 * 3^(e_k) t` (T^k is affine with slope `3^(e_k)/2^k` on the class), so the
post-flip value `b = (T^k(n) - 1)/2 = b_0 + 3^(e_k) t` runs through all residues mod `2^l` as `t` runs mod `2^l`. So
the next `l` parities after a flip are uniform and independent of the prefix (checked for `k+1 <= 7`, `l <= 8`).
This is why one high flip usually suffices. For a flip at a *stopping time*, the continuation is a fresh Collatz
word (up to the sparse other flips), with drift `ln(3/4)/2` per step. It must fall by about `ln(w*/3)`, which is of
order `L^(1/3)` for the words that carry `rho^peak`. It does so within the remaining time with probability tending
to 1, provided the flip leaves `>> L^(1/3)` steps. A proof of `delta_L <= poly(L) rho^peak_L` would still have to
control three things: flips that are not at stopping times, partner interactions, and frozen-pair blocking. **OPEN.**

## 8. Failures and caveats

- **Pairing construction, two failed versions.** Without partner isolation the greedy peak catch stranded partners
  `2i` of earlier flips (even sources whose own pair was already flipped; e.g. `n = 890154` at `L = 12`, `50346` and
  `827946` at `L = 14`). A partner-safe but non-isolated version still stranded `n = 735217` at `L = 16`: its path
  entered a flipped partner `1102826 = 2i` from its odd preimage `(4i-1)/3`, and every pair on the path was frozen.
  Partner isolation `v = 2 (mod 3)` removed all failures in the tested range. That the final greedy never gets stuck
  for larger `n` is not proved.
- The edit set of the peak construction has density `0.59–0.78 rho^peak`; `eps_L` may be below `rho^peak` by more than
  this sharing factor. Only the polynomial window `[rho^peak/M*_L, rho^peak]` is proved; at `L = 128` it spans a
  factor `~5.8e4`.
- The elementary constants of Theorem 3 are poor (`0.36` / `31.8` vs `2.108`); only the order `L^(1/3)` is elementary.
  The sharp constant needs the cited Mogul'skii lemma, which I did not re-derive.
- The fits cannot separate `kappa` from the prefactor exponent `beta` well: the free q=3 fit gives `2.158`
  (+2.4%), and the q=3 normalised quantity oscillates with `frac(cL)`, so only the ratio is fitted for q=3.
- The digest subagent's claims are confirmed: `eps_L = rho^peak_L` up to `poly(L)` (Theorem 1), and for q=3
  `rho^peak/rho_L = exp(-Theta(L^(1/3)))` by strip confinement (Theorems 3–4). Its numbers are reproduced (§6). For
  `q >= 5` the Mogul'skii constant has `z = qc/(1-c)`, not `q`.

## 9. Contrast with the strategy cube (THM-4479)

In the cube the price is polynomially sharp: `N_k <= delta_k <= |Bad_k|` with `N_k >= 2^(hk)/(3k^2) >= |Bad_k|/(3k^2)`
(THM-4479, Theorems 1–2), both `2^(hk) k^(-3/2)` up to slowly varying factors. For arbitrary edits the price is
`rho_L exp(-Theta(L^(1/3)))` (q=3).

The difference is height. A periodic strategy edits residue classes, and the Haar cost of a class is `2^-k` whatever
the height of the orbit points it catches. An arbitrary edit at `v` costs density weight `~1/v`, so catching an orbit
at height `w* n` costs `1/w*` of catching it at `n`. Residue classes cannot see height; the necklace lower bound of
THM-4479 is a purely 2-adic obstruction. The same archimedean-versus-2-adic split separated catch-high from periodic
edits in the orchestrator's §2 for `q >= 5`. The present result shows it also changes the second-order term for
`q = 3`, where the exponents agree.

## 10. Reproduction

```
python3 04-computation/experiments/procgen_peak_20260926_run.py > 05-knowledge/results/procgen_peak_20260926.out
```

- Python 3.10, numpy 2.2.6, scipy 1.15.3, mpmath 1.3.0. Deterministic; the only run-dependent lines of the output
  are those starting with `[time]`.
- Final run (macOS, 8 GB machine shared with other lanes): wall time 252.9 s, peak RSS 394.9 MB
  (`/usr/bin/time -l`; the runner's own report is 377 MiB). A trial run gave output identical apart from `[time]`
  lines and five lines added later (four strata prints and the Lemma F check).
- The output has 220 lines; all `ok:` lines passed and it ends with `ALL CHECKS PASSED`.

| file | sha256 (raw bytes) |
|---|---|
| `04-computation/experiments/procgen_peak_20260926_lib.py` | `ecc5f8a27c13ca40fad95ae8f6486a88a982d4cfa386d5b3808ac307ec778070` |
| `04-computation/experiments/procgen_peak_20260926_integers.py` | `b443fac49dc41f114e9ba622d4a0d61aba289461db46553b6768c5c486108b50` |
| `04-computation/experiments/procgen_peak_20260926_pairing.py` | `5b8915c08a89f908d09e6cfdbb824d4d5f3334b3fbae65b989626a1bed28811c` |
| `04-computation/experiments/procgen_peak_20260926_run.py` | `ccc11e68718529ca9336b4b64a449234895f99b4da177a96c8ffe255009e818e` |
| `05-knowledge/results/procgen_peak_20260926.out` | `9144dfecbbb30e9ee74a550db827393323ba04ead6c02b5bb8650029805ff679` |
| same output without the `[time]` lines (`grep -v '^\[time\]'`) | `a255e1ab6ed06d8b85850e22f5bdb45273ba2b5115b9e695c310f9f115801614` |

What each script does:

- `lib`:
  - exact floors `floor(mc)` by integer comparison;
  - enumeration of all words;
  - the vectorised strip DP, exact (Python integers) or float, which gives `F(s)` at every breakpoint;
  - the layer-cake sum with a certified truncation bracket;
  - `M*_L`, `M_L(K)`;
  - the Chernoff majorant and the block minorant;
  - the elementary Lemma C bounds.
- `integers`: the peak construction on `n <= 10^6` and its verification.
- `pairing`: the THM-4475 rule and the partner-isolated peak catch, with final verification of every
  `3 <= n <= 10^6`.
- `run`: all checks.

Temporary files lived only under `scratch/procgen_peak/` (not part of the deliverable).
