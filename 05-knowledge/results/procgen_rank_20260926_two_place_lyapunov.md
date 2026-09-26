# What a Collatz rank must look like: every point of every backward orbit of an expanding cycle carries a forced 2-adic charge; infinite banks, adaptive centers, heights

**Status.**
- **PROVED** (full hand proofs below; each statement is also exercised by the runner):
  - **Lemma L (local expansion).** For a signed bank of valuation counters at rational centers with `sum |c_beta| (1 + log2 H(beta)) < infinity` ("height-summable"), the potential near any rational point `zeta` is `c_zeta D + Phi*(zeta) + o(1)` at integers of 2-adic depth `D` and size `<= 2^(D+B)`. The error bound is explicit.
  - **Theorem A (forced charges).** Let `R(n) = a log n + h(n) + sum_beta c_beta v2(n - beta)` with `a > 0`, `h` bounded, the bank height-summable. Suppose `R` never increases along Collatz steps from some `n_1` on, or descends with bounded lookahead `L`. Then every point `z` of the backward orbit of every expanding periodic orbit `x` carries the charge `c_z >= a chi(x)`, where `chi(x) = log(3^a/2^p)/p`. Stepwise, the charge is also non-increasing along forward orbits and constant on the cycle, and the rest of the bank is a tension on the cycle.
  - **Corollary A1.** No height-summable bank works for Collatz (signed, with any bounded `h`, stepwise or with bounded lookahead). The backward orbit of `-1` alone has `2^(j-1)` points at depth `j`, each needing charge `>= a log(3/2)`. Explicit violating integers are given.
  - **Theorem R+ (strategy cube).** For a level-`k` sign strategy, a height-summable bank rank exists iff the strategy is class (i). So infinite height-summable banks add nothing to THM-4482's periodic and bounded corrections.
  - **Proposition S (the shadow is paid, the seam is not).** A charge `c > a chi` on the points of one expanding cycle, plus a periodic `h`, makes `R` descend with margin at every step inside the cycle's shadow (checked for `-1`, `-5`, `-17`, the 3n−1 cycle `5, 7, 10`, and a level-5 sign strategy). The entry step from any uncharged preimage fails, with slope exactly `c`.
  - **Theorem C (nonnegative banks are nowhere height-summable).** If `c >= 0` and `R` is finite at one odd and one even integer, then the total mass is finite. If `R` also descends, every nonempty open set of `Z_2` carries centers with infinite height moment `sum c_beta (1 + log2 H(beta))`.
  - **Theorem D (Collatz-completeness).** A nonnegative bank at rational centers, with `R(Tn) <= R(n) - eta` for all `n >= 2`, exists **iff every Collatz orbit reaches 1**. The construction puts one "private" atom `z_m = m + 2^(D_m)/3` next to each integer and stores the stopping time of `m` in the height `2^(D_m)`.
  - **Propositions Q2.1–Q2.3 (adaptive centers).** Consider ranks `a log n + c Lam(n)`, where `Lam(n)` is the 2-adic closeness of `n` to the nearest expanding periodic point.
    - With the period bounded by `P`, they fail at `2^(V+1) - 2`, even with lookahead (Q2.1).
    - With preperiodic centers of bounded preperiod `Q`, they fail at `2^(Q+1)(2^V - 1)` (Q2.2).
    - With unbounded period, `Lam` is an excursion time in disguise: `l* <= Lam < 2 l*`, where `l*` is the last time the multiplier `3^(a_l)/2^l` exceeds 1. This rank fails by unbounded amounts on explicit "hovering" integers built from lower Christoffel words of the upper best approximations `2/3, 7/11, 12/19, 53/84, ...` of `log_3 2` (Q2.3).
  - **Proposition Q2.5.** The windowed future envelope `max_(i <= theta log2 m) log T^i m` never increases iff every orbit point's future peak comes within `theta log2 m` steps.
  - **Q3 identities.**
    - The Liouville bound `v2(n - beta) <= log2((n+1) H(beta))`.
    - The shadow identity: inside a shadow, the numerator of `n - x` loses one factor 2 per step and gains one factor 3 per odd step. Its prime-to-6 part is invariant, and the local changes at `infinity`, 2 and 3 sum to 0.
    - The forced charge is the ratio of the archimedean to the 2-adic expansion rate: `chi = log 2 · log|lam|_inf / log|lam|_2`.
- **FINITE-EXACT.**
  - The forced mass per period `M(p)` for `p <= 64` (`M(64) = 8.2·10^15`).
  - The strategy-cube census:
    - no strategy of levels 2–5 has only isolated expanding cycles (exhaustive, 65,812 strategies);
    - the loop at `-1` can be an isolated SCC from level 6 on (CP-SAT witness);
    - "loop at `-1` isolated and every other cycle contracting" is INFEASIBLE for levels 4–10 (CP-SAT).
  - Theorem D's bank, verified on the window of all orbits of `2..1000`.
  - The hovering integers for `2/3` (`r <= 5`), `7/11` (`r <= 2`) and `12/19` (`r = 1`); the largest has 946,509 bits.
- **EMPIRICAL.** Violation statistics of the adaptive candidates, measured on:
  - every step of the orbits of `3..5000` (254,283 steps);
  - the path and delay records up to `10^12` (49,771 steps);
  - 600 random starts in `[10^11, 10^12]` (110,295 steps).

  Every violation of a current-center rank at a point `>= 11` is a reset, never a shadow step. The largest unbounded-period violation is `126.4`, at `m = 414028311267890` on the orbit of `881715740415`. Future-peak horizons `theta* = 9.46, 12.93, 8.48` suffice on the three sets.
- **CITED.**
  - OEIS A006884 (path records) and A006877 (delay records), terms `<= 10^12`, used only as test inputs.
  - 2 is a primitive root modulo every power of 3 (standard).
  - THM-4474, THM-4479, THM-4482 (in-repo).
- **ANALOGY.**
  - Banks as two-place local heights (Néron/Green functions at `infinity` and 2) of an infinite divisor.
  - Counters as 2-adic local canonical heights of repelling cycles.
- **OPEN.**
  - Nonnegative banks that are summable but not height-summable. Theorem D shows that deciding whether one works is equivalent to Collatz. For the natural period-only profiles on periodic centers, the band `sum_p f(p) P(p) < infinity = sum_p p f(p) P(p)` is not excluded.
  - Whether a strategy with finitely many (but at least one) expanding periodic points exists at any level.
  - The windowed-envelope horizon: this is a peak-time conjecture.
  - Collatz itself.
- No HYP or THM file was created. Nothing here proves or disproves Collatz. These are statements about which *kinds* of Lyapunov functions can exist.

Session `collatz-procgen-20260922`, rank lane, 2026-09-26.
- **Scripts:** `04-computation/experiments/procgen_rank_20260926_{lib,q1,cube,q2,q3,run}.py`.
- **Output:** [procgen_rank_20260926.out](procgen_rank_20260926.out).
- **Parents:**
  - [THM-4482](../../01-canon/theorems/THM-4482-ranks-are-tensions-strategy-cube.md) and its [note](procgen_tension_20260926_ranks_christoffel_duality.md) §2.1–2.4 (ranks = tensions; Corollary K on finite banks);
  - the [Kuratowski reframe](kuratowski_reframe_20260925.md) §3–5 and §7 (repetition lemma, seam identity, the reset family `n_H`, the centers `alpha_d`, the open "resource that pays for resets");
  - the [Bernoulli boundary note](bernoulli_boundary_20260925.md) §3–4;
  - [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md) (two places) and [THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md) (height beats residues);
  - [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) and [THM-4479](../../01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md).

## 0. The answers in brief

| question | answer | status |
|---|---|---|
| Q1: weights forced by one expanding cycle | `c_x >= a chi(x)` at every point `x` of the cycle. Inside the shadow the counter drops `p` per period while `log n` grows `log lambda`. The bound is sharp: `c > a chi` plus a periodic `h` pays the whole shadow (Prop. S). | PROVED |
| Q1: the weights forced overall | The same charge at **every point of the backward orbit** of the cycle (`2^(j-1)` points at depth `j` for `-1`). This is the *entry* (seam) obstruction; the cycle gives only the circulation obstruction. | PROVED (Thm A) |
| Q1: summable? is `R` finite? | No. The forced mass per period is `M(p) = 2^(h p + O(log p))`, and it is already infinite on a single backward orbit. With the forced profile, `R = +infinity` at every positive integer. | PROVED + FINITE-EXACT table |
| Q1: impossibility class | Every signed bank with `sum |c|(1 + log2 H) < infinity`, plus any bounded `h`, stepwise or with bounded lookahead. Nonnegative banks must be height-non-summable in every open set. In the cube: bank-provable ⟺ class (i). | PROVED (Thms A, C, R+) |
| Q1: non-vacuity | Locally yes: every shadow is payable (Prop. S), including the 3n−1 cycle `5, 7, 10` and a level-5 sign strategy. Globally, a nonnegative bank exists **iff Collatz holds** (Thm D). In the cube, no strategy of levels 2–5 has finitely many expanding periodic points. | PROVED + FINITE-EXACT |
| Q2: adaptive centers | Current center with period `<= P`: fails at `2^(V+1) - 2`. Preperiod `<= Q`: fails at `2^(Q+1)(2^V - 1)`. All periods: `Lam` is an excursion time (`l* <= Lam < 2 l*`) and fails on hovering integers, by at least `c p floor(log 2/log lambda) - log 2` (6.0, 48.4, 431.5, 12400 for `2/3, 7/11, 12/19, 53/84`). Every observed violation is a reset. | PROVED + EMPIRICAL |
| Q2: what does work on data | Windowed future envelopes with horizon `theta log2 m`. They never increase iff future peaks come within the window, a peak-time conjecture. On the tested sets `theta = 13` suffices. | PROVED (equivalence) + EMPIRICAL |
| Q3: heights | The product formula bounds 2-adic closeness by height, which is why height-summable banks are tame. In shadows the dynamics moves valuation from the place 2 to the place 3, and the forced charge is the ratio of the local expansion rates at `infinity` and 2. No height-type distance to the rational cycles decreases along shadows. | PROVED identities + ANALOGY |
| Q4: minimal feature | A rank must grow like `a chi(x) v2(n - z)` near **every** point `z` of the (dense) backward orbits of **all** expanding cycles. This is impossible for summable linear aggregation with bounded height. A successful rank must either store archimedean data in the heights of its centers (Thm D; nowhere height-summable), or aggregate max-plus over the future (unbounded, height-dependent lookahead). | PROVED dichotomy |

## 1. Setting and notation

* `T` is the shortcut Collatz map, `T(x) = x/2` (`x` even) and `(3x+1)/2` (`x` odd), on `Z_(2)`, the rationals with odd denominator (a subring of `Z_2`; parity is the parity of the numerator). For a level-`k` sign strategy `sigma` (THM-4474), `T_sigma(x) = (3x + sigma(x mod 2^k))/2` on odd `x`.
* `v2` is the 2-adic valuation. For `u/D` in lowest terms, `H(u/D) = max(|u|, |D|)` is the naive height and `h(x) = log H(x)` the logarithmic Weil height. Write `kappa = log(3/2)`, `c_* = log_3 2` and `h* = h_bin(c_*) = 0.94996`.
* **Periodic points.** A parity word `w` of length `p` with `a` ones has the unique periodic point `x_w = C_w/(2^p - 3^a)` (THM-4474 Lemma 2). It is **expanding** iff `lambda = 3^a/2^p > 1`, with Lyapunov exponent `chi(x) = log(lambda)/p`. Expanding periodic points of `T` are negative (`C_w > 0`). Their heights satisfy `H(x_w) <= 3^(a-1) 2^p` (checked for all words of length `<= 14`: `log2 H <= 1.543 p`).
* **Backward orbit.** `T^(-inf)(x) = {z : T^j z = x for some j >= 0}`. Every `y` in `Z_(2)` has exactly the two preimages `2y` and `(2y-1)/3`; the only self-preimage is `-1`. So the backward tree of `-1` has `2^(j-1)` points at depth `j >= 1` (checked to depth 14). It is dense in `Z_2`: `T^j` maps each class mod `2^j` affinely onto `Z_2` (Terras).
* **Banks.** A bank is a countable family of centers `beta` in `Z_(2)`, none a positive integer, with real weights `c_beta` (atoms at the same point are merged, so `c_x` is the total weight at `x`, and `c_x = 0` if `x` is not a center). Its potential is `Phi(n) = sum_beta c_beta v2(n - beta)`.
  * **Height-summable** means `S_H := sum_beta |c_beta| (1 + log2 H(beta)) < infinity`.
  * A **rank** is `R(n) = a log n + h(n) + Phi(n)` with `a > 0` and `h : Z_(>0) -> R` bounded.
  * **(S)** means `R(Tn) <= R(n)` for all `n >= n_1`. **(L)** means that every `n >= n_1` has some `1 <= j <= L` with `R(T^j n) <= R(n)`.
* **Potential form.** `v2(x - y) = sum_(s>=1) 1[x = y mod 2^s]`, so `Phi(x) = sum_(s>=1) mu(x + 2^s Z_2)` for the (signed) measure `mu = sum c_beta delta_beta`. A bank is the 2-adic logarithmic potential of a charge distribution, the non-archimedean analogue of `log|x - y|^(-1)`.
* **Good integers.** For `z` in `Z_(2)` and `D >= 1`, `y_D(z)` denotes the unique integer in `[2^(D+1), 2^(D+2))` with `v2(y - z) = D` exactly.
* **Shadows (THM-4474 Lemma 1).** If `v2(y - z) = D >= k + i`, then `T^i y - T^i z = (3^(a_i)/2^i)(y - z)` exactly, with `a_i` the number of odd points among `z, ..., T^(i-1) z`. So `v2(T^i y - T^i z) = D - i`.

## 2. Q1: banks of valuation counters

### 2.1 The product-formula bound and Lemma L

**Liouville bound (PROVED).** For an integer `n >= 1` and `beta = u/D` in `Z_(2)` with `beta != n`:
`v2(n - beta) = v2(nD - u) <= log2 |nD - u| <= log2((n+1) H(beta))`.
For two rationals, `v2(zeta - beta) <= log2(2 H(zeta) H(beta))`. This is the product formula `prod_v |y|_v = 1` in its weakest form: 2-adic closeness costs height. Checked on 23,994 pairs, 12,000 of them deliberately 2-adically close.

**Lemma L (local expansion; PROVED).** Let the bank be height-summable, let `zeta` be in `Z_(2)` and `B >= 2`. Put
`Phi*(zeta) = sum_(beta != zeta) c_beta v2(zeta - beta)` (absolutely convergent) and
`delta_zeta(D) = sum_(beta != zeta, v2(beta - zeta) >= D) |c_beta| (1 + log2 H(beta))`, which tends to 0.
Then for every integer `y >= 1` with `v2(y - zeta) = D` and `y + 1 <= 2^(D+B)`:

    | Phi(y) - c_zeta D - Phi*(zeta) | <= (B + log2 H(zeta)) delta_zeta(D).

*Proof.* The difference is `sum_(beta != zeta) c_beta [v2(y - beta) - v2(zeta - beta)]`. Split the sum:
* If `v2(zeta - beta) < D`, the term is `0`.
* If `v2(zeta - beta) > D`, then `v2(y - beta) = D`, and the bracket lies in `[-v2(zeta - beta), 0]`, bounded by `log2(2 H(zeta) H(beta))`.
* If `v2(zeta - beta) = D`, the bracket is `v2(y - beta) - D`, which lies in `[0, B + log2 H(beta)]` by the Liouville bound.

In both nonzero cases the bracket is at most `(B + log2 H(zeta))(1 + log2 H(beta))`. ∎

*Checks.* Signed banks with 240 near atoms `zeta + 2^m/3` of weights `~ m^(-3)` (heights `~ 2^m`) at seven points `zeta`: `|error| <= 0.031 × bound`, and the bound is `~2·10^-3` at `D = 160`.

**Remark (truncated counters).** A truncated counter `min(v2(n - beta), N)` equals the potential of the Haar probability on `beta + 2^N Z_2`, minus the indicator of that ball. Everything below holds with `log2 H(beta)` replaced by `min(N_beta, log2 H(beta))`. A truncated family with `sum |c_beta| N_beta < infinity` is just a bounded `h`, which covers the Bernoulli-boundary periodic corrections.

### 2.2 Theorem A: forced charges

**Theorem A.** Let `R` be a rank whose bank is height-summable (and `h` bounded), and suppose `R` satisfies (S) or (L). Let `x` be an expanding periodic point of `T`, with orbit `x_0 = x, ..., x_(p-1)` and exponent `chi = chi(x)`. Let `z` be in `T^(-inf)(x)`, with `T^j z = x_0`. Then:
* **(i)** `c_z >= a chi`.
* **(ii)** Under (S), `c` is non-increasing along the forward orbit `z, Tz, ..., x_0` and constant on the cycle, `c_(x_i) = c_x >= a chi`.
* **(iii)** Under (S) with `h = 0`, the bank is a tension on the cycle: `a w(x_i) - c_x + Phi*(x_(i+1)) - Phi*(x_i) <= 0` for every `i`, where `w = log(3/2)` at odd and `-log 2` at even points.

The same holds for every level-`k` strategy `T_sigma` in place of `T`.

*Proof.* Write `z_i = T^i z`, so `z_(j+i) = x_(i mod p)`, and `y = y_D(z)`.
* **Setup.** For each fixed `i`, as `D -> infinity`:
  * `T^i y` is at exact depth `D - i` from `z_i` (shadows);
  * `T^i y + 1 <= 2^(D - i + 4 + 2i)`;
  * `T^i y >= 2^(D+1-i) -> infinity`, so all these points eventually exceed `n_1`;
  * `log(T^(i+1) y / T^i y) -> w(z_i)`.

  By Lemma L (with `B = 4 + 2i`), `R(T^i y) = a log T^i y + h(T^i y) + c_(z_i)(D - i) + Phi*(z_i) + o(1)`.
* **(ii), monotonicity.** (S) at `T^i y` reads
  `a log(T^(i+1)y/T^i y) + [h(T^(i+1)y) - h(T^i y)] + c_(z_(i+1))(D - i - 1) - c_(z_i)(D - i) + O(1) <= 0`.
  The `h`-difference is at most `2||h||`. Dividing by `D` and letting `D -> infinity` gives `c_(z_(i+1)) <= c_(z_i)`. On the cycle `c` is periodic and non-increasing, hence constant.
* **(i) under (S).** Fix `m >= 1` and take `z = x_0` (depth `D`). Chaining the `mp` steps and using `log(T^(mp) y / y) -> m log lambda` gives
  `a m log lambda - c_x m p <= 2||h|| + o(1)`.
  Let `D -> infinity`, then `m -> infinity`: `c_x >= a chi`. With (ii), `c_z >= c_x`.
* **(iii).** With `h = 0`, the one-step inequality inside the cycle has limit `a w(x_i) - c_x + Phi*(x_(i+1)) - Phi*(x_i) <= 0`.
* **(i) under (L).** Fix `N` and `T_max = j + (Np+1)L + L`.
  * For each `D`, follow the lookahead jumps from `y`: times `0 = t_0 < t_1 < ...` with gaps `<= L` and `R(T^(t_(s+1)) y) <= R(T^(t_s) y)`, up to time `T_max`. There are finitely many jump patterns, so one pattern occurs for infinitely many `D`; work along those `D`.
  * Each jump gives `c_(z_(t_(s+1))) <= c_(z_(t_s))` (divide by `D`, as above).
  * At least `Np + 1` jump times lie in `[j, j + (Np+1)L]`. So some residue class mod `p` contains `N + 1` of them, `tau_0 < ... < tau_N`, all at the same cycle point `x'`, with `tau_N - tau_0 = Mp` for some `M >= N`.
  * Chaining the jumps from `tau_0` to `tau_N` gives `a M log lambda - c_(x') M p <= 2||h|| + o(1)`, hence `c_(x') >= a chi - 2||h||/(Np)`.
  * Since `c` is non-increasing along jumps from `t_0 = 0`, `c_z >= c_(x') >= a chi - 2||h||/(Np)`. Let `N -> infinity`. ∎

*Where the hypothesis enters.* Only through Lemma L at the finitely many points `z_0, ..., z_(T_max)`, and through `Phi*` there. So local height-summability near these points, together with `sum |c| < infinity`, suffices.

*Sharpness.* For a negative expanding `x`, `T^p n = lambda n - (lambda - 1) x > lambda n`. So `c = a chi` exactly still fails inside the shadow, while any `c > a chi` works there (Proposition S). The forced charge is `a chi`, with strict inequality required.

### 2.3 Consequences for Collatz; the forced mass

**Corollary A1 (PROVED).** For Collatz, no height-summable signed bank, with any bounded `h`, satisfies (S) or (L).
* The backward tree of `-1` has `2^(j-1)` points at depth `j`, each with `c_z >= a log(3/2)` by Theorem A. So `sum |c| = infinity`, contradicting `S_H < infinity`.
* **Explicit witnesses.** Pick `z` in `T^(-inf)(-1)` with `c_z < a kappa`; one exists since `|c_beta| -> 0`. Walk the chain `z, Tz, ..., -1`. Either the charge increases at some step `z_i -> z_(i+1)`, or it is non-increasing and `c_(-1) < a kappa`.
  * In the first case the step `T^i y -> T^(i+1) y` of the good integer at depth `D` violates (S), by `(c_(z_(i+1)) - c_(z_i)) D + O(1)`.
  * In the second case every step inside the shadow of `-1` violates (S), by about `a kappa - c_(-1)`.
* **Checked** (output §C) on four banks: `{-1: 0.9 kappa}`; `{-2^j: 1.2 kappa 2^-j}`; `{-2^j: 1.2 kappa, j <= 20}` (it fails at the uncharged odd preimage `-5/3` of `-2`); and a random nonnegative bank on all expanding periodic points of period `<= 8` plus the tree of `-1` to depth 6.
  * The measured increments at depths 60, 120, 240 are exactly linear, with the predicted slopes `0.24328`, `0.48656` and `0.19323`.
  * Lookahead `L = 5, 20, 50` fails at the good integers near `-8`.

**The forced mass (PROVED + FINITE-EXACT).** On the periodic points alone, the forced charge per period is
`M(p) = sum over expanding primitive necklaces of length p of log(3^a/2^p)`.
* **Bounds.** `M(p) <= log 3 · sum_(a > c_* p) C(p,a) <= log 3 · 2^(h* p)` by the entropy bound, since `c_* > 1/2`. Conversely, the necklaces with `a = ceil(c_* p) + 1` give `M(p) >= 2^(h* p)/poly(p)`. So `M(p) = 2^(h* p + O(log p))`.
* **Table.** `|log2 M(p) - h* p| <= 2 log2 p + 2` for `p = 16..64`.

| `p` | expanding points of period `p` | `M(p)` |
|---|---|---|
| 1 | 1 | 0.405 |
| 3 | 3 | 0.118 |
| 8 | 32 | 5.28 |
| 12 | 768 | 64.2 |
| 20 | 137,800 | 7,810 |
| 32 | 236,611,808 | 1.25·10^7 |
| 64 | 3.02·10^17 | 8.16·10^15 |

(`M(2) = 0`: no primitive expanding word of length 2.)
* **Consequences.**
  * The orbits `O_L` of `1^(L-1)0` alone need total charge `log lambda_L = (L-1) log 3 - L log 2 -> infinity`.
  * **With the forced profile (nonnegative), `R(n) = +infinity` at every positive integer.** Indeed `Phi(n) >= sum over centers in the parity class of n`, and both parity classes carry infinite forced mass: every expanding cycle other than `{-1}` has odd and even points, and half the tree of `-1` is even.
* **Weight profiles.**
  * A profile depending only on the period, `c_x = f(p(x))` with `f` non-increasing, fails. Either `f` stays above `a chi(O_L)`, which tends to `a log(3/2)`, and the mass is infinite; or `f(L) < a chi(O_L)` for some large `L`, and Theorem A applies whenever `sum_p p f(p) P(p) < infinity`.
  * "Finitely many weights above each threshold" is implied by summability, hence covered.
  * The only unexcluded period-only profiles have `sum_p f(p) P(p) < infinity = sum_p p f(p) P(p)` (OPEN; see the Status block).

**Theorem R+ (the strategy cube; PROVED).** For a level-`k` strategy `sigma`, the following are equivalent:
* (a) `sigma` is class (i);
* (d) some rank `a log n + h(n) + sum c_beta v2(n - beta)`, with `h` bounded and a height-summable bank, satisfies (S) or (L) for `T_sigma`.

*Proof.* (a) ⇒ (d): THM-4482 (b), with `c = 0`.
(d) ⇒ (a): suppose `G_sigma` has an expanding cycle `gamma`. Its periodic point `x_gamma` is a rational `!= 0`. The even preimages `2^j x_gamma` (for all large `j`, outside the cycle) are distinct points of `T_sigma^(-inf)(x_gamma)`. By Theorem A for `T_sigma`, each has charge `>= a chi(gamma)`, so `sum |c| = infinity`. If `x_gamma` is a positive integer it cannot be a center, and `c_(x_gamma) = 0 < a chi` already contradicts (i). ∎

So the four obstructions (THM-4482's bounded corrections, finite banks and bounded lookahead, plus height-summable infinite banks) are all the same fact, §5.

### 2.4 Proposition S: the shadow is paid, the seam is not (local non-vacuity)

**Proposition S (PROVED).** Let `x_0, ..., x_(p-1)` be the orbit of an expanding periodic point (of `T`, of the 3n−1 map, or of any `T_sigma`). Fix `c > a chi` and set `mubar = chi a - c < 0`. Define:
* `tau_i = sum_(l != i) v2(x_i - x_l)`, and `K = 2 + max_(i != l) v2(x_i - x_l)`;
* `g_0 = 0` and `g_(i+1) = g_i + mubar - a w(x_i) + c` (consistent around the cycle);
* `h(n) = -c tau_i + g_i` if `n = x_i (mod 2^K)`, and `h(n) = 0` otherwise.

Then `R = a log n + c sum_i v2(n - x_i) + h` satisfies `R(Tn) - R(n) = mubar + a eta(n)` whenever `v2(n - x_i) >= K + 1` for some `i`. Here `eta(n) = log(1 + 1/(3n))` at odd `n` and `0` at even `n`, so the step descends by at least `|mubar|/2` for `n >= 2a/(3|mubar|)`.

*Proof.* At depth `D >= K + 1` from `x_i`, `Phi(n) = cD + c tau_i` exactly: the other counters sit at their mutual valuations, which are smaller than `D`. The step goes to depth `D - 1` from `x_(i+1)`, still a deep node, so

    R(Tn) - R(n) = a log(Tn/n) - c + c(tau_(i+1) - tau_i) + h_(i+1) - h_i
                 = a log(Tn/n) - a w(x_i) + mubar.  ∎

*Checks* (output §D and cube.b). With `c = chi + 0.02`:
* **Shadows.** Every shadow step (317 per cycle) has `R(Tn) - R(n) <= -0.02` for:
  * `-1` (`chi = 0.40547`);
  * the `-5` cycle (`chi = 0.03926`);
  * the `-17` cycle (`chi = 0.00597`, `K = 6`);
  * the **3n−1 cycle `5 -> 7 -> 10`**, with centers at the positive integers themselves: the rank is evaluated only off the cycle;
  * the level-5 sign strategy #23876, whose isolated expanding 2-cycle `{13, 19}` carries the periodic orbit `{1/5, -1/5}` (multiplier `9/4`).
* **Seams.** The entry step from an uncharged preimage fails, with slope exactly `c`. The preimages used are `-2` (of `-1`), `-11/3` (of `-5`), `-35/3` (of `-17`), `11/3` (of `5`, on the 3n−1 sheet) and `2/5` (of `1/5`). The increments are exactly linear in the depth.

This is the "3n−1 map restricted near its cycles" asked for: near each cycle, the natural infinite-depth counter is a consistent rank. What no summable bank can supply is the charge at the dense set of seams.

### 2.5 Theorem C: nonnegative banks are nowhere height-summable

**Theorem C (PROVED).** Let `c >= 0`. Suppose `Phi(n) < infinity` at one odd and one even integer, and `R = a log n + h + Phi` (`h` bounded) satisfies (S) or (L). Then:
* **(i)** `sum_beta c_beta < infinity`;
* **(ii)** for every nonempty open `U` of `Z_2`, `sum_(beta in U) c_beta (1 + log2 H(beta)) = infinity`;
* **(iii)** for every expanding periodic `x`, every `z` in `T^(-inf)(x)` with `T^j z = x`, and every good integer `y` at depth `D = V + j` from `z` (`V` large), `Phi(y) >= a chi(x) V - C(z, L, ||h||)`.

*Proof.*
* **(i).** `Phi(n) >= sum over beta = n (mod 2) of c_beta`, since every such counter is at least 1.
* **(iii).** Follow the orbit of `y` through the shadow of `x` for `V - 1 - L` further steps (under (L), stop at a jump time in the last window of length `L`). The rank does not increase along the chain, and `Phi >= 0`, `|h| <= ||h||`. So `Phi(y) >= a log(T^t y / y) - 2||h||`, and `log(T^t y/y) >= log(T^j y/y) + (t - j - p) chi(x) - O(1)` inside the shadow.
* **(ii).** Suppose `z` in `T^(-inf)(-1)` has a neighbourhood with finite height moment.
  * Then `Phi(y) = sum_(s <= D) mu(z + 2^s Z_2) + tail(y)`.
  * The first sum is `c_z D + o(D)`, because `mu(z + 2^s Z_2 \ {z}) -> 0` for a finite measure.
  * The tail satisfies `tail <= sum over beta in z + 2^D Z_2, beta != z of c_beta (B + log2 H(beta)) -> 0`, by the Liouville bound and the local moment.
  * With (iii) this gives `c_z >= a log(3/2)`.

  By (i), only finitely many `z` can do this. So the closed set `{x : every neighbourhood of x has infinite height moment}` contains `T^(-inf)(-1)` minus a finite set. That set is dense, so the closed set is all of `Z_2`. ∎

Compared with Theorem A, (ii) needs **no** summability hypothesis: for nonnegative banks the shadow lower bound (iii) is unconditional. The escape it leaves is exactly the one Theorem D uses.

### 2.6 Theorem D: nonnegative banks are exactly as hard as Collatz

**Theorem D (PROVED).** The following are equivalent:
* **(a)** every Collatz orbit of a positive integer reaches 1;
* **(b)** there are `a > 0`, `eta > 0` and a nonnegative bank at rational centers (odd denominators, no positive integer), finite at every positive integer, with `R(Tn) <= R(n) - eta` for all `n >= 2`.

*Proof.*
* **(b) ⇒ (a).** `R >= a log 2 > 0` for `n >= 2`, and `R` drops by `eta` at every step while the orbit stays `>= 2`. So the orbit reaches 1 within `R(n)/eta` steps.
* **(a) ⇒ (b).** Let `sigma(n)` be the stopping time (`T^(sigma(n)) n = 1`), so `sigma(Tn) = sigma(n) - 1` for `n >= 2` and `sigma(n) >= log2 n`.
  * **Parameters.** Take `C = 2a`, `F(n) = C sigma(n) - a ln n + K`, `delta <= a/2`, `K >= 5.64 delta`, `eps_m = delta/(m(m+1))`, and atoms `z_m = m + 2^(D_m)/3` with weight `eps_m`. The integers `D_m >= log2 m + 2` are chosen by induction on `m`.
  * **Leak bound.** For `n != m`:
    * if `v2(n - m) < D_m`, then `v2(n - z_m) = v2(n - m) <= log2(n + m)`;
    * otherwise `|n - m| >= 2^(D_m) >= 4m`, so `n >= 2^(D_m)` and `v2(n - z_m) = v2(3(n - m) - 2^(D_m)) <= log2(4n)`.

    Hence the leak `L(n) = sum_(m != n) eps_m v2(n - z_m) <= delta (3 + S_1 + log2 n)`, with `S_1 = sum_m log2 m/(m(m+1)) = 1.1376`.
  * **Well-defined induction.** `L(n)` depends on `D_m` only for `m < n` (the case `n >= 2^(D_m)`). For `m > n` it involves only `v2(n - m) < log2 m < D_m`.
  * **The choice of `D_n`.** Set `D_n = floor((F(n) - L(n))/eps_n)`. Then `Phi(n) = eps_n D_n + L(n)` lies in `(F(n) - eps_n, F(n)]`.
  * **Admissibility.** `D_n >= (F - L)/eps_n - 1 >= log2 n + 2` holds because `F - L - eps_n (log2 n + 3) >= (1.307a - 1.5 delta) log2 n + K - 5.64 delta >= 0`, using `2a sigma(n) >= 2a log2 n`, `a ln n = 0.693 a log2 n`, `L(n) <= delta(4.14 + log2 n)` and `eps_n <= delta/2`.
  * **Descent.** `R(n) = a ln n + Phi(n)` lies in `(C sigma(n) + K - eps_n, C sigma(n) + K]`, so `R(Tn) - R(n) < -C + eps_n <= -(C - delta/2)`. ∎

*Remarks.*
* The centers can be taken inside `T^(-inf)(-1)` (a dense set of "natural" centers): use the preperiodic point whose word agrees with `m` for `D_m` letters, differs at the next letter, then continues `1^inf`. The leak estimates change only by constants (the cancellation case is bounded by `O(log n)` via the Liouville bound).
* The height of `z_m` is about `2^(D_m)`, with `D_m ~ F(m)/eps_m`: **the rank stores the stopping time of `m` in the height of its private center.** By Theorem C this is unavoidable: any working nonnegative bank has infinite height moment in every open set.

*Window check* (output §F). With `a = 1`, `delta = 0.05`, `K = 2.5`, take the finite bank of atoms `m <= M = 125252` (the largest value on the orbits of `2..1000`). Then `min D_m = 99 > log2 M + 2`, so the bank's potential is exact on the window. All 39,889 steps of those orbits satisfy `R(Tn) - R(n) <= -1.9991 <= -C + delta/2 = -1.975`.
* The height moment `sum eps_m D_m` is `7.6·10^4, 1.1·10^6, 1.4·10^7, 1.7·10^7` for `M = 10^3, 10^4, 10^5, 125252`.
* The shadow lower bound (iii) holds at all 107 window points `2^j(2^V - 1)`.

A finite bank that works on a finite window does not contradict THM-4482's Corollary K, which exhibits failures only at larger integers.

### 2.7 The cube side: are there strategies with finitely many expanding periodic points?

A strategy has finitely many expanding periodic points iff every SCC of `G_sigma` that contains an expanding cycle **is** that single cycle ("class (i')").
* **Why it would matter.** Such a strategy would be provable, with unbounded lookahead (sketch, not used: each orbit crosses each isolated expanding cycle at most once, every other cycle contracts, and a Lemma P potential on the contracting SCCs does the rest). By R+ it would still admit **no** height-summable bank, since `T^(-inf)` of its cycle points is infinite. So it would be the cleanest modified map separating "provable" from "bank-provable".
* **Census (FINITE-EXACT).**

| level | class (i) | class (i') | other | isolated expanding SCCs |
|---|---|---|---|---|
| 2 | 1 | 0 | 3 | 0 |
| 3 | 1 | 0 | 15 | 0 |
| 4 | 16 | 0 | 240 | 0 |
| 5 | 1,052 | 0 | 64,484 | 4 |

* **Details.**
  * The class-(i) counts reproduce THM-4474.
  * The four level-5 isolated expanding SCCs are all the 2-cycle `{13, 19}` (periodic orbit `{1/5, -1/5}`), inside strategies that also have non-isolated expanding cycles.
  * CP-SAT (2 workers) finds a level-6 strategy whose loop at `-1` is an isolated SCC (re-verified by BFS). But "loop at `-1` isolated **and** every other cycle contracting" is INFEASIBLE at every level 4–10 (4.6 s in total). The model: a closed reachable set from the exit `2^(k-1) - 1` avoiding `-1`, plus Lemma P potentials at `F_(2^k)` on all edges avoiding `-1`.
  * Whether class (i') is empty at every level is OPEN.

## 3. Q2: adaptive centers

### 3.1 The candidates

The parity-vector map is a 2-adic isometry (Lagarias; re-checked on 1,750 random cases), so `v2(m - x_w)` is the length of the common prefix of the parity word `Q(m)` and `w^inf`. Hence, for the expanding periodic points of period `<= P`,

    Lam_P(m) := max_x v2(m - x) = max over p <= P with Q(m)[0:p] expanding of (length of the p-periodic prefix of Q(m)).

The candidates are (with `a = 1`, `c = 1.1 log(3/2) = 0.446`):
* **A_P** `= log m + c Lam_P(m)`: the "current center", i.e. the nearest expanding periodic point of period `<= P`, as a **max-plus** bank.
* **A_inf**: the same with every period.
* **B_(P,Q)** `= log m + c max_(i<=Q) [i + Lam_P(T^i m)]`: closeness `v2(m - z)` to the expanding **pre**periodic points `z` of preperiod `<= Q`.
* **D_theta** `= max_(0 <= i <= theta log2 m) log T^i m` and **C_theta** `= max_(i <= theta log2 m) A_20(T^i m)`: windowed future envelopes.

### 3.2 Bounded period, bounded preperiod: the seam families (PROVED)

**Proposition Q2.1.** For every `P >= 1`, `c > 0`, `a > 0` and `L`:

    min_(1<=j<=L) [A_P(T^j y_V) - A_P(y_V)] >= c(V - L + 1 - P) - a log 2  ->  +infinity,   y_V = 2^(V+1) - 2.

*Proof.* The parity word of `y_V` is `0 1^V 0 ...`.
* The expanding blocks `Q[0:p]`, `p <= P`, are `0 1^(p-1)` for `p >= 3`, and each has periodic run exactly `p` (the letter at position `p` is 1, not 0). So `Lam_P(y_V) = P` for `P >= 3` and `0` for `P <= 2`, for all `V >= P`. This is `Lam_P(-2)`.
* `T^j y_V = 3^(j-1) 2^(V-j+1) - 1`, so `Lam_P(T^j y_V) >= v2(T^j y_V + 1) = V - j + 1`, and `log(T^j y_V / y_V) >= -log 2`. ∎

The Kuratowski reset `n_H = (2^(H+3) - 13)/9` (`H = 6j + 5`) is the same phenomenon three steps earlier:
* `n_H -> -13/9` 2-adically, and `-13/9` is not periodic, so `Lam_20(n_H) = Lam_20(-13/9) = 22` for `H >= 35`;
* `Lam_20(T^3 n_H) >= H`;
* `A_20(T^3 n_H) - A_20(n_H) = 5.9, 19.3, 46.1, 99.6` for `H = 35, 65, 125, 245`.

**Proposition Q2.2.** With preperiodic centers of preperiod `<= Q` (period `<= P`), take `y = 2^(Q+1)(2^V - 1)`. Then `B(y) <= Q + P` while `B(Ty) >= Q + V`, so `B_(P,Q)(Ty) - B_(P,Q)(y) >= c(V - P) - a log 2`.

*Proof.* `T^i y = 2^(Q+1-i)(2^V - 1)` for `i <= Q`. Its word `0^g 1^V ...` (`g >= 1`) has expanding blocks of length `<= P` only of the form `0^min(g,p) 1^...`, each with run exactly `p`. And `T^(Q+1) y = 2^V - 1` has `Lam >= V`. ∎

Checked: `Q = 5, 10`, `V = 40, 80, 160`, with slope exactly `c`. The preimage `-2^(Q+1)` of `-1` sits just outside the window.

Both are the seam obstruction of Theorem A, seen by max-plus banks: a center set that is **not** backward-closed misses an entry point.

### 3.3 Unbounded period: the counter is an excursion time, and it fails on hovering orbits (PROVED)

Let `M_l(m) = 3^(a_l(m))/2^l` be the prefix multiplier and `l*(m) = max{l >= 1 : M_l(m) > 1}` (0 if none). By the affine form of `T^l`, `M_l > 1` implies `T^l m > m`, so `l*(m)` is at most the last passage of the orbit above `m`. (Equality holds for 2,893 of the 2,997 starts `3 <= n < 3000`; `l*` never exceeds the last passage.)

**Proposition Q2.3.**
* **(a)** `l*(m) <= Lam_inf(m) < 2 l*(m)` whenever `l*(m) > 0`, and `Lam_inf(m) = 0` iff `l*(m) = 0`.
* **(b)** At odd `m`, `l*(Tm) <= l*(m) - 1`. At even `m`, `l*(m) = 1 + max{l : M_l(m/2) > 2}` (0 if no such `l`).
* **(c)** Let `a/p > log_3 2` be in lowest terms, `w` the lower Christoffel word of slope `a/p` (prefix counts `floor(aj/p)`), `lambda = 3^a/2^p`, and `1 <= r <= log 2/log lambda`. There is a positive integer `k` whose parity word is `w^r 0^s (1 0)^inf` for some `s >= 1`. For it, `A_inf(k) - A_inf(2k) >= c r p - a log 2`. Along the upper best approximations `a/p` of `log_3 2` (`1/1, 2/3, 7/11, 12/19, 53/84, ...`), `lambda -> 1` and `p floor(log 2/log lambda) -> infinity`. Hence `sup_m [A_inf(Tm) - A_inf(m)] = +infinity` for every `c > 0`.

*Proof.*
* **(a).** Take `p = l*`: the block `Q[0:l*]` is expanding with run `>= l*`, so `Lam_inf >= l*`. Conversely, let an expanding block `B = Q[0:p]` have run `rr`. Then `Q[0:kp] = B^k` for `k = floor(rr/p)`, and `B^k` is expanding, so `kp <= l*`; also `p <= l*`. Hence `rr < (k+1)p <= 2 l*`.
* **(b).** If `Q(Tm)[0:l]` is expanding, then prefixing a 1 multiplies the multiplier by `3/2`, so `Q(m)[0:l+1]` is expanding. Prefixing a 0 halves it.
* **(c), the integer `k`.** Let `T^(rp)(x) = (3^(ra) x + C)/2^(rp)` on the cylinder of `w^r`. The last 1 of `w^r` contributes `2^(position)` to `C`, so `3` does not divide `C`. Since 2 is a primitive root mod `3^(ra)`, the solutions `s` of `2^(rp+s) = C (mod 3^(ra))` form a residue class mod `2·3^(ra-1)`; take one with `s >= 1` and `2^(rp+s) > C`. Then `k = (2^(rp+s) - C)/3^(ra)` is a positive integer; the inverse branches map `Z_2` into the cylinder of `w^r`, and `T^(rp) k = 2^s`, after which the orbit is `2^(s-1), ..., 1, 2, 1, ...`.
* **(c), the multipliers.** The prefix counts of `w^inf` are `floor(aj/p)`, so `M_j(k) <= lambda^(j/p) <= lambda^r <= 2` for `j <= rp`, and `M_(rp)(k) = lambda^r > 1`. Afterwards it is `lambda^r 2^(-i) <= 2` during the `s` halvings, at most `lambda^r 2^(-s) (3/2) <= 3/2` at the first tail letter, and it keeps falling (factor `3/4` per tail pair). So `M_l(k) <= 2` for all `l`: by (b), `l*(2k) = 0` and `Lam_inf(2k) = 0`. Also `l*(k) >= rp`, so by (a) `Lam_inf(k) >= rp`. The step `2k -> k` changes `log` by `-log 2`. ∎

*Checks* (output §Q2.3). (a) and (b) hold at all 16,911 relevant orbit points (max ratio `Lam/l* = 1.818`). The hovering integers were built for `2/3` (`r = 1, 3, 5`; `s = 3, 267, 22137`), `7/11` (`r = 1, 2`; `k` has 946,509 bits) and `12/19` (`r = 1`). Each follows its word exactly, reaches `2^s`, and has `l*(k) = rp` and `l*(2k) = 0`. The proved lower bound `c p floor(log 2/log lambda) - log 2` is 6.0, 48.4, 431.5 and 12,400 for `2/3, 7/11, 12/19, 53/84`.

*Meaning.* The "nearest expanding periodic point of any period" is dense enough to see the future: `Lam_inf` is, up to a factor 2, the time until the orbit's multiplier drops below 1 for good. So `A_inf` is no longer a static 2-adic invariant but a disguised excursion time. It fails where the multiplier **hovers** between 1 and 2 for a long time, i.e. along near-critical expanding words. This is the unbounded-period face of the seam obstruction: in the step `2k -> k`, the preimage `2k` sees nothing expanding, while `k` sees a long near-critical shadow.

### 3.4 Where the candidates fail on data (EMPIRICAL)

All steps `m -> Tm` with `m >= 3`, before the orbit reaches 1. Entries are "violations / largest increment".

| set | steps | A_1 | A_3 | A_20 | A_inf | B_(3,10) | D_2 | D_8 | D_16 | C_4 | C_8 | C_16 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| orbits of 3..5000 | 254,283 | 34,396 / 5.11 | 53,759 / 4.42 | 45,121 / 10.90 | 39,297 / 36.33 | 71,308 / 4.87 | 9,677 / 1.22 | 757 / 0.41 | 0 | 38 / 2.66 | 0 | 0 |
| path + delay records `<= 10^12` | 49,771 | 6,893 / 9.12 | 11,291 / 7.78 | 10,317 / 11.80 | 9,118 / 126.42 | 14,233 / 7.78 | 2,100 / 1.22 | 108 / 2.43 | 0 | 198 / 4.53 | 23 / 0.65 | 0 |
| 600 random in `[10^11, 10^12]` | 110,295 | 14,132 / 6.44 | 22,623 / 6.20 | 19,810 / 13.58 | 18,881 / 100.11 | 31,091 / 6.20 | 1,358 / 1.22 | 85 / 0.24 | 0 | 29 / 2.00 | 0 | 0 |

* **Every violation of `A_1, A_3, A_20, A_inf` at a point `m >= 11` is a reset** (`Delta Lam >= 0`), on all three sets. This is forced: if `Delta Lam <= -1`, then `Delta R <= log((3m+1)/(2m)) - c < 0` at odd `m >= 11` and `Delta R <= -log 2 - c` at even `m`. Shadow steps, where the counter falls by one, are always paid. This is "the shadow is paid" of §2.4 in max-plus form.
* **The largest violations** (records set):
  * `A_20`: `11.80` at the even step from `m = 9314882` (orbit of 106239), where `Lam_20` jumps from 0 to 28.
  * `A_inf`: `126.42` at the even step from `m = 414028311267890` (orbit of 881715740415), where `Lam_inf` and `l*` jump from 0 to 285. A natural hovering orbit: the multiplier of `m/2` stays in `(1, 2]` for 285 steps.

### 3.5 Windowed future envelopes (PROVED equivalence + EMPIRICAL horizon)

**Proposition Q2.5.** Consider convergent orbits. `D_theta` never increases on the steps from points `>= n_1` iff every orbit point `m >= n_1` (before the orbit reaches 1) attains the maximum of its remaining orbit within `theta log2 m` steps.

*Proof.* The "if" direction holds because the window then always contains the future maximum, which cannot increase. For "only if", let `m` violate the condition. The values before reaching 1 are distinct, so at the first time `t'` at which the future maximum of `m` enters the window, `D_theta` jumps strictly. ∎

A divergent orbit would make `D_theta` increase infinitely often. So "`D_theta` works" is a peak-time conjecture of Collatz strength (it excludes divergence, not cycles).

The measured horizon `theta* = max over orbit points of (time to future peak)/log2 m` is:

| set | `theta*` | violations at `theta*/2` |
|---|---|---|
| orbits of 3..5000 | 9.464 | 9,919 |
| records `<= 10^12` | 12.932 | 376 |
| random `[10^11, 10^12]` | 8.478 | 640 |

`D_theta` has **no** violation at `theta = theta* + 0.01` on each set. `C_16` (the envelope of `A_20`) has none on any set.

### 3.6 Cycle controls

On the contracting 5n+1 cycles through 13 and 17 (period 7, `5^3 < 2^7`), `A_20` (with the 5n+1 expanding words) is finite and periodic. Its increments sum to 0 with a positive maximum, so the harness flags them. On the 3n−1 cycles through 5 and 17, the cycle points **are** positive expanding periodic points, so `Lam = infinity` there (the run fills the whole computed sequence, to 881 letters) and the candidate is not even finite.

## 4. Q3: heights and the product formula (short, typed)

* **PROVED (Liouville).** `v2(n - beta) log 2 <= log(n+1) + h(beta)`. So 2-adic closeness to a rational center is paid for in height, which is exactly why height-summable banks are tame (Lemma L) and why Theorem D must use centers of huge height.
* **Heights of cycle points (PROVED/FINITE).** `H(x_w) <= 3^(a-1) 2^p`, so `h(x_w) <= a log 3 + p log 2`. For near-critical words (`a ~ c_* p`) this is about `2p log 2`. The observed maximum of `log2 H/p` over lengths `<= 14` is `1.543`, attained at `p = 14` by `0 1^13`, a rotation of the dense word `1^13 0` (rate close to `log2 3`). The backward tree of `-1` has `log2 H(z) <= 1.472 j` at depth `j <= 14`.
* **PROVED (shadow identity).** If `v2(n - x) = D` for a periodic `x`, then for `j < D`, `y_j := T^j n - T^j x = (3^(a_j)/2^j) y_0` exactly. The numerator of `y_j` has `v2 = D - j` and `v3` = `v3(numerator of y_0) + a_j`; its prime-to-6 part and its denominator are invariant. So `log|y_j|_inf - log|y_0|_inf = -(log|y_j|_2 - log|y_0|_2) - (log|y_j|_3 - log|y_0|_3)`. Checked exactly on 1,685 shadow steps (`-1`, `-5`, `O_5`, `O_7`, `-17`).
* **PROVED (the forced charge as a ratio of local rates).** For the multiplier `lambda = 3^a/2^p`, `log|lambda|_inf + log|lambda|_2 + log|lambda|_3 = 0`, and `chi = log 2 · log|lambda|_inf / log|lambda|_2`. The rank balance over one period is `a log|lambda|_inf - (c/log 2) log|lambda|_2 <= 0`: a weighted product-formula inequality at the two places the rank uses. The third place, 3, absorbs the difference and is invisible to the rank. Checked for all 1,767 expanding words of length `<= 12`.
* **Answer to "a height-based distance to the nearest rational cycle".** Candidate: `d(n, x) = h(n - x) - v2(n - x) log 2`.
  * For large `n` this is the logarithm of the prime-to-2 part of the numerator of `n - x`, i.e. the odd places.
  * By the shadow identity, it **increases** by `log 3` at every odd step inside a shadow, and only its prime-to-6 part is invariant.
  * So no function of this type decreases along shadows; its decrease is not implied by Collatz; and the only shadow-invariant piece (the prime-to-6 part) carries no descent.

  By the product formula, every height-type distance splits into archimedean size (`log n`), the 2-adic counter (the bank) and odd-place parts, and the dynamics only moves mass among them. Height-based distances collapse to banks. (PROVED for the identities; the conclusion "no useful height distance" is a statement about this natural family only.)
* **ANALOGY (two-place local heights).**
  * `a log n + sum c_beta v2(n - beta) log 2` is the sum of the local Weil (Néron) functions at `infinity` and at 2 of the divisor `a(infinity) + sum (c_beta)(beta)` on `P^1`.
  * The global height (all places) of a finite divisor is `deg · h + O(1)`, blind to Collatz. The rank deliberately drops the odd places, where the dynamics stores the 3-adic gain.
  * The forced charges make the divisor infinite (infinite degree), so no global height exists. Theorem D's bank is an infinite divisor whose points carry the archimedean information in their heights.
  * The counters `v2(x - alpha)` behave like 2-adic local canonical heights of repelling cycles (they drop by exactly `p` per period), in the spirit of Call–Silverman canonical heights. No theorem from that theory is used.

## 5. Q4: which rank classes are ruled out, and what a successful rank must have

**One fact, now in two halves.** THM-4482 identified three obstructions (bounded corrections, finite banks, bounded lookahead) with one fact: an expanding cycle is a positive circulation, so no bounded tension exists. Theorem A splits this into two halves.
* **The circulation half (cycle).** On each expanding cycle the rank needs a logarithmic singularity of strength `>= a chi` in the 2-adic depth, and the rest of the rank must be a tension there (Thm A(ii)(iii)). This is THM-4482's obstruction localized. Proposition S shows it is also sufficient locally.
* **The entry half (backward orbit; new).** In-edges (seams) force the same singularity at **every** point of `T^(-inf)(x)`: a dense, exponentially branching set (`2^(j-1)` points at depth `j`). Bounded, finite, periodic and summable-with-bounded-height potentials cannot carry a dense set of equal singularities. This is the obstruction behind the Kuratowski reset family `n_H`, the centers `alpha_d`, and the even preimages `2^(V+1) - 2`.

**Ruled out (PROVED):**

| rank class | why | reference |
|---|---|---|
| `a log n + h`, `h` bounded (any period) | no singularity | THM-4482 (b'), Bernoulli §4 |
| finite banks + periodic `h` | finitely many singularities | THM-4482 Cor. K |
| bounded lookahead (any of the above) | chains of jumps | THM-4474, THM-4482; Thm A (L) |
| **signed banks with `sum |c|(1 + log2 H) < infinity`, any bounded `h`, stepwise or bounded lookahead** | forced charge on a whole backward orbit | **Thm A, Cor. A1** |
| **nonnegative banks with finite height moment in some open set** | same, via the unconditional shadow bound | **Thm C** |
| **truncated counters with `sum |c| min(N, log2 H) < infinity`** | Remark §2.1 | **Thm A** |
| **max-plus banks at periodic centers of bounded period (current center)** | misses the seam `-2` | **Q2.1** |
| **max-plus banks at preperiodic centers of bounded preperiod** | misses `-2^(Q+1)` | **Q2.2** |
| **max-plus banks at all expanding periodic points (unbounded period)** | an excursion time; hovering integers | **Q2.3** |
| **in the strategy cube: any of the linear classes beyond class (i)** | Thm R+ | **R+** |

**Not ruled out, and exactly characterized:**
* **Nonnegative banks with infinite height moment in every open set.** They exist iff Collatz (Thm D). Any such bank stores archimedean data (stopping times) in the heights of its centers (Thm C forces this).
* **Max-plus / future envelopes with height-dependent horizon `theta log2 m`.** They never increase iff future peaks come within the window (Prop. Q2.5); on tested data, `theta ≈ 13` suffices (EMPIRICAL).

**The minimal feature (PROVED as a necessary condition).**
* **Necessary.** A successful rank must, near every point `z` of `T^(-inf)(x)` for every expanding cycle `x`, grow at least like `a chi(x) v2(n - z) - O_z(1)` along integers (Thm C(iii) for nonnegative corrections; Thm A for height-summable signed ones). Its 2-adic part therefore has **unbounded resolution at a dense set of centers**.
* **Two routes.** Resolution must be paid for either
  * in **height**: centers of unbounded height whose height moment diverges in every open set, i.e. archimedean information disguised as 2-adic data (Thm D is the extreme case: height = stopping time); or
  * by **nonlinear, max-plus aggregation over the future**, i.e. unbounded lookahead with a height-dependent horizon (Prop. Q2.5).
* **Bounded height fails.** "Unbounded 2-adic resolution with bounded height" (height-summable banks, any period or preperiod bound) is impossible.

These two routes are THM-4476's two places and THM-4480's "height beats residues" once more: the obstruction is archimedean, and 2-adic data can only certify descent after being coupled to the archimedean size.

## 6. What this does and does not do; failures

**Does.**
* It answers Q1 with a single mechanism (forced charges on backward orbits), which contains the Bernoulli, Kuratowski §5 and Corollary K obstructions and extends them to all height-summable infinite banks, signed, with bounded lookahead.
* It shows that height-summability is the exact dividing line for nonnegative banks: summable-with-bounded-height is impossible; nowhere-height-summable is Collatz-complete.
* It shows that the natural adaptive center, with unbounded period, is a disguised excursion time, and it gives explicit, proved failure families.

**Does not.** Nothing here bears on the truth of Collatz. Theorem D is an equivalence, not progress on either side. No HYP or THM file is proposed.

**Failures and corrected guesses** (all checked):
* **"An isolated expanding loop gives a map where a finite bank works."** False: even for an isolated loop, the entry from the preimage `2x` needs its own charge (Prop. S seam checks; Thm A). Separately, an isolated `-1` loop with a contracting rest does not exist for levels 4–10 (CP-SAT), and no class-(i') strategy exists for levels `<= 5`.
* **"Random sampling shows the `-1` loop is never isolated."** It found none in 2,000 samples per level for `k = 6..11`, but CP-SAT finds one at level 6 at once. Isolation is rare, not impossible; the impossibility is only with a contracting rest.
* **"`Lam_inf` fails at `2^(V+1) - 2` like `Lam_P`."** False: long periodic blocks `0 1^(p-1)` approximate `-2`, so `Lam_inf(2^(V+1)-2) ≈ V`. The failure moves to the hovering integers (Q2.3).
* **"The unbounded-period counter can be analysed with naive tails."** Q2.3(c) needs the multiplier of the whole future controlled. That is why the witnesses are built to land exactly on `2^s`: 2 is a primitive root mod `3^m`, the same trick as the Kuratowski `alpha_d` family.
* **"`M(p)/(2^(h p) p^(-3/2))` converges."** Not established: the ratio drifts from 1.11 to 2.09 over `p = 16..64`. Only `M(p) = 2^(h p + O(log p))` is claimed.
* **Cycle controls.** A first version asked `A_20` to have summing increments on the 3n−1 cycles. It failed, because those cycle points are expanding **periodic** points, where the counter is infinite. The control was corrected (§3.6).

## 7. Reproduction

```bash
python3 -u 04-computation/experiments/procgen_rank_20260926_run.py > 05-knowledge/results/procgen_rank_20260926.out
```

* **Requirements.** `numpy`, `ortools` (CP-SAT, 2 workers). No compiled code, no caches, no network (the OEIS record lists are embedded as inputs).
* **Cost.** One process. Final run: wall time 92.0 s (Q1 0.5 s, cube 27.6 s, Q2 63.7 s, Q3 0.2 s); peak RSS 381 MB by the runner's own `getrusage` (399 MB by `/usr/bin/time -l`), checked in the output to be below 700 MB. The output contains 80 `ok:` checks.
* **Checks.** Every `ok:` line is a `check(...)` that aborts the run on failure. Unlabelled indented lines print the tables quoted above.
* **Scratch.** Exploration scripts are in `scratch/procgen_rank/` (not committed).
* **SHA-256** (raw bytes):
  * `procgen_rank_20260926_lib.py` `7fd831c01a1f67282bc7f90f1c0ceebded8c11fc98349a7b953d91f965b59cbf`
  * `procgen_rank_20260926_q1.py` `49c14b77d400071dbc41fda68307cc69e875f81f51a603a750e3e53bbb195793`
  * `procgen_rank_20260926_cube.py` `062f804603395fe09abbff2d537ca6c60596f98d700f4195dfd5713792dbb8e4`
  * `procgen_rank_20260926_q2.py` `6fd71ec7844bfac7b0dc802c27d16bca81ad9347ce4854dbe03d0b1e0fa6cce2`
  * `procgen_rank_20260926_q3.py` `d80e102d7e555a68488ad60fc249c4be4ecb575344531b50cbcac6db23e7d602`
  * `procgen_rank_20260926_run.py` `7fb1df408904fdc90120340f4f71eab00783b159e2111c21daa23551f964b519`
  * `procgen_rank_20260926.out` `ab32013335a0348243729640aa13d3d334c2ec31bca5325755db758490552bf2`
