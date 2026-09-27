# Conjecture R holds with a constant: N_m(L) <= 1.0039 A_(m+15)(L) for all m and L

Lane `robin2`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_robin2_20260926_{lib,run}.py`.
Output: [procgen_robin2_20260926.out](procgen_robin2_20260926.out) (one runner, 28 `[OK]` checks, ends with `ALL CHECKS PASSED`).

## Status

**Summary.**
- **Theorem R_15 (PROVED: hand proofs plus computer-assisted finite parts; pre-audited, see below).** For all `m >= 2` and `L >= 1`,

      N_m(L) <= 1.0039 · A_(m+15)(L).

  This is HYP-9142 (Conjecture R) with the fixed shift `c0 = 15` and the constant `K_0 = 1.0039`. It removes the polynomial factors `(m+2)^17` (robin note, Corollary 3) and `4e(m+3)L` (THM-4488 (1)). For `m = 1` the statement is trivial: `N_1(L) = 0` for `L >= 1`.
- **Consequence (given pairpeak Theorems B-C).** Pairpeak Theorem C with `K = 15`, `K_0 = 1.0039` gives
  `π_L <= (1.0039 · 3^15 (27L + 18) + 4L 2^(-L)) ρ^peak_L`, i.e. `π_L = O(L) ρ^peak_L`. This is the `O(L)` factor in HYP-9142's title; it was `O(L^18)` (robin note) and `O(L^3)` (THM-4488). The constant `1.0039·3^15·27 ≈ 3.9·10^8` is large.
- **Bounded `m`, small shifts (Theorem R_small, FINITE-EXACT via the same Proposition 1).** For `c0 ∈ {1,2,3}` and `2 <= m <= 24`, `N_m(L) <= Γ* A_(m+c0)(L)` for every `L`, with `Γ* <= 1.1562 / 1.0474 / 1.0148`.
- **Shift `c0 = 1` for all `m` stays OPEN.** Exact data: `max N_m/A_(m+1) = 1.0032546` at `(m,L) = (7,50)`; `max N_m/A_(m+2) = 1.0000889`; `max N_m/A_(m+3) = 1 + 1.99·10^-6`; `max N_m/A_(m+15) = 1 + 6.8·10^-26` (runner §7, `m <= 60`). The method below needs the top wall of the hard strip far enough (`K` about 15; see §4) to beat its union bounds.

**The mechanism in four sentences.**
1. In the integer coordinate `V_t = U_t - floor(ct)` both walks are lazy nearest-neighbour walks in the same Sturmian environment; the Robin walk differs from the hard-wall walk only at the zone site `m`.
2. An exact telescoping identity writes `N - A` as a sum over zone visits of the gradient `V_k(m-1) - V_k(m)` of the hard-wall count, `k` = remaining time. Once this gradient is `<= 0` for all remaining times `k >= k1(m)` ("eventual monotonicity", EM), only the last `k1(m)` steps can create excess.
3. EM follows from one inequality between two bridge counts (a TP2 argument), which is proved by a leading-zero injection, the FKG inequality on the lattice of `e`-subsets, a real-valued cycle lemma and Hoeffding's bound; `k1(m) ≈ (m+1)/μ` is the natural descent time.
4. On the last `k1(m)` steps the Robin/hard-wall ratio is at most the half-line/hard-wall ratio at the top site `m`, i.e. `1/(1 - q)` with `q` = P(touch the top | survive the bottom). `q <= 3.85·10^-3` by four elementary regimes (Kolmogorov, exact survival, binomial tails, and a bridge decomposition for large `m`).

**Labels.**
- PROVED (hand proofs in §2): Lemmas 1-6, Propositions 1-3, Theorem R_15. The finite parts are:
  - FINITE-EXACT (exact integer computation): EM for `2 <= m <= 25` over all Sturmian factors (1025 landing factors, completeness certified); the exact worst-phase survival for `m <= 60` (§2.7, regime B); Theorem R_small (§2.9, `c0 <= 3`, `m <= 24`).
  - Rigorous floating evaluation of explicit closed forms with large margins (float error `< 10^-12`, margins `>= 1.5·10^-3` relative): the EM criterion at `m = 26`, the Hoeffding sum `G`, the binomial-tail bounds for `60 < m < 2550`, the bridge bound for `2550 <= m <= 10^6`, and a monotone majorant beyond `10^6`.
- CITED: Hoeffding (1963) Theorems 1 and 4; Fortuin-Kasteleyn-Ginibre (1971); Robbins (1955); Kolmogorov's and Ville's maximal inequalities; Morse-Hedlund (1940) (a Sturmian word has exactly `n+1` factors of length `n`).
- FINITE-EXACT (sanity, not used by the proof): the telescoping identity, Lemma 2, TP2, Lemma 3, Lemma 4 (exhaustive on 763 small bridges), Lemma 5 (4000 random sequences), Lemmas 5-6 on exact bridge probabilities; `N_m(L) <= 1.0039 A_(m+15)(L)` exactly for `m <= 20, L <= 2500`; `21 <= m <= 40, L <= 1500`; `41 <= m <= 60, L <= 400`.
- EMPIRICAL: the ratio profiles for `c0 = 1, 2, 3` (§3); the true top-hit ratio is `1 + 1.4·10^-7` on sampled phases, against the proved `1 + 3.9·10^-3`.

**Pre-audit.** A fresh agent audited the note independently, writing its own code without importing the lane's library (`scratch/procgen_robin2/audit/`, not part of the deliverable).
- It rates every lemma, proposition and the assembly SOUND.
- It found no counterexample.
- It checked Proposition 1 exactly (true EM onset, `Γ*` over all factors) for seven `(m,K)` pairs including shift 1, and Lemma 4 on 320 large random bridges.
- It checked regime D on exact instances (`m ∈ {30,60,120}`, `K ∈ {2,5,10,15}`; the exact `q` at `K = 15` is at most `1.3·10^-7`).
- It extended the exact EM check to `m = 26-45, 50, 60, 80`. The minimum ratio falls from 1.0775 to 1.0275, consistent with the analytic part.
- It reproduced all constants.
- It asked for one wording fix in regime D, now made.
- It noted that the constants are evaluated in floating point with margins far above rounding error, not in interval arithmetic.
- §2.9 (Theorem R_small) was added after the audit. It uses only the audited Proposition 1 and Lemma 3 together with exact computation. The auditor's own value `Γ* = 1.0897` for shift 1 at `m = 7` matches §2.9.

Nothing here bears on the Collatz conjecture. These are finite-horizon counting statements. No novelty or priority claim.

## 0. Setting: V-coordinates

- `c = log_3 2`, `d = 1 - c`, `μ = c - 1/2 = 0.130930`. `F_t = floor(ct)` (exact: largest `n` with `3^n <= 2^t`), `δ_t = F_(t+1) - F_t ∈ {0,1}`. The word `δ` is the characteristic Sturmian word of slope `c`; since `c > 1/2` it has no factor `00`.
- A binary word `ε_0 ε_1 ...` drives `V_(t+1) = V_t + ε_t - δ_t`, `V_0 = 0`. With `U_t = ε_0 + ... + ε_(t-1)`: `V_t = U_t - F_t`, and the slope level of the pairpeak/robin notes is `S_t = U_t - ct = V_t - {ct}`. For `t >= 1`, `{ct} ∈ (0,1)`, hence

      S_t > 0  <=>  V_t >= 1,        S_t < M  <=>  V_t <= M.

- **Hard wall.** `A_M(L)` = number of words with `V_t ∈ [1, M]` for `1 <= t <= L-1` and `V_L >= 1`.
- **Robin barrier.** The zone `Z_m = [m-1, m-c)` of the modified walk `S'` is exactly the site `V' = m`. Proof: `S' ∈ Z_m` iff the integer `V'` lies in `[m-1+{ct}, m-c+{ct})`, which contains an integer (namely `m`) iff `{ct} > c`; and `V' = m` is reached only from `m-1` through `δ_(t-1) = 0`, which forces `{ct} = {c(t-1)} + c ∈ (c,1)`. At the zone `δ_t = 1`, so both letters send `m` to `m-1`. Hence `N_m(L)` counts the paths of this lazy walk on `{1, ..., m}` surviving `V >= 1` at times `1..L`, with weight 2 per visit to `m`. Runner §1 checks equality with the pairpeak library for `m <= 14`, `M <= 29`, `L <= 160`.
- **Continuation counts** from site `x` at time `s`, horizon `k`:
  - `V^(s)_k(x)`: `V ∈ [1,M]` at times `s+1..s+k-1`, `V >= 1` at time `s+k` (the hard wall, relaxed at the last step as in `A_M`);
  - `H^(s)_k(x)`: only `V >= 1` at times `s+1..s+k` (half-line);
  - `W^(s)_k(x)`: the Robin weighted count.
- **Kernel.** `P_(s,s+n)(x,y)` = number of words from `(s,x)` to `(s+n,y)` with `V ∈ [1,M]` at times `s+1..s+n`.
- A **zone time** is a time `t` at which the Robin walk can sit at `m`; then `δ_(t-1)δ_t = 01`. The following time `s = t+1` is a **landing time**; landing times satisfy `δ_(s-2)δ_(s-1) = 01`.
- Throughout `M = m + K` with `K = 15`.

## 1. The theorem

**Theorem R_15.** For all integers `m >= 2` and `L >= 1`: `N_m(L) <= Γ · A_(m+15)(L)` with `Γ = 1.0039`.

Parameters used in the proof: `K = 15`, `w = 10^-3`,

    n1(m) = ceil((m+1)/(μ - w)),        k1(m) = n1(m) + 1.

`n1(m)` is the natural descent time `(m+1)/μ` inflated by `0.8%`.

The proof has two halves.
- **Long remaining times (Proposition 2, EM).** For every landing time `s` and every `k >= k1(m)`: `V^(s)_k(m) >= V^(s)_k(m-1)`.
- **Short remaining times (Proposition 3, TB).** For every time `u` and every `1 <= k <= k1(m)`: `H^(u)_k(m) <= Γ V^(u)_k(m)`.

Proposition 1 turns these into the theorem.

## 2. Proofs

### 2.1 Hybrid telescoping and the phase decomposition

For `0 <= t <= L` let `R_t(x)` be the Robin weighted count of paths from `(0,0)` to `(t,x)`, and

    Z_t = Σ_x R_t(x) V^(t)_(L-t)(x)      (Robin for t steps, hard wall for the remaining L - t steps).

**Lemma 1.** `Z_0 = A_M(L)`, `Z_L = N_m(L)`, and

    Z_(t+1) - Z_t = R_t(m) · E^(t+1)_(L-t-1),        E^(s)_k := V^(s)_k(m-1) - V^(s)_k(m).

Hence `N_m(L) - A_M(L) = Σ_(t<L) R_t(m) E^(t+1)_(L-t-1)`.

*Proof.* Expand `R_(t+1)` by the last Robin step and `V^(t)_(L-t)` by the first free step.
- From a non-zone site `x ∈ [0, m-1]` both steps go to `x - δ_t` or `x + 1 - δ_t <= m <= M`, with the same bottom killing, so the contributions agree. (The top wall is inactive there; for `L - t = 1` the relaxed count `V_1` agrees with the Robin step too.)
- From the zone site `m` (then `δ_t = 1`) the Robin walk contributes `2 V^(t+1)(m-1)` and the hard-wall walk `V^(t+1)(m-1) + V^(t+1)(m)`; the difference is `E^(t+1)_(L-t-1)`. ∎

(Runner §2: the identity holds exactly for `(m,K) ∈ {(4,1),(5,2),(7,3)}`, `L <= 70`.)

**Proposition 1 (phase decomposition).** Let `k1 >= 1`. Assume `E^(s)_k <= 0` for every landing time `s` and every `k >= k1`. Put

    Γ* = sup { W^(u)_k(x) / V^(u)_k(x) :  u >= 1,  x ∈ [1,m] reachable by the Robin walk at time u,  0 <= k <= k1 }.

Then `N_m(L) <= Γ* A_M(L)` for all `L >= 1`.

*Proof.*
- **`L > k1`.** For `t <= L - 1 - k1` the remaining time `L - t - 1` is `>= k1`, and `R_t(m) > 0` only at zone times, so every increment of Lemma 1 is `<= 0`: `Z_(L-k1) <= Z_0 = A_M(L)`. By the Markov property of the Robin weights, `N_m(L) = Z_L = Σ_x R_(L-k1)(x) W^(L-k1)_(k1)(x) <= Γ* Σ_x R_(L-k1)(x) V^(L-k1)_(k1)(x) = Γ* Z_(L-k1)`.
- **`L <= k1`.** The first step from `(0,0)` is forced to site 1, so `N_m(L) = W^(1)_(L-1)(1)` and `A_M(L) = V^(1)_(L-1)(1)`. ∎

(Runner §8d: at small shifts `K ∈ {3,4,5}`, with `k1` = the observed EM onset and `Γ*` sampled over start times, `N_m(L) <= Γ* A_(m+K)(L)` holds on `L <= 700`.)

### 2.2 Robin below the half-line; the hard wall relative to the half-line

**Lemma 2.**
- (a) `W^(u)_k(x) <= H^(u)_k(x)` for every reachable `(u,x)`.
- (b) For fixed `u, k`, the ratio `V^(u)_k(x) / H^(u)_k(x)` is nonincreasing in `x ∈ [1, M]`.
- Consequently `Γ* <= sup_(u, k <= k1) H^(u)_k(m) / V^(u)_k(m)`.

*Proof.*
- (a) The Robin path is the free path of the same word minus the number of flips so far. So whenever the Robin path survives, the free path survives the bottom.
- (b) All counts are positive: the word `ε_t = δ_t` keeps `V` constant. Put `r^V(x) = V(x+1)/V(x)` and `r^H` likewise, for `x ∈ [1, M-1]`. We show `r^V_k <= r^H_k` by induction on `k`.
  - For `k <= 1` the two counts coincide, since the top is relaxed at the last step.
  - For the step, write `V^(u)_(k+1)(x) = V(x-δ)1[x-δ >= 1] + V(x+1-δ)1[x+1-δ <= M]` with `V = V^(u+1)_k`, `δ = δ_u`. The same formula without the top indicator holds for `H`.
  - With `α = V(x-δ)1[·]`, `β = V(x+1-δ)`, `γ = V(x+2-δ)1[x+2-δ <= M]`: `r^V_(k+1)(x) = (1 + γ/β)/(1 + α/β)`. Here `γ/β <= γ'/β'` and `α/β >= α'/β'` by the induction hypothesis (or trivially at the walls), where the primes denote the `H`-quantities. So `r^V_(k+1)(x) <= r^H_(k+1)(x)`.
- The consequence: for `x <= m`, `W/V <= H/V = (H(x)/V(x)) <= H(m)/V(m)` by (a) and (b). ∎

(Runner §2: (a) and (b) are checked exactly on samples, `k <= 120`.)

### 2.3 Eventual monotonicity from one bridge inequality

**Lemma 3.** Let `s` be a landing time and `n >= 1`. If `P_(s,s+n)(m,1) >= P_(s,s+n)(m-1,1) > 0`, then `E^(s)_k <= 0` for every `k >= n+1`.

*Proof.*
- **TP2.** For `x < x'` and `y < y'` in `[1,M]`, `P(x,y)P(x',y') >= P(x,y')P(x',y)`. A path from `x` to `y'` and a path from `x'` to `y` live on the same integer grid and move by `ε - δ` with the same `δ`. Their difference changes by at most 1 per step and changes sign, so they meet. Swapping the tails after the first meeting injects such pairs into pairs `(x→y, x'→y')`, keeping the constraint `[1,M]`.
- With `y = 1` and `ρ = P(m,1)/P(m-1,1) >= 1`, TP2 gives `P_(s,s+n)(m,y) >= ρ P_(s,s+n)(m-1,y) >= P_(s,s+n)(m-1,y)` for every `y`.
- By Chapman-Kolmogorov with nonnegative kernels the entrywise inequality persists to every `n' >= n`.
- Finally `V^(s)_k(x) = Σ_y P_(s,s+k-1)(x,y) V^(s+k-1)_1(y)`, so `V^(s)_k(m) >= V^(s)_k(m-1)` once `k - 1 >= n`. ∎

(Runner §2: TP2 on adjacent pairs, and the conclusion on 222 landing times up to `k = 200`. Scratch: the ratio `P_(s,s+n)(m,1)/P_(s,s+n)(m-1,1)` is nondecreasing in `n`; it crosses 1 near `n ≈ (m-1)/μ`.)

### 2.4 Leading zeros and FKG

Fix a landing time `s`, `n >= 1`, and let `e = 1 - m + F_(s+n) - F_s`: the number of ones of a word taking site `m` at time `s` to site 1 at time `s+n`.
- Identify such words with `e`-subsets `A ⊂ {0, ..., n-1}` (positions of ones), taken uniformly.
- `BOK` = {`V >= 1` at times `s+1..s+n`} and `TH` = {`V >= M+1` at some time in `s+1..s+n-1`}.
- `C = BOK ∖ TH` (confined bridges), so `|C| = P_(s,s+n)(m,1)`.
- `f(A) = min A` is the number of leading zeros.

**Lemma 4.** If `(n - e)/(e + 1) <= 1 - P(TH)/P(BOK)` and `P(BOK) > 0`, then `P_(s,s+n)(m,1) >= P_(s,s+n)(m-1,1)`.

*Proof.*
1. **Injection.** A confined bridge from `m-1` has `e+1` ones. Removing its first one gives a confined bridge from `m`.
   - Before that position both words have only zeros, and the new path is the old one plus 1, hence in `[2, m]`.
   - At that position the two paths merge, and they coincide afterwards.
   - A given `A'` has at most `f(A')` preimages (the removed one sat before `min A'`).
   - So `P_(s,s+n)(m-1,1) <= Σ_(A'∈C) f(A')`.
2. **FKG.** Order `e`-subsets by `A ≼ B` iff `b_i <= a_i` for all `i` (sorted positions: `B` has its ones earlier). This is a distributive lattice (componentwise max/min of increasing sequences), and the uniform measure satisfies the FKG lattice condition with equality.
   - The path `V` is pointwise increasing in this order, so `1_BOK` is increasing.
   - `f` is decreasing.
   - FKG (Fortuin-Kasteleyn-Ginibre 1971) gives `E[f 1_BOK] <= E[f] P(BOK)`.
3. **Hockey stick.** `E[f] = Σ_(r>=1) C(n-r, e)/C(n, e) = C(n, e+1)/C(n, e) = (n-e)/(e+1)`.
4. **Conclusion.** `Σ_C f = C(n,e) E[f 1_C] <= C(n,e) E[f 1_BOK] <= C(n,e) E[f] P(BOK) <= C(n,e)(P(BOK) - P(TH)) <= C(n,e) P(C) = |C|`. ∎

(Runner §2: the injection bound, the FKG inequality and the formula for `E[f]` are checked exhaustively on 763 bridges.)

### 2.5 Two bridge estimates

For the bridge of §2.4 put
- `a = m - {cs}`, the S-height at time `s`; `b = 1 - {c(s+n)} ∈ (0,1]`;
- `T = a - b = cn - e`, and `v = T/n` (the descent speed of the chord).

Since `F_(s+n) - F_s ∈ {floor(cn), floor(cn)+1}`, `T ∈ {m-2+{cn}, m-1+{cn}}`; in particular `m - 2 < T < m`.

**Lemma 5 (real cycle lemma).** Let `x_1, ..., x_n` be reals with `max x_i <= x_max` and total `T > 0`. At least `ceil(T/x_max)` of the `n` cyclic rotations have all partial sums `> 0`. Consequently `P(BOK) >= ceil(T/c)/n >= v/c`.

*Proof.*
- Extend the partial sums periodically, `P_(i+n) = P_i + T`. Rotation `j` is good iff `P_i > P_j` for all `i ∈ (j, j+n]`, iff (by periodicity) `j` is a strict future minimum: `P_i > P_j` for all `i > j`.
- The set `F` of future-minimum times is `n`-periodic and nonempty. For `j ∈ F` the next element of `F` is the last minimiser of `P` on `(j, ∞)`, whose value is at most `P_(j+1) <= P_j + x_max`. So consecutive values of `P` on `F` increase by at most `x_max`, and one period of `F` has at least `T/x_max` elements.
- **Application.** Read the bridge backwards from `b`. Its increments are `c - ε ∈ {c, -d}` (max `c`) with total `T`. A good rotation keeps the reversed path above `b > 0`, so `BOK` holds. Rotations preserve the uniform measure, so `P(good) = E[#good rotations]/n`. ∎

(Runner §2: 4000 random sequences, exhaustive over rotations; on the runner's exact bridges `P(BOK)` is at least `1.30` times the bound.)

**Lemma 6 (Hoeffding).** With `G(K, v) = Σ_(t>=1) exp(-2(K + vt)^2/t)`:

    P(TH) <= Σ_(t=1)^(n-1) exp(-2 (M - a + vt)^2 / min(t, n-t)) <= 2 G(K, v).

*Proof.*
- `TH` needs `S > M` at some time `s + t` with `1 <= t <= n-1`, i.e. `U_t - te/n > M - a + vt`.
- `U_t` is a sample of size `t` without replacement from `e` ones and `n - e` zeros. Hoeffding (1963, Theorem 4 with Theorem 1) gives `exp(-2D^2/t)`. The complementary sample gives `exp(-2D^2/(n-t))`.
- Split at `t = n/2`, use `M - a >= K`, and use `vt >= v(n-t)` for `t > n/2`. ∎

### 2.6 Proposition 2 (eventual monotonicity)

**Proposition 2.** For every `m >= 2` and every landing time `s`: `P_(s,s+n1(m))(m,1) >= P_(s,s+n1(m))(m-1,1) > 0`. Hence, by Lemma 3, `E^(s)_k <= 0` for all `k >= k1(m)`.

*Proof.*
- **Positivity.** `P(m-1, 1) > 0`: descend with `ε = 0` at the `δ = 1` steps and then stay with `ε = δ`. The window contains `F_(s+n1) - F_s >= m-2` such steps.
- **Exact part, `2 <= m <= 25` (FINITE-EXACT, runner §3).**
  - `P_(s,s+n)` depends only on the factor `δ_s ... δ_(s+n-1)`.
  - The landing factors are the suffixes of the factors of length `n+2` that begin with `01`.
  - Scanning `δ` finds `n+3` distinct factors of length `n+2`. By Morse-Hedlund this is all of them.
  - The inequality holds for all 1025 landing factors. The smallest ratio `P(m,1)/P(m-1,1)` is `1.0849` (`m = 25`); see the table in §3.
- **Analytic part, `m >= 26`.**
  - Since `(μ - w) n1 >= m + 1` and `T < m`: `e > cn1 - m >= (1/2 + w) n1 + 1`. Hence `(n1 - e)/(e + 1) < (1/2 - w)/(1/2 + w)`.
  - Since `T > m - 2` and `n1 <= (m+1)/(μ-w) + 1`: `v > v_*(m) := (μ - w)(m-2)/(m+2)`.
  - Lemmas 5-6 give `P(TH)/P(BOK) <= 2cG(K, v)/v <= 2cG(K, v_*)/v_*`. This is nonincreasing in `m`, because `G` decreases in `v` and `v_*` increases in `m`.
  - At `m = 26`: `v_* = 0.111368`, `G(15, v_*) <= 2.1613·10^-4`, `(1/2-w)/(1/2+w) = 0.996008`, `2cG/v_* = 2.4489·10^-3`, sum `= 0.998457 < 1` (runner §4).
  - `G` is evaluated as a partial sum to `t = 40000` plus the tail bound `e^(-4Kv - 2v^2(T1+1))/(1 - e^(-2v^2))` (from `(K+vt)^2/t >= 2Kv + v^2 t`), with a factor `1 + 10^-9`.
  - Lemma 4 applies. ∎

### 2.7 Proposition 3 (the top-hit bound on the last k1 steps)

Fix `m`, a time `u` and `1 <= k <= k1(m)`. Consider the free walk (fair letters) from site `m` at time `u`, at S-height `a = m - {cu} ∈ (m-1, m]`. Let
- `BOK_k` = {`V >= 1` at times `u+1..u+k`};
- `TH'` = {`V >= M+1` at some time in `u+1..u+k-1`};
- `q = P(TH' | BOK_k)`.

Then `H^(u)_k(m) - V^(u)_k(m) = 2^k P(TH' ∩ BOK_k)`, so `H/V = 1/(1-q)`.

**Proposition 3.** `q <= q_max := 3.8517·10^-3` for all `m >= 2`, `u`, and `k <= k1(m)`.

*Proof.* Two facts are used throughout.
- `S_t = U_t - ct` has i.i.d. increments `ε - c` with `E 3^(ε - c) = (3^d + 3^-c)/2 = (3/2 + 1/2)/2 = 1`. So `3^(S_t)` is a martingale, and Ville's inequality gives `P(sup_t S_t >= x) <= 3^-x`. (The positive root of `(e^(θd) + e^(-θc))/2 = 1` is `θ* = ln 3` exactly.) In particular `P(TH') <= 3^-(M-a) <= 3^-K`.
- `P(BOK_k)` is nondecreasing in `a` and nonincreasing in `k`.

Let `k_b(m)` be the largest `k` with `m - 1 - μk >= sqrt(k)`.

- **A. `k <= k_b(m)`, all `m`.** `M_t = S_t + μt` is a martingale with increment variance `1/4`. `BOK_k` can only fail if `M_t <= -a + μt <= -(a - μk) < -sqrt(k)` for some `t <= k`. Kolmogorov's inequality bounds this by `1/4`. So `q <= P(TH')/P(BOK) <= (4/3) 3^-15 = 9.3·10^-8`.
- **B. `k_b < k <= k1`, `m <= 60` (FINITE-EXACT).** Letting `a → (m-1)+` gives `P(BOK_k) >= β_m(k) := H^(0)_k(m-1)/2^k`, the survival from S-height exactly `m-1`. It is computed exactly. `q <= 3^-15/β_m(k) <= 8.6·10^-7`; the worst case is `(m,k) = (2,24)`, `β = 0.0812`.
- **C. `k_b < k <= k1`, `60 < m < 2550`.** Here `P(BOK_k) >= P(BOK_(k1))`. Stop at the first `τ` with `a + S_τ <= 0` and use the strong Markov property with Ville's inequality:
  `P(BOK_(k1)) >= P(S_(k1) > D - m + 1) - 3^-D = P(U_(k1) >= F_(k1) + D - m + 2) - 3^-D`.
  - The binomial tail is bounded below by `J` consecutive terms, `J = ceil(sqrt(k1))`. Each term is at least `0.6753 k^(-1/2) e^(-k KL(p||1/2))` (Robbins: `sqrt(2/π) e^(-1/6) = 0.675395`). Here `KL(1/2+y||1/2) <= 2y^2/(1-4y^2)`, from the series `Σ (2y)^(2j)/(2j(2j-1))`.
  - Optimising `D ∈ [3,29]`, the worst case over `60 < m < 2550` is `P_lb = 0.01197` at `m = 68`. So `q <= 5.8·10^-6`.
- **D. `k_b < k <= k1`, `m >= 2550` (bridge decomposition).**
  - Condition on the number of ones `e`. The endpoint is `z_e = a + e - ck`, the chord speed `v_e = (a - z_e)/k`, and `P_e` is the uniform arrangement.
  - Call `e` *good* if `0 < z_e <= a - v_min k`, with `v_min = 0.11`.
  - For good `e`, the proofs of Lemmas 5-6 apply verbatim to the bridge from `a` to `z_e`: `P_e(TH') <= 2G(K, v_e)` and `P_e(BOK) >= v_e/c`.
  - Bad endpoints number at most `B = #{S_k > -v_min k} <= 2^k e^(-k KL(c - v_min||1/2))` (Chernoff).
  - Hence

        q <= 2cG(K, v_min)/v_min + (c/v_min) · B / Σ_(good) C(k,e).

  - The denominator is at least the single term `C(k, e*)`, `e* = floor(ck - a) + 1`. Then `z_(e*) ∈ (0,1]`, `1 <= e* <= k-1`, and `e*` is good because `m - 2 >= v_min k1`.
  - By Robbins, `C(k,e*) 2^-k >= 0.6753 k^(-1/2) e^(-k KL(p*||1/2))` with `|p* - 1/2| <= y_max = max(|μ - (m-2)/k1|, |m/(k_b+1) - μ|)`.
  - So the second term is at most `(c/v_min) R(m)`, where

        R(m) = (sqrt(k1)/0.6753) exp(-(k_b+1)(KL(c - v_min||1/2) - 2y_max^2/(1 - 4y_max^2))).

  - With `h = 2cG(15, 0.11)/0.11 = 2.978·10^-3` and `KL(c - 0.11||1/2) = 8.8·10^-4`:
    - evaluating every `m ∈ [2550, 10^6]` gives `max (h + (c/v_min)R(m)) = 3.8517·10^-3`, at `m = 2550`;
    - for `m > 10^6` use `k1 <= (m+1)/(μ-w) + 2`, `k_b + 1 >= (m - 1 - sqrt((m-1)/μ))/μ`, and `y_max <= Y = 1.0004·10^-3`. The two pieces of `y_max` themselves jitter with the floors and ceilings, but they are bounded by the explicit majorants `μ - (m-2)/((m+1)/(μ-w)+2)` and `mμ/(m-1-sqrt((m-1)/μ)) - μ`, and these do decrease in `m`. This gives a majorant `R_maj(m)` whose logarithm has negative derivative for `m >= 10^6`. Also `m - 2 >= v_min k1` persists because `v_min/(μ - w) = 0.847 < 1`. ∎

(Runner §5; §8 samples the true `q`: `H/V - 1 <= 1.4·10^-7`.)

### 2.8 Assembly

- By Proposition 2 the hypothesis of Proposition 1 holds with `k1 = k1(m)`.
- By Lemma 2 and Proposition 3, `Γ* <= 1/(1 - q_max) = 1.0038665 <= 1.0039`.
- Proposition 1 gives `N_m(L) <= 1.0039 A_(m+15)(L)`. ∎

### 2.9 Small shifts for bounded `m` (FINITE-EXACT, via Proposition 1)

Proposition 1, Lemma 2 and Lemma 3 hold for every shift `K`; only Propositions 2-3 used `K = 15`. For a single `m`, both hypotheses of Proposition 1 can be verified exactly.
- **EM.** Take every landing factor (complete Sturmian factor set of length `N_max = 6(m-1)/μ + 60`) and find the first `n` with `P(m,1) >= P(m-1,1) > 0`. By Lemma 3, EM holds from `k1 = 1 + max` of these onsets.
- **`Γ*`.** The supremum in Proposition 1 depends on the start time only through the factor `δ_(u-1) ... δ_(u+k1-1)`; the first letter decides whether `x = m` is reachable. So `Γ*` is a maximum over the complete set of factors of length `k1 + 1`, computed exactly with rational arithmetic.

**Theorem R_small.** For `c0 ∈ {1, 2, 3}`, every `2 <= m <= 24` and **every** `L >= 1`:

    N_m(L) <= Γ*_(c0)(m) · A_(m+c0)(L),   max_(m<=24) Γ*_1 = 1.156186,   max Γ*_2 = 1.047346,   max Γ*_3 = 1.014780.

(Runner §9, 28th check.) For shift 1 and `4 <= m <= 24` this is, as far as we know, the first bound uniform in `L`. The robin note's Theorems 1-2 need shift 2, and THM-4488 needs shift 5. For `m <= 3` the constant is exactly 1.

| `m` | 4 | 6 | 8 | 10 | 12 | 16 | 20 | 24 |
|---|---|---|---|---|---|---|---|---|
| `c0 = 1`: `Γ*` (EM onset `n1`) | 1.0164 (18) | 1.0750 (35) | 1.1095 (59) | 1.1249 (87) | 1.1371 (119) | 1.1489 (193) | 1.1540 (277) | 1.1562 (368) |
| `c0 = 2` | 1 (16) | 1.0055 (30) | 1.0173 (46) | 1.0286 (62) | 1.0351 (78) | 1.0424 (111) | 1.0459 (144) | 1.0473 (179) |
| `c0 = 3` | 1 (16) | 1 (30) | 1.0012 (44) | 1.0049 (60) | 1.0072 (75) | 1.0118 (106) | 1.0138 (136) | 1.0148 (168) |

- **Beyond `m = 24` (EMPIRICAL, scratch, same code, not re-checked by the runner).**
  - `c0 = 1`: `Γ*(30) = 1.157532`, `Γ*(31) = 1.157642`.
  - `c0 = 2`: `Γ*(30) = 1.048396`.
  - `c0 = 3`: `Γ*(30) = 1.015468`.
  - All increase slowly and appear to converge.
- **EM onset.**
  - For `c0 = 2, 3` it stays near the descent time: `n1 ≈ (0.93-1.04)(m-1)/μ`.
  - For `c0 = 1` it drifts beyond it: `n1 = 2.09 (m-1)/μ` at `m = 24`, `2.34 (m-1)/μ` at `m = 31`.
  - This is why this route does not obviously extend to all `m` at shift 1.
- **Where the maximum sits.** In every case the maximiser of `W/V` is the zone start `x = m` with `k ≈ 3m`.
- **Loss.** For shift 1, `Γ* ≈ 1.157` against the observed `max_L N_m/A_(m+1) = 1.0033`.

## 3. Tables

**EM, exact part (runner §3).** Minimum over all landing factors of `P_(s,s+n1)(m,1)/P_(s,s+n1)(m-1,1)`:

| `m` | 2 | 3 | 5 | 8 | 10 | 12 | 15 | 18 | 20 | 22 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `n1(m)` | 24 | 31 | 47 | 70 | 85 | 101 | 124 | 147 | 162 | 178 | 201 |
| factors | 10 | 13 | 19 | 27 | 33 | 39 | 47 | 55 | 61 | 67 | 75 |
| min ratio | 2.573 | 1.844 | 1.441 | 1.265 | 1.202 | 1.169 | 1.137 | 1.118 | 1.101 | 1.093 | 1.085 |

**Constants of the proof.**

| quantity | value |
|---|---|
| EM criterion at `m = 26` | `0.996008 + 0.002449 = 0.998457 < 1` |
| TB regime A (`k <= k_b`) | `q <= 9.3·10^-8` |
| TB regime B (`m <= 60`, exact) | `q <= 8.6·10^-7` |
| TB regime C (`60 < m < 2550`) | `q <= 5.8·10^-6` |
| TB regime D (`m >= 2550`) | `q <= 3.8517·10^-3` (attained at `m = 2550`) |
| `Γ = 1/(1 - q_max)` | `1.0038665`, stated as `1.0039` |

**Observed ratio profile (EMPIRICAL, exact counts; runner §7).** `max_L N_m(L)/A_(m+c0)(L)` over the verified range (`L <= 2500` for `m <= 20`):

| `m` | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `c0 = 1` | 1.000773 | 1.001629 | 1.003155 | **1.003255** | 1.003204 | 1.002996 | 1.002530 | 1.001516 | 1.000878 | 1.000498 | 1.000274 | 1.000150 |
| `c0 = 2` | 1 | 1 | 1.0000059 | 1.0000250 | 1.0000562 | **1.0000889** | 1.0000705 | 1.0000523 | 1.0000276 | 1.0000117 | 1.0000053 | 1.0000022 |

- For `c0 = 3` the maximum is `1 + 1.99·10^-6` (`m = 12`).
- For `c0 = 15` the ratio exceeds 1 only by `<= 6.8·10^-26` (`28 <= m <= 60`).
- The excess over 1 exists for every shift tested (`c0 <= 6`, scratch, `m <= 40`, `L <= 1500`). It decays by a factor of about 30-50 per unit of `c0`, and in `m` roughly like `e^(-0.3 m)` for `c0 = 1`.
- So `N_m <= A_(m+c0)` with constant exactly 1 is false for every tested `c0`. A constant `> 1` is needed.

## 4. Remarks

- **How tight is `Γ`.**
  - The proved `Γ - 1 = 3.9·10^-3` comes entirely from regime D at `m = 2550`: the Hoeffding union bound (`h = 3.0·10^-3`) plus the Chernoff/Robbins ratio.
  - The sampled true top-hit excess is `<= 1.4·10^-7` (runner §8c), and the observed `N/A_(m+15)` excess is `<= 7·10^-26`.
- **Why `K = 15`.** The shift enters twice, both times through Hoeffding union bounds whose exponent is about `8vK` with `v ≈ μ` (i.e. `e^(-1.05K)` up to a factor `sqrt(K/v^3)`):
  - the EM criterion at `m = 26` (margin `1.5·10^-3`);
  - the regime-D term `h`.

  A sharper route should reach `K ≈ 8-10`: bound the density of the first half of a bridge against i.i.d. letters by an explicit local-limit factor (about 1.8), then use Ville's inequality for the i.i.d. walk. At the descent speed the exponent is then `θ* = ln 3` with no union factor. This was not carried out.
- **The EM onset is the descent time.** Scratch data (40-120 phases per case, `K = 3..40`, `m <= 80`): the ratio `P_(s,s+n)(m,1)/P_(s,s+n)(m-1,1)` is nondecreasing in `n` and crosses 1 at `n ≈ (0.87-1.01)(m-1)/μ`. With exact bridge probabilities the FKG criterion of Lemma 4 already holds at `≈ (0.99-1.00)(m-1)/μ`.
- **Why a constant needs the far top.**
  - The Robin top absorbs `26%` of the tilted mass per zone visit. The hard wall at `m + c0 + 1` absorbs everything, but only rarely.
  - When the top is barely reached (the Gaussian-tail regime `L ~ m^2/2`), the soft Robin top loses slightly less than a hard wall two or three units higher. This is the source of the excess `N_m/A_(m+c0) > 1`.
  - For `c0 = 1` the gradient dip is deep (`min_k V_k(m)/V_k(m-1) = 0.857` at `m = 30`) and lasts until `k ≈ 16m`, well past the descent time `≈ 7.6m`.
  - For `c0 = 3` the dip is `0.984` and ends at `≈ 6.5m` (`m = 120`), approaching `7.6m` from below as `m` grows.
  - Heuristically, for `m >> e^(2 ln3 · c0)` the dip ends at `(1 + O(3^-c0))` times the descent time, where survival is exponentially rare. This is why Proposition 3 needs the bridge decomposition, not just `P(TH)/P(BOK)`.
- **The uniform dip bound** (scratch sketch; not used by the proof and not audited). The ratio map `r_(k+1)(x) = f(r_k(x-δ), r_k(x+1-δ))`, `f(a,b) = a(1+b)/(1+a)`, is monotone. The free eigenfunction `1 - 3^(S-M)` (eigenvalue 2) satisfies the sign conditions of a stationary subsolution. This gives `V_k(m)/V_k(m-1) >= 1 - O(3^-c0)` for all `k`: `0.972` for `c0 = 3`, against the observed `0.984`. Free eigenfunctions cannot give eventual monotonicity: every positive combination with bulk ratio `<= 1` contains eigenvalue-2 components, which dominate forever.

## 5. Dead ends and what replaced them

1. **Suggested route 1** (Doob transform with the Lemma E sine functions and constant loss).
   - Time-independent sub/supersolutions bound `N` and `A` only up to the endpoint factor. `λ^S` is not an eigenfunction, and the survivors end near the bottom where the sine profiles are small.
   - A constant needs matching prefactors in both the half-line regime (`L << m^2`) and the spectral regime. No such pair of explicit functions was found.
   - Replaced by the hybrid telescoping (Lemma 1), which compares the two walks path by path.
2. **Ratio (Riccati) subsolutions built from free eigenfunctions.** They prove the uniform dip bound but cannot see the bottom-driven transient (§4). Replaced by the TP2 envelope (Lemma 3), whose monotonicity in `n` is exact.
3. **The robin note's Theorem 6** (RM with the `F`-slack). It needs eventual monotonicity plus a lower bound on `F_k/V_k` on the window before it. That is again a survival-conditioned estimate. Proposition 1 needs no `F`.
4. **Crude top-hit bound `q <= P(TH')/P(BOK)` on the whole window.** It fails for large `m`, because `k1(m)` exceeds the descent time by a fixed fraction and `P(BOK)` is then exponentially small in `m`. Replaced by the bridge decomposition (regime D).
5. **Ballot lower bounds for bridges.**
   - "Last `ℓ` letters are zeros" plus a Hoeffding union bound gives `c_B ~ 2^-ℓ` and needs `T_0 ~ 200` exact steps plus a without-replacement correction, which forces `m >~ 2000`.
   - The robin note's cycle lemma gives only `1/n`.
   - Both are replaced by the real cycle lemma (Lemma 5): `P(BOK) >= v/c ≈ 0.2` at the descent speed, with no loss in `m`.
6. **Injections from Robin words into hard-wall words.** They are impossible with constant 1 for every tested shift (§3).

## 6. OPEN

- **Conjecture R with shift `c0 <= 14` for all `m`, in particular `c0 = 1`.** Bounded `m` is settled for `c0 <= 3`, `m <= 24` (Theorem R_small). Numerically the sharp constants are `1.0032546` (`c0 = 1`) and `1.0000889` (`c0 = 2`), attained at moderate `m`.
  - Proposition 1 holds for every shift.
  - What is missing for small `c0` is (i) an explicit EM onset (it moves beyond the descent time), and (ii) a conditional top-hit bound on windows where survival is exponentially rare.
- **A smaller shift with this method** (`K ≈ 8-10`), via the local-limit/Ville sharpening of Lemma 6 (§4).
- **Whether `sup_m` of the excess `N_m/A_(m+c0) - 1` is attained at bounded `m`** for every `c0` (observed for `c0 <= 6`).
- **HYP-9140 for the consistent price `δ_L`** is untouched. So are Collatz and the globally consistent pairing.

## 7. Reproduction

```
python3 -u 04-computation/experiments/procgen_robin2_20260926_run.py > 05-knowledge/results/procgen_robin2_20260926.out
```

- **Environment.** Python 3.10.0, standard library only.
- **Imports.** The runner imports `procgen_pairpeak_20260926_lib.py` read-only, for the §1 cross-check of the model.
- **Run.** Wall time 188 s, of which §9 takes about 170 s. Peak RSS 27 MB (`ru_maxrss`). 98 output lines, 28 checks, `ALL CHECKS PASSED`.
- **Runner sections.**
  - 0 constants;
  - 1 model equality;
  - 2 lemma sanity checks;
  - 3 EM exact;
  - 4 EM analytic;
  - 5 top-hit regimes A-D;
  - 6 assembly;
  - 7 exact verification and profiles;
  - 8 consistency checks (8a worst phase; 8b EM directly at `K = 15`; 8c sampled `q`; 8d Proposition 1 at small shifts);
  - 9 Theorem R_small (`c0 ∈ {1,2,3}`, `m <= 24`).
- **Parameters** (top of the runner): `K = 15`, `W = 10^-3`, `M_E = 26`, `M_X = 60`, `M_2 = 2550`, `V_MIN = 0.11`, `M_BIG = 10^6`.

| file | sha256 |
|---|---|
| `04-computation/experiments/procgen_robin2_20260926_lib.py` | `842df233831b308936399b1d563df72ca2e4812e05c3ff8b9f6e38610c20a318` |
| `04-computation/experiments/procgen_robin2_20260926_run.py` | `4243144d3554661db8ac32e185095dcf53b936265350b2e32a9f5340c4dcdd18` |
| `05-knowledge/results/procgen_robin2_20260926.out` | `02ab8aeada4b518de8defba7bfb253a9a1033955451a4400ccdcfe8fc58bfbfe` |
| same output without timing lines (`grep -v '(t='`, `wall time`, `peak RSS`) | `bd13719bcd3d3dd177c422053916d764436d3d4f672ddb0e5b9b995f14812030` |
| `04-computation/experiments/procgen_pairpeak_20260926_lib.py` (imported, unchanged) | `aa915a333d913f2ada76a90ee1dc6602a17179f38b24b7f23267053f39b4ec56` |

**References.**
- W. Hoeffding, Probability inequalities for sums of bounded random variables, JASA 58 (1963) 13-30: Theorem 1 (bounded summands), Theorem 4 and §6 (sampling without replacement).
- C. M. Fortuin, P. W. Kasteleyn, J. Ginibre, Correlation inequalities on some partially ordered sets, Comm. Math. Phys. 22 (1971) 89-103.
- H. Robbins, A remark on Stirling's formula, Amer. Math. Monthly 62 (1955) 26-29.
- M. Morse, G. A. Hedlund, Symbolic dynamics II. Sturmian trajectories, Amer. J. Math. 62 (1940) 1-42.
- Kolmogorov's maximal inequality and Ville's inequality for nonnegative (super)martingales (standard).

Lemma 5 generalises the Dvoretzky-Motzkin cycle lemma to real steps bounded above; it is proved here. The TP2 path-switching argument is the one of THM-4488 §2, re-proved here for the two-sided strip.

## 8. Hypothesis bookkeeping (no files created or edited)

- **HYP-9142.** Suggested status:
  - PROVED with shift 15 and constant 1.0039 (this note, Theorem R_15; pre-audited by a fresh agent);
  - for shifts 1, 2, 3 PROVED for `m <= 24` with constants 1.1562, 1.0474, 1.0148 (Theorem R_small);
  - OPEN for shifts `c0 <= 14` for all `m`, in particular the stated `c0 = 1`.
- **HYP-9140, private price.** Given pairpeak Theorems B-C: `π_L <= (1.0039·3^15(27L+18) + 4L 2^(-L)) ρ^peak_L = O(L) ρ^peak_L`. This improves `O(L^3)` (THM-4488) in the polynomial degree only; the constant is huge. The consistent price stays OPEN.
- **Promotion candidates after an independent audit:**
  - Theorem R_15 with Lemmas 1-6;
  - Lemma 5 (real cycle lemma) and Lemma 3 (EM from one bridge inequality) as reusable tools.
