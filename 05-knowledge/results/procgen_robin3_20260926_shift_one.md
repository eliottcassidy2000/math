# Shift 2 for all m: N_m(L) <= 2.1285 A_(m+2)(L), and N_m(L) <= 1.1632 A_(m+3)(L); shift 1 stays open

Lane `robin3`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_robin3_20260926_{lib,run}.py`. The robin2 library is imported read-only, and the orchestrator's counters are executed read-only in runner §1.
Output: [procgen_robin3_20260926.out](procgen_robin3_20260926.out) (one runner, 30 `[OK]` checks, ends with `ALL CHECKS PASSED`).

## Status

**Summary.**
- **Theorem R_2 (PROVED: hand proofs plus computer-assisted finite parts; see labels).** For all `m >= 2` and `L >= 1`,

      N_m(L) <= 2.1285 · A_(m+2)(L).

  This is HYP-9142 with the uniform shift `c0 = 2`. THM-4513 had `c0 = 15`.
- **Theorem R_3 (same status).** For all `m >= 2` and `L >= 1`,

      N_m(L) <= 1.1632 · A_(m+3)(L).

- **HYP-9142 as conjectured (shift 1, all `m`) stays OPEN.** Shift 1 is known for `m <= 24` (THM-4513, Theorem R_small). Section 4 shows, with data, why this method stops at shift 2.
- **Consequence (pairpeak Theorem C with `K = 2`, `K_0 = 2.1285`).** `π_L <= (517.3 L + 344.9 + 4L 2^(−L)) ρ^peak_L`. With THM-4513 the coefficient of `L` was `3.9·10^8`.

**What is new** (§2):
1. **An EM criterion without `1/P(BOK)` (Lemma A, Proposition EM\*).**
   - Every event `{l_t <= V_t <= u_t}` is a sublattice of the lattice of `e`-subsets, because join and meet act as pointwise max and min of paths.
   - FKG on the top-avoiding sublattice turns robin2's criterion `(n−e)/(e+1) <= 1 − P(TH)/P(BOK)` into `(n−e)/(e+1) <= 1 − P(TH)`.
2. **Two-sided bridge top-hit bounds (Lemmas LR, TOP, LR-, TOP-LOW).**
   - At the hitting time, the likelihood ratio of a bridge prefix against i.i.d. letters lies between explicit bounds.
   - With Ville's inequality, or with the exact i.i.d. hitting law, this gives `P_e(TH) <= e^(1/(6n)) e^(−ηy) + (Hoeffding tail)` and a matching lower bound.
3. **A refined EM criterion (Proposition EM\*\*).** It keeps the top-avoidance gain after `j <= 12` leading zeros. It proves EM at shift 2 for `m >= 150`, where EM\* provably fails.
4. **Decoupled and room-weighted top-hit bounds (Lemma B', Proposition D').**
   - Lemma B' compares survival from `M+1` with survival from `m−1` through the martingale `e^(θS) g(θ)^(−t)`, using two exact survival sequences per `m`.
   - Proposition D' averages the cycle-lemma loss `c/v` over the endpoint heights of the bridge decomposition. Bridges with room `z > j` survive with probability close to 1.

**Constants** (all upper bounds rounded up):

| shift | `m <= 24` | middle range (Lemma B') | beyond (regimes A + D') | overall |
|---|---|---|---|---|
| 2 | 1.047346 (R_small) | `25..149`: 1.640368; `150..1000`: 1.780765 | `m > 1000`: 2.128393 | **2.1285** |
| 3 | 1.014780 (R_small) | `25..500`: 1.163123 | `m > 500`: 1.155568 | **1.1632** |

Without Theorem R_small, shift 3 holds with `2.8251` for all `m >= 2`. At shift 2, Lemma B' fails for `m <= 4`.

**Labels.**
- PROVED (hand proofs, §2): Lemma A, Propositions EM\* and EM\*\*, Lemmas LR, LR-, TOP, TOP-LOW, DIP, B', Propositions EM3, EM2, TB3, TB2 and D', Theorems R_3 and R_2.
  - The robin2 results used are audited parts of THM-4513: Lemmas 1-5, Proposition 1, regime A of Proposition 3, Theorem R_small.
- Computer-assisted parts:
  - FINITE-EXACT, float-rigorous: float64 dynamic programs with nonnegative sums and exact halving. The relative error is `<= 1.01·(steps + size)·2^(−53)`, compared with margins `>= 2.1 n 2^(−53)` and `10^(−8)`. These cover:
    - EM at shift 3 for `2 <= m <= 79` and at shift 2 for `25 <= m <= 149`, over all landing Sturmian factors;
    - Lemma B' for every `2 <= m <= 500` (shift 3) and every `25 <= m <= 1000` (shift 2).
  - Interval arithmetic (mpmath.iv, 40 digits): the analytic EM constant at `m_E = 80` (shift 3), `h_0` in EM\*\*, the D' constants at `m_X = 500` (shift 3) and `1000` (shift 2), and all exponents `η_+`.
  - Float with safety factors `1 ± 10^(−9)`, errors `< 10^(−11)`: the partial sums of EM\*\* (the exact i.i.d. hitting-law dynamic program and the Lemma LR- bounds) and of the DIP series.
  - Monotonicity in `m` carries each analytic constant from its threshold to all larger `m` (§2.6, §2.7, §2.12, §2.13).
- CITED: Hoeffding (1963) Theorems 1 and 4; Fortuin-Kasteleyn-Ginibre (1971); Robbins (1955); Ville's and Kolmogorov's maximal inequalities; Morse-Hedlund (1940).
- FINITE-EXACT sanity checks (not used by the proof; runner §2, §7, §8h):
  - Lemma A on 384 exhaustive bridges;
  - Lemma LR on 3980 exact hitting configurations and Lemma LR- on 2838;
  - Lemma TOP (111 bridges, exact/bound `<= 0.786`) and Lemma TOP-LOW (18 cases, exact/bound `>= 1.093`) against exact bridge probabilities;
  - Lemma B' against the exact `q`; Lemma DIP on 200 exact bridges;
  - both theorems exactly for `m <= 40`, `L <= 800/500`.
- EMPIRICAL: the shift-1 onset drift and the thin worst-phase margin (§4); the growth of `N_3/A_3` (shift 0 is false).
- Heuristic, labelled in the runner: the i.i.d. values of the refined criterion (runner §9e).

**Audit.**
- A fresh agent audited the note with its own code. Every proof item of the shift-3 part is SOUND and every constant was reproduced; see §7.
- The shift-2 extension (Lemmas LR-, TOP-LOW, EM\*\*, D', Theorem R_2) was audited separately; see §7.

Nothing here bears on the Collatz conjecture. These are finite-horizon counting statements. No novelty or priority claim.

## 0. Setting

Everything is as in robin2 §0 ([procgen_robin2_20260926_constant_robin.md](procgen_robin2_20260926_constant_robin.md)). In brief:
- `c = log_3 2`, `d = 1 − c`, `μ = c − 1/2`; `F_t = floor(ct)`, `δ_t = F_(t+1) − F_t`.
- `V_(t+1) = V_t + ε_t − δ_t`, `S_t = V_t − {ct}`. For `t >= 1`: `V_t >= 1 ⟺ S_t > 0` and `V_t <= M ⟺ S_t < M`.
- In `S`-coordinates the free walk has i.i.d. steps `+d` (letter 1) and `−c` (letter 0).
- `A_M(L)`: hard wall, `V ∈ [1, M]` at times `1..L−1` and `V_L >= 1`. `N_m(L)`: Robin barrier with zone site `m`.
- `V^(s)_k(x)`, `H^(s)_k(x)`, `W^(s)_k(x)`: hard-wall, half-line and Robin continuation counts. `P_(s,s+n)(x,y)`: the confined kernel.
- A **landing time** `s` has `δ_(s−2)δ_(s−1) = 01`, so `{cs} ∈ [2c−1, c)`. The closed end is attained only at `s = 2`.
- Throughout `M = m + c0`.
- `n1(m) = ceil((m+1)/v0)` is the EM bridge length and `k1(m) = n1(m) + 1` the top-hit window.

**Parameters** (runner §0, §8).

| | `v0` | `m_E` | `m_X` | `v_min` | `v_small` | `J` (D') |
|---|---|---|---|---|---|---|
| shift 3 | `11/100` | 80 | 500 | `21/200` | `2/25` | 19 |
| shift 2 | `1/10` for `25 <= m <= 149`; `2/25` for `m >= 150` | 150 (EM\*\*, `J = 12`, `t0 = 600`) | 1000 | `3/40` | `1/20` | 40 |

The window parameter may depend on `m`: robin2's Proposition 1 is applied for each `m` separately.

**The exponent `η_+`.** For `0 < p < c` put `φ_p(η) = p e^(ηd) + (1−p) e^(−ηc)`.
- `φ_p` is convex, `φ_p(0) = 1` and `φ_p'(0) = p − c < 0`.
- For `κ >= 0` with `min φ_p < e^(−κ)`, let `η_+(p, κ)` be the larger root of `φ_p = e^(−κ)`.
- `{φ_p <= e^(−κ)} = [η_−, η_+]`. Since `φ_p(η)` increases in `p` for `η > 0`, `η_+` decreases in `p`; it also decreases in `κ`.
- **Certificate.** If `φ_(p_hi)(η_lo) <= e^(−κ)` in interval arithmetic, then `η_+(p, κ)` exists and is `>= η_lo` for every `p <= p_hi`.
- Under i.i.d. letters `P_p`, `exp(η_+(S_t − S_0) + κt)` is a martingale.

## 1. The theorems and the architecture

**Theorem R_2.** For all `m >= 2` and `L >= 1`: `N_m(L) <= 2.1285 A_(m+2)(L)`.

**Theorem R_3.** For all `m >= 2` and `L >= 1`: `N_m(L) <= 1.1632 A_(m+3)(L)`.

Both follow from robin2's architecture:
- **Proposition 1 (robin2).** If `E^(s)_k := V^(s)_k(m−1) − V^(s)_k(m) <= 0` for every landing time `s` and every `k >= k1`, then `N_m(L) <= Γ* A_M(L)`, where `Γ*` is a supremum of `W/V` over the last `k1` steps.
- **Lemma 2 (robin2).** `Γ* <= sup_(u>=1, 0<=k<=k1) 1/(1 − q_(u,k))`, with `q = P(TH' | BOK_k)` for the free walk from site `m` at time `u`:
  - `TH' = {S_t > M` for some `1 <= t <= k−1}`;
  - `BOK_k = {S_t > 0` for `1 <= t <= k}`.
- **Lemma 3 (robin2).** EM for all `k >= n+1` follows from the single bridge inequality `P_(s,s+n)(m,1) >= P_(s,s+n)(m−1,1) > 0`.

For each shift, the proof supplies the two hypotheses:
- **EM.** The bridge inequality at `n = n1(m)` for every `m` and every landing time:
  - shift 3: Proposition EM3, exact for `m < 80`, analytic via EM\* for `m >= 80`;
  - shift 2: Proposition EM2, exact for `m < 150`, analytic via EM\*\* for `m >= 150`.
- **TB.** A bound `q <= q(m) < 1` for all `u >= 1` and `k <= k1(m)`:
  - shift 3: Proposition TB3 (Lemma B' for `m <= 500`, regimes A and D' beyond);
  - shift 2: Proposition TB2 (Lemma B' for `m <= 1000`, regimes A and D' beyond).
- For `m <= 24` both theorems use robin2's Theorem R_small instead.

## 2. Proofs

### 2.1 Lemma A (sublattice FKG) and the EM criterion EM\*

Fix a landing time `s`, `n >= 1`, a start site `x` and `e >= 1`.
- Identify the words of length `n` with `e` ones with `e`-subsets `A = {a_1 < ... < a_e} ⊂ {0..n−1}`, the positions of the ones.
- Order them by `A ≼ B ⟺ b_i <= a_i` for all `i` (in `B` the ones come earlier). This is a distributive lattice, a sublattice of `Z^e` under componentwise min/max.
- The path is `V_t = x + U_t − (F_(s+t) − F_s)`, with `U_t = #{i : a_i < t}`.

**Lemma A.**
- (a) The path of `A ∨ B` (positions `min(a_i, b_i)`) is the pointwise maximum of the paths of `A` and `B`. The path of `A ∧ B` is the pointwise minimum. Hence every event `{l_t <= V_t <= u_t for t ∈ I}` is a sublattice.
- (b) Let `T = {V_t <= M for 1 <= t <= n−1}`, `BOK = {V_t >= 1 for 1 <= t <= n}`, `C = T ∩ BOK`, and `f(A) = min A`. Under the uniform measure,

      E[f 1_C] <= E[f 1_T] P(C)/P(T) <= E[f] P(C)/P(T).

*Proof.*
- (a) `{i : a_i < t}` is an initial segment `{1..U^A_t}`. Also `min(a_i,b_i) < t ⟺ a_i < t` or `b_i < t`. So `U^(A∨B)_t = max(U^A_t, U^B_t)`. Likewise `U^(A∧B)_t = min(U^A_t, U^B_t)`.
- (b) By (a), `T` is a sublattice. So the uniform measure `μ_T` on `T` satisfies the FKG lattice condition, with equality on `T × T` and a zero right side otherwise.
  - `1_BOK` is increasing and `f` is decreasing.
  - FKG (Fortuin-Kasteleyn-Ginibre 1971) gives `E_(μ_T)[f 1_BOK] <= E_(μ_T)[f] μ_T(BOK)`.
  - Multiplying by `P(T)` gives the first inequality; `f 1_T <= f` gives the second. ∎

**Proposition EM\*.** Consider the uniform bridge from site `m` at landing time `s` to site 1 at time `s+n`, with `e` ones, and let `TH = T^c`. If `(n−e)/(e+1) <= 1 − P_e(TH)` and `P_(s,s+n)(m−1,1) > 0`, then `P_(s,s+n)(m,1) >= P_(s,s+n)(m−1,1)`.

*Proof.*
- `P_(s,s+n)(x,1) = |C_x|`, the number of confined bridges.
- robin2 Lemma 4 step 1 (remove the first one of a confined bridge from `m−1`) gives `|C_(m−1)| <= Σ_(C_m) f = C(n,e) E[f 1_(C_m)]`.
- By Lemma A this is at most `|C_m| E[f]/P(T)`.
- `E[f] = (n−e)/(e+1)` by the hockey-stick identity. ∎

(robin2's Lemma 4 needed `(n−e)/(e+1) <= 1 − P(TH)/P(BOK)`, with `P(BOK) ≈ 0.1-0.2`. That is what forced its shift `15`.)

### 2.2 Lemma LR (bridge prefix against i.i.d. letters: upper bound)

For `1 <= e <= n−1` put `p = e/n`, `q = 1 − p`, and `LR_t(i) = C(n−t, e−i) / (C(n,e) p^i q^(t−i))`. This is the ratio of the probability of a fixed prefix of length `t` with `i` ones under the uniform bridge to that under i.i.d. Bernoulli(`p`) letters.

**Lemma LR.** Assume:
- `p ∈ [0.4, 0.6]` and `1 <= t <= n/2`;
- `δ := i − pt >= 2`;
- `e − i >= 1` and `(n−e) − (t−i) >= 1`;
- `w := δ/(n−t) <= 1/4`.

Then `LR_t(i) <= exp(t/n + 1/(6n))`.

*Proof.* Put `f = n − e`, `j = t − i`, `z = t/n`, `ψ(x) = (1−x) log(1−x) + x`, and let `r_N ∈ (1/(12N+1), 1/(12N))` be Robbins' remainders.
1. **Representation.** `LR = [e!/((e−i)! e^i)] [f!/((f−j)! f^j)] [(n−t)! n^t / n!]`. With Robbins' form of Stirling:

       log LR = ½ log[e f (n−t) / ((e−i)(f−j) n)] − eψ(i/e) − fψ(j/f) + nψ(z) + Δr,
       Δr = (r_e − r_(e−i)) + (r_f − r_(f−j)) + (r_(n−t) − r_n).

2. **Robbins.** `r_N` is decreasing: `r_(N') > 1/(12N'+1) >= 1/(12N−11) > r_N` for `N' < N`. So the first two brackets are `<= 0` and `Δr < 1/(12(n−t)) <= 1/(6n)`.
3. **Convexity.** Put `G(δ) = pnψ(z + δ/(pn)) + qnψ(z − δ/(qn))`, so `G(0) = nψ(z)` and `G'(0) = 0`. Since `ψ'' = 1/(1−x) >= 1`, `G'' >= 1/(pqn)`. As `i/e = z + δ/(pn)` and `j/f = z − δ/(qn)`, this gives `−eψ(i/e) − fψ(j/f) + nψ(z) <= −δ²/(2pqn)`.
4. **Square-root term.** It equals `−½ log(1−z) − ½ log[(1 − w/p)(1 + w/q)]`, since `(e−i)/e = (1−z)(1−w/p)` and `(f−j)/f = (1−z)(1+w/q)`.
   - `(1−w/p)(1+w/q) = 1 + w(p−q−w)/(pq) >= 1 − x`, with `x := w(w + |p−q|)/(pq) <= 0.25·0.45/0.24 < 0.47`.
   - So the second part is `<= x/(2(1−x))`.
   - This is `<= δ²/(2pqn) = w²(n−t)²/(2pqn)` iff `(w + |p−q|)/(1−x) <= δ(n−t)/n`. The left side is `< 0.85`; the right side is `>= δ/2 >= 1`.
5. **Assembly.** `log LR <= −½ log(1−z) + 1/(6n) <= z + 1/(6n)`, since `z + ½ log(1−z) >= 0` for `z ∈ [0, ½]`. ∎

The bound is asymptotically sharp as `t/n → 0`. The auditor found ratios up to `0.992`.

### 2.3 Lemma TOP (top hits of a bridge: upper bound)

**Lemma TOP.** Consider a uniform bridge `(n, e)` from `S`-height `a`, with `p = e/n = c − v`, `v > 0`, top level `a + y`, and `TH = {S_t > a+y` for some `1 <= t <= n−1}`. Assume:
- `p ∈ [0.4, 0.6]` and `y >= 2`;
- `2(y+1)/n + v <= 1/4` and `n(c/2 − v) >= y + 2`;
- `η := η_+(p, 1/n)` exists.

Then

    P_e(TH) <= e^(1/(6n)) e^(−ηy) + e^(−v²n − 4vy) / (1 − e^(−4v²)).

*Proof.* Let `τ` be the first `t` with `S_t > a + y`.
- **Hitting configuration.** `U_τ = i_τ` is the unique integer in `(y + cτ, y + cτ + d]`, so `δ = i_τ − pτ ∈ (y + vτ, y + vτ + d]`.
- **Change of measure.** `{τ = t}` is determined by the prefix of length `t`, and both prefix probabilities depend only on its number of ones. Hence `P_e(τ = t) = LR_t(i_t) P_p(τ = t)`.
- **`t <= n/2`.** Here `δ > 2` and `w < 2(y+1)/n + v <= 1/4`. Also `e − i_t > n(c/2 − v) − y − 1 >= 1` and `(n−e) − (t−i_t) > dn/2 >= 1`. So Lemma LR applies, and

      Σ_(t<=n/2) P_e(τ = t) <= e^(1/(6n)) E_p[e^(τ/n); τ < ∞] <= e^(1/(6n)) e^(−ηy).

  The last step: `exp(η(S_t − a) + t/n)` is a nonnegative martingale, and its value at `τ` is `>= e^(ηy + τ/n)` (optional stopping and Fatou).
- **`t > n/2`.** Put `s = n − t`. The last `s` letters are a sample without replacement, and `U_t >= i_t` means their number of ones is `<= ps − δ`. Hoeffding (1963, Theorems 4 and 1) gives `exp(−2(y + v(n−s))²/s)`.
- **Summing the tail.** With `A = y + vn/2` and `x = n/2 − s > 0`, `2(A + vx)²/s >= 4A²/n + 4v²x`. So the tail is at most `e^(−v²n − 4vy)/(1 − e^(−4v²))`. ∎

### 2.4 Lemma DIP (bottom dips of a high bridge)

**Lemma DIP.** For a uniform bridge of length `k` from `S`-height `a` to `z > 0`, with chord speed `v = (a−z)/k > 0`:

    P_e(S_t < 0 for some 1 <= t <= k−1) <= Σ_(s>=1) exp(−2(z + vs)²/s).

*Proof.* `S_t < 0 ⟺ U_t − pt < −(z + vs)` with `s = k − t`. Apply Hoeffding to the complementary sample of size `s`, then sum over `t`. ∎

The series is decreasing in `z` and `v`. It is evaluated as a partial sum plus the tail bound `e^(−4zv − 2v²(S+1))/(1 − e^(−2v²))`.

### 2.5 Lemma B' (a decoupled top-hit bound)

Let `β_h(k)` be the survival probability (`V >= 1` at times `1..k`) of the free walk from `V = h` at time 0, i.e. from `S`-height exactly `h`. For `θ > 0` put `g(θ) = (e^(θd) + e^(−θc))/2`.

**Lemma B'.** For every start time `u >= 1`, every `k >= 1` and every `θ > 0`:

    q_(u,k) <= e^(−θ c0) g(θ)^k max_(1<=r<=k−1) [g(θ)^(−r) β_(M+1)(r)] / β_(m−1)(k).

*Proof.* Site `m` at time `u` has `S`-height `a ∈ (m−1, m)`, and survival is monotone in the `S`-height. Hence `P(BOK_k) >= β_(m−1)(k)`.

For the numerator, let `τ` be the first time with `V >= M+1`.
- `V` is skip-free upward, so `V_τ = M+1`, at `S`-height `<= M+1`.
- By the strong Markov property, `P(TH' ∩ BOK_k) <= Σ_(t<k) P(τ = t) β_(M+1)(k−t)`.
- Write `β_(M+1)(k−t) = g^(−t) · g^t β_(M+1)(k−t)` and use `E[g^(−τ); τ < ∞] <= e^(−θ(M−a)) <= e^(−θc0)`. This holds because `e^(θ(S_t − a)) g^(−t)` is a nonnegative martingale for the fair walk. ∎

In the runner:
- `β_(m−1)` is computed with an extra kill above `m − 1 + 80` (a lower bound);
- `β_(M+1)` is computed with all mass above `M + 81` counted as surviving (an upper bound);
- the bound is minimised over 180 values of `θ ∈ [0.2, 1.095]`. Optimising `θ` separately for each `k` gives the same values.

### 2.6 Proposition EM3 (eventual monotonicity at shift 3)

**Proposition EM3.** For every `m >= 2` and every landing time `s`: `P_(s,s+n1(m))(m,1) >= P_(s,s+n1(m))(m−1,1) > 0`, with `M = m+3` and `v0 = 11/100`.

*Proof.*
- **Exact part, `2 <= m <= 79` (FINITE-EXACT, runner §3).**
  - The landing factors of length `n` are the suffixes of the factors of length `n+2` beginning with `01`. All `n+3` of the latter are found (Morse-Hedlund), so the list is complete.
  - Float64 dynamic programs in probabilities (each step a sum of two nonnegative numbers, then an exact halving) have relative error `<= 1.01 n 2^(−53)`.
  - For every factor, `P(m,1) >= P(m−1,1)(1 + 2.1 n 2^(−53))` and `P(m−1,1) > 2^(−n−1)`. Exact values are multiples of `2^(−n)`.
  - The smallest ratio is `1.0823` (`m = 79`). The auditor confirmed it with exact integers.
- **Analytic part, `m >= 80` (runner §4).** Let `n = n1(m)`.
  - Parameter ranges: `a = m − {cs} ∈ (m−c, m−2c+1]`, `b = 1 − {c(s+n)} ∈ (0,1)`, `T = a − b ∈ (m−1−c, m+1−2c]`, `v = T/n`, `y = M − a ∈ [c0+2c−1, c0+c)`.
  - (i) `(n−e)/(e+1) = (dn+T)/(cn−T+1) < (d+v0)/(c−v0)`, since `T/n < v0`.
  - (ii) `v > v_*(m) := v0(m−1−c)/(m+1+v0)`, which increases in `m`. So `η_+(c−v, 1/n) >= η_lo := η_+(c − v_*(80), 1/n1(80))`.
  - (iii) Lemma TOP applies. Its side conditions are checked at `m = 80` and only improve with `m`: `p ∈ [0.5209, 0.5247]`, `2(y+1)/n + v <= 0.1226`, `y >= 3.26`, `n(c/2 − v0) − (y+2) >= 145`.
  - Every term of `E_f + Top + Late` is nonincreasing in `m`. At `m = 80` (interval arithmetic) the certified total is `<= 0.977998 < 1`, with components `E_f <= 0.919645`, `Top <= 0.056983`, `Late <= 0.001371` (each rounded up). On a grid the value keeps decreasing: `0.972848` at `m = 150`, `0.969417` at `m = 1000`.
  - Proposition EM\* applies. Positivity is robin2's. ∎

### 2.7 Proposition TB3 and regimes A, D, D' (the top-hit bound)

**Proposition TB3.** For every `u >= 1` and `0 <= k <= k1(m)` (shift 3):
- `q <= 0.64603` for `m <= 24`;
- `q <= 0.14025` for `25 <= m <= 500`;
- `q <= 0.13463` for `m > 500`.

*Proof.*
- **`m <= 500`.** Lemma B', float-rigorous, for each `m` (runner §5a). Values, rounded up: `m = 2`: `0.6461`; `5`: `0.3039`; `25`: `0.1403`; `100`: `0.0956`; `500`: `0.0804`.
- **`m > 500`, `k <= k_b(m)`.** robin2 regime A: `P(BOK) >= 3/4` (Kolmogorov) and `P(TH') <= 3^(−3)` (Ville), so `q <= 4/81`.
- **`m > 500`, `k_b(m) < k <= k1(m)`.** Regime D' below, at `m_X = 500`: `q <= 0.13463`. ∎

**Regime D (bridge decomposition).**
- Condition on the number `e` of ones among the `k` letters, with weights `w_e = C(k,e)2^(−k)`. Given `e`, the letters form a uniform bridge from `a` to `z_e = a + e − ck`, with speed `v_e = c − e/k` and `p = e/k = c − v_e`.
- Then `q = Σ_e w_e P_e(TH' ∩ BOK) / Σ_e w_e P_e(BOK)`. Split the `e` with `z_e > 0` into three classes:
  - L: `v_e >= v_min`;
  - H: `v_small <= v_e < v_min`;
  - B: `v_e < v_small`.
- By the mediant inequality, `q <= max(A_L/B_L, A_H/B_H) + W_B/B_L`, where `A_X, B_X` are the class sums and `W_B` is the class-B weight.
- **Lemma TOP for classes L and H.** Its hypotheses hold for every `k > k_b(m) >= k_b(m_X)`:
  - `v_e <= m/(k_b+1) <= μ + (1 + sqrt(k_b+1))/(k_b+1) =: v_max` (by the definition of `k_b`), so `p ∈ [c − v_max, c − v_small] ⊂ [0.4, 0.6]`;
  - `2(y+1)/k + v_e <= 2(c0+2)/k_min + v_max <= 1/4`, and `k(c/2 − v_e) >= k_min(c/2 − v_max) >= c0 + 3`;
  - `η_+` exists and is bounded below by the interval certificate at `p = c − v_min` (class L) or `c − v_small` (class H), with `κ = 1/k_min`.
  - At `m_X = 500`: `p >= 0.482`, `2(c0+2)/k_min + v_max = 0.152`, `k_min(c/2 − v_max) = 562`.
- **Class H.**
  - `z_e > z_*(m) := m − 1 − v_min((m+1)/v0 + 2)`, which increases in `m`.
  - Lemma DIP gives `P_e(BOK) >= 1 − Σ_s exp(−2(z_* + v_small s)²/s)`.
  - Lemma TOP bounds `P_e(TH')`.
  - Together: `A_H/B_H <= Q_H`.
- **Class B.**
  - `W_B <= exp(−2k(c − v_small − ½)²)` (Chernoff, Pinsker).
  - `B_L >= w_(e*) v_(e*)/c` with `e* = floor(ck − a) + 1` (so `z_(e*) ∈ (0,1]`), `v_(e*) >= (m−2)/K1 >= v_min`, and `w_(e*) >= 0.6753 k^(−1/2) e^(−k KL(e*/k||½))` (Robbins).
  - This gives `R_B <= ρ(m) = (sqrt(K1)/0.6753)(cK1/(m−2)) exp(−(k_b+1)Δ)`, with `K1 = (m+1)/v0 + 2`, `Δ = 2(c − v_small − ½)² − 2y²/(1−4y²)` and `y = max(μ − (m−2)/K1, (1+sqrt(k_b+1))/(k_b+1)) >= |e*/k − ½|`.
  - `y` is nonincreasing in `m`, so `Δ` is taken at `m_X`. Since `k_b(m+1) >= k_b(m) + 1`, `ρ(m+1)/ρ(m) <= (1+1/(m+1))^(1.5) e^(−Δ) < 1` once `1.5/(m+1) < Δ` (checked).
- **Class L, plain version.** `P_e(BOK) >= v_e/c` (robin2 Lemma 5) and `P_e(TH') <= Top_L := Top(k_min, c − v_min, c0)`, so `A_L/B_L <= Top_L c/v_min`. This gave `q_D <= 0.43059` in the first version of this note.

**Proposition D' (room-weighted class L).** Write `e = e* + j`, so `z_e ∈ (j, j+1]`; class L is `{j <= j_L}` with `j_L >= J` (side condition `J + 1 <= z_*(m_X)`).
- Let `D_j` be nondecreasing lower bounds for `P_e(BOK)` on class L:
  - `D_j = v_min/c` for `j <= 2` (cycle lemma);
  - `D_j = max(v_min/c, 1 − e^(1/(6k)) e^(−η_L j) − e^(−v_min² k − 4v_min j)/(1 − e^(−4v_min²)))` for `j >= 3`.
- The second form is Lemma TOP for the time-reversed bridge. It dips below 0 iff its partial sums `Σ(ε − c)` exceed `z_e > j`, which is the same problem with `y = z_e`. The side conditions `2(J+2)/k_min + v_max <= 1/4` and `k_min(c/2 − v_max) >= J + 4` are checked.
- Let `r_lo <= w_(e+1)/w_e` for `e* <= e < e* + J`.

Then

    A_L / B_L  <=  Top_L · Φ,        Φ = Σ_(j<=J) r_lo^j / Σ_(j<=J) r_lo^j D_j.

*Proof.*
- `A_L/B_L <= Top_L Σ_L w_e / Σ_L w_e D_(j(e))`.
- **Truncation to `j <= J`.** `D` is nondecreasing, so dropping the terms with `j > J` can only lower the weighted average of `D`.
- **Monotone likelihood ratio.** On `{0..J}` the weights `π_j ∝ w_(e*+j)` and `π'_j ∝ r_lo^j` satisfy `π_(j+1)/π_j >= π'_(j+1)/π'_j`, so `π >=_st π'`. Since `D` is nondecreasing, `E_π[D] >= E_(π')[D]`.
- **The decay bound.** `w_(e+1)/w_e = (k−e)/(e+1)` decreases in `e`. With `e*/k <= p*_hi := c − (m−2)/K1` and `k >= k_min`, `r_lo = (1 − p*_hi − J/k_min)/(p*_hi + (J+1)/k_min)` is a valid lower bound that increases with `m`. ∎

**Shift-3 constants (runner §5b), at `m_X = 500`:**
- `Top_L = 0.07166`, `Φ <= 1.37515` (`r_lo >= 0.89602`), so `Q_L' <= 0.09854`;
- `Q_H <= 0.13428` (dip `<= 5.2·10^(−4)`, `z_* = 20.56`), `R_B <= 3.51·10^(−4)` (`Δ = 0.00425`);
- `q_D' <= 0.13463`.

Every piece is nonincreasing in `m`. On a grid the bound is `0.13452, 0.13394, 0.13338, 0.13286, 0.13269` at `m = 507, 600, 1000, 3000, 10000` (runner §5c; scratch).

### 2.8 Assembly: Theorem R_3

- By Proposition EM3 and robin2 Lemma 3, EM holds with `k1 = n1(m) + 1`.
- By robin2 Lemma 2 and Proposition TB3, `Γ* <= 1/(1 − q(m))`:
  - `1/(1 − 0.140245) <= 1.163123` for `25 <= m <= 500`;
  - `<= 1.155568` for `m > 500`.
- For `m <= 24`, robin2 Theorem R_small gives `1.014780`. (Proposition TB3 alone gives `<= 2.8251`, from `m = 2`.)
- Hence `N_m(L) <= 1.1632 A_(m+3)(L)`. ∎

### 2.9 Lemma LR- (a lower bound for the likelihood ratio)

**Lemma LR-.** Let `1 <= e <= n'−1`, `p' = e/n'`, `q' = 1 − p'`, and `1 <= t < n'`. Let `i` be an integer with `1 <= i <= e−1` and `0 <= j' := t − i <= n'−e−1`. Put `δ = i − p't > 0` and `w = δ/(n'−t) < p'`. Then

    log LR_t(i) >= −(δ²/(2(n'−t)))(1/(p'−w) + 1/q') − w(p'−q')⁺/(2p'q') − 1/(12(e−i)) − 1/(12(n'−e−j')).

*Proof.* Use the representation of Lemma LR, step 1.
- **Square-root term.** It is `>= −½ log[1 + w(p'−q'−w)/(p'q')] >= −w(p'−q')⁺/(2p'q')`.
- **Convexity.** `nψ(z) − G(δ) >= −(δ²/2) sup_(s<=δ) G''(s)`, and `G''(s) <= (1/(p'−w) + 1/q')/(n'−t)`, because `p'n'(1 − z − δ/(p'n')) = (n'−t)(p'−w)`.
- **Robbins.** `r_(n'−t) − r_(n') > 0` and `0 < r_N < 1/(12N)`. ∎

### 2.10 Lemma TOP-LOW (a lower bound for bridge top hits)

**Lemma TOP-LOW.** Take a uniform bridge `(n', e)` with `p' = e/n' < c`, top distance `y' > 0`, and let `τ` be the i.i.d. (`P_(p')`) hitting time. Let `L(t)`, `1 <= t <= t0 < n'`, be nonincreasing with `L(t) <= LR_t(i_t)` whenever `P_(p')(τ = t) > 0`, and put `L(t0+1) = 0`. Then

    P_e(TH) >= Σ_(t=1)^(t0) (L(t) − L(t+1)) P_(p')(τ <= t),

and `P_(p')(τ <= t) >= P_(p_lo)(τ_(y_hi) <= t)` for `p' >= p_lo` and `y' <= y_hi`.

*Proof.*
- The change of measure of Lemma TOP gives `P_e(TH) >= Σ_(t<=t0) L(t) P_(p')(τ = t)`.
- Abel summation, with `P(τ <= 0) = 0`, gives the displayed sum.
- For the monotonicity, couple the letters so that the `p'`-walk dominates the `p_lo`-walk pathwise, and note that the barrier `y' + ct` lies below `y_hi + ct`. ∎

In the runner, `P_(p_lo)(τ_(y_hi) <= t)` is computed by an exact dynamic program; its barrier is integral (`U_t > c0 + c(j+1+t) ⟺ U_t >= c0 + F_(j+1+t) + 1`).

### 2.11 Proposition EM\*\* (the refined criterion)

**Proposition EM\*\*.** Take the landing bridge `(n, e)` from site `m`, with `n >= (m+1)/v0` (so `(n−e)/n < d + v0`) and `P_(s,s+n)(m−1,1) > 0`. Let `h_0 >= P_e(TH)`. For `1 <= j <= J`, let `h_j` be a lower bound for the top-hit probability of the bridge that remains after `j` leading zeros. That bridge has length `n−j`, `e` ones, start `a − cj` and top distance `y + cj`. If

    [ Σ_(j=1)^J (d+v0)^j (1 − h_j) + (d+v0)^(J+1)/(1 − d − v0) ] / (1 − h_0)  <=  1,

then `P_(s,s+n)(m,1) >= P_(s,s+n)(m−1,1)`.

*Proof.*
- By Lemma A, `E_(C_m)[f] <= E_T[f] = Σ_(j>=1) P(f >= j) P(T | f >= j)/P(T)`.
- `P(f >= j) = C(n−j,e)/C(n,e) <= ((n−e)/n)^j <= (d+v0)^j`.
- Given `f >= j` (the first `j` letters are zeros), the path descends during those steps. The remaining letters are a uniform `e`-subset, so `P(T | f >= j)` is the top-avoidance probability of the remaining bridge: `<= 1 − h_j` for `j <= J`, and `<= 1` otherwise.
- Hence `E_(C_m)[f] <= 1`, and robin2 Lemma 4 step 1 gives the claim. ∎

EM\* is the case `J = 0`. At shift 2 it cannot work: the best asymptotic value, over `v0`, of `(d+v0)/(c−v0) + e^(−η(v0) y_lo)` is `1.0293` (runner §9b).

### 2.12 Proposition EM2 (eventual monotonicity at shift 2)

**Proposition EM2.** With `M = m+2`:
- for `25 <= m <= 149`, the bridge inequality holds at `n1(m) = 10(m+1)`;
- for `m >= 150`, it holds at `n1(m) = ceil((m+1)/(2/25))`.

This holds for every landing time.

*Proof.*
- **`25 <= m <= 149` (FINITE-EXACT, runner §8a).** As in EM3; the smallest ratio is `1.0535` (`m = 149`).
- **`m >= 150` (runner §8b).** Proposition EM\*\* with `v0 = 2/25`, `J = 12`, `t0 = 600`:
  - `h_0 <= 0.22787`, from Lemma TOP with `y >= y_lo = 2c + 1` in interval arithmetic, as in EM3;
  - `h_j` from Lemma TOP-LOW, with the i.i.d. law at `(p_lo, y_hi) = (c − v0, 2 + c(j+1))`;
  - `L(t) = exp(`Lemma LR- at `n' = n − j`, `p' ∈ [c − v0, (c − v_*)(1 + J/(n−J))]`, `δ <= 2 + c(j+1) + v0 t + 1)`, made nonincreasing;
  - result: `h_1, ..., h_4 >= 0.0923, 0.0585, 0.0371, 0.0234`, and the criterion value is `<= 0.98059 < 1` at `m = 150`.
- **Monotonicity in `m`.** `h_0` is nonincreasing (as in EM3). For the `h_j`: `n' = n1(m) − j` grows with `m`; the box `p' ∈ [c − v0, (c − v_*(m))(1 + J/(n−J))]` shrinks, since `v_*` grows and `J/(n−J)` falls; `δ_max(t)` and the i.i.d. laws at `(p_lo, y_hi)` do not depend on `m`; and each term of the Lemma LR- bound is nondecreasing in `n'` and in the box. So the `h_j` are nondecreasing and the criterion value is nonincreasing in `m`. On a grid it decreases steadily: `0.98017` at 153, `0.97569` at 198, `0.96834` at 400, `0.96310` at 1500, `0.96177` at 5000 (scratch; the auditor scanned 511 values of `m`). ∎

### 2.13 Proposition TB2 and Theorem R_2

**Proposition TB2.** For every `u >= 1` and `0 <= k <= k1(m)` (shift 2, `k1 = n1 + 1` with the `n1` of EM2):
- `q <= 0.39039` for `25 <= m <= 149`;
- `q <= 0.43845` for `150 <= m <= 1000`;
- `q <= max(4/27, 0.53017)` for `m > 1000`.

*Proof.*
- **`m <= 1000`: Lemma B' (runner §8d).**
- **`m > 1000`.** Regime A (`q <= (4/3)3^(−2) = 0.14815`), and Proposition D' at `m_X = 1000` (runner §8e): `Top_L = 0.28380`, `Φ <= 1.86815` (`r_lo >= 0.79544`), `Q_L' <= 0.53017`, `Q_H <= 0.43160`, `R_B <= 2.7·10^(−21)`. The bound is monotone in `m`: `0.52970, 0.52689, 0.51835, 0.50705, 0.49948` at `m = 1013, 1100, 1500, 3000, 10000`. ∎

**Theorem R_2 (proof).**
- robin2's Proposition 1 with EM2 and TB2 gives `Γ* <= 1/(1−q)`:
  - `<= 1.640368` for `25 <= m <= 149`;
  - `<= 1.780765` for `150 <= m <= 1000`;
  - `<= 2.128393` for `m > 1000`.
- robin2 Theorem R_small gives `1.047346` for `m <= 24`.
- Hence `N_m(L) <= 2.1285 A_(m+2)(L)`. ∎

## 3. Tables

**EM, exact parts.** Minimum over all landing factors of `P(m,1)/P(m−1,1)` at `n = n1(m)`:

| shift 3, `m` | 2 | 3 | 5 | 10 | 20 | 40 | 60 | 79 |
|---|---|---|---|---|---|---|---|---|
| `n1(m)` | 28 | 37 | 55 | 100 | 191 | 373 | 555 | 728 |
| min ratio | 2.6225 | 1.9108 | 1.4910 | 1.2574 | 1.1527 | 1.1051 | 1.0899 | 1.0823 |

At shift 2 (`n1 = 10(m+1)`, `25 <= m <= 149`) the smallest ratio is `1.0535`.

**EM, analytic parts.**
- Shift 3, EM\*: `<= 0.977998` at `m = 80`, then `0.972587`, `0.969564`, `0.968908` at `m = 160, 800, 8000`.
- Shift 2, EM\*\*: `<= 0.98059` at `m = 150`, then `0.97073`, `0.96310` at `m = 300, 1500`.

**Lemma B'** (values rounded up).

| shift 3, `m` | 2 | 3 | 5 | 12 | 24 | 25 | 50 | 100 | 200 | 500 |
|---|---|---|---|---|---|---|---|---|---|---|
| `q_B'` | 0.6461 | 0.4084 | 0.3039 | 0.1904 | 0.1386 | 0.1403 | 0.1111 | 0.0956 | 0.0873 | 0.0804 |

| shift 2, `m` | 25 | 50 | 100 | 149 | 150 | 300 | 500 | 1000 |
|---|---|---|---|---|---|---|---|---|
| `q_B'` | 0.3904 | 0.3199 | 0.2848 | 0.2768 | 0.4367 | 0.4232 | 0.4157 | 0.4129 |

The jump at `m = 150` is the switch from `v0 = 1/10` to `v0 = 2/25`, i.e. to longer windows.

**Exact truth** (auditor; exact over all phases and all `k <= k1`, shift 3): `sup q = 0.0843, 0.0699, 0.0649, 0.0619` at `m = 25, 60, 100, 150`. The largest ratio of exact to `q_B'` is `0.841`.

**Regime D' constants.**

| | `m_X` | `Top_L` | `Φ` | `Q_L'` | `Q_H` | `R_B` | `q_D'` |
|---|---|---|---|---|---|---|---|
| shift 3 | 500 | 0.07166 | 1.37515 | 0.09854 | 0.13428 | `3.51·10^(−4)` | 0.13463 |
| shift 2 | 1000 | 0.28380 | 1.86815 | 0.53017 | 0.43160 | `2.7·10^(−21)` | 0.53017 |

**Observed.**
- `max N_m/A_(m+2) = 1 + 8.894·10^(−5)` at `(m,L) = (9,58)`.
- `max N_m/A_(m+3) = 1 + 1.985·10^(−6)` at `(12,77)` (runner §7, §8h; `m <= 40`).
- These are unchanged from robin2 §3.

## 4. Why the method stops at shift 2

**Shift 1: eventual monotonicity is late and thin (EMPIRICAL, runner §9a, §9d; scratch).**
- The bridge crossing `P_n(m,1) >= P_n(m−1,1)`, maximised over landing phases (including those closest to `{cs} = 2c−1`), happens at `n/((m−1)/μ)` equal to:
  - `1.27, 1.91, 2.58, 3.08` for `m = 10, 20, 40, 80`;
  - scratch: `3.27` (`m = 120`) and `3.37` (`m = 160`); a fit `C_∞ − b/m` gives `C_∞ ≈ 3.65`.
- At the worst phase (the zone top, zone site at distance `→ 1.262` below the wall) the eventual gradient ratio `V_k(m)/V_k(m−1)` is only `1.062, 1.027, 1.017, 1.014`. The margin is about 1%.
- Heuristically (not proved) this ratio tends to `e^(θ_0) u(y)/u(y+1)`, with `e^(θ_0) = c/d`, where `u(y) = y + 1 − E[{cτ}]` is the harmonic function of the zero-drift tilted walk killed at the top. A proof along these lines would need the law of the hitting phase to about `0.1%`.
- On EM-sized windows (`~3.6` descent times) the top-hit tools fail:
  - Lemma B' gives `q_B' = 3.05, 2.17, 1.91, 1.80` (`m = 10, 20, 40, 80`);
  - in D' the weights decay fast (`r_lo ≈ 0.67`) while the needed room is large (`η ≈ 0.3`), so `Φ` blows up (scratch estimate `Q_L' ≈ 3`);
  - the i.i.d. value of the refined criterion at `v ≈ 0.035` is `1.14-1.28 > 1` (scratch).

**Shift 0 is false (FINITE-EXACT values, EMPIRICAL growth, runner §9c).** `N_3(L)/A_3(L) = 381.9` at `L = 1000` and `1.593·10^5` at `L = 2000`.

## 5. Dead ends and what replaced them

1. **robin2's FKG step with `P(C) >= P(BOK) − P(TH)`.** Replaced by FKG on the top-avoiding sublattice (Lemma A).
2. **Hoeffding union bounds for top hits.** Replaced by the likelihood ratio at the hitting time plus Ville's inequality (Lemmas LR, TOP), which have the exact exponential rate and no union factor.
3. **Proposition EM\* at shift 2.** It fails asymptotically (`>= 1.0293`). Replaced by EM\*\*, which needs lower bounds on hitting probabilities (Lemmas LR-, TOP-LOW).
4. **The top-hit bound `q <= P(TH')/β` (robin2 regimes B/C).** It fails once survival is rare. Replaced by the decoupling at the hitting time (Lemma B').
5. **The cycle-lemma loss `c/v` in class L.** Plain regime D fails at shift 2 and gives `0.43` at shift 3. Replaced by averaging over the endpoint heights (Proposition D').
6. **An upper bound on bridge survival (a reverse cycle lemma).** Not found, and not needed after item 5.
7. **A single long-window parameter at shift 2.** With `v0 = 2/25` for all `m >= 25`, Lemma B' gives `q = 0.587` at `m = 25`. Using `v0 = 1/10` below `m = 150`, where exact EM allows the shorter window, lowers it to `0.390`.
8. **A bug fixed during the work, in sanity checks only.** The first version of the exact bridge top-hit routine (runner §2c; also scratch explorations) truncated the all-paths array at `top + 2`, so it underestimated long-bridge top-hit probabilities. The proof never used it. The corrected routine covers every reachable site; against it Lemma TOP holds with ratio `<= 0.786`, and Lemma TOP-LOW with ratio `>= 1.093`.

## 6. OPEN

- **HYP-9142 as conjectured (shift 1).**
  - The data give the constant `1.0033`. Theorem R_small covers `m <= 24`.
  - Through Proposition 1, EM at shift 1 needs an onset of `~3.6` descent times and resolves a ~1% effect at the zone-top phase. The top-hit tools of this note fail on such windows (§4).
  - Two routes to try: phase-resolved EM\*\* with exact i.i.d. hitting laws in both directions and the overshoot included; or a Robin-adapted comparison function in Proposition 1.
- **Better constants.**
  - At shift 2 the bottleneck is `Q_L' = 0.530` (`m > 1000`), then `q_B' = 0.4385` (`m = 161`); at shift 3 it is `q_B' = 0.1403` (`m = 25`).
  - The exact `q` is about `0.06-0.08`.
- HYP-9140 for the consistent price `δ_L`, Collatz, and the globally consistent pairing are untouched.

## 7. Reproduction and audit

```
python3 -u 04-computation/experiments/procgen_robin3_20260926_run.py > 05-knowledge/results/procgen_robin3_20260926.out
```

- **Environment.** Python 3.10.0, numpy 2.2.6, mpmath 1.3.0.
- **Run.** Wall time 110.5 s, peak RSS 96 MB (`ru_maxrss`), 70 output lines, 30 checks, `ALL CHECKS PASSED`.
- **Runner sections.**
  - 0 constants;
  - 1 counters equal the orchestrator's;
  - 2 lemma sanity checks;
  - 3-4 EM at shift 3;
  - 5 top-hit bound at shift 3;
  - 6 assembly;
  - 7 exact verification;
  - 8 Theorem R_2;
  - 9 diagnostics (shift 1, shift 2 criterion, shift 0).

**Audit 1** (fresh agent, own code, no import of this lane's library; `scratch/procgen_robin3/audit/`, scripts `a01`-`a09`). It covered the shift-3 part (§2.1-§2.8).
- **Verdicts.** Every item SOUND: Lemmas A, LR, TOP, DIP, B', Propositions EM\*, EM3, TB3, and the assembly with robin2's Proposition 1 and Lemma 2.
  - It noted that Proposition 1's `L <= k1` branch needs TB for all `0 <= k <= k1`, which the proof supplies.
  - It checked that the landing-phase range `[2c−1, c)` is right, with its closed end attained at `s = 2`, the worst phase.
- **Independent numbers.**
  - EM exact integers: `1.257393`, `1.122006`, `1.082319` at `m = 10, 30, 79`.
  - Every `m = 2..79` passes, and so do `m = 80, 100, 150, 200`.
  - Lemma B' values (`0.646023, 0.140245, 0.095600, 0.080361`).
  - The EM total `0.9779979` and all regime-D constants.
  - The exact supremum of `q` over all phases.
  - `max N_m/A_(m+3) = 1 + 1.9851·10^(−6)` at `(12,77)`, `m <= 20`, `L <= 500`.
  - Shift-2 exact EM for `m = 25..150` (min `1.053371`) and the shift-2 Lemma B' maximum `0.390381`.
- **Issues raised, all cosmetic, now fixed.**
  - Upper bounds rounded to nearest instead of up.
  - A missing list of side conditions in §2.7.
  - An unused hypothesis in Lemma LR.
  - The `.out`/hash bookkeeping.

**Audit 2** (same agent, own code `b01`-`b06`, no lane imports; the shift-2 extension, §2.9-§2.13, and the use of D' in Theorem R_3). It audited lib sha `7f877f6c` and run sha `fc00ab31`, which are the versions below.
- **Verdicts.** Every item SOUND: Lemmas LR- and TOP-LOW, Propositions EM\*\*, D', EM2, TB2, and Theorems R_2 and R_3 with D'.
  - Choosing the window `k1(m)` separately for each `m` is legitimate, because Proposition 1 and Lemmas 2-3 are per `m`.
- **Independent numbers.**
  - EM\*\* at `m = 150`: `h_0 = 0.2278647`, `h_1..h_4 = 0.09234, 0.05856, 0.03711, 0.02343`, value `0.9805898`; `0.9707342` at 300, `0.9631009` at 1500.
  - The EM\*\* value never increases over 511 values of `m` in `[150, 3000]`.
  - D': shift 3 `q = 0.1346244`; shift 2 `q = 0.5301619`. All side conditions hold.
  - Lemma B' at shift 2: max `0.390381` over `25..150` and `0.4384433` over `150..1000` (at `m = 161`).
  - `Γ_2 = 2.128393`, `Γ_3 = 1.1631222`.
- **Exact cross-checks.**
  - Lemma LR- on 6000 exact configurations, with no violation.
  - Lemma TOP-LOW on 8 landing bridges at `m = 150`, `j <= 12`, using an untruncated top-only DP; the bound/exact ratio is `<= 0.93`.
  - Exact `E_T[f] <= 0.907`.
  - P_e(BOK) >= D_j and the ratio condition on the weights hold on exact bridges. The exact class-L ratio is `<= 0.051` (bound `0.0985`, shift 3) and `<= 0.207` (bound `0.530`, shift 2).
  - Exact EM with the `2/25` windows at `m = 150, 161, 200`.
- **Issues raised, now fixed.**
  - Required: restate Lemmas LR and TOP with `δ >= 2`, `y >= 2`. The proof covers it, since the right side is `>= δ/2 >= 1 > 0.85`.
  - Wording of TOP-LOW: `L(t) <= LR_t(i_t)` only where `P(τ = t) > 0`.
  - State why the EM\*\* bound is monotone in `m`: `n'` grows, the `p'`-box shrinks, `δ_max` and the i.i.d. laws do not depend on `m`, and each Lemma LR- term is nondecreasing in `n'`.
  - Round five upper bounds up.

| file | sha256 |
|---|---|
| `04-computation/experiments/procgen_robin3_20260926_lib.py` | `7f877f6ccc480824fe89c552184b7ffcd7758b3da510d8bee548d8a060d1b336` |
| `04-computation/experiments/procgen_robin3_20260926_run.py` | `fc00ab311d4bed835cff2bf6d6763e68c4bdf866c1cce3f23ae108a357d0175b` |
| `05-knowledge/results/procgen_robin3_20260926.out` | `42f498334881c06497c68410ebc3864dc0357c2fa26559d03b09db7a92079472` |
| same output without timing fields (see below) | `a799795a83d470683f1c41128a878728236b3c397477ac2336aa4cdcb011e559` |
| `04-computation/experiments/procgen_robin2_20260926_lib.py` (imported, unchanged) | `842df233831b308936399b1d563df72ca2e4812e05c3ff8b9f6e38610c20a318` |

The timing-free hash is taken over `sed -E 's/\(t=[0-9.]+s\)//g; s/\([0-9.]+s(, rounded up)?\)//g; s/ \([0-9.]+s\)$//' out | grep -v "wall time\|peak RSS"`.

**References.**
- W. Hoeffding, Probability inequalities for sums of bounded random variables, JASA 58 (1963) 13-30.
- C. M. Fortuin, P. W. Kasteleyn, J. Ginibre, Correlation inequalities on some partially ordered sets, Comm. Math. Phys. 22 (1971) 89-103.
- H. Robbins, A remark on Stirling's formula, Amer. Math. Monthly 62 (1955) 26-29.
- M. Morse, G. A. Hedlund, Symbolic dynamics II. Sturmian trajectories, Amer. J. Math. 62 (1940) 1-42.
- Ville's and Kolmogorov's maximal inequalities (standard).

## 8. Hypothesis bookkeeping (no files created or edited)

- **HYP-9142.** Suggested status line:
  - PROVED WITH A CONSTANT at shift 2: `N_m(L) <= 2.1285 A_(m+2)(L)` for all `m, L` (Theorem R_2);
  - at shift 3 with constant `1.1632` (Theorem R_3);
  - shift 1 PROVED for `m <= 24` (THM-4513 R_small);
  - shift 1 for all `m` OPEN.
- **HYP-9140, private price.** Given pairpeak Theorems B-C with `K = 2`: `π_L <= (517.3 L + 344.9 + 4L 2^(−L)) ρ^peak_L`.
- **Promotion candidates after audit:** Theorems R_2 and R_3, with Lemma A / EM\* / EM\*\*, Lemmas LR, LR-, TOP, TOP-LOW, Lemma B' and Proposition D' as reusable tools.
