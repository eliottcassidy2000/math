# Positive drift and provability in the strategy cube: an entropy law for sign flips, no provable strategy for q >= 23, and the Haar question for 5n+1 reduced to concentration

**Status.**
- **PROVED** (hand proofs in §§2–4; the runner checks their numerical content):
  - **Theorem 1 (stationary gain identity).** Take any odd `q`, level `k`, flip set `R`, and any stationary law `pi` of the uniform-lift chain on `G_sigma`. Then
    `pi(odd) = 1/2 - sum_{s in R} pi(s) g(s)`.
    - `g(s)` is the number of odd steps that the flip at `s` removes from the next `k-1` forced steps.
    - The same holds for an adversary appending Bernoulli(`theta`) parities, with `theta` in place of `1/2`.
    - It sharpens Proposition 4 of the cube-distance lane: every closed class of a class-(i) strategy has `pi(R) > (1/2 - log_q 2)/(k-1)`.
  - **Theorem 2 (entropy / merge law).** For every strategy of every `q n ± 1` cube,
    `1 - h(pi(odd)) <= sum_{merge pairs} (pi(s)+pi(s*)) h(pi(s)/(pi(s)+pi(s*))) <= pi(R ∪ R*) <= pi(odd)`.
    - `h` is the binary entropy.
    - `s* = s - 2q^(-1) mod 2^k` is the unflipped *partner* that shares the target pair of the flip `s`.
    - The uniform-lift chain produces exactly 1 bit per step. Its parity sequence can carry at most `h(pi(odd))` bits. The deficit can only be destroyed at merges, where a flip and its partner send the orbit to the same place.
  - **Corollary 2.1 (universal floor; no provable strategy for q >= 23).**
    - Every sign strategy of every `q n ± 1` cube, at every level, has `rho_max >= pi(odd) >= p0 = 0.2270922`, where `p0` is the root of `p + h(p) = 1`.
    - Hence **class (i) is empty at every level for every odd q >= 23**.
  - **Corollary 2.2 (constant flip mass; answers "is pi(R) >= const on closed classes?" — yes).** Every closed class of every class-(i) strategy of the 5n±1 cube has
    - `pi(R ∪ R*) > 1 - h(log_5 2) = 0.01391` and
    - `pi(R) > 0.00126`.
    
    These constants replace the `0.0693/k` of Proposition 4. For `q = 7, 9` the constants are `0.0605` and `0.1006`.
  - **Corollary 2.3 (the Haar question reduces to concentration).**
    - The Haar distance of any class-(i) strategy is at least `0.01391/λ`, where `λ` is the maximum invariant density `dpi/dU` on `R ∪ R*`.
    - It is also at least `0.01391^2/||dpi/dU||_2^2`.
    - So **5n+1 is in the Haar closure of class (i) only if the invariant densities of the approximating strategies blow up on their flips**.
    - A uniform bound on those densities would give a positive distance.
  - **Lemma 3 (carry lemma for 5n+1).** A flip replaces `T(r)` by `T(r) - 1`. That `-1` travels unchanged through the `10` pairs of the parity word, turning each into `01`.
    - At the first `11` pair it removes two odd steps (`11 -> 00`) and couples the orbits as `z` against `25z+8`.
    - At a `0x` pair it adds one odd step.
    - This explains where efficient flips sit (§6).
- **FINITE-EXACT.**
  - **5n+1.** Class (i) is empty for `k <= 6`; `delta_7(5) = 29`, re-deriving the cube-distance lane with independent code. **`44 <= delta_8(5) <= 50`**, which is new: lower bound from 10 implicit-hitting-set rounds, upper bound from a certified set.
  - **Min-max cycle density `rho*(q,k)`** (the least `rho_max` over all level-`k` strategies), each value SAT-certified. Class (i) at level `k` is nonempty iff `rho*(q,k) < log_q 2`.

    | `q` | `rho*(q,k)` | levels | class (i) |
    |---|---|---|---|
    | 3 | `1/2` | every `k` | nonempty |
    | 5 | `1/2` for `k <= 6`, `3/7` for `k = 7, 8` | `k <= 8` | nonempty from `k = 7` |
    | 7 | same values as `q = 5` | `k <= 9` | **empty** |
    | 9 | same values as `q = 5` | `k <= 8` | **empty** |
    | 11 | same values as `q = 5` | `k <= 8` | **empty** |
    | 13 | `1/2` | `k <= 8` | **empty** |

- **VERIFIED.** Explicit class-(i) flip sets for 5n+1 at every `k <= 19`, each with an integer potential checked on every edge. They come from lifting and greedily pruning the `k = 7` optimum.
  - Haar fraction: `0.391` (`k = 8`) down to `0.180` (`k = 19`).
  - `k × Haar` rises slowly from `3.09` (`k = 10`) to `3.41` (`k = 19`).
  - In the final sets every flip is critical.
- **EMPIRICAL.**
  - Over `10 <= k <= 19` the upper bounds fit `≍ 1/k` with a rising constant, `C k^(-0.84)`, and `a + b/k` with `a ≈ 0.035 > 0` about equally well.
  - The certified sets show no concentration: `max dpi/dU` on `R ∪ R*` is at most 2.4, and `||dpi/dU||_2^2` at most 1.46 (`k <= 14`).
- **OPEN.**
  - **Is 5n+1 in the Haar closure of class (i)?** There is neither an `o(1)` construction with proof (Q1) nor a positive Haar floor (Q2). By Corollary 2.3 the question is whether provable strategies can concentrate their invariant measure on their flips.
  - Is class (i) ever nonempty for `q = 7, 9, ..., 21`?
- **REFUTED.** Nothing. A wrong guess of mine, that 7n±1 is provable at level 2, was caught by the engine tests before it was recorded (§8).
- No HYP or THM file was created. Collatz itself is untouched (§9).

Session `collatz-procgen-20260922`, drift lane, 2026-09-26.
- Scripts: `04-computation/experiments/procgen_drift_20260926_{lib,ihs,local,minmax,run}.py` and the C engine `procgen_drift_20260926_engine.c`.
- Output: [procgen_drift_20260926.out](procgen_drift_20260926.out).
- Parents: [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) (the strategy cube; Theorems A, C) and the cube-distance note [procgen_cubedist_20260925_distance_to_provability.md](procgen_cubedist_20260925_distance_to_provability.md) (its Theorem 1 for 3n+1, the necklace Theorem 2, Proposition 4, the DRIFT data).

## 0. The answer in brief

**Setting** (THM-4474).
- A level-`k` sign strategy `sigma` sets `T(n) = n/2` for even `n` and `(q n + sigma(n mod 2^k))/2` for odd `n`, where `q` is odd.
- `R` is the set of residues with `sigma = -1`, the *flips* relative to `q n + 1`.
- **Class (i)** (bounded-lookahead provability) means that every cycle of the parity graph `G_sigma` has odd density `< c_q = log_q 2`.
- `delta_k(q)` is the least `|R|` of a class-(i) strategy, and `delta_k(q)/2^(k-1)` its **Haar distance**.
- For `q = 3` the cube-distance lane proved that this distance tends to 0 at the rate `2^(-(1-h)k)` (up to polynomial factors). The question here is `q = 5`, whose drift `(1/2)log(5/4)` is positive.

**One inequality decides the picture.** Run the uniform-lift chain on `G_sigma`: an adversary that appends a fair random parity bit at every step.
- Its node sequence carries exactly 1 bit per step.
- Class (i) forces the long-run frequency of odd steps below `c_q`. For `q >= 5`, `c_q < 1/2`, so the parity sequence carries at most `h(c_q) < 1` bits per step.
- The missing `1 - h(c_q)` bits per step must be destroyed. The only place where information is destroyed is a **merge**: a flipped residue `s` and its partner `s* = s - 2/q` send the orbit to the same target pair.
- Hence **flips and partners carry a constant stationary mass** (Theorem 2): `0.0139` for 5n+1, independent of `k`.

**Consequences.**
- **`q >= 23`.** Even all the odd mass cannot pay the deficit: `pi(odd) + h(pi(odd)) >= 1` forces `pi(odd) >= 0.227 > log_q 2`. No sign strategy of `q n ± 1` is provable at any level (Corollary 2.1).
- **`q = 3`.** The drift is negative, `c_3 > 1/2`, and the inequality is void. This is why flips can be exponentially sparse there.
- **5n+1.** The Haar distance can tend to 0 only if the flipped chain's invariant density **concentrates** on `R ∪ R*`, to at least `0.0139/Haar` somewhere (Corollary 2.3).
  - The best explicit sets, certified up to `k = 19`, show no concentration (density `<= 2.4` where it was computed, `k <= 14`), and their Haar fraction decreases like `3.1/k–3.4/k`.
  - Whether it tends to 0 or levels off is **OPEN**. It is the concentration question.

| | `q = 3` | `q = 5` | `q = 7` | `q = 9` | `q >= 23` |
|---|---|---|---|---|---|
| drift `(1/2) log(q/4)` | `-0.144` | `+0.112` | `+0.280` | `+0.405` | `>= +0.874` |
| `c_q = log_q 2` | `0.631` | `0.431` | `0.356` | `0.315` | `<= 0.221` |
| stationary constraint (Thm 2) | none | `pi(R ∪ R*) > 0.0139` | `> 0.0605` | `> 0.1006` | contradiction |
| class (i) | nonempty for `k >= 2` | nonempty for `k >= 7` | empty for `k <= 9` | empty for `k <= 8` | **empty at every level** |
| `rho*(q,k)` | `1/2` (every `k`) | `3/7` (`k = 7, 8`) | `3/7` (`k = 7..9`) | `3/7` (`k = 7, 8`) | `>= 0.227` |
| Haar distance | `2^(-(1-h)k + O(log k))` (cube-distance lane) | necklace `≈ 1.2/k–1.6/k` `<=` distance `<=` `3.1/k–3.4/k` (`k <= 19`); limit OPEN | undefined (`k <= 9`) | undefined (`k <= 8`) | undefined |

## 1. Setting and notation

Fix an odd `q >= 3` and `k >= 2`; put `N = 2^k` and `H = 2^(k-1)`.

**Parity words (Terras).**
- `Phi_m(x)` is the word of the first `m` parities of the all-plus orbit of `x` (`T_0 n = (qn+1)/2` or `n/2`).
- `Phi_m` is a bijection `Z/2^m -> {0,1}^m`, and `Phi_(m-1)(T_0 x)` is `Phi_m(x)` shifted left by one letter.
- Node `s` of `Z/N` carries the word `w(s) = Phi_k(s)`; it is odd iff `w(s)_0 = 1`.

**Edges of `G_sigma`.**
- Even `s`, or odd `s` outside `R`: `s -> w(s)_1…w(s)_(k-1) b` for `b in {0,1}` (the de Bruijn edges).
- Flipped `s` in `R`: `s -> u(s) b`, where `u(s) = Phi_(k-1)(y_-(s))` and `y_±(s) = (qs ± 1)/2`.
- The two lifts of any `t' mod H` have the words `Phi_(k-1)(t')0` and `Phi_(k-1)(t')1`.

**Chains.**
- `Q^theta_sigma` moves to the successor whose last letter is 1 with probability `theta`: the adversary appends a Bernoulli(`theta`) parity.
- `P_sigma = Q^(1/2)_sigma` is the **uniform-lift chain**.
- Stationary laws live on closed classes.

**Sandwich** (THM-4474 Theorem C). The proof uses only that a stationary edge flow is a circulation, so it applies to every `Q^theta_sigma`. On every closed class `pi(odd) <= rho_max`; for class (i), `pi(odd) < c_q`.

**Gain.** For odd `s`, `g(s) = |w(s)_(1..k-1)| - |u(s)|`, where `|.|` counts ones. It is the number of odd steps among the next `k-1` forced steps that the flip removes; it is negative if the flip adds some. Always `|g| <= k-1`.

**Partner.** `s* = s - 2q^(-1) mod N`. It is odd, and `y_+(s*) ≡ y_-(s) mod H`.

## 2. Theorem 1: the stationary gain identity

**Theorem 1.** For every strategy `sigma` (flip set `R`), every `theta` in `(0,1)` and every stationary law `pi` of `Q^theta_sigma`,

`pi(odd) = theta - sum_(s in R) pi(s) g(s)`.

*Proof.*
1. **The base chain mixes in `k` steps.** Let `Q_0 = Q^theta_∅`. In word coordinates it shifts left and appends a Bernoulli(`theta`) letter. After `k` steps all letters are fresh, so `Q_0^k(s,·) = B_theta` (the Bernoulli product law) for every `s`.
2. **The perturbation.** Write `Q_sigma = Q_0 + Delta`.
   - `Delta(s,·) = 0` for `s` outside `R`.
   - For `s` in `R`, `Delta(s,·) = nu_s^- - nu_s^+`. Here `nu_s^+` is the law of the word `w(s)_(1..k-1) b` and `nu_s^-` that of `u(s) b`, with `b ~ Bernoulli(theta)`.
3. **Expansion.** Iterating `pi = pi Q_0 + pi Delta` gives `pi = pi Q_0^k + sum_(j<k) (pi Delta) Q_0^j`, and `pi Q_0^k = B_theta`.
4. **Evaluation on the odd nodes `O`.** `B_theta(O) = theta`. Let `nu` be the law of a fixed `(k-1)`-word `v` followed by a Bernoulli letter.
   - For `j <= k-2`, `(nu Q_0^j)(O) = [v_j = 1]`: after `j` shifts the first letter is `v_j`.
   - For `j = k-1`, it equals `theta`.
   
   Hence `sum_(j<k) (nu_s^- - nu_s^+) Q_0^j (O) = |u(s)| - |w(s)_(1..k-1)| = -g(s)`. ∎

**Corollary 1.1 (Proposition 4 sharpened).** If `sigma` is in class (i), then for every `theta` and every stationary law of `Q^theta_sigma`,
`sum_(s in R) pi(s) g(s) > theta - c_q`, hence `pi(R) > (theta - c_q)/(k-1)`.
- For 5n+1 with `theta = 1/2` this gives `0.0693/(k-1)`.
- *Proof:* the sandwich gives `pi(odd) < c_q`, and `g <= k-1`. ∎

**Remark (gain tail).**
- For a uniform odd `s`, `y_+(s)` and `y_-(s)` are each uniform mod `H`. So `|Phi_(k-1)(y_±)|` is `Binomial(k-1, 1/2)`, and Hoeffding gives `#{s : g(s) >= x} <= 2^(k-1) · 2 exp(-x^2/(2(k-1)))`.
- Numerically `Var g ≈ 0.52–0.53 k`, and `max g = 4, 6, 9, 12` at `k = 8, 12, 16, 20` (runner §B3).
- **Consequence.** If `pi <= λ U` on `R`, then Theorem 1 gives `Haar(R) ≳ (1/2 - c_q)/(1.5 λ sqrt(2k ln k))`, up to lower-order terms. Theorem 2 below replaces this by a constant.

## 3. Theorem 2: the entropy / merge law

**Theorem 2.** Let `sigma` be any strategy of any odd `q` at any level `k`, and let `pi` be a stationary law of `P_sigma` with `p = pi(odd)`. Then

`1 - h(p) <= sum_(s in R, s* not in R) (pi(s)+pi(s*)) h( pi(s)/(pi(s)+pi(s*)) ) <= pi(R ∪ R*) <= p`,

where `R* = {s* : s in R, s* not in R}`, `h` is the binary entropy, and `0·h(0/0) = 0`.

*Proof.*
- **(a) Predecessors.** Fix `t'` mod `H`.
  - Its odd predecessors are `s_+(t') = q^(-1)(2t'-1)` if that residue is unflipped, and `s_-(t') = q^(-1)(2t'+1)` if that residue is flipped. They solve `qs ± 1 ≡ 2t' mod N`, and `s_+(t') = s_-(t')*`.
  - So the pair `{t', t'+H}` has two odd predecessors iff `s_-(t')` is in `R` and its partner is not. Call it a **merge pair**.
  - The even predecessor `2t' mod N` is unique.
- **(b) Entropy rates.** Let `(Y_m)` be the stationary chain and `X_m = Y_m mod 2`.
  - Each node has two distinct successors, each taken with probability `1/2`, so `H(Y_(m+1)|Y_m) = 1` and `H(Y_0..Y_n) = H(Y_0) + n`.
  - `H(X_0..X_n) <= (n+1) h(p)`.
- **(c) Chain rule backwards.** `H(Y_0^n | X_0^n) <= H(Y_n) + sum_(m<n) H(Y_m | Y_(m+1), X_m) <= k + n H(Y_0 | Y_1, X_0)`, by conditioning and stationarity.
  - With `H(Y_0^n) = H(X_0^n) + H(Y_0^n | X_0^n)`, divide by `n` and let `n -> ∞`: `1 <= h(p) + H(Y_0|Y_1,X_0)`.
- **(d) The conditional entropy.**
  - Given `Y_1 = t` and `X_0 = 0`, `Y_0` is the even predecessor.
  - Given `X_0 = 1`, `Y_0` is ambiguous only at a merge pair with flip `s`.
  - There each predecessor sends half its mass to each lift. So for each lift, `P(Y_1 = t, X_0 = 1) = (pi(s)+pi(s*))/2`, and the conditional law of `Y_0` is `(pi(s), pi(s*))/(pi(s)+pi(s*))`.
  - Summing over both lifts gives the middle expression.
- **(e) The last two inequalities.** The second follows from `h <= 1`. For the third: the nodes `s`, `s*` over all merge pairs are distinct odd nodes, since `R ∩ R* = ∅` and `s -> s*` is injective. ∎

**Corollary 2.1 (universal floor).**
- `p + h(p) >= 1` for every stationary law of every strategy (every odd `q`, every level).
- `p + h(p)` is increasing on `(0, 1/2]`, so `p >= p0 = 0.227092195…`, its root there.
- The sandwich then gives `rho_max >= p0` for every strategy.
- **Class (i) is empty at every level for every odd `q >= 23`**, because `log_q 2 <= log_23 2 = 0.2211 < p0`.
- The bound does not reach `q = 21`: `log_21 2 = 0.2277` is just above `p0`.
- A separate run of `..._minmax.py 23 2 3 4 5 6 7` (not part of the runner) gives `rho*(23,k) = 1/2` for `k <= 7`, consistent with this.

**Corollary 2.2 (constant flip mass).** If `sigma` is in class (i), then on every closed class
- `pi(R ∪ R*) > 1 - h(c_q)`, and
- `pi(R) > A_q`, where `A log_2(e/A) = 1 - h(c_q)`.

| `q` | `c_q` | `1 - h(c_q)` | `A_q` |
|---|---|---|---|
| 5 | 0.43068 | 0.01391 | 0.00126 |
| 7 | 0.35621 | 0.06051 | 0.00704 |
| 9 | 0.31546 | 0.10062 | 0.01307 |
| 11 | 0.28906 | 0.13249 | 0.01838 |
| 21 | 0.22767 | 0.22607 | 0.03631 |

*Proof.*
- Class (i) gives `p < c_q < 1/2`, hence `1 - h(p) > 1 - h(c_q)`.
- For the second bound, write `a = pi(s)` and `b = pi(s*)`. Then `(a+b)h(a/(a+b)) = a log((a+b)/a) + b log(1 + a/b) <= a log((a+b)/a) + a log_2 e`.
- Summing with Jensen (weights `a_s/A`) and using `A + B <= 1` gives at most `A log_2(e/A)`, which is increasing on `(0,1)`. ∎

**Corollary 2.3 (the Haar reduction).** Let `f = dpi/dU` on a closed class of a class-(i) strategy. Then

`Haar(R) = |R|/2^(k-1) >= U(R ∪ R*) > (1 - h(c_q)) / max_(R ∪ R*) f`, and `Haar(R) > (1 - h(c_q))^2 / ||f||_2^2` (Cauchy–Schwarz).

- So any sequence of class-(i) strategies of 5n+1 whose Haar distance tends to 0 has `max_(R∪R*) f -> ∞` and `||f||_2 -> ∞`.
- **Conversely, a uniform bound on the invariant densities of provable strategies would keep 5n+1 out of the Haar closure**, at distance `> 0.0139/sup f`.

**Remarks.**
- **Why `q = 3` escapes.**
  - Theorem 2 holds for `q = 3`, but `c_3 > 1/2` allows `pi(odd) = 1/2`, where `1 - h = 0`.
  - The flipped undecided strategies of the cube-distance lane show this concretely (runner §H1): `pi(odd) ≈ 0.42–0.43` and `pi(R) = 0.041, 0.034, 0.029` at `k = 8, 10, 12`. There is no floor.
  - Theorem 1 with `theta > c_3` does constrain `q = 3`. Its reference law `B_theta` puts its mass on words with about a `theta`-fraction of ones, a set of Haar measure about `2^(-(1-h(theta))k)`, so it imposes no Haar bound.
  - The sign of the drift, i.e. `c_q` against `1/2`, is exactly what makes the Haar-reference chain constrained.
- **Slack.**
  - For the certified 5n+1 sets the merge entropy is `0.19–0.25`, against a requirement of `0.0139`. Their stationary odd frequency is `≈ 0.31`, far below `c_5 = 0.43`. The binding constraints are worst-case cycles, not averages (§5).
  - The universal floor `p0 = 0.227` lies below every stationary odd frequency the runner computes. The greedy max-halving strategy gives exactly `1/3` (to `1e-6`) for `q = 5, 7, 9`, and the certified sets give `0.308–0.313`.

## 4. Lemma 3: the carry lemma for 5n+1

**Lemma 3.** Let `T` be 5n+1 and `y ∈ Z_2`. The orbits of `y` and `y - 1` compare as follows.

| `y mod 4` | parities of `(y, Ty)` | parities of `(y-1, T(y-1))` | relation after two steps |
|---|---|---|---|
| 1 | `11` | `00` | `T^2 y = 25 T^2(y-1) + 8`; the next 3 parities agree |
| 3 | `10` | `01` | `T^2(y-1) = T^2 y - 1` |
| 2 | `01` | `11` | `T^2(y-1) = 5 T^2 y - 7` |
| 0 | `00` | `10` | `T^2(y-1) = 5 T^2 y - 1` |

*Proof.* Direct computation.
- `y ≡ 1 (mod 4)`: `Ty = (5y+1)/2` is odd and `y-1 ≡ 0 (mod 4)`. Also `T^2 y = (25y+7)/4 = 25(y-1)/4 + 8`, and `25z + 8 ≡ z (mod 8)`, which gives the three agreeing parities (Terras).
- `y ≡ 3 (mod 4)`: `T^2(y-1) = (5y-3)/4 = T^2 y - 1`.
- `y ≡ 2 (mod 4)`: `T(y-1) = 5Ty - 2` is odd, and `T^2(y-1) = 5T^2 y - 7`.
- `y ≡ 0 (mod 4)`: `T(y-1) = 5Ty - 2` is even, and `T^2(y-1) = 5T^2 y - 1`.

The identities are polynomial; the runner checks them on 400,001 integers. ∎

**Reading.**
- A flip at `r` replaces `y = T(r)` by `y - 1`. By the `y ≡ 3 (mod 4)` case, the `-1` passes unchanged through every pair `10` of the parity word of `y`, turning it into `01`.
- At the first pair that is not `10`:
  - on a `11` pair the flipped orbit reads `00`, a **gain of two odd steps**, and continues as `z` against `25z + 8`;
  - on a `0x` pair it gains an odd step, so the flip is locally harmful.
- This is the local mechanism of the certified sets (runner §J, `k = 14`):
  - flips sit on windows `1(10)^j 11…`;
  - `1110`-windows are flipped 84% of the time and `1111`-windows 25%;
  - `101…`-windows (`r ≡ 7 (mod 8)`), where the flip is locally harmful, are flipped only 10% of the time, presumably for cycle breaking;
  - `100…`-windows (`r ≡ 3 (mod 8)`) are never flipped.
- The sign of the exact gain agrees with this two-letter prediction for 80% of odd residues at `k = 12` (runner §D2). Later disagreements between the coupled orbits account for the rest.

## 5. 5n+1: exact values, certified upper bounds, fits, diagnostics

**Exact values** (`..._ihs.py`; runner §E).
- *Method.* An implicit hitting set:
  - RC2 (core-guided MaxSAT) solves the hitting-set problem exactly.
  - Its no-goods are expanding cycles recorded with the signs on their odd nodes. Every no-good is re-derived as a genuine closed walk with `5^a > 2^p`.
  - The seeds are all expanding de Bruijn cycles of length `<= k+2`.
- *`k <= 6`:* class (i) is empty. The no-goods are UNSAT (Glucose4).
- *`k = 7`:* **`delta_7 = 29`** (Haar `0.4531`).
  - No hitting set of size 28 exists over 1215 no-goods.
  - The optimum has `rho_max = 3/7`, and an independent pure-Python Karp agrees.
- *`k = 8`:* **`44 <= delta_8 <= 50`** (Haar `0.344–0.391`).
  - The lower bound is the exact RC2 optimum after 10 IHS rounds (3116 no-goods, 1102 s).
  - The upper bound is the best of 21 seeded prunings of the lifted `k = 7` optimum.
  - Longer runs (not reproduced by the runner) did not move either bound: an 11th RC2 round ran for over 40 minutes without finishing, and local search found no class-(i) set below 50.

**Certified upper bounds** (runner §G).
- *Construction.* Lift the `k = 7` optimum level by level (the map is unchanged). At each level, greedily restore `sigma = +` wherever class (i) survives.
- *Exactness of each test.* It uses the potential at the exact threshold `F`, the best lower approximation of `log_5 2` with denominator `<= 2^k`. A Farey-neighbour certificate checks `F` exactly, so every simple cycle above `F` is expanding.
- *Certificate.* Every final set carries an integer potential, checked on all `2^(k+1)` edges.

| `k` | `\|R\|` | Haar | `k` × Haar | necklace lower bound `N_k(5)` | its Haar | threshold `F` |
|---|---|---|---|---|---|---|
| 7 | 29 (exact, §E2) | 0.4531 | 3.172 | 10 | 0.1563 | `31/72` |
| 8 | 50 | 0.3906 | 3.125 | 23 | 0.1797 | `59/137` |
| 9 | 94 | 0.3672 | 3.305 | 44 | 0.1719 | `205/476` |
| 10 | 158 | 0.3086 | 3.086 | 67 | 0.1309 | `351/815` |
| 11 | 295 | 0.2881 | 3.169 | 136 | 0.1328 | `643/1493` |
| 12 | 529 | 0.2583 | 3.100 | 216 | 0.1055 | `643/1493` |
| 13 | 1000 | 0.2441 | 3.174 | 448 | 0.1094 | `643/1493` |
| 14 | 1890 | 0.2307 | 3.230 | 714 | 0.0872 | `4647/10790` |
| 15 | 3546 | 0.2164 | 3.246 | 1525 | 0.0931 | `8651/20087` |
| 16 | 6754 | 0.2061 | 3.298 | 3178 | 0.0970 | `21306/49471` |
| 17 | 12837 | 0.1959 | 3.330 | 5286 | 0.0807 | `21306/49471` |
| 18 | 24762 | 0.1889 | 3.401 | 11091 | 0.0846 | `97879/227268` |
| 19 | 47081 | 0.1796 | 3.412 | 18660 | 0.0712 | `97879/227268` |

**Fits** (EMPIRICAL; these are upper bounds from a heuristic, not `delta_k`), over `10 <= k <= 19`:
- `log_2 Haar = 1.838 + 0.0328 k - 1.161 log_2 k` (rms `0.0110`);
- `log_2 Haar = 1.091 - 0.842 log_2 k`, a power law `k^(-0.842)` (rms `0.0139`);
- `Haar = 0.0355 + 2.729/k` (rms `0.0022`);
- `Haar = 3.203/k` (rms `0.0076`).
- The data cannot separate three shapes:
  - `≍ 1/k` with a slowly rising constant: `k` × Haar goes from `3.09` to `3.41`;
  - a power `k^(-0.84)`;
  - `a + b/k` with a small positive `a`.
- The three-parameter fit has a positive linear coefficient. So there is no exponential decay here: the local decay is polynomial and slower than `1/k`.
- By Corollary 2.3, a positive floor is what bounded invariant densities would force, and decay to 0 is what concentration would allow.
- The necklace bound (cube-distance Theorem 2, which does not depend on the multiplier) gives `≈ 1.2/k–1.6/k` here and tends to `2/k`. It remains the best unconditional lower bound.

**Concentration diagnostics** (runner §C4: the closed class of each certified set, `f = dpi/dU`).

| `k` | Haar | `pi(odd)` | `1 - h(pi(odd))` | merge entropy | `pi(R)` | `pi(R*)` | `U(R ∪ R*)` | `max f` on `R ∪ R*` | `\|\|f\|\|_2^2` |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 0.4531 | 0.3095 | 0.1074 | 0.2515 | 0.1577 | 0.1139 | 0.4141 | 1.310 | 1.265 |
| 8 | 0.3906 | 0.3080 | 0.1092 | 0.2352 | 0.1472 | 0.1043 | 0.3672 | 1.534 | 1.320 |
| 9 | 0.3672 | 0.3082 | 0.1090 | 0.2320 | 0.1397 | 0.1047 | 0.3516 | 1.576 | 1.344 |
| 10 | 0.3086 | 0.3092 | 0.1077 | 0.2131 | 0.1290 | 0.0953 | 0.2998 | 1.857 | 1.389 |
| 11 | 0.2881 | 0.3096 | 0.1072 | 0.2075 | 0.1214 | 0.0947 | 0.2817 | 1.825 | 1.404 |
| 12 | 0.2583 | 0.3117 | 0.1049 | 0.1967 | 0.1121 | 0.0909 | 0.2554 | 1.968 | 1.412 |
| 13 | 0.2441 | 0.3120 | 0.1045 | 0.1917 | 0.1096 | 0.0878 | 0.2418 | 2.204 | 1.438 |
| 14 | 0.2307 | 0.3130 | 0.1034 | 0.1864 | 0.1047 | 0.0862 | 0.2291 | 2.393 | 1.456 |

- The flips carry a stationary mass of `0.7–0.9 ×` their Haar mass (`pi(R)` against `U(R) = Haar/2`): there is no concentration.
- `max f` on `R ∪ R*` creeps up slowly (1.3 to 2.4).
- An exploratory 4-minute random-toggle hill-climb found little more concentration even when maximizing it. It maximized `||f||_2^2` over class-(i) strategies near these sets, with any number of flips. It reached only `1.39` at `k = 8` and `1.59` at `k = 10` (`max f <= 2`). This run is not part of the runner.

## 6. Q1: what fails, and what the certified sets look like

**Simple explicit rules all fail** (runner §I, `k = 8, 10, 12, 14`). Each leaves a verified expanding cycle, mostly of length 2 to 6:
- all 111-windows (`r ≡ 5 (mod 8)`);
- the max-halving strategy (`r ≡ 1 (mod 4)`);
- the undecided set `Bad_k(5)`, i.e. the 3n+1 construction; also `Bad_k ∩ 111`;
- one flip per expanding necklace, at its lowest or its highest ballot rotation (Haar `N_k/2^(k-1)`);
- all positive-gain residues `g >= 1`.

**Why the 3n+1 construction does not transfer.**
- (i) `Bad_k(5)` has positive density, because 5n+1 drifts upward.
- (ii) By Lemma 3, a flip is a descent only on `1(10)^j 11` windows. On `10…` windows it adds an odd step.
- (iii) Deeper, by Theorem 2: every provable strategy must merge at a constant stationary rate. Blocks that each descend on their own therefore cannot suffice. The flips must be revisited with positive frequency, yet occupy a vanishing Haar fraction.

**Structure of the certified sets** (runner §J).
- Every flip is critical: removing any one leaves a verified expanding cycle, of median length `11, 15, 17` at `k = 10, 12, 14`.
- `rho_max = 3/7`, attained by cycles of length 7 or 14.
- Flipped fraction per prefix class at `k = 14`: `111`: 0.55, `110`: 0.27, `101`: 0.10, `100`: 0.
- Flipped residues are visited about twice as often as unflipped residues of the same class: `f ≈ 1.08` against `0.54` for the 111 class. This is mild concentration that does not grow with `k`.

**An `o(1)` construction is OPEN.**
- Corollary 2.3 says what such a construction must do: its orbits must return to the flips at a positive frequency while the flips' Haar measure vanishes. That is, it needs a *funnel*.
- The obvious funnels are too expensive. Forcing an orbit into a flip `j` steps after the previous flip requires flipping all `2^(j+1)` completions, because the adversary's free bits can avoid any single one.

## 7. Q4: other multipliers

**The min-max density** (runner §F; `..._minmax.py`).
- *Method.* A descent on the value.
  - The current strategy's exact `rho_max` comes from Dinkelbach iteration, with a verified cycle and a verified potential.
  - A SAT search then looks for a strategy all of whose cycles are strictly below that value.
  - The final value is certified by an UNSAT set of re-derived cycle no-goods.
- *Properties.* `rho*` is non-increasing in `k` (lifting keeps the map) and is at least `p0 = 0.227` (Corollary 2.1).
- *Values* (Status table).
  - For `q = 5, 7, 9, 11` it drops from `1/2` to `3/7` at `k = 7`.
  - `3/7 < c_5`, but `3/7 > c_7, c_9, c_11`. So 7n+1, 9n+1 and 11n+1 **have no provable sign strategy at the levels computed**.
  - `q = 13` stays at `1/2` through `k = 8`.
  - `q = 3` is pinned at `1/2` for every `k`, a two-line proof. With `sigma(1) = +` the cycle `1 -> 2 -> 1` has density `1/2`; with `sigma(1) = -` the point 1 is fixed with density 1. The level-2 strategy that flips `r ≡ 3 (mod 4)` makes every odd step even-followed, so `rho_max <= 1/2`.
- That `q = 5, 7, 9, 11` share the same values is observed, not explained.

**Does the rate depend on the drift?** Yes. There are three regimes, set by `c_q = log_q 2`, equivalently by the drift `(1/2) log(q/4)`.

1. **Negative drift (`q = 3`).** There is no stationary constraint (`c_3 > 1/2`). The Haar distance tends to 0 exponentially (cube-distance Theorem 1).
2. **Positive drift with `c_q > p0` (`5 <= q <= 21`).** The flips carry a constant stationary mass (Theorem 2).
   - For 5n+1 the Haar distance decreases polynomially on the data. Its limit is OPEN and equivalent to a concentration question.
   - For `q = 7, 9, 11` no provable strategy exists at the computed levels.
3. **Large drift with `c_q <= p0` (`q >= 23`, drift `>= 0.874`).** There is no provable strategy at any level (Corollary 2.1).

For `q = 7, 9` the necklace lower bounds (runner §H2) are `N_k(q)/2^(k-1) ≈ 1.65/k–1.9/k` for `k = 8..16`, tending to `2/k`. They would apply if class (i) became nonempty.

## 8. Controls and consistency checks

- **`q = 3`** (runner §H1). Flipping `Bad_k` is class (i) at `k = 8, 10, 12, 14`, with Haar `0.148 → 0.090`.
  - Its stationary odd frequency is `0.42–0.43`: below `1/2`, and far below `c_3 = 0.631`.
  - `pi(R)` falls with `k` (`0.041, 0.034, 0.029`). For `q = 3` there is no constant floor, unlike `q = 5`.
- **Engine** (runner §A).
  - It agrees with an independent pure-Python Karp on 300 random strategies (`q = 3, 5, 7, 9`, `k = 3..5`). This covers Karp, the potential certificate and Dinkelbach.
  - Incremental toggles agree with Karp on 200 random toggles from the `k = 7` optimum, and the maintained potential is valid after each.
- **Theorem 1.** Checked on 360 (strategy, `theta`, closed class) triples: error `3e-16`, float64.
- **Theorem 2.** Checked on 150 closed classes of random strategies (`q = 3..11`, `k = 4..9`) and on every certified 5n+1 set.
- **A caught mistake.** At first I thought the max-halving strategy of 7n±1 is provable at level 2, on the grounds that every odd step is followed by two halvings. In fact `(qn±1)/2` already contains one halving, so only blocks `10` are guaranteed, with multiplier `7/4 > 1`. The engine test returned `rho_max = 1/2`, and the claim was never recorded.
- **Implementation notes.**
  - The C engine keeps one global state. The Python wrapper refuses to use a stale instance; an early test mixed two instances, and this guard was added after that test was fixed.
  - An experimental Dijkstra pre-check for toggles was slower on least potentials, which have many tight edges, and was removed.

## 9. What this does and does not say

**It does:**
- It proves an exact information-theoretic law for sign flips (Theorem 2) that separates the three drift regimes.
- It settles an infinite family: no `q n ± 1` with odd `q >= 23` has a bounded-lookahead-provable sign strategy at any level.
- It answers the lane's question "is `pi(R) >= const`?": yes, `pi(R) > 0.00126`.
- It reduces "is 5n+1 in the Haar closure of class (i)?" to the concentration of invariant densities.

**It does not:**
- It does not decide that question: a decay like `1/k` and a positive floor both remain consistent with everything proved and computed.
- Nothing here bears on Collatz itself (3n+1 is in regime 1, where Theorem 2 is void), nor on actual 5n+1 or 7n+1 orbits.

**Relation to the other lanes of 2026-09-26.**
- **The same constant, in a different role.** [THM-4480](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md) concerns arbitrary, non-periodic, fixed-horizon edits. There the price is `2^(-(1-H(log_q 2))L + o(L))` for every odd `q`, including `0.013911` for 5n+1.
  - That is the constant `1 - h(log_5 2) = 0.01391` of Corollary 2.2.
  - There it is the *exponent of a vanishing price*: arbitrary edits can see height and catch each undecided orbit at its peak. Here it is a *constant stationary mass*: residue classes cannot see height.
  - Both are the entropy deficit of the critical tube, the parity words of density `log_q 2`.
- **Periodic deletions versus flips.** The deletion lane (`procgen_mykk_20260926`, in progress when this note was written) treats periodic edits that send the residues of `R` to a fixed value. It reports a price `≍ 1/k` for `q = 5` and the chain `delta_k >= FVS^odd >= FVS >= N_k`. Theorem 2 explains why sign flips can be qualitatively harder:
  - an edit to a fixed value merges *all* edited orbits into one point, so it destroys many bits per visit;
  - a flip merges *exactly two* orbits (`s` and `s*`), so it destroys at most one bit per visit.
- **The periodic 3n+1 cube.** The periodic 3n+1 cube is settled by [THM-4479](../../01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md), the cube-distance theorem, now promoted. There Theorem 2 is void.

## 10. Reproduction

```bash
DRIFT_KMAX=19 python3 04-computation/experiments/procgen_drift_20260926_run.py > 05-knowledge/results/procgen_drift_20260926.out
```

- **Requirements.** `numpy`, `scipy`, `mpmath`, `pysat` (RC2, Glucose4) and a C compiler.
- **Engine build.** The engine is compiled into `scratch/procgen_drift/engine_<source-hash>.so`, which is not committed.
- **Environment variables.**
  - `DRIFT_KMAX` (script default 18; the output above used 19) is the top level of §G.
  - `DRIFT_K8_ITERS` (default 10) is the number of IHS rounds at `k = 8`.
- **Cost.** One process; wall time `1493 s` on a shared 8-core machine whose load average from other agents was about 20.
  - Of this, the ten RC2 rounds at `k = 8` took `1102 s` and the `k = 19` pruning `262 s`.
  - Peak RSS of the runner: `284 MB` (self-reported); `/usr/bin/time` reports `298 MB`.
  - RC2 and Glucose4 are single-threaded.
- **Checks.** Every claim in the output is a `check(...)` that raises on failure. The output ends with `ALL CHECKS PASSED`.
- **SHA-256** (raw bytes):
  - `procgen_drift_20260926_engine.c` `d26eb8df59efb184fafcba497f18d8b60fce4cbe30ee390421fe2e3974d32a43`
  - `procgen_drift_20260926_lib.py` `7d7fd5c88892a731966c2d1555ccdaaf33ff8c283a8c00da0591ff759046d896`
  - `procgen_drift_20260926_ihs.py` `b781c768364530ea6657cfb77d71c3f068b3da501ead209872912d3dab445d23`
  - `procgen_drift_20260926_local.py` `42a6f934030c9ec8a8bb31601f648c0cb2fb740675f66186773e02200ea123eb`
  - `procgen_drift_20260926_minmax.py` `0f38db32e7e4c15ef92f99dbe0a7b20bff7d7b9e0f50e5502813598b73d5d85b`
  - `procgen_drift_20260926_run.py` `c9d6e124d22172eaff4aed5c2a4d90d1051af3cfde00c0f401363c5465edf6d4`
  - `procgen_drift_20260926.out` `b72eaa28c95dbb642ba25dfed9eac09cb5ce3139ccd1d67ef397397e9fd59f83`
  - The runner prints the six script hashes itself.
