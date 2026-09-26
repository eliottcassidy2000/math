# The Robin inequality holds up to a polynomial: N_m(L) <= (m+2)^17 A_(m+2)(L), so the private pairing price is poly(L) rho^peak_L; the second-order term of rho^peak_L and pi_L is -kappa_3 L^(1/3) + O(log L), elementary

Lane `robin`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_robin_20260926_{lib,run}.py`.
Output: [procgen_robin_20260926.out](procgen_robin_20260926.out) (one runner; every `ok:` line is a `check(...)` that aborts on failure; it ends with `ALL CHECKS PASSED`).

## Status

**Summary.**
- Conjecture R (HYP-9142) with a *constant* `K_0` is **not** proved.
- It **is** proved with shift `K = 2` and a polynomial factor, which is the fallback the brief named:
  `N_m(L) <= (m+2)^17 A_(m+2)(L)` for all `m >= 2`, `L >= 1`.
- This suffices for HYP-9140 for the private price (which asks for `poly(L)`) and for the sharp constant of `pi_L`. It gives `pi_L = O(L^18) rho^peak_L`, not the `O(L)` in HYP-9142's title; that factor would need a constant `K_0`.
- Two consequences follow:
  - HYP-9140 holds for the *private* pairing price, `pi_L <= poly(L) rho^peak_L`, given Theorems B–C of the pairpeak note.
  - The second-order term of `rho^peak_L` is `-kappa_3 L^(1/3) + O(log L)`, by an elementary proof. The same holds for `pi_L`, given pairpeak Theorem B.
- **The mechanism.** The free letter walk (steps `+(1-c)` w.p. `c`, `-c` w.p. `1-c`) has *exact* eigenfunctions `e^(beta s) sin(theta s + phi)` on the whole line (Lemma E). Two of them are used:
  - a shifted one is a supersolution for the reflected (Robin) walk at width `m + 2` (Theorem 1);
  - an unshifted one is a subsolution for the hard-wall walk at width `M` (Theorem 2).
  - At `M = m + 2` both have the same eigenvalue `kappa(pi/(m+2))`.
- **The only lossy step** is the endpoint. Counting measure ends paths near the bottom, while one-sided bounds see the top. A two-phase bridge lemma (cycle lemma + Hoeffding without replacement) repairs this at polynomial cost.

**PROVED (full hand proofs below).**
- **Lemma E (free eigenfunctions).** For `0 < theta < pi/c` put
  - `beta(theta) = ln[(1-c) sin(c theta)/(c sin((1-c) theta))] < 0`;
  - `kappa(theta) = c e^(beta(1-c)) cos((1-c)theta) + (1-c) e^(-beta c) cos(c theta)`.

  Then `c f(s+1-c) + (1-c) f(s-c) = kappa(theta) f(s)` for every real `s` and every `f(s) = e^(beta s) sin(theta s + phi)`. Moreover `ln kappa(theta) <= -0.1 theta^2` on `(0, pi/2]`, and `kappa(theta) = 1 - sigma^2 theta^2/2 + O(theta^4)` with explicit constants on `(0, 0.15]`.
- **Lemma Z (zone condition).** `F(theta) = kappa sin(2 theta) - 2(1-c) e^(-beta c) sin((2+c)theta) > 0` for `0 < theta <= pi/4`.
  - Analytic on `(0, 0.15]`.
  - **Computer-assisted** on `[0.15, pi/4]`: interval arithmetic (`mpmath.iv`), 129 boxes.
- **Theorem 1 (reflected barrier, sharp rate).**
  `N_m(L) <= e^(|beta(theta_m)|(m-c)) 2^(HL) kappa(theta_m)^L <= e^(0.063) 2^(HL) kappa(pi/(m+2))^L` for all `m >= 2`, `L >= 0`, where `theta_m = pi/(m+2)`.
  - Per step this is `1 - 1.149/(m+2)^2 + ...`.
  - It replaces Theorem A of the pairpeak note, whose per-step factor is `1 - 0.0056/(m+1)^2`.
- **Theorem 2 (hard wall, lower bound).** `A_M(L) >= 2^(HL) kappa(pi/M)^L / P_2(M)` for every real `M >= 4` and every `L >= 1`.
  - `P_2` is explicit (§3.4) and `ln P_2(M) <= 16.83 ln M + 1.31` for `M >= 3400`.
- **Corollary 3 (HYP-9142 up to a polynomial, `K = 2`).**
  `N_m(L) <= P(m) A_(m+2)(L)` with `P(m) = e^(|beta(theta_m)|(m-c)) P_2(m+2) <= (m+2)^17`, for all `m >= 2`, `L >= 1`.
- **Corollary 4 (HYP-9140 for the private price).** `pi_L <= (9(27L+18)(L+2)^17 + 4L 2^(-L)) rho^peak_L` for all `L >= 3`.
  - This is Theorem C of the pairpeak note with `K = 2`, `K_0 = (L+2)^17`.
  - It rests on that note's Theorems B–C (lane-proved in the same session, not yet independently audited).
- **Corollary 5 (sharp second-order term, `q = 3`, elementary).** Let `Lambda(L) = max_(a>=2) [L ln kappa(pi/a) - a ln 3]`. Then `Lambda(L) = -kappa_3 L^(1/3) + O(L^(-1/3))`, and
  - `|ln(2^((1-H)L) rho^peak_L) - Lambda(L)| = O(log L)`;
  - `ln 2 + ln(2^((1-H)L) rho^peak_L) <= ln(2^((1-H)L) pi_L) <= Lambda(L) + O(log L)`.
  - Hence `ln rho^peak_L` and `ln pi_L` both equal `-(1-H)L ln 2 - kappa_3 L^(1/3) + O(log L)`, with `kappa_3 = 2.107580`.
  - The `rho^peak` statement uses only Theorem 2 and a Dirichlet supersolution. The `pi_L` statement also uses Theorem B of the pairpeak note.
  - For `q = 3` this sharpens THM-4480 Theorem 4 (`(1+o(1))`, modulo Mogul'skii) to an elementary `O(log L)` error.
  - It also settles the sharp constant of `pi_L`, which the pairpeak note had conditional on Conjecture R and Mogul'skii.
- **Theorem 6 (reduction of the sharp Conjecture R).** `N_m(L) <= C A_(m+K)(L) - (C-1) F_L(0)` follows from a *one-walk* inequality for the hard-wall walk alone, RM(K, C):
  `V_k(z-c) - V_k(z+1-c) <= (2 - 2/C) F_k(z-c)` at the zone points `z` along the time diagonal.

**Pre-audit.** A fresh agent audited the note independently, writing its own code before reading the lane's scripts (`scratch/procgen_robin/audit/`, not part of the deliverable).
- It found no mathematical error, and marks every hand proof SOUND.
- It reproduced every finite-exact number and table.
- It re-derived all constants.
- It confirmed the float checks of Theorem 1 at `L <= 4096` with exact integers.
- It flagged three presentation gaps, all fixed here:
  - the `B`-free `theta^4` term in Lemma E(iii);
  - interval endpoints at the float `pi/4` and `pi/2`, now outward enclosures;
  - an overstated application claim.

**CITED (classical).**
- Hoeffding (1963, Thm 4 and §6): the Chernoff–Hoeffding bound holds for sampling without replacement.
- Robbins (1955): Stirling bounds.

The cycle lemma is proved here.

**FINITE-EXACT** (exact integer DPs).
- Theorem 1 against `N_m(L)` for all `1 <= L <= 400`, `2 <= m <= 401`: `max ln(N/bound) = -0.622`.
- Theorem 2 against `A_M(L)` for `4 <= M <= 59`, `L <= 400`: minimum log-margin `2.59`.
- `max N_m(L)/A_(m+2)(L) = 1.0000889` at `(L,m) = (58,9)` over `L <= 400` and all `m`. For comparison, `max N_m/A_(m+1) = 1.003255` at `(50,7)`.
- RM(1, 1.25) and RM(2, 1.075) hold at the first 60 zone times, for `k <= 400` and `m in {4,...,20}`.
- `N_m(L) <= 1.25 A_(m+1)(L) - 0.25 F_L(0)` for `m < 40`, `L <= 400`.

**VERIFIED (float DP).** Theorem 1 at `L = 256, 512, 1024, 2048, 4096` for `m <= 60` and `m = 65, 70, ..., 120`, with margins `e^-4.2` to `e^-6.6`.

**EMPIRICAL.**
- Per-step growth rates are ordered `Robin(m) < H ln2 + ln kappa(pi/(m+2)) < Dirichlet(m+2)` (`L = 20000`).
- Effective widths from the growth rates at `L = 20000`: Robin `m + 0.53` (`m = 3`) rising to `m + 1.15` (`m = 16`); Dirichlet `M + 0.43`.
- `ln kappa(theta) = -sigma^2 theta^2/2 - 0.00496 theta^4 + ...`.
- `A_M(L) >= 2^(HL) kappa(pi/M)^L / M^3` is observed (runner: `M^2/2 <= L <= 400`; the audit: all `L <= 1024`, `M <= 59`). So the proved `P_2` is crude.

**OPEN.**
- Conjecture R with a constant `K_0` (for any `K`).
- The polynomial version with `K = 1`. It needs boundary-layer-corrected test functions (§7); pure trigonometric ones cannot reach it.
- RM for all `k` and all zone phases.
- HYP-9140 for the *consistent* price `delta_L` (interference; pairpeak note §§6–8).

Nothing here bears on the Collatz conjecture. These are fixed-horizon modification prices. No novelty or priority claim.

## 0. Setting

Notation follows the pairpeak note §0 and §3:
- `c = log_3 2`, `sigma^2 = c(1-c) = 0.232857`, `lambda = (1-c)/c = 0.584963`, `H = H_2(c) = 0.949956` bits;
- `Q` = iid Bernoulli(`c`) letters;
- `S_t = e_t - ct` is the letter walk of a word (`e_t` = number of ones among the first `t` letters).

**Counting via the tilt** (peak note §4). For every word `w` of length `L`, `Q(w) 2^(HL) lambda^(e(w) - cL) = 1`. Hence every set `E` of words satisfies `|E| = 2^(HL) E_Q[lambda^(e - cL); E]`.

**Hard wall.** For real `M > 0`:

    A_M(L) = #{u in {0,1}^L : S_t > 0 (1<=t<=L), S_t < M (0<=t<=L-1)} = 2^(HL) E_Q[lambda^(S_L); same event].

**Reflected barrier** `Pi^(m)` (pairpeak note §3).
- The modified word `omega` drives `S'_t = U_t - ct`.
- `Z_m = [m-1, m-c)` is the zone. At a zone state both letters move `S'` down by `c`, and the letter 1 is a flip (`J += 1`).
- `F_L = {S'_t > 0, 1 <= t <= L}` and `N_m(L) = |F_L|`.
- Since `e(omega) - cL = S'_L + J_L` (pairpeak Lemma 3), `N_m(L) = 2^(HL) E_Q[lambda^(S'_L + J_L); F_L]`.
- `S'` never reaches `m - c` (pairpeak Lemma 4(a)). The Robin state space is therefore `D_m = [0, m-c)`.

## 1. Free eigenfunctions and the zone condition

**Lemma E.** Let `0 < theta < pi/c` and define `beta(theta)`, `kappa(theta)` as in the Status.

**(i)** For every `phi` and every real `s`, `f(s) = e^(beta s) sin(theta s + phi)` satisfies `c f(s+1-c) + (1-c) f(s-c) = kappa(theta) f(s)`.

*Proof.*
- Let `g(z) = c e^(z(1-c)) + (1-c) e^(-zc)` (the `Q`-mgf of one step) and `z = beta + i theta`.
- `Im g(z) = c e^(beta(1-c)) sin((1-c)theta) - (1-c) e^(-beta c) sin(c theta)`. This vanishes iff `e^beta = (1-c) sin(c theta)/(c sin((1-c)theta))`, which is the definition of `beta`.
- So `g(z) = Re g(z) = kappa(theta)`, and `c f(s+1-c) + (1-c) f(s-c) = Im[e^(i phi) g(z) e^(zs)] = kappa f(s)`. ∎

(Runner §1: 300 random triples, relative error `4.6e-49` at 50 digits.)

**(ii)** `beta(theta) < 0`.

*Proof.* `h(x) = ln(sin x/x)` is strictly decreasing on `(0, pi)` (`h' = cot x - 1/x < 0`), and `beta = h(c theta) - h((1-c) theta)` with `c theta > (1-c) theta`. ∎

Since `h <= 0`, also `|beta(theta)| <= -h(c theta) <= x^2/6 + x^4/150` with `x = c theta <= 1`. (The series of `-h` has positive coefficients `1/6, 1/180, 1/2835, ...`.)

**(iii) Bounds on `kappa`.** Let `0 < theta <= 0.15`, `B = |beta| <= b_1 theta^2` with `b_1 = c^2 (1/6 + (0.15c)^2/150) = 0.066369`, `x = c theta`, `x' = (1-c) theta`.

*Lower bound.*
- Use `e^(-y) >= 1-y`, `e^y >= 1+y` and `cos u >= 1 - u^2/2`; all factors are nonnegative.
- `kappa >= c(1 - B(1-c))(1 - x'^2/2) + (1-c)(1 + Bc)(1 - x^2/2) = 1 - sigma^2 theta^2/2 - sigma^2 B (2c-1) theta^2/2`.
- Hence `kappa >= 1 - sigma^2 theta^2/2 - eps_4 theta^4`, with `eps_4 = sigma^2 (2c-1) b_1/2 = 0.00202`.

*Upper bound.*
- Use `e^(-y) <= 1 - y + y^2/2`, `e^y <= 1 + y + y^2 e^y/2` and `cos u <= 1 - u^2/2 + u^4/24`.
- The `B`-free terms give `1 - sigma^2 theta^2/2 + (c(1-c)^4 + (1-c)c^4) theta^4/24`, and the last term is `0.00292 theta^4`.
- The terms linear in `B` combine to `sigma^2 B[-(2c-1)theta^2/2 + (c^4 - (1-c)^4) theta^4/24]`, which is at most `sigma^2 B (c^4-(1-c)^4) theta^4/24`.
- The `B^2` terms are at most `sigma^2 B^2 ((1-c) + c e^(Bc))/2`.
- Hence `kappa <= 1 - sigma^2 theta^2/2 + C_4 theta^4`, with `C_4 = 0.00344`, and so `ln kappa <= kappa - 1 <= -(sigma^2/2 - C_4 (0.15)^2) theta^2 <= -0.1163 theta^2`.

*On `[0.15, pi/2]`.* Interval arithmetic gives `ln kappa(theta) < -0.1 theta^2` (1831 boxes).

So `0 < kappa(theta) < 1` and `ln kappa(theta) <= -0.1 theta^2` on `(0, pi/2]`. ∎ (Runner §1 re-checks both Taylor bounds on a grid.)

**Lemma Z.** `F(theta) = kappa(theta) sin 2theta - 2(1-c) e^(-beta c) sin((2+c)theta) > 0` for `0 < theta <= pi/4`.

*Proof.*
- **Analytic part, `theta <= 0.15`.** Here `e^(-beta c) = e^(Bc) <= 1 + b_1 c e_b theta^2` with `e_b = e^(b_1 c (0.15)^2)`. Also `sin 2theta >= 2theta - 4theta^3/3 > 0` and `0 <= sin((2+c)theta) <= (2+c)theta`.
- With (iii):

      F/theta >= (1 - sigma^2 theta^2/2 - eps_4 theta^4)(2 - 4theta^2/3) - 2(1-c)(2+c)(1 + b_1 c e_b theta^2) >= G_0 - G_1 theta^2,

  where `G_0 = 2(c + c^2 - 1) = 0.058004` and `G_1 = 4/3 + sigma^2 + 2(1-c)(2+c) b_1 c e_b = 1.647587`.
  - The dropped `theta^4` coefficient `(sigma^2/2)(4/3) - 2 eps_4 = 0.151` is positive.
  - `G_0 - G_1 (0.15)^2 = 0.0209 > 0`.
- **On `[0.15, pi/4]`.** Adaptive interval bisection with `mpmath.iv` at 30 digits proves `F > 0` with 129 boxes (runner §2). ∎

*Margin.* The leading coefficient is `2(1 - (1-c)(2+c)) = 0.058`. So the zone condition at width `m + 2` holds with 2.9% to spare as `m -> infinity`, and with much more for small `m` (`F/theta = 0.34` at `theta = pi/4`).

## 2. Theorem 1: the reflected barrier at the sharp rate

**Theorem 1.** For all integers `m >= 2` and `L >= 0`, with `theta = pi/(m+2)`,

    N_m(L) <= e^(|beta(theta)|(m-c)) 2^(HL) kappa(theta)^L,     |beta(theta)|(m-c) <= 0.063.

*Proof.*
1. **Tilt.** `N_m(L) = 2^(HL) E_Q[lambda^(S'_L + J_L); F_L] <= 2^(HL) E_Q[lambda^(J_L); F_L]`, since `lambda < 1` and `S'_L > 0` on `F_L`.
2. **Kernel.** Under `Q` the letters are iid Bernoulli(`c`) and drive `S'`, so `E_Q[lambda^(J_L); F_L] = (P_R^L 1)(0)`, where `P_R` acts on functions on `D_m = [0, m-c)`:
   - `P_R f(s) = c f(s+1-c) + (1-c) f(s-c) 1[s-c > 0]` for `s < m-1`;
   - `P_R f(s) = (c lambda + 1 - c) f(s-c) = 2(1-c) f(s-c)` for `s in Z_m`.

   From `s < m-1` the up-step lands below `m - c`. On `Z_m`, `s - c > 0` since `m >= 2`.
3. **Supersolution.** Let `g(s) = e^(beta s) sin(theta(s+c))`. On `D_m`, `theta(s+c) in [theta c, theta m) ⊂ (0, pi)`, so `g > 0`. We show `P_R g <= kappa g` on `D_m`.
   - **`s < m-1`.** By Lemma E, `c g(s+1-c) + (1-c) g(s-c) = kappa g(s)`. So `P_R g(s) = kappa g(s) - (1-c) g(s-c) 1[s-c <= 0]`. For `s - c in (-c, 0]`, `g(s-c) = e^(beta(s-c)) sin(theta s) >= 0`.
   - **`s in Z_m`.** Put `v = m + 2 - c - s in (2, 3-c]`. Then `theta(s+c) = pi - theta v` and `theta s = pi - theta(v+c)`, so the claim reads `2(1-c) e^(-beta c) sin(theta(v+c)) <= kappa sin(theta v)`.
     - `v -> sin(theta(v+c))/sin(theta v)` is decreasing on `(0, pi/theta - c)`: its log-derivative is `theta[cot theta(v+c) - cot theta v] < 0`.
     - The range is admissible, since `3 theta <= 3pi/4 < pi`.
     - So the worst case is `v -> 2`, which is `F(theta) >= 0`: Lemma Z, as `theta <= pi/4`.
4. **Iteration.** `P_R` is positive, so `P_R^L g <= kappa^L g` on `D_m`. With `1 <= g/inf_(D_m) g` this gives `(P_R^L 1)(0) <= kappa^L g(0)/inf g`.
5. **Constant.**
   - `g(0) = sin theta c`.
   - On `D_m`, `e^(beta s) >= e^(-|beta|(m-c))`.
   - By concavity, `sin theta(s+c) >= min(sin theta c, sin theta m)`. Here `sin theta m = sin 2theta > sin theta c`, because `theta c < 2theta <= pi/2`.
   - Hence `g(0)/inf g <= e^(|beta|(m-c))`.
   - By Lemma E(ii), `|beta|(m-c) <= ((c pi)^2/6 + (c pi)^4/(150(m+2)^2)) (m-c)/(m+2)^2 <= 0.0625` (the maximum `0.06244` is at `m = 3`). Numerically `|beta|(m-c) <= 0.0411`. ∎

**Remarks.**
- **Tightness of the width.** The width `m + 2` can be lowered to `m + delta_R(m)`, where `delta_R(m)` is the least `delta` with the zone condition at `theta = pi/(m + delta)`:

  | `m` | 2 | 3 | 5 | 10 | 32 | 128 | `infinity` |
  |---|---|---|---|---|---|---|---|
  | `delta_R` | 0.993 | 1.212 | 1.447 | 1.650 | 1.761 | 1.777 | `2c(1-c)/(2c-1) = 1.7785` |

  The limit is forced: the zone condition says the sine's top zero lies at least `2c(1-c)/(2c-1)` above the top `m - c` of `D_m`, and the bottom zero is at `-c`. So `2` is the best integer shift for this family.
- **Comparison with Theorem A** (pairpeak note), all quantities `/2^L`:

  | `L` | `m` | exact | Theorem 1 | Theorem A |
  |---|---|---|---|---|
  | 180 | 4 | 2.00e-8 | 6.04e-6 | 1.94e-3 |
  | 180 | 8 | 4.72e-6 | 2.51e-4 | 1.94e-3 |
  | 400 | 4 | 8.11e-17 | 2.40e-12 | 8.98e-7 |
  | 400 | 8 | 9.67e-11 | 9.62e-9 | 9.42e-7 |
  | 400 | 16 | 1.10e-9 | 2.32e-7 | 9.42e-7 |

  Theorem 1 has the right exponential rate up to the width shift. Theorem A loses a factor about 200 in the exponent.
- **Per-step rates at `L = 20000`, EMPIRICAL.** The ordering `Robin(m) < H ln 2 + ln kappa(pi/(m+2)) < ln-rate of A_(m+2)` holds; for example, at `m = 8`: `0.644217 < 0.646919 < 0.647848`.

## 3. Theorem 2: the hard wall from below

**Theorem 2.** For every real `M >= 4` and integer `L >= 1`, `A_M(L) >= 2^(HL) kappa(pi/M)^L / P_2(M)` with `P_2` as in §3.4.

### 3.1 Subsolution

Let `psi(s) = e^(beta s) sin(pi s/M)` with `beta = beta(pi/M) < 0`, and let `P_M` be the `Q`-kernel killed outside `(0, M)`.

**Claim.** `P_M psi >= kappa psi` on `(0, M)`, where `kappa = kappa(pi/M)`.
- By Lemma E, killing only drops terms at landing points `(-c, 0]` and `[M, M+1-c)`, where `sin(pi s'/M) <= 0` (since `M > 1`).
- Also `0 <= psi <= 1` on `(0, M)`.

Hence, for `x in (0, M)`, with `C_k = {S_t in (0, M), 1 <= t <= k}`:

    Q_x(C_k) = (P_M^k 1)(x) >= (P_M^k psi)(x) >= kappa^k psi(x).        (3.1)

Also `psi(1-c) >= e^(-|beta|) sin(pi(1-c)/M) >= 2(1-c) e^(-|beta|)/M`.

### 3.2 Three bridge lemmas

Given reals `x, x'` and `n >= 1` with `u = x' - x + cn` an integer in `[0, n]`, the bridge `x -> x'` is the uniform random arrangement of `u` steps `+(1-c)` and `n - u` steps `-c`. Then:
- `Q_x(S_n = x') = b(n,u) = C(n,u) c^u (1-c)^(n-u)`;
- conditionally on `S_n = x'`, the increments form the bridge.

**B1 (binomial).** For `1 <= u <= n-1`, `b(n,u) >= e^(-n D(u/n || c))/(2 sqrt n)`, and `D(p||c) <= (p-c)^2/(2v)` with `v = min x(1-x)` over the segment between `p` and `c`.
- Robbins gives `C(n,u) >= e^(-1/6) sqrt(n/(2 pi u(n-u))) e^(n H_e(u/n)) >= 0.675 e^(n H_e)/sqrt n`.
- `c^u (1-c)^(n-u) = e^(-n(H_e(p) + D(p||c)))`.
- `D(p||c) = int_c^p (p-x)/(x(1-x)) dx`. ∎
- For `p in [1/2, 2c-1/2]`: `v = sigma^2` if `p < c`, and `v >= v_up = (2c-1/2)(3/2-2c) = 0.181430` if `p > c`.

**B2 (cycle lemma, proved).** Let `eta_1, ..., eta_n` be reals with partial sums `P_k` and total `P_n > -b` for some `b >= 0`. Then some cyclic rotation has all partial sums `> -b`; if `b = 0` they are all `> 0`. Hence a uniformly random arrangement of a fixed multiset with total `> -b` has all partial sums `> -b` with probability `>= 1/n`.

*Proof.*
- Let `j` be the **last** index in `{0, ..., n-1}` with `P_j = min_(0<=i<=n-1) P_i`, and rotate to start after `j`.
- For `1 <= k <= n-1-j` the new partial sums are `P_(j+k) - P_j > 0`, since `j` is the last minimiser.
- For `k >= n - j` they are `P_n - P_j + P_(k-n+j) >= P_n > -b`, because `P_i >= P_j` and `P_j <= P_0 = 0`.
- **Probability.** Rotations permute arrangements. Summing over all arrangements `a` the number of good rotations of `a` gives `n` times the number of good arrangements, and each `a` has at least one good rotation. ∎

(Runner §4: exhaustive for `n <= 10`.)

**B3 (Hoeffding 1963, Thm 4, CITED).** For a uniform arrangement of a multiset of reals in an interval of length 1 with mean `mu`, `P(P_k - k mu >= t) <= e^(-2t^2/k)`, and the same with `<= -t`. (Runner §4: exhaustive check for `n <= 12`.)

### 3.3 The end lemma

Let `n(M) = max{n : 4 n ln(2n) <= (M-2)^2}` and `n_min(M) = ceil((M/2+1)/(c-1/2))`.

For real `M >= M_1 = 110`, `n(M) >= n_min(M)`:
- On `[110, 3000]` this is an interval-safe grid check. Both functions are nondecreasing, and the check is `n(M_a) >= n_min(M_a + 1/4)`.
- For `M >= 3000`: `n(M) >= floor((M-2)^2/(8 ln(M-2))) >= 4M >= n_min(M)`. The first inequality holds because `n' = floor((M-2)^2/(4 ln((M-2)^2)))` has `2n' <= (M-2)^2`, so `4n' ln(2n') <= 4n' ln((M-2)^2) <= (M-2)^2`.

**Lemma B4 (end lemma).** Let `M >= 110` be real and `n = n(M)`. For every real `y in (0, M)`:

    Q_y(C_(2n), S_(2n) in (0,1]) >= q(M) := (16 n^3)^(-1) exp(-(M/2+1)^2 (1/v_up + 1/sigma^2)/(2n)).

*Proof.*
- **Targets.**
  - `x_1` = the point of `y + Z - cn` in `[M/2, M/2+1)`;
  - `e` = the point of `x_1 + Z - cn` in `(0, 1]`.

  Each unit interval contains exactly one point of a coset of `Z`. By the Markov property, `Q_y(C_(2n), S_(2n) = e) >= Q_y(C_n, S_n = x_1) · Q_(x_1)(C_n, S_n = e)`. Each factor is a binomial term times a bridge confinement fraction.
- **Binomial parameters.**
  - `|x_1 - y| < M/2 + 1` and `|e - x_1| < M/2 + 1`.
  - With `n >= n_min`, `u/n in (1/2, 2c-1/2)` and `1 <= u <= n-1`.
  - By B1 each binomial term is at least `(2 sqrt n)^(-1) e^(-(M/2+1)^2/(2vn))`:
    - phase A upward (`y < x_1`): `v = v_up`;
    - phase A downward: `v = sigma^2`;
    - phase B, always downward: `v = sigma^2`.
- **Confinement, margin `t = M/2 - 1`.** The choice of `n` gives `n e^(-2t^2/n) <= 1/(4n)`, i.e. `4n^2 e^(-(M-2)^2/(2n)) <= 1`.
  - **(A1) `y <= M/2`.** The total `x_1 - y = u - cn` is `>= 0` and nonzero (`cn` is irrational), so by B2 at least `1/n` of the arrangements have all partial sums `> 0` and stay above `y > 0`. Leaving through `M` needs a deviation `>= M - max(y, x_1) > t` above the mean line, which by B3 and a union bound has probability `<= n e^(-2t^2/n)`.
  - **(A2) `M/2 < y < x_1`.** Both walls are more than `t` from the mean line; B3 twice.
  - **(A3) `y > x_1`.** Apply B2 to the negated increments (all partial sums `< 0`, hence the path stays below `y < M`). The bottom is more than `t` below the mean line; B3.
  - **(B) `x_1 -> e`.** Reverse and negate: `eta_j = -xi_(n+1-j)` is again a uniform arrangement, with total `x_1 - e > 0`. By B2, at least `1/n` of the arrangements have all partial sums `> 0`, i.e. `x_1 + P_k = e + Pi_(n-k) > 0`. The top is more than `t` from the mean line; B3.

  In every case the confinement fraction is `>= 1/n - 1/(4n) >= 1/(2n)`.
- **Multiply.** `(1/(2n))^2 (2 sqrt n)^(-2) = 1/(16 n^3)`, and the worst exponent is `(M/2+1)^2 (1/v_up + 1/sigma^2)/(2n)`. ∎

Runner §4 checks the exact bridge fractions and the conclusion at `M = 110, 150`. The true minimum of `Q_y(...)` over 60 starting points is `5.4e-16` and `2.7e-17`, against `q = 2.5e-25` and `5.4e-27`. The lemma is crude by about `e^21`.

### 3.4 Assembly

`A_M(L) >= 2^(HL) E_0[lambda^(S_L); C_L]`, and the first step must be `+(1-c)`, so `E_0[lambda^(S_L); C_L] = c E_(1-c)[lambda^(S_(L-1)); C_(L-1)]`.

- **(crude, all `M >= 4`).** On `C_L`, `lambda^(S_L) > lambda^M`, so by (3.1)
  `E_0[...] >= lambda^M c kappa^(L-1) psi(1-c) >= kappa^L lambda^M c 2(1-c) e^(-|beta|)/M`.
  This gives `P_crude(M) = M e^(|beta|)/(2c(1-c) lambda^M)`.
- **(long, `M >= 110`, `L - 1 >= 2n`).**
  - Keep only paths with `S_(L-1) in (0,1]`, where `lambda^(S_(L-1)) >= lambda`.
  - Apply the Markov property at time `L - 1 - 2n` and Lemma B4 (valid for *every* position in `(0, M)`), then (3.1):
    `E_(1-c)[lambda^(S_(L-1)); C_(L-1)] >= lambda q(M) Q_(1-c)(C_(L-1-2n)) >= lambda q(M) kappa^(L-1-2n) psi(1-c)`.
  - Since `kappa < 1`, this gives `P_long(M) = M e^(|beta|)/(2c(1-c) lambda q(M))`.
- **(short, `M >= 110`, `8 <= L-1 < 2n`).**
  - There is one bridge `1-c -> e` of length `k = L-1`, with `|u - ck| < 1`.
  - By B1 the binomial term is `>= e^(-1/(16 v_up))/(2 sqrt k)`.
  - B2 with `b = 1-c` keeps the path above 0 for a fraction `>= 1/k`.
  - B3 for the top (margin `M - 1`): `k e^(-2(M-1)^2/k) <= 1/(2k)`, using `k < 2n` and `4n ln 2n <= (M-2)^2`.
  - This gives `P_short = 4 (2n)^(3/2) e^(1/(16 v_up))/(c lambda)`.
- **(tiny, `L <= 8`).** The all-ones word is counted (`7(1-c) < 4 <= M`). So `A_M(L) >= 1 >= 2^(HL) kappa^L 2^(-8H)`.
- **Result.**
  - `P_2(M) = P_crude(M)` for `M < 110`;
  - `P_2(M) = min(P_crude, max(P_long, P_short, 2^(8H)))` for `M >= 110`.

**Values of `P_2`.**

| `M` | 4 | 16 | 64 | 108 | 200 | `10^3` | `10^4` | `10^6` |
|---|---|---|---|---|---|---|---|---|
| `ln P_2(M)` | 4.3 | 12.1 | 39.2 | 63.4 | 70.7 | 94.2 | 130.0 | 203.7 |

**Growth of `P_2`.** Write `ln P_long = ln M + |beta| - ln(2c(1-c) lambda) + ln(16 n^3) + (M/2+1)^2 (1/v_up + 1/sigma^2)/(2n)`, with `-ln(2c(1-c)lambda) = 1.3004` and `1/v_up + 1/sigma^2 = 9.80625`.
- From `n <= (M-2)^2/(4 ln 2)`: `ln(16 n^3) <= 6 ln M - 0.2867`.
- From `n >= floor((M-2)^2/(8 ln(M-2)))`: `(M/2+1)^2/(2n) <= ((M+2)/(M-2))^2 ln(M-2)(1 + 10^-5) <= 1.00237 ln M` for `M >= 3400`.
- Hence `ln P_long(M) <= 16.8295 ln M + 1.02` for `M >= 3400`. `P_short` and `2^(8H)` are much smaller.
- So `ln P_2(M) <= 16.83 ln M + 1.31` for `M >= 3400`. The exact asymptotic is `16.81 ln M - 7.9 ln ln M + O(1)`. ∎

## 4. Corollaries

**Corollary 3.** For all `m >= 2`, `L >= 1`: `N_m(L) <= P(m) A_(m+2)(L)` with `P(m) = e^(|beta(pi/(m+2))|(m-c)) P_2(m+2) <= (m+2)^17`.

*Proof.* Theorems 1 and 2 at `M = m + 2` have the same `kappa(pi/(m+2))`. The bound `P(m) <= (m+2)^17` holds by computation for `m < 10^5` (`max ln P(m) - 17 ln(m+2) = -15.85`), and by `ln P(m) <= 16.83 ln(m+2) + 1.38` for `m + 2 >= 3400`. ∎

This is HYP-9142 with `K = 2` and a polynomial `K_0`. The observed constant is `1.0000889`, and `1.003255` for `K = 1`.

**How the true ratio behaves in `m` and `L`** (exact, `L <= 400`, runner §5):

| `m` | `K=1`: max over `L` (at `L`) | `K=1`: at `L=400` | `K=2`: max over `L` (at `L`) | `K=2`: at `L=400` |
|---|---|---|---|---|
| 4 | 1.000773 (15) | 6.8e-3 | 1 (all `L <= 14`) | 1.0e-4 |
| 6 | 1.003155 (39) | 0.33 | 1.0000059 (23) | 0.067 |
| 7 | 1.003255 (50) | 0.54 | 1.0000250 (31) | 0.19 |
| 9 | 1.002996 (88) | 0.81 | 1.0000889 (58) | 0.50 |
| 12 | 1.001516 (153) | 0.95 | 1.0000523 (104) | 0.82 |
| 15 | 1.000667 (226) | 0.99 | 1.0000186 (142) | 0.95 |

- The excess over 1 appears only in the crossover `L ~ m^2`.
- For fixed `m` the ratio then decays. EMPIRICAL: the Robin growth rate is strictly below the Dirichlet one (§2, Remarks). Theorems 1–2 prove only the non-strict comparison, and only for `K = 2`.
- For `m >= 20` the maximum lies near or beyond `L = 400`.

**Corollary 4 (HYP-9140 for the private price).** For `L >= 3`,

    pi_L <= (9 (27L + 18)(L+2)^17 + 4L 2^(-L)) rho^peak_L.

*Proof.* Apply Theorem C of the pairpeak note with `K = 2` and `K_0 = max_(2<=m<=L) P(m) <= (L+2)^17`. It uses Theorem B(2) with `m_top = L >= 3`. ∎

It depends on the pairpeak note's Theorem B (private certificates from the barrier flips, Lemma 5) and Theorem C.

**Corollary 5 (the second-order term).** Let `Lambda(L) = max_(a>=2) [L ln kappa(pi/a) - a ln 3]`. Then:
1. `Lambda(L) = -kappa_3 L^(1/3) + O(L^(-1/3))`, where `kappa_3 = min_theta [theta ln 3 + pi^2 sigma^2/(2theta^2)] = 2.107580`.
2. `ln(2^((1-H)L) rho^peak_L) >= Lambda(L) - O(log L)` (Theorem 2).
3. `ln(2^((1-H)L) rho^peak_L) <= Lambda(L) + O(log L)` (Dirichlet supersolution, below).
4. `ln(2^((1-H)L) pi_L) <= Lambda(L) + O(log L)` (Theorem 1 and pairpeak Theorem B(2)).

Since `2 rho^peak_L <= pi_L` (pairpeak Theorem B(1)), both quantities are `-kappa_3 L^(1/3) + O(log L)`.

*Proof.*
1. **Asymptotics of `Lambda`.** Put `phi(a) = a ln 3 - L ln kappa(pi/a)`, so `Lambda(L) = -min_a phi(a)`.
   - For `a >= pi/0.15`, Lemma E(iii) gives `-ln kappa(pi/a) >= sigma^2 pi^2/(2a^2) - C_4 pi^4/a^4`.
   - If `a >= eps L^(1/3)`: `phi(a) >= kappa_3 L^(1/3) - C_4 pi^4/(eps^4 L^(1/3))`.
   - If `a < eps L^(1/3)`: `phi(a) >= 0.1 pi^2 L/a^2 >= kappa_3 L^(1/3)` for `eps` small, using `ln kappa <= -0.1 theta^2`.
   - At `a* = theta* L^(1/3)`, `theta* = (pi^2 sigma^2/ln 3)^(1/3) = 1.278935`: `phi(a*) <= kappa_3 L^(1/3) + O(L^(-1/3))` by the lower Taylor bound.
   - Numerically `(Lambda(L) + kappa_3 L^(1/3)) L^(1/3) -> -0.1807`. This matches `-0.00496 pi^4/theta*^4`.
2. **Lower bound.**
   - A word in the event of `A_a(L)` is bad and has peak `M(u) < a`.
   - So `2^((1-H)L) rho^peak_L >= 3^(-a) 2^(-HL) A_a(L) >= 3^(-a) kappa(pi/a)^L/P_2(a)` by Theorem 2 at `a = a*` (`a* >= 4` for `L >= 31`).
   - `ln P_2(a*) = O(log L)`. This is asymptotic: `P_2(a*)` is polynomial only once `a* >= 110`, i.e. `L >~ 6·10^5`. Below that the crude `P_crude(a*) = e^(O(a*))` applies, which is what the table uses.
3. **Upper bound for `rho^peak`.** By the displayed strata bound below, `2^((1-H)L) rho^peak_L <= (c/(1-c)) sum_(m>=0) 3^(-m) e^(|beta|(m+1)) kappa(pi/(m+2))^(L-1)`.
   - `|beta(pi/(m+2))|(m+1) <= 0.171` for all `m >= 0`, by Lemma E(ii); numerically the maximum is `0.113`.
   - With `a = m + 2`, each term with `m <= L` is at most `(c/(1-c)) 9 e^0.171 · 3^(-a) kappa(pi/a)^(L-1) <= (c/(1-c)) 9 e^0.171 e^(Lambda(L-1))`.
   - `Lambda(L-1) = -kappa_3 L^(1/3) + O(L^(-1/3))` by step 1.
   - The terms with `m > L` sum to at most `(c/(1-c)) e^0.171 · 3^(-L) · 3/2`, which is negligible.
   - So the bound is `O(L) e^(-kappa_3 L^(1/3))`.
4. **Upper bound for `pi_L`.**
   - Theorem B(2) of the pairpeak note with `m_top = L` bounds `pi_L` by `2L[3^(1-L) rho_L + sum_(m=2)^(L-1) 3^(1-m) N_(m+1)(L)/2^L] + 2 N_2(L)/2^L`.
   - By Theorem 1, `3^(1-m) N_(m+1)(L)/2^L <= e^0.063 2^(-(1-H)L) 81 · 3^(-(m+3)) kappa(pi/(m+3))^L <= 81 e^0.063 2^(-(1-H)L) e^(Lambda(L))`. The term `N_2` is treated the same way.
   - `3^(1-L) rho_L <= 3^(1-L) 2^(-(1-H)L)` is negligible.
   - Hence `2^((1-H)L) pi_L <= O(L^2) e^(Lambda(L))`. ∎

**The Dirichlet supersolution and the strata bound used in step 3.**
- *Dirichlet supersolution.* For real `W >= 1`, `f(s) = e^(beta s) sin(pi(s+c)/(W+1))` satisfies `P_W f <= kappa(pi/(W+1)) f` on `[0, W)`: it is positive on `(-c, W+1-c)`, which contains every killed landing point. So `Q_0(0 < S_j < W, 1<=j<=n) <= (c/(1-c)) e^(|beta|W) kappa(pi/(W+1))^n`, using `sin(theta c)/sin(theta(1-c)) < c/(1-c)`, i.e. `beta < 0`.
- *Summing over the peak strata `M(u) in [m, m+1)`* gives

      2^((1-H)L) rho^peak_L <= (c/(1-c)) sum_(m>=0) 3^(-m) e^(|beta|(m+1)) kappa(pi/(m+2))^(L-1).

**Rigorous bracket (runner §6).** Values of `ln(2^((1-H)L) rho^peak_L)`:

| `L` | `Lambda(L)` | `-kappa_3 L^(1/3)` | lower (Thm 2) | exact | upper (supersolution) | upper for `ln(2^((1-H)L) pi_L)` |
|---|---|---|---|---|---|---|
| 100 | -9.822 | -9.783 | -15.24 | -10.67 | -5.75 | 1.12 |
| 200 | -12.356 | -12.325 | -18.76 | -13.58 | -8.20 | -0.58 |
| 400 | -15.553 | -15.529 | -23.16 | -17.19 | -11.30 | -2.97 |
| 1000 | -21.094 | -21.076 | -30.68 | -23.33 | -16.71 | -7.46 |
| 2000 | -26.568 | -26.554 | -38.04 | -29.31 | -22.08 | -12.13 |

At `L = 1000` the elementary bracket of the peak note was `[-319.9, -2.09]`; it is now `[-30.7, -16.7]`.

## 5. The sharp Conjecture R: reduction to one walk (Theorem 6)

**Definitions.** Fix `m >= 2`, `K >= 1`, `M = m + K`.
- For real `x`, `V_k(x)` (resp. `F_k(x)`) counts the words `w in {0,1}^k` with `x + S^w_j in (0, M)` (resp. `(0, m-1)`) for `1 <= j <= k-1` and `x + S^w_k > 0`.
- The zone contains a lattice point at time `t` iff `{ct} > c`, and then that point is `z_t = m - {ct}`.

**Theorem 6.** Let `C >= 1` and `L >= 1`. Suppose that for every zone time `1 <= t <= L-1`,

    V_k(z_t - c) - V_k(z_t + 1 - c) <= (2 - 2/C) F_k(z_t - c)   with k = L - 1 - t.        (RM)

Then `N_m(L) <= C A_M(L) - (C-1) F_L(0)`.

*Proof.* Let `W_n(x)` be the reflected-walk count of the `n`-step continuations from a point `x` below the zone. We prove by strong induction on `n` the stronger statement `W_n(x) <= C V_n(x) - (C-1) F_n(x)`, for all `x` below the zone and all `n` that occur.
- **First zone visit.** Split `W_n(x)` and `V_n(x)` at the first zone visit, at step `j` and point `z` (the walk cannot jump over the zone, which has width `1-c`):
  - `W_n(x) = F_n(x) + sum_j E_j(x) · 2 W_(n-j-1)(z - c)`;
  - `V_n(x) = F_n(x) + sum_j E_j(x) [V_(n-j-1)(z-c) + V_(n-j-1)(z+1-c)]`.
- **Induction step.** By the induction hypothesis and (RM):
  - `2W <= 2CV(z-c) - 2(C-1)F(z-c) <= C[V(z-c) + V(z+1-c)]`;
  - summing gives `W_n <= F_n + C(V_n - F_n)`.
- **Which `k` are needed.** The needed `k` is the remaining time `L - 1 - t`, because the horizon is fixed. ∎

**Why `F`?** Without it, (RM) with `C = 1` would read `V_k(z-c) <= V_k(z+1-c)`. That fails for intermediate `k`, for every `K`: before the bottom is felt, the walk from `z + 1 - c` is simply closer to the top wall. So `K_0 = 1` is not provable this way, and indeed it is false for `K = 1, 2`.

**Numbers.** Maximum of `(V_k(y) - V_k(y+1))/F_k(y)`, `y = z_t - c`, over the first 60 zone times (`t <= 163`) and `k <= 400`:

| `K` \ `m` | 4 | 6 | 8 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|
| 1 | 0.077 | 0.255 | 0.315 | 0.361 | 0.374 | 0.379 |
| 2 | 0.000 | 0.023 | 0.064 | 0.107 | 0.123 | 0.128 |

So RM holds on this set with `C = 1.25` (`K = 1`) and `C = 1.075` (`K = 2`).

The ratio peaks around `k ≈ 2-3m` (scratch exploration; the runner checks only the maxima). Two regimes can be seen:
- **Half-line regime** (the bottom not yet reached). The ratio rises toward a limit near 0.39. This regime is a half-line problem for the drifting counting walk, and a Lundberg-type bound on up-crossing probabilities should control it (not carried out).
- **Long-time regime.** The ratio turns negative, because `V_k(y+1)/V_k(y)` tends to a limit above 1 (about 1.05 for `K = 1`, `m = 12`). This limit is the ratio of the principal eigenfunction, which carries the counting-measure factor `e^(alpha s)`.

A proof of RM for all `k` needs this long-time eigenfunction monotonicity with quantitative convergence. It is OPEN.

## 6. What is now proved about HYP-9142, HYP-9140 and THM-4480

- **HYP-9142 (Conjecture R).** Proved up to a polynomial with `K = 2`: `N_m(L) <= (m+2)^17 A_(m+2)(L)` (Corollary 3).
  - The constant version stays OPEN. It is reduced to the one-walk inequality RM (Theorem 6).
  - Finite-exact evidence: `sup N_m/A_(m+1) = 1.003255`, `sup N_m/A_(m+2) = 1.0000889`.
- **HYP-9140, private price.** PROVED, given Theorems B–C of the pairpeak note: `pi_L <= O(L^18) rho^peak_L` (Corollary 4).
  - The polynomial is crude; on the tested range the private cost is `6-9 x 2 rho^peak` (pairpeak §4).
  - The *consistent* price `delta_L` is untouched. Blocking by frozen certificates remains the obstruction.
- **Sharp constant of `pi_L`.** PROVED, elementary, given pairpeak Theorem B: `ln pi_L = -(1-H)L ln 2 - kappa_3 L^(1/3) + O(log L)`. Neither Conjecture R nor Mogul'skii is needed.
- **THM-4480 Theorem 4 (`q = 3`).** Its statement `ln rho^peak_L = -(1-H)L ln 2 - kappa_3 L^(1/3)(1+o(1))`, modulo Mogul'skii, can be strengthened to an elementary `+ O(log L)` (Corollary 5).
  - Since `rho^peak_L/M*_L <= eps_L(3) <= rho^peak_L` with `M*_L = O(L^3)` (THM-4480 Theorem 1), the same holds for the arbitrary-edit price `eps_L(3)`.
  - With `rho_L = 2^(-(1-H)L)/poly(L)` (peak note Consequence (b)), `pi_L/rho_L` and `rho^peak_L/rho_L` are `exp(-kappa_3 L^(1/3) + O(log L))`. This sharpens pairpeak Theorem B(4), which had `exp(-Theta(L^(1/3)))`.
  - The method uses `lambda < 1`, so `q >= 5` is not covered here.

## 7. Failures and caveats

- **No injection.** A direct injection from barrier survivors into hard-wall words is impossible for `K = 1, 2`, since the ratios exceed 1. Swap-based encodings of the zone bits also fail to decode: the multiplicities are up to `2^(#zone visits)`.
- **Monotonicity alone fails.** The one-walk inequality `V_k(z-c) <= V_k(z+1-c)` would give `K_0 = 1`. It fails at intermediate `k`; at `m = 12` the minimum of `V_k(y+1)/V_k(y)` is `0.8644` for `K = 1` and `0.998366` for `K = 4` (runner §7). Hence the `F`-slack of Theorem 6.
- **One-sided bounds lose the endpoint.** In the counting picture the eigenfunctions carry a factor `e^(alpha s)`, `alpha = ln(1/lambda)`. So a bound of the form `sum_paths 1 >= sum psi(end)/max psi` loses `e^(alpha M)`: paths end at the bottom, `psi` peaks at the top. Lemma B4 is the repair, at polynomial cost.
- **A block argument was superseded.** Blocks of length `A m^2` with a per-block Dirichlet supersolution give the sharp constant only with `exp(o(L^(1/3)))` accuracy. Theorem 1 gives it directly.
- **`K = 1` is out of reach of pure trigonometric test functions.**
  - The Robin supersolution needs width `>= m + 1.7785` asymptotically.
  - The Dirichlet subsolution has width exactly `M`.
  - The true effective widths (EMPIRICAL, runner §3) are about `m + 1.15` and `M + 0.43` at `m = 16`. So `K = 1` holds with a width margin of about 0.28 there; a scratch computation gives about 0.23 at `m = 40`.
  - Reaching it needs boundary-layer-corrected functions on both sides.
- **Crude constants.** The proved `P(m)` is `(m+2)^17`, against an observed ratio `<= 1.0001`. Lemma B4 is crude by about `e^21` at `M = 110`. The rigorous bracket of Corollary 5 has width 9.5 at `L = 100` and 16.0 at `L = 2000`.
- **Computer-assisted steps.**
  - Lemma Z on `[0.15, pi/4]` (129 interval boxes).
  - The bound `ln kappa < -0.1 theta^2` on `[0.15, pi/2]` (1831 boxes).
  - The grid check of `n(M) >= n_min(M)` on `[110, 3000]`.
  - All three are re-run by the runner.
- **Dependencies.** Corollary 4 and the `pi_L` part of Corollary 5 rely on Theorems B–C of the pairpeak note (lane-proved, same session). Theorems 1–2, Corollary 3 and the `rho^peak` part of Corollary 5 do not.
- **Floating point.** The float DP checks (Theorem 1 at `L >= 256`, the rates) are floating point. The exact checks use Python integers. An earlier float comparison in §7 of the runner produced a false alarm; it now uses the integer form `4N <= 5A - F`.

## 8. Reproduction

```
python3 04-computation/experiments/procgen_robin_20260926_run.py > 05-knowledge/results/procgen_robin_20260926.out
```

- **Environment.** Python 3.10.0, numpy 2.2.6, mpmath 1.3.0 (interval arithmetic `mpmath.iv`).
- **Dependencies.** The runner imports read-only:
  - `procgen_pairpeak_20260926_lib.py` (exact floors, `N_m(L)` and `A_M(L)` DPs);
  - `procgen_peak_20260926_lib.py` (`H_2`, `kappa_3`, `rho_L` and `rho^peak_L` float DPs).
- **Final run.**
  - Wall time 47.5 s on the shared 8 GB machine.
  - Peak RSS 183 MB by `/usr/bin/time -l`; the runner itself reports 175 MB.
  - 178 output lines, 55 checks, `ALL CHECKS PASSED`.
- **Determinism.** A second run gave identical output apart from the `[time]` lines (47.8 s, 195 MB).

| file | sha256 (raw bytes) |
|---|---|
| `04-computation/experiments/procgen_robin_20260926_lib.py` | `6bcb072646218a5a9653e2f0b57e81228b13009d7b3337bdca6442482c7a1831` |
| `04-computation/experiments/procgen_robin_20260926_run.py` | `96d13fab3ef6acf17cca132665fadfb3683af28acc08fb1c233a5bd2a0874852` |
| `05-knowledge/results/procgen_robin_20260926.out` | `cfea2a012e96f5c71cfb6834caf0e0be7c7280a243a3b7491f7bcd51ac577d80` |
| same output without `[time]` lines (`grep -v '^\[time\]'`) | `55a42cc0c3ac150aa1cecee0cdeec6a588c89912e901a8c80cd6886a40449042` |
| `04-computation/experiments/procgen_pairpeak_20260926_lib.py` (imported, unchanged) | `aa915a333d913f2ada76a90ee1dc6602a17179f38b24b7f23267053f39b4ec56` |
| `04-computation/experiments/procgen_peak_20260926_lib.py` (imported, unchanged) | `ecc5f8a27c13ca40fad95ae8f6486a88a982d4cfa386d5b3808ac307ec778070` |

**Runner sections.**
- 0 constants;
- 1 Lemma E and the `kappa` bounds;
- 2 Lemma Z;
- 3 Theorem 1;
- 4 Theorem 2 and its lemmas;
- 5 Corollary 3;
- 6 Corollaries 4–5;
- 7 Theorem 6.

Temporary files lived only under `scratch/procgen_robin/` (not part of the deliverable).

## 9. Hypothesis bookkeeping (no files created)

- **HYP-9142.** Status should become:
  - PROVED up to a polynomial with `K = 2` (Corollary 3);
  - OPEN for a constant `K_0`;
  - reduced to RM (Theorem 6).
- **HYP-9140.** PROVED for the private price `pi_L` (Corollary 4, given pairpeak Theorems B–C). OPEN for `delta_L`.
- **Promotion candidates after audit:**
  - Theorems 1–2 with Corollaries 3 and 5;
  - an addendum to THM-4480 Theorem 4 for `q = 3` (elementary `O(log L)`).
