# Theta round 2: theta_3(2^10/3^9) is not quadratic (HYP-9132, partial), and a no-go for cube determinants (HYP-9127)

**Status.**
* **AUDIT (orchestrator, 2026-09-23, partial).** The place-exchange dictionary `mu_bar = gamma/(gamma-1)` is consistent across the lane's results: `gamma = 7/3 -> 7/4` (H1), `3 -> 3/2`, `2 -> 2` (Stihl/Choulet for Euler E1/E2), and `3.27694 -> 1.43918` (NQ). The value `mu_bar(2^10/3^9) = 1.42647 < 1.43918` was re-computed. The transcribed KRVZ inputs (Propositions 2 and 4, Lemma 3 / (4.17)) were NOT independently re-derived. Theorem NQ is therefore PROVED modulo those cited propositions. Proposition C's exact `e3(n)` law and Theorem NG are elementary and as stated.
* **PROVED, given [KRVZ] Propositions 2 and 4 and Lemma 3 / (4.17)** (the same cited inputs as round-2 Lemma D and Theorem H3):
  * **Theorem NQ.** Let `rho = 2^L/M_0` with `M_0 > 2^L` odd and `mu_bar = log2(M_0)/L`. If `mu_bar < c_NQ = 126 pi^2/(79 pi^2 + 72 sqrt3 Im Li_2(e^(2 pi i/3))) = 1.43918`, then `Theta(rho) = sum_(k>=0) rho^(k^2)` lies in no quadratic field.
  * **HYP-9132, partial.** `rho = 2^10/3^9` has `mu_bar = 1.42647 < c_NQ`. So `theta_3(2^10/3^9) = 2 Theta - 1` has degree `>= 3` over `Q`, or is transcendental.
  * **The transcription does not break.** Every step of KRVZ's proof of their Theorem 1 goes through 2-adically; section 1.3 lists the steps.
  * **Corollary NQ-Y.** No 2-adic integer of degree `<= 2` over `Q` has an eventually-`Y` parity vector, where `Y` is the 3x+1 square-swap word. The same holds for every square-swap block word under `mx + r` with `m^A > 2^L` and `mu_bar < 1.43918`.
* **PROVED (elementary; section 2), for the cubes:**
  * **Lemma P (parabola).** The tails of `X = sum rho^(k^3)` form the parabola `{(3s^2, 3s)}` in the exponent plane of `G(u,v) = sum u^k v^(k^2) rho^(k^3)`. It contains no nondegenerate parallelogram and meets every line in at most 2 points.
  * **Theorem NG (two-place no-go).** Take a determinant of monomially normalized tails with any index pattern. Suppose there is no cancellation at the place 2 or at the place `M_0`. Then it certifies nothing that a single tail does not.
  * **Proposition L.** For Hankel-type patterns `s_ij = x_i + y_j` and normalizations `s^3 + c s^2`, there is no cancellation at either place, except in two cases. At `c = 0` (round-2 K7) the leading layer at `M_0` has rank one. At `c = 3` the leading layer at 2 has rank one.
  * **Proposition C (peeling).** Both rank-one cases are computed exactly by row differences.
    * At `c = 0` the `xi`-order is `e3(n) = (n-1)(3n^2 - 15n + 19)` for all `n >= 3`. This **proves** the law that round 2 found FINITE-EXACT.
    * At `c = 3`, `v_2 = L(n-1)(3n^2 - 3n + 1)`.
  * **Proposition V (the second variable).** `Z_n = G(1, rho^n)` is an exponential sum in `n`, and `v_2 det(Z_(i+j)) = L n(n-1)^2(5n-4)/12` exactly, which is quartic. But every entry except `Z_0 = X` lies off the `X`-orbit.
  * **Corollary H.** A Hankel certificate at size `n`, with any normalization `s^3 + c s^2` (`c` an integer), needs forced cyclotomic factors of total degree `>= 2(1 - 1/mu_bar) n^4 - O(n^3) = 0.598 n^4 - O(n^3)`.
* **FINITE-EXACT** (`.out`):
  * **Theorem NQ checks.** The algebra holds in `Z[sqrt 17]` and `Z[sqrt -7]` for `n <= 5`. For `Theta`-approximants in those fields, `v_2(iota_1 N_n)` equals the predicted Hankel valuation exactly.
  * **Margin at `rho = 2^10/3^9`.** `A_NQ(n)` is positive for `73 <= n <= 8000`. Quadratic equations with coefficients up to `2^10`, `2^100`, `2^1000` and `2^10000` are excluded at `n = 75, 87, 143, 351`.
  * **Cubes:**
    * the Proposition L classification: 16 patterns, 7 normalizations, 2 places;
    * Proposition C: `n <= 12` and `n <= 9`;
    * Hankel certificate margins for `n <= 10`. Every normalization reaches less far than a single tail: best 368.9 bits against 376.9.
    * The cyclotomic content grows like `<= 0.234 n^3` for `n <= 12`, and `m_1(n) = e_1(n)` there.
* **CITED:** [KRVZ] (local copy from round 2). This round read, `[P]`: Theorem 1, (2.1)–(2.6), Proposition 4 with Remark 6, and section 5 ((5.4)–(5.12) and the proof of Theorems 1 and 2). There was no web access this round.
* **OPEN:**
  * **HYP-9127** (cubes). No natural determinant family works (section 2.8). The one loophole left is a quartic cyclotomic content (HC-R3-2).
  * **HYP-9132 in full.** Degree `>= 3` is proved; transcendence is not. The KRVZ method stops at degree 2, because degree 3 would need `mu_bar < 0.9595`.
  * **Non-quadraticity for `mu_bar >= 1.43918`.** For example the HYP-9131 map, `mu_bar = 1.61831`.
  * **Non-quadraticity of Euler's function at `2^10/3^9`.** This method needs `mu_bar < 1.11860`.

Session `collatz-procgen-20260922`, theta round-2 lane (machine `mac-mini`, 2026-09-23). It follows the round-2 note [theta beyond phi](collatz_procgen_20260923_theta_beyond_phi.md), cited below as TBP, whose Theorem H1, Lemma D and Proposition K7 it uses.
* **No priority is claimed.** Theorem NQ is KRVZ's Theorem 1 with the places exchanged, so that the 2-adic place is analytic and the archimedean places are arithmetic.
* **Files.** No HYP or THM file was created.
  * Scripts: `04-computation/experiments/collatz_procgen_20260923_theta2_{nq,cube}.py` and `..._theta2_run.sh`.
  * Output: [`collatz_procgen_20260923_theta2.out`](collatz_procgen_20260923_theta2.out).

## 0. Answers in brief

1. **HYP-9132 (partial): yes.**
   * KRVZ's non-quadraticity theorem transcribes to `Q_2` step by step, with threshold `mu_bar < 1.43918`. This is the dictionary image `gamma/(gamma - 1)` of their `gamma > 3.27694`.
   * The one new ingredient is that `K = Q(Theta)` sits inside `Q_2`, so 2 splits in `K`. Then one 2-adic embedding carries the whole Hankel valuation, and the other contributes `>= 0`.
   * The margin at `2^10/3^9` is small but positive: `0.206 n^3` bits.
   * `theta_3(2^10/3^9)` has degree at least 3.
   * Euler's function and the HYP-9131 map lie outside the threshold.
2. **HYP-9127 (cubes): nothing works, and here is why.**
   * **What made squares work.** Square tails are exponential sums in the shift `s`, so every level of a Hankel matrix has rank one. That gives the Vandermonde cancellation.
   * **What cubes have instead.** Cube tails are the points of a parabola (Lemma P), and each level's exponent is a strictly convex quadratic in `s`.
   * **Theorem NG.** Determinants without cancellation are no better than their best single entry.
   * **Proposition L.** A normalization can make at most one level rank one, at one place, and Proposition C computes that exactly. The gain is cubic against a quartic clearing (Corollary H).
   * **The second variable (Proposition V).** It does have exponential-sum structure, with quartic valuation, but its entries are new unknowns.
   * **Weighted and residue-class companions** stay Gaussian, and mod 2 they collapse to polynomials in `F`.
   * **Loophole.** The only one left is forced cyclotomic factors of quartic total degree. The observed degree is `~0.23 n^3`.

## 1. Theorem NQ: `theta_3(2^10/3^9)` is not quadratic

Notation follows TBP section 1:
* `S_m = sum_(k<=m) rho^(k^2)`;
* `t_s = rho^(-s^2)(Theta - S_s)`;
* `T_n = det(t_(i+j))_(i,j<n)`, which has `v_2(T_n) = L n(n+1)(2n-1)/2` and so is nonzero;
* `K_n = 4 sum_(i<n) i^2`.

### 1.1 Statement and proof

**Theorem NQ.** If `M_0 > 2^L` is odd and `mu_bar < c_NQ = 14/(11 - 24 C_cyc) = 126 pi^2/(79 pi^2 + 72 sqrt3 Im Li_2(e^(2 pi i/3))) = 1.43918`, then `[Q(Theta):Q] >= 3` (or `Theta` is transcendental). Here `C_cyc = 5/54 - Im Li_2(e^(2 pi i/3))/(pi^2 sqrt3) = 0.05301135`.

*Proof.* Suppose `[Q(Theta):Q] <= 2`. The case `K = Q` is TBP Theorem H3 (`mu_bar < 2.87837`), so let `K = Q(Theta)` be quadratic.
1. **2 splits.** `K` is a subfield of `Q_2`, via the inclusion `iota_1`. So 2 has a place of degree 1 and ramification index 1 in `K`, that is, 2 splits. Let `iota_2` be the other embedding `K -> Q_2`.
   * Choose an integer `w >= 1` with `alpha := w Theta` in `O_K`.
2. **The determinant.** Put `G_m = alpha M_0^(m^2) - w sum_(k<=m) 2^(L k^2) M_0^(m^2-k^2)`, which lies in `O_K`, and `N_n = det(G_(i+j) 2^(L(i-j)^2))`, which lies in `O_K`.
   * Under `iota_1` we have `iota_1(G_m) = w M_0^(m^2)(Theta - S_m)`.
   * So TBP H1 step 3 gives `iota_1(N_n) = 2^(L K_n) w^n T_n`.
   * Hence `N_n != 0` and `v_2(iota_1 N_n) = L K_n + n v_2(w) + L n(n+1)(2n-1)/2`.
3. **The other 2-adic embedding.** `iota_2(O_K)` lies in `Z_2`, so `v_2(iota_2 N_n) >= 0`.
   * The norm `Norm(N_n) = iota_1(N_n) iota_2(N_n)` is a nonzero rational integer.
   * `v_2(Norm N_n) >= L K_n + L n(n+1)(2n-1)/2`.
4. **The forced odd divisor.** TBP Lemma D is a divisibility of polynomials in `Z[xi, a, b]`, transferred from [KRVZ] Propositions 2 and 4 in `Z[q, alpha, mu]`. Specializing `a = alpha`, `b = w` gives `G_n | N_n` in `O_K`, where
   * `G_n = M_0^(e'(n)) prod_(1<=l<n/2) Phi_l^hom(M_0^2, 2^(2L))^(e_l(n))`, an odd integer;
   * `e'(n) = 2 e_0(n) - n(n-1)`;
   * `e_l(n) = sum_(i<n)(floor((i+l)/3l) + floor(i/3l))`.

   Hence `G_n^2 | Norm(N_n)`.
5. **Archimedean embeddings.** TBP H1 step 4 gives, for each complex embedding `sigma` of `K`:
   * `|sigma(N_n)| <= n! C_sigma^n M_0^(K_n)`, where `C_sigma = |sigma(alpha)| + (2n-1) w`.
6. **Count.** Steps 3 to 5 give

   `A_NQ(n) := L K_n + L n(n+1)(2n-1)/2 + 2 log2 G_n - 2 log2 n! - 2 K_n log2 M_0 <= n log2(C_1 C_2)`.

   * By TBP Lemma D(iii),(iv): `A_NQ(n) = [(7/3) L - 2(11/12 - 2 C_cyc) log2 M_0] n^3 + O(n^2 log^2 n)`.
   * The bracket is positive exactly when `mu_bar < 14/(11 - 24 C_cyc)`.
   * The right side is `O(n log n)`. This is a contradiction for large `n`. ∎

**Remarks.**
* **Degree 3 is out of reach.** Degree 3 would need `(7/3) L > 3(11/12 - 2 C_cyc) log2 M_0`, that is, `mu_bar < 0.9595 < 1`. So, as in KRVZ (5.11), only degrees 1 and 2 are reachable.
* **Heights.** A root of `A x^2 + B x + C` with `max|A|, |B|, |C| <= H` has `alpha = A Theta` integral, `w = |A|` and `|sigma(alpha)| <= 2H`. The exclusion therefore needs `A_NQ(n) > 2n log2((2n+1) H)`.
* **Dictionary.** KRVZ (5.12) with `d = 2` and `A = 1/2`, `C = 2/3`, `B = 65/216 - Im Li_2/(pi^2 sqrt3) = 5/24 + C_cyc` gives `gamma > 3.27694`. Our constant is `3.27694/2.27694 = 1.43918` (`.out` Q0). Their `B` is exactly our forced-divisor exponent, measured in `q`-degree.

### 1.2 The numbers at `rho = 2^10/3^9` (`.out` Q1–Q3)

* **Leading coefficient.** `(7/3) 10 - 2(11/12 - 2 C_cyc) log2 3^9 = +0.20621` bits per `n^3`.
* **Sign of the margin.** `A_NQ(n) <= 0` for `n <= 72`, and `A_NQ(n) > 0` for `73 <= n <= 8000`, where every `n` was checked.
* **Growth.** `A_NQ(n)/n^3` equals 0.057, 0.131, 0.197 and 0.204 at `n = 100`, 200, 1600 and 6400; it tends to 0.206.
  * `e_l(n)` is used in closed form, checked against its definition for `n < 120`.
  * `log2 G_n` agrees with the exact integer to `1.5e-11` bits.
* **Explicit exclusions.** Quadratic equations with coefficients `<= 2^10`, `2^100`, `2^1000` and `2^10000` are excluded at `n = 75`, 87, 143 and 351.
* **Algebra in `O_K`** (Q1). For random `alpha` in `Z[sqrt 17]` and `Z[sqrt -7]` (2 splits in both, since `17 = -7 = 1 mod 8`), with 12 cases each and `n <= 5`:
  * `G_n^2 | Norm(N_n)` holds;
  * the archimedean bound holds with 3.8 to 4.2 bits to spare.
* **The one-embedding valuation** (Q2). LLL at precision `2^4000` gives `alpha = u + v sqrt d` and `w`, all of size `2^1332`, with `iota_1(alpha)/w = Theta mod 2^4000`. For `n = 2, ..., 5`:
  * `v_2(iota_1 N_n) = 130, 500, 1260, 2550`, exactly the prediction `L K_n + v_2(T_n)`;
  * `v_2(iota_2 N_n)` is between 2 and 10;
  * `G_n^2 | Norm` holds, `Norm != 0`, and the archimedean bound holds.

### 1.3 The transcription step by step (where it could have broken)

| KRVZ, section 5 | 2-adic version | status |
|---|---|---|
| `K = Q(alpha, lambda, mu)`, `d = [K:Q]` | `K = Q(Theta)`, a subfield of `Q_2`, so 2 splits | new, elementary |
| (5.4): `V_n` small at the analytic place | `v_2(iota_1 N_n)` exact (TBP H1) | proved |
| (5.5): forced divisor `Delta_n` | TBP Lemma D, specialized in `O_K` | cited ([KRVZ] Props. 2, 4, Lemma 3) |
| Lemma 5: conjugates at most `|q|^(C n^3)` | archimedean embeddings: permutation bound; second 2-adic embedding: `v_2 >= 0` | proved |
| Lemma 4: `V_n != 0` infinitely often | `N_n != 0` for **every** `n` (unique minimal Cauchy–Binet term) | proved, and stronger |
| (5.9): `sigma^(...) Omega(n) V~_n` in `Z_K` | `N_n` in `O_K` by construction | proved |
| (5.10)–(5.12) with `d = 2` | step 6 | proved |

No step breaks. The 2-adic version is slightly easier at Lemma 4: non-vanishing holds for every `n`, and no density argument is needed.

### 1.4 Corollaries and limits

* **Corollary NQ-Y.** Let `B != B'` be blocks of length `L` with `A` ones, with `m^A > 2^L` and `mu_bar = A log2(m)/L < 1.43918`.
  * The block identity of the [hard-class note](collatz_procgen_20260922_hard_class.md) (section 2.4) gives `Phi_T(Y_(B,B')) = -m^(-A)[R_B/(1 - rho) + (R_(B') - R_B)(Theta(rho) - 1)]` with `R_(B') != R_B`. This is a non-constant rational affine image of `Theta`.
  * Suppose a 2-adic integer `x` of degree `<= 2` had parity vector `u Y_(B,B')`. Then `T^|u|(x) = Phi_T(Y_(B,B'))` would lie in `Q(x)`, and `Theta(rho)` would have degree `<= 2`.
  * So no such `x` exists. The rational case is Theorem Y.
  * Examples: `Y` itself (3x+1, `2^10/3^9`), 3x+1 with `2^8/3^7` (1.38684), 5x+1 with `2^5/5^3` (1.39316), 7x+1 with `2^2/7` (1.40368).
* **Not reached.**
  * The HYP-9131 map (`2^33/5^23`, 1.61831).
  * Euler's function at `2^10/3^9`. The same count with KRVZ's `lambda != 0` exponents has leading coefficient `-2.75` bits; the threshold is `mu_bar < 1/(1 - 2 C_cyc) = 1.11860` (Q4).
* **What HYP-9132 still needs.** Transcendence, which is beyond any Hankel/KRVZ count (degree 3 would need `mu_bar < 0.96`). A 2-adic Nesterenko-type theorem is still the natural route; p-adically only transcendence degree `>= 2` is known (see HYP-9132).

## 2. Cubes: the natural determinant families and a no-go

### 2.1 Setting: two places and the parabola

**Tails and single tails.**
* `X = sum_(k>=0) rho^(k^3)` and `R_(s+1) = X - S_s = sum_(k>s) rho^(k^3)`.
* If `X = a/b`, then `G_s = M_0^(s^3)(a - b S_s)` is an integer with `v_2(G_s) = L(s+1)^3 + v_2(b)` and `|G_s| <= C_s M_0^(s^3)`, where `C_s = |a| + (s+1)|b|`.
* **Single-tail certificate.** `G_s != 0`, so `a/b` is excluded as soon as `L[(s+1)^3 - mu_bar s^3] > log2 C_s - v_2(b)`.
* At `mu_bar = 1.42647` the best `s` is 5, and the single-tail reach is `L B(mu_bar) = 376.9` bits.

**Determinant certificates.** A determinant certificate combines three things for the integer `N`, and compares them with the permutation bound:
* its 2-adic valuation;
* its forced `M_0`-power, which is the `xi`-order of the polynomial identity, with `xi = M_0/2^L`;
* its forced cyclotomic factors `Phi_d^hom(M_0, 2^L)`.

**Lemma P (PROVED; `.out` C1).**
* **The functional equation.** `G(u,v) = sum_k u^k v^(k^2) rho^(k^3)` satisfies `G(u,v) = 1 + u v rho G(u v^2 rho^3, v rho^3)`.
* **The tails lie on a parabola.** `R_s = rho^(s^3) G(rho^(3s^2), rho^(3s))`. So the `X`-orbit is the parabola `P(s) = (3s^2, 3s)` in the plane of exponents `(a, b)` of `G(rho^a, rho^b)`.
* **No parallelograms, no three collinear points.** The parabola has no nondegenerate parallelogram, and meets every line in at most 2 points.

*Proof.* `P(s1) + P(s2) = P(s3) + P(s4)` gives equal sums and equal sums of squares, hence equal pairs. Equivalently, for a sum pattern:

`(x1+y1)^2 + (x2+y2)^2 - (x1+y2)^2 - (x2+y1)^2 = 2(x1-x2)(y1-y2) != 0`. ∎

**Why this is the obstruction.**
* **Squares.** The level-`m` exponent is `m^2 + 2ms`, affine in `s` for every `m`. So every level of a Hankel matrix has rank one: that is the Vandermonde cancellation of TBP H1.
* **Cubes.** The level-`m` exponent is `m^3 + 3m^2 s + 3m s^2`, the linear functional `<(m, m^2), P(s)> + m^3` on the parabola. Its quadratic coefficient `3m` depends on `m`, so no normalization `rho^(-psi(s))` makes two levels affine at once.

### 2.2 Theorem NG: without cancellation, a determinant is no better than its best tail

**Setting.**
* **Entries.** Let `A_ij = lambda_ij b R_(s_ij + 1)`, with any index pattern `s_ij` and any monomials `lambda_ij = +-2^(alpha_ij) M_0^(beta_ij)`.
* **Clearing.** Clear the entries to integers `n_ij = 2^(u_i + w_j) M_0^(p_i + r_j) A_ij` by row and column factors.
* **Per-entry quantities.** For each entry:
  * `v_ij` is its 2-adic valuation under `X = a/b`;
  * `w_ij = p_i + r_j + beta_ij - s_ij^3` is its forced `M_0`-exponent;
  * `h_ij` is `log2` of its permutation-bound size;
  * the single-tail margin is

    `mu_ij = v_ij + w_ij log2 M_0 - h_ij = L[(s+1)^3 - mu_bar s^3] - log2 C_s + v_2(b)`, with `s = s_ij`.

  The margin `mu_ij` does not depend on the normalization.
* **Quantities for `N = det(n_ij)`.**
  * `V_2 = v_2(N)`;
  * `V_M` is the forced `M_0`-exponent;
  * `U = log2 n! + max_pi sum_i h_(i pi(i))`.

**Theorem NG (PROVED).** Suppose `V_2 = min_pi sum_i v_(i pi(i))`, attained at `pi_0`, and `V_M = min_pi sum_i w_(i pi(i))`. That is, `N` has no cancellation at the place 2 or at `M_0`. Then

`V_2 + V_M log2 M_0 - U = sum_i mu_(i pi_0(i)) - log2 n! - Gap`,

where `Gap = [max_pi sum h - sum_(pi_0) h] + [sum_(pi_0) w - V_M] log2 M_0 >= 0`.

So if the determinant excludes `a/b` and has no forced cyclotomic factor, then a single tail `R_(s+1)`, with `s = s_(i pi_0(i))` for some `i`, already excludes `a/b`. Forced cyclotomic factors add their size to the right side.

*Proof.* `V_2 = sum_(pi_0) v = sum_(pi_0)(mu + h - w log2 M_0)`. Add `V_M log2 M_0` and subtract `U`. ∎

*Remark (sketch).* The permutation bound is essentially attained when `|b| << |a|`: the `a^n` coefficient is dominated by one permutation (for a sum pattern with convex `psi`, the identity). So no better archimedean bound holds uniformly over rationals of height `<= H`.

### 2.3 Proposition L: which normalizations cancel (PROVED; `.out` C2)

Take a sum pattern `s_ij = x_i + y_j`, with `x` and `y` strictly increasing, and the normalization `rho^(-psi(s_ij))` times row and column monomials.

* **(a) Leading terms.**
  * At 2, the leading term of every entry is level `m = 1`, with exponent `L c2(s)`, where `c2(s) = (s+1)^3 - psi(s)`. Higher levels add at least `7L`.
  * At `M_0` (the `xi`-adic order, as a polynomial identity in `a` and `b`), the leading term is `k = s`, with exponent `cM(s) = psi(s) - s^3`. Other terms add at least 1. At `s = 0` the `a`-term merges with it, with coefficient `a - b`.
* **(b) No cancellation when the cost is strictly convex or concave.** Suppose `c2` is strictly convex (resp. concave) on the range of `s`.
  * Then the leading assignment problem has the unique optimum `pi(i) = n - 1 - i` (resp. the identity), by the strict Monge exchange.
  * So there is no cancellation at 2, and `V_2 = L sum_i c2(s_(i pi(i)))`, plus the separable row and column terms.
  * The same holds at `M_0` with `cM`.
* **(c) One place stays strict.** `c2 + cM = (s+1)^3 - s^3 = 3s^2 + 3s + 1` is strictly convex. So on every rectangle of the pattern at least one place is strictly Monge.
* **(d) The two degenerate cases.** For `psi = s^3 + c s^2` (plus affine terms; `c` an integer), `c2 = (3-c)s^2 + 3s + 1` and `cM = c s^2`.
  * Both places are free of cancellation unless `c = 3` or `c = 0`.
  * At `c = 3`, `c2` is affine: the leading layer at 2 is rank one and all `n!` matchings tie.
  * At `c = 0`, `cM = 0`: this is TBP K7.
* **Check.** 16 patterns (Hankel and random, `n <= 6`), 7 normalizations (raw, `c = -1, ..., 4`), 2 places. In every non-degenerate case the valuation equals the unique minimal matching. There is cancellation in exactly the two degenerate cases, 16 patterns each.
* **Leading order for `psi ~ lambda s^3` with `lambda != 1`** (a computation, not a full proof).
  * The quartic part of `V_2 + V_M log2 M_0 - U` is `L n^4 [(1-lambda)(1 - 2 mu_bar) - 2 lambda(mu_bar - 1)]` for `lambda <= 1`, and `L n^4 [2 - mu_bar(1 + lambda)]` for `lambda >= 1`.
  * Both are negative for `mu_bar > 1`.

### 2.4 Proposition C: the two rank-one normalizations, exactly (PROVED; `.out` C3)

Both parts use the Hankel pattern `s = i + j`.

**(i) `psi = s^3` (TBP K7).**
1. **The entries.** `F_ij(xi) = b t_(i+j) = a xi^(s^3) - b sum_(k<=s) xi^(s^3 - k^3)`.
2. **Row differences.** Replace row `i` by row `i` minus row `i-1`, for `i = n-1` down to 1. The new entries are `b(t_(sigma+1) - t_sigma)`, with `sigma = i + j - 1`.
3. **Their `xi`-orders.**
   * Row 0: order 0.
   * Entry `(1, 0)`: order 0, with coefficient `-a`.
   * Otherwise: order `E(sigma) = 3 sigma^2 - 3 sigma + 1`, with coefficient `+b` (or `-(a - b)` at `sigma = 1`).
4. **The optimal assignment.** Put `E'(0) = 0` and `E'(sigma) = E(sigma)` otherwise; `E'` is strictly increasing and strictly convex. The unique optimal assignment sends row 0 to column `n-1` and row `i >= 1` to column `n-1-i`, so that all `sigma = n - 2`.
5. **The order.** Hence `ord_xi det F = (n-1) E(n-2) = e3(n) = (n-1)(3n^2 - 15n + 19)` for `n >= 3`, and 0 for `n = 2`.
   * The leading coefficient is `+-b^n` for `n >= 4`, so this is also the exact forced order.
   * This **proves** TBP's law, which was FINITE-EXACT for `3 <= n <= 12`.

**(ii) `psi = (s+1)^3`, equivalently `c = 3`.**
1. **The entries.** `t'_s = rho^(-(s+1)^3) R_(s+1) = 1 + rho^(3s^2 + 9s + 7) + ...`, so the leading 2-adic layer is the all-ones matrix.
2. **The same row differences** leave row 0 of order 0, and rows `i >= 1` of order `L E_2(i+j-1)`, with `E_2(sigma) = 3 sigma^2 + 9 sigma + 7`.
3. **The same assignment** gives `v_2 det(t'_(i+j)) = L(n-1) E_2(n-2) = L(n-1)(3n^2 - 3n + 1)`, with leading coefficient `+-1`.
4. **At `M_0` there is no cancellation**, and `V_M = n(3n^2 - 3n + 1)`.

**Checks.**
* (i): `n <= 12`, from exact forced divisors; uniqueness of the assignment enumerated for `n <= 8`.
* (ii): `n <= 9`.

**Conclusion.** One layer peels off, and the gain is cubic (about `3n^3`) at one place. The next layer is again a strictly convex Gaussian, as Lemma P predicts.

### 2.5 The Hankel family: margins and the cyclotomic loophole (`.out` C4, C4b)

**Margin formula.** For the normalizations `psi_c`:

`margin_n = V_2 + V_M log2 M_0 + sum_d m_d log2|Phi_d^hom(M_0, 2^L)| - (mu_bar - 1) L sum_(i<n) psi(2i) - log2 n!`.

A rational with `|a| + (2n-1)|b| <= C` is excluded at size `n` if `margin_n > n log2 C`.
* `V_2` and `V_M` are exact for every integer `c` (Propositions L and C), and cubic.
* The permutation term is `-2(mu_bar - 1) L n^4 (1 + O(1/n))` for every `c`.

**Corollary H (PROVED).** A Hankel certificate for `X` at size `n`, with any normalization `s^3 + c s^2`, needs the forced cyclotomic part of `det F` to have degree `>= 2(1 - 1/mu_bar) n^4 - O(n^3) = 0.598 n^4 - O(n^3)`.

**FINITE-EXACT data.**

| normalization | best reach (bits) | at `n` | margin negative from |
|---|---|---|---|
| `c = 0` (K7) | 200.1 | 5 | `n = 8` |
| `c = 1` | 84.9 | 3 | `n = 5` |
| `c = 2` | 73.6 | 3 | `n = 5` |
| `c = 3` | 368.9 | 5 | `n = 8` |
| single tail | 376.9 | `s = 5` | – |

* **Cyclotomic content, `c = 0`.** The degree is 6, 13, 20, 37, 67, 106, 152, 221, 307, 404 for `n = 3, ..., 12`. That is at most `0.234 n^3`, against the roughly `0.598 n^4` needed (12399 at `n = 12`).
  * Only `Phi_d` with `d <= 30` occur.
  * `m_1(n) = e_1(n)` exactly, with `e_1` from [KRVZ] Remark 6, for all `n <= 12`.
  * `m_2 = m_3 = m_6` throughout.
* **Cyclotomic content, `c = 3`.** It agrees with `c = 0` to within 2 for `n <= 10`.

In range, no determinant certificate beats a single tail. Asymptotically they all fail, unless the cyclotomic content turns quartic (HC-R3-2).

### 2.6 The second variable (Proposition V; PROVED; `.out` C5)

* **An exponential sum in `n`.** `Z_n = G(1, rho^n) = sum_k rho^(k^3) (rho^(k^2))^n`.
* **Cauchy–Binet.** It gives `det(Z_(i+j)) = sum_(k_1<...<k_n) prod rho^(k_l^3) prod_(l<l') (rho^(k_l'^2) - rho^(k_l^2))^2`.
  * The term for `{k_l}` has valuation `L[sum k_l^3 + 2 sum_(l<l') k_l^2]`.
  * It is uniquely minimal at `{0, ..., n-1}`.
  * So `v_2 = L n(n-1)^2(5n-4)/12`, checked for `n <= 8`.
* **So the answer to "does the second variable give an exponential sum?" is yes**, with quartic valuation `(5/12) L n^4`.
* **But the entries are new unknowns.** The points `(0, n)` lie on the line `a = 0`, which meets the parabola only at `X = Z_0`.
  * `Z_1, Z_2, ...` are values of `G` off the `X`-orbit. No relation with `Q + Q X` is known or expected.
  * Each point `(0, n)` generates its own orbit `{(3s^2 + 2ns, 3s + n)}` under the functional equation, so a Hankel matrix of size `n` in `v` uses `2n - 1` orbits.
  * Under "`X` rational" the determinant is not an integer.
* **What it does prove.** Only a joint statement: `Z_0, ..., Z_(2n-2)` are not all rationals with a common denominator and numerators `<= H`, unless `n log2 H + log2 n! >= L n(n-1)^2(5n-4)/12`.

### 2.7 Hermite–Padé companions (`.out` C6)

* **Weighted tails.** The tails `sum_(k>s) k^j rho^(k^3)` have the same Gaussian leading exponent; the weights shift valuations only by `j v_2(s+1)`.
  * Two-sequence block-Hankel determinants (weights `k^0` and `k^1`, `n = 2, ..., 8`) have unique minimal matchings, and their valuation equals it.
  * So Theorem NG applies.
* **Mod 2 the companions collapse** (Frobenius).
  * `sum_k k^j x^(k^3) = sum_(k odd) x^(k^3) = F + F^8` in `F_2[[x]]` for `j >= 1`.
  * The even part is `F_even = F(x^8) = F^8`.
  * So every weighted or mod-2 residue companion is a polynomial in `F` over `F_2`. The `F_2` parity half of a zero estimate for `(F, F_1, ...)` gains no new linear information.
* **Residue classes mod 3.** `F_(0 mod 3)(x) = F(x^27)` brings in `X(rho^27)`, the same problem at another `rho` with the same `mu_bar`. The classes `k = +-1 mod 3` give new Gaussian sums.
* **Dimension count** (a Siegel's-lemma heuristic). Type-I forms in `r + 1` companions with common support `{k^3}` need `mu_bar < (r+2)/(r+1)`. At 1.4265 only `r <= 1` survives, and `r = 1` is a harder variant of HYP-9130.

### 2.8 Verdict on HYP-9127

**HYP-9127 stays OPEN.** The natural determinant families all fail:
* **Hankel and sum-pattern determinants**, with any cubic normalization: Theorem NG, Propositions L and C, Corollary H.
* **Weighted and residue-class block-Hankel determinants**: no cancellation, so Theorem NG applies.
* **The second-variable Hankel**: exponential-sum structure and quartic valuation, but off the orbit.

None certifies more than a single tail.

**On "quartic forced valuation".** Quartic 2-adic valuation alone is easy and useless.
* The raw Hankel determinant has `V_2 = L n^4` exactly, against a clearing of `2 mu_bar L n^4`.
* A certificate needs cancellation beyond the minimal matching at some place. Lemma P limits that to one cubic layer.

**The one loophole** is forced cyclotomic content of quartic degree. The observed degree is about `0.23 n^3`.

## 3. Hypothesis candidates (for `INDEX.md`; no HYP files created)

* **HYP-9132 → PARTIAL.** Suggested status line: "degree `>= 3` PROVED (theta_round2 Theorem NQ, given [KRVZ] Props. 2, 4); transcendence OPEN; no quadratic 2-adic integer has an eventually-`Y` parity vector."
* **HC-R3-1 (non-quadraticity up to the irrationality threshold; OPEN).** `Theta(rho)` is not quadratic whenever `mu_bar < 2.87837`. The first case is the HYP-9131 map. It needs either a forced divisor about twice KRVZ's, or conjugate control better than the trivial bound.
* **HC-R3-2 (cyclotomic content of cube Hankel determinants; OPEN, FINITE-EXACT for `n <= 12`).** The forced cyclotomic part of the K7 determinant has degree `O(n^3)`, and `m_1(n) = e_1(n)`. A proof would make the Hankel no-go for HYP-9127 unconditional (Corollary H).
* **HC-R3-3 (Euler non-quadraticity; OPEN).** `(rho;rho)_inf` at `2^10/3^9` is not quadratic. This method reaches only `mu_bar < 1.11860`.

## 4. Sources

* **[KRVZ]** C. Krattenthaler, I. Rochev, K. Väänänen, W. Zudilin, *On the non-quadraticity of values of the q-exponential function and related q-series*, Acta Arith. 136 (2009) 243–269, arXiv 0812.2921v1.
  * `[P]` this round, from the local copy downloaded in round 2: Theorem 1, (2.1)–(2.6), Proposition 4 with Remark 6, and section 5 ((5.4)–(5.12), the proof of Theorems 1 and 2, Lemmas 4 and 5 as used there).
  * `[P]` in round 2: Propositions 1 and 2 with proofs, Lemma 3 and (4.17).
* **Internal.**
  * The [TBP note](collatz_procgen_20260923_theta_beyond_phi.md): Theorem H1, Lemma D, Proposition K7, the D6 exact forced divisors.
  * The [hard-class note](collatz_procgen_20260922_hard_class.md): block identity, Theorem Y.
  * The files `HYP-9127*`, `HYP-9130*` and `HYP-9132*`, read and not modified.
* No web requests this round.

## 5. Reproduction

```bash
bash 04-computation/experiments/collatz_procgen_20260923_theta2_run.sh          # full: ~100 s, peak ~100 MB
bash 04-computation/experiments/collatz_procgen_20260923_theta2_run.sh --quick  # smoke: ~10 s
```

* The full run writes `05-knowledge/results/collatz_procgen_20260923_theta2.out`: the `nq` sections Q0–Q4 and the `cube` sections C1–C6.
* Requirements: python3 with gmpy2, python-flint 0.9, sympy, mpmath and numpy.
* Timing probes: `scratch/procgen_theta2/timing_{det,gcd}.py`.
