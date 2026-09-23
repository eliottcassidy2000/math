# Theta values beyond phi: a 2-adic Hankel-determinant proof of HYP-9131, every square-swap word under 5x+1 and 7x+1, Euler's function, and why cubes stay out of reach

**Status.**
* **AUDIT (orchestrator, 2026-09-23).** Theorem H1 was re-derived in full, with no gap found:
  * the tails are exponential sums;
  * Cauchy–Binet gives a unique minimal-valuation term at `{1..n}`, so `v_2(T_n) = L n(n+1)(2n-1)/2`;
  * the row and column scaling makes `N_n` an integer with `v_2 = L n(2n-1)(7n-1)/6`;
  * the permutation bound gives `|N_n| <= n! C^n M_0^(K_n)`;
  * the leading terms compare `7/3` against `4/3`, giving `mu_bar < 7/4`.

  HYP-9131 (`mu_bar = 1.61831`) is therefore PROVED by the elementary H1 alone. Theorem E's exponent arithmetic re-checks (`v_2(N_n) = L n^3`), but its integrality input is the transcribed [KRVZ] (2.6), which is cited, not re-derived. H2 and H3 rest on the transcribed [KRVZ] Propositions 2 and 4, also cited.
* **PROVED** (this note; exact checks in the `.out`):
  * **Theorem H1** (elementary and self-contained). `Theta(rho) = sum_(k>=0) rho^(k^2)` is irrational in `Q_2` whenever `rho = 2^L/M_0` (`M_0` odd) has `mu_bar < 7/4`.
  * **HYP-9131 is PROVED.** For the 5x+1 square-swap word with 23 ones in 33 letters, `mu_bar = 1.61831 < 7/4`.
  * **Theorem H2**, using [KRVZ] Proposition 2: irrational for `mu_bar < 28/11 = 2.54545`. So **every square-swap block word under 5x+1** is settled.
  * **Theorem H3**, using [KRVZ] Propositions 2 and 4: irrational for `mu_bar < 252 pi^2/(79 pi^2 + 72 sqrt3 Im Li_2(e^(2 pi i/3))) = 2.87837`. So **every square-swap word under 7x+1**, and under 9x+1 when `A/L < 0.908`, is settled.
  * **Theorem Q+.** The same thresholds hold for every single quadratic family with signs: triangular numbers, either pentagonal branch, `2k^2 + k`, and so on. H1 is proved for all families. H2 and H3 transfer to the listed families by the same specialization, but that step is a sketch; only their `v_2` laws are checked.
  * **Theorem E (Euler's function).** `(rho;rho)_inf = prod(1 - rho^n) = sum_(k in Z) (-1)^k rho^(k(3k-1)/2)` is irrational in `Q_2` for `mu_bar < 3/2` (elementary), `< 2` (with [KRVZ] Proposition 1) and `< 2.23719` (with Proposition 4).
    * This covers **`rho = 2^10/3^9`** (`mu_bar = 1.42647`, elementary level).
    * It also covers every signed bilateral pentagonal swap word under `3x+r`. Such words exist: (L, A) = (9, 6), (10, 7), (11, 7), ....
  * **Unsigned bilateral pentagonal.** `sum_(k in Z) rho^(k(3k-1)/2)` is irrational for `mu_bar < 21/16` (elementary).
  * **Proposition K7 (the method's limit for cubes).** The cube Hankel determinant has `v_2 = L(3n^3 - 3n^2 + n)` exactly, against a quartic clearing, so the method yields only finite certificates for cubes. **HYP-9127 stays OPEN.**
  * **Lemma F2 and Lemma E** (HYP-9130 in restricted form). The parity half of the zero estimate follows from `F_2`-independence of the Padé constraint rows. `F_2`-independence holds, provably, on the echelon family `E = N^3 - 1` (`w -> 1`).
* **FINITE-EXACT:**
  * `v_2` laws of all Hankel determinants: 8 maps, `n <= 12`, and families, `n <= 9`.
  * Integer identities, archimedean bounds, and the divisibilities from [KRVZ] Propositions 1, 2 and 4, transcribed: `n <= 10`.
  * Margin tables. At the elementary level, heights `<= 2^(10^9)` are excluded for HYP-9131 at `n = 13305`.
  * The exact forced polynomial divisor of the Hankel determinants: squares `n <= 14`, cubes `n <= 7`.
  * The cube `X`-order law `e3(n) = 3n^3 - 18n^2 + 34n - 19` for `3 <= n <= 12`.
  * HYP-9130 data:
    * the `F_2`-rank frontier satisfies `W_2(E) >= 1.4265 E` for all 446 tested `E <= 1500`;
    * a Berlekamp–Massey profile to `N = 40000`;
    * odd-leading witnesses for every tested `E <= 600`.
  * The "1089" tests.
* **CITED:**
  * `[P]`: [KRVZ] (Theorems 1 and 2, the construction of section 2, Lemma 1, Propositions 1 and 2 with proofs, Proposition 4, Lemma 3, (4.17), Proposition 5).
  * `[R via KRVZ]`: Bézivin 1998; Choulet 2001; Stihl 1983; Tschakaloff and Bundschuh.
* **OPEN:**
  * HYP-9127 (cubes); HYP-9130 in full.
  * Square swaps with `mu_bar >= 2.87837`. The cheapest is 9x+1 with 10 ones in 11 letters, `mu_bar = 2.88175`.
  * The unsigned pentagonal sum for `mu_bar >= 21/16`.
  * HYP-9132 (transcendence). The non-quadraticity transcription is not done (section 8, HC-TH2).
* **REFUTED:**
  * The previous note's framing that `phi` is the frontier for 2-adic theta values ("cheapest open instance #1", HC-CT4 as the route). The literature had passed `phi` since Bézivin 1998; that lane's search missed it.
  * "Structured cancellation mod 1089 or mod 11 can close the 0.299 cube gap" (section 5).
* **NUMEROLOGY:** "33^2 = 1089" (section 6).

Session `collatz-procgen-20260922`, theta-beyond-phi lane (machine `mac-mini`, 2026-09-23). It follows the [cube-theta note](collatz_procgen_20260923_cube_theta.md) and the [hard-class note](collatz_procgen_20260922_hard_class.md).
* **No priority is claimed.** The method is Bézivin's (1998), refined by Choulet (2001) and by Krattenthaler–Rochev–Väänänen–Zudilin [KRVZ]. What is new here is the exchange of places: the 2-adic place becomes the analytic one and the archimedean place the arithmetic one. The rest is the application to Collatz swap words.
* Whether Bézivin or Choulet also state 2-adic versions is UNVERIFIED; their papers were not read.
* **Files.** No HYP or THM file was created.
  * Scripts: `04-computation/experiments/collatz_procgen_20260923_theta_{hankel,euler,defect,zero,1089}.py` and `..._theta_run.sh`.
  * Output: [`collatz_procgen_20260923_theta.out`](collatz_procgen_20260923_theta.out).

## 0. Answers in brief

1. **HYP-9131: yes, past `phi`, and far past it.**
   * Bézivin's Hankel determinant of the tails, transcribed to `Q_2`, proves `Theta(rho)` irrational for `mu_bar < 7/4` with a one-page elementary proof (H1). The 23/33 word (`1.61831`) is therefore settled.
   * Adding the proven `q`-order and cyclotomic divisibilities of [KRVZ] gives `28/11` (H2) and `2.87837` (H3). These settle every square-swap word under 5x+1 and under 7x+1.
   * **The gain does not come from the theta-null tie.** The method works for every Tschakaloff value. The tie does add forced cyclotomic factors beyond KRVZ's, but in range they grow more slowly than the proven part (section 5).
   * The exact missing estimates of the previous note (the square zero estimate, HC-CT4) are bypassed. The ultrametric Cauchy–Binet expansion makes the determinant nonzero automatically.
2. **HYP-9130 (cubes).**
   * Restricted results only:
     * Lemma F2 (parity from `F_2` rank);
     * Lemma E (an echelon family, `w -> 1`);
     * finite data for every tested `E` (`F_2` frontier about `2E`, witnesses of height `<= 2^92` at `E = 600`, `w = 1.9`).
   * The Hankel route that proves squares fails for cubes, for a proved reason (K7): cube tails are Gaussian sums, not exponential sums. Their Hankel determinants have cubic 2-adic valuation and a cubic forced divisor against a quartic clearing.
3. **"33^2 = 1089": NUMEROLOGY.**
   * The arithmetic facts are real: `3^5 = 2*11^2 + 1`, so 11 is a base-3 Wieferich prime.
   * No periodicity, identity or functional equation of `X` or of the 33-letter word was found modulo 1089, 11-adically, or through block lengths 33 or 1089. The forced divisors carry no factor 11.
   * The proof of HYP-9131 does not use 33.
4. **Pentagonal.**
   * Euler's function is `F_q(-1; 1)` in KRVZ's `q`-exponential family. Its Hankel tails are exponential sums, so `(rho;rho)_inf` is irrational for `mu_bar < 3/2` (elementary), 2 and 2.23719. That includes `rho = 2^10/3^9`.
   * The genuine "pentagonal theorem": every signed bilateral pentagonal swap word under `3x+r` has an irrational Bernstein number.
5. **The owner's principle.**
   * **Cubes.** The defect of the Liouville bookkeeping is zero: the heights of partial sums are exact up to `log2(K+1)` bits. The large gcds of periodic approximants are pure non-canonicity. So no cancellation closes the 0.299 cube gap.
   * **Squares.** The same principle is exactly what drives H2 and H3: the calculable defect "a nonzero integer is at least its 2-part times its forced odd divisor" is the lever. It moves the height exponent per unit of 2-adic valuation along

     `1.4265 -> 0.951 -> 0.882 -> 0.815 -> 0.560 -> 0.496`

     (periodic approximants, Zudilin with `m = 0`, Zudilin at the optimal `m`, H1, H2, H3; section 5).

## 1. Theorem H1: the 2-adic Hankel theorem (elementary)

Let `L >= 1` and `M_0 > 2^L` be odd. Put `rho = 2^L/M_0`, `mu_bar = log2(M_0)/L`, `Theta = sum_(k>=0) rho^(k^2)` in `Z_2`, and `S_m = sum_(k<=m) rho^(k^2)`. The case `M_0 < 2^L` is Liouville-easy, and the same proof works there with `2^L` in place of `M_0`.

**Theorem H1.** If `mu_bar < 7/4`, then `Theta` is irrational. More precisely, if `Theta = a/b` with `b` odd, then for every `n >= 1`

`L n(2n-1)(7n-1)/6 <= log2(n!) + n log2(|a| + (2n-1)|b|) + K_n log2 M_0`, where `K_n = 4 sum_(i<n) i^2`.

*Proof.*
1. **Tails are exponential sums.**

   `t_s := rho^(-s^2)(Theta - S_s) = sum_(m>=1) rho^(m^2) (rho^(2m))^s`.
2. **Cauchy–Binet with an ultrametric minimum.** Let `T_n = det(t_(i+j))_(0<=i,j<n)`. Writing `(t_(i+j)) = Y C Y^T` with `Y_(i,m) = rho^(2mi)` and `C = diag(rho^(m^2))`, we get

   `T_n = sum_(1<=m_1<...<m_n) prod_k rho^(m_k^2) prod_(k<l) (rho^(2m_l) - rho^(2m_k))^2`

   (a 2-adically convergent sum).
   * The term for `M = {m_k}` has exact valuation `L[sum_k m_k^2 + 4 sum_(k<l) m_k]`, because `v_2(rho^(2m_l) - rho^(2m_k)) = 2L m_k`.
   * Since `m_k >= k`, the unique minimum is at `M = {1, ..., n}`.
   * Hence **`v_2(T_n) = L n(n+1)(2n-1)/2`**, and in particular `T_n != 0`.
3. **Integrality.** Suppose `Theta = a/b`. Then `b t_m = G_m/2^(L m^2)` with

   `G_m = a M_0^(m^2) - b sum_(k<=m) 2^(L k^2) M_0^(m^2-k^2)`, an integer.

   Scale row `i` and column `j` by `2^(2L i^2)` and `2^(2L j^2)`. Since `2i^2 + 2j^2 - (i+j)^2 = (i-j)^2`, this gives

   `N_n := det(G_(i+j) 2^(L(i-j)^2)) = 2^(L K_n) b^n T_n`, an integer.

   So `N_n != 0` and `v_2(N_n) = L K_n + v_2(T_n) = L n(2n-1)(7n-1)/6`.
4. **Archimedean bound.**
   * `|G_m| <= (|a| + (m+1)|b|) M_0^(m^2)`.
   * So entry `(i,j)` is at most `C M_0^((i+j)^2 + (i-j)^2) = C M_0^(2i^2+2j^2)`, with `C = |a| + (2n-1)|b|`.
   * Pulling the row and column factors out gives `|N_n| <= n! C^n M_0^(K_n)`.
5. **Conclusion.** `2^(v_2(N_n)) <= |N_n|` gives the displayed inequality. Its leading terms are `(7/3) L n^3` against `(4/3) n^3 log2 M_0`. So it fails for large `n` when `7L > 4 log2 M_0`. ∎

**Checks** (`.out` K1, K2).
* `v_2(T_n)` equals the formula exactly for `n <= 12` on 8 maps, `mu_bar` from 1.43 to 2.85.
* `N_n = 2^(L K_n) b^n T_n(a/b)` holds exactly for random fake rationals.
* The archimedean bound is tight to about 1 bit.

**HYP-9131** (`L = 33`, `M_0 = 5^23`, `mu_bar = 1.61831`).
* The margin is `A1(n) = 5.794 n^3 + 57.3 n^2 - O(n log n)`, positive for all `n`.
* A rational of height `<= 2^h` is excluded at `n` = 2 (`h = 100`), 37 (`h = 10^4`), 416 (`h = 10^6`) and 13305 (`h = 10^9`).
* **HYP-9131 is PROVED.**

## 2. Refinements H2 and H3: the forced odd divisors ([KRVZ] transcribed)

**Dictionary.**
* **KRVZ's setup.** They study `F_q(z; lambda) = sum_n z^n / prod_(j<=n)(q^j - lambda)` with `|q| > 1`, the tails `v_n` (their (2.4)), and the Hankel determinant `V_n = det(v_(i+j))`, a polynomial in `q`, `alpha` and `mu` with integer coefficients.
* **Our specialization.** Put `xi = 1/rho = M_0/2^L`, `q = xi^2`, `alpha = xi`, `lambda = 0`. Then `T_q(alpha) = Theta` and `v_m = xi^(m(m+1)) (mu - S_m)`, as a polynomial identity.
* **Identity (checked exactly).** `b^n V_n = xi^(n(n-1)) det F(xi)`, with `F_ij(xi) = a xi^((i+j)^2) - b sum_(k<=i+j) xi^((i+j)^2 - k^2)` in `Z[xi, a, b]`, and `N_n = Y^(K_n) det F(X/Y)` with `X = M_0`, `Y = 2^L`.
* **Transfer of divisibility.** A divisibility of `V_n` by a polynomial in `q` specializes to a divisibility of `det F`; this uses Gauss's lemma, since cyclotomic polynomials are monic. Homogenizing then turns it into divisibility of the integer `N_n` by odd integers.

**Lemma D (forced divisors; PROVED given [KRVZ]).**
* (i) [KRVZ, Prop. 2 (`lambda = 0`)]: `V_n` lies in `q^(e_0(n)) Z[q, alpha, mu]`, with `e_0 = n(n-2)(5n-2)/24` for even `n` and `n(n-1)(5n-7)/24` for odd `n`. Hence `M_0^(e'(n))` divides `N_n`, with `e'(n) = max(0, 2e_0(n) - n(n-1)) ~ (5/12) n^3`.
* (ii) [KRVZ, Prop. 4 (any polynomial `b`)]: `Phi_l(q)^(e_l(n))` divides `V_n` for `1 <= l < n/2`, with `e_l(n) = sum_(i<n)(floor((i+l)/3l) + floor(i/3l))`. The polynomials `q`, `Phi_1(q)`, `Phi_2(q)`, ... are pairwise coprime, so their product divides too. Hence `M_0^(e'(n)) prod_l Phi_l^hom(M_0^2, 2^(2L))^(e_l(n))` divides `N_n`. These are odd integers, since `M_0` is odd.
* (iii) **Size of the cyclotomic values.** With `t = 2^(2L)/M_0^2 < 1`:

  `log|Phi_l^hom(M_0^2, 2^(2L))| = 2 phi(l) log M_0 + log|Phi_l(t)|`, and `|log|Phi_l(t)|| <= sum_(d>=1) |log(1 - t^d)| < infinity`.

  Since `e_l(n) <= n^2/(3l) + n`, the correction totals `O(n^2 log n)`.
* (iv) By [KRVZ, Lemma 3 and (4.17)], `sum_l phi(l) e_l(n) = C_cyc n^3 + O(n^2 log^2 n)` with

  `C_cyc = 5/54 - Im Li_2(e^(2 pi i/3))/(pi^2 sqrt3) = 0.05301135`.

**Theorems H2 and H3.** Combine `|N_n| >= 2^(v_2(N_n)) * (forced odd divisor)` with the bound of step 4.
* **H2:** `(7/3) L + (5/12) log2 M_0 > (4/3) log2 M_0`, i.e. `mu_bar < 28/11`.
* **H3:** `(7/3) L > (11/12 - 2 C_cyc) log2 M_0`, i.e. `mu_bar < 252 pi^2/(79 pi^2 + 72 sqrt3 Im Li_2(e^(2 pi i/3))) = 2.87837`.
* These are exactly `gamma/(gamma - 1)` for Choulet's `gamma = 28/17` and KRVZ's `gamma = 1.53237645`. The real-case constants map one-to-one through `mu_bar = gamma/(gamma - 1)`; see the K0 table.

**Checks** (`.out` K3, K4).
* **`X`-order.** With `X = 10007` prime, the `X`-order of `N_n` is `e(n) = 0, 0, 2, 9, 23, 47, ...` and always at least `e'(n)`.
* **Why `e(n)` exceeds `e'(n)`.** The difference is the `alpha`-degree term of KRVZ's leading monomial, which is `O(n^2)`. For `n <= 12` the tied `xi`-order equals exactly `2 e_0(n) + deg_alpha - n(n-1)`.
* **Divisibility.** `M_0^(e') prod Phi_l^hom(...)^(e_l)` divides `N_n` for random `(a, b)` on four maps, `n <= 10`.
* **Leading `n^3` coefficients of the margin** (positive means proved at that level):

  | map | `mu_bar` | H1 | H2 | H3 |
  |---|---|---|---|---|
  | Y (3x+1) | 1.426 | +4.31 | +10.26 | +11.77 |
  | HYP-9131 | 1.618 | +5.79 | +28.05 | +33.71 |
  | 5x+1 with `A = L` | 2.322 | –0.76 | +0.21 | +0.45 |
  | 7x+1 with `A = L` | 2.807 | –1.41 | –0.24 | +0.058 |

  For 7x+1 with `A = L`, the exact margin `A3` is positive from `n = 27` on (+1170 at `n = 40`). For 5x+1 with `A = L`, `A2` is positive from `n = 7` on.

## 3. Consequences

**Square-swap words.**
* Setup: blocks `B != B'` of length `L` with the same number `A` of 1s, under `x/2, (mx+r)/2`. By the block identity, rationality of the Bernstein number is equivalent to rationality of `Theta(rho)`, with `rho = 2^L/m^A`.
* Coverage (`.out` K5):

  | map | H1 covers `A/L <` | H2 | H3 | all blocks covered at level |
  |---|---|---|---|---|
  | 3x+r (and Mahler's map) | 1 | 1 | 1 | **H1** |
  | 5x+1 | 0.7537 | 1 | 1 | **H2** |
  | 7x+1 | 0.6234 | 0.9067 | 1 | **H3** |
  | 9x+1 | 0.5521 | 0.8030 | 0.9080 | partial |
  | 11x+1 | 0.5059 | 0.7358 | 0.8320 | partial |

* Nearest open square instances, `mu_bar` just above 2.87837: 9x+1 with 10 ones in 11 letters (2.88175), 11x+1 with 5 in 6 (2.88286), 9x+1 with 31 in 34 (2.89023).

**Theorem Q+ (quadratic families).**
* **Setting.** `P(m) = a m(m-1)/2 + s m`, normalized so that `a, s >= 1`, and `mu = sum_k eps^k rho^(P(k))`.
* **Tails.** `t_n = sum_(m>=1) eps^m rho^(P(m)) (rho^(am))^n` is again an exponential sum.
* **Valuation.** `v_2(T_n) = L[sum_(m<=n) P(m) + a n(n+1)(n-1)/3]`, verified for squares, triangular numbers, both pentagonal branches and `2k^2+k` (`.out` K6).
* **Clearing.** Use `R_i = a i^2 + s i`. Then `R_i + R_j - P(i+j) = (a/2)((i-j)^2 + i + j)`, which is a nonnegative integer.
* **Conclusion.** The same count gives `(7/6) a L > (2/3) a log2 M_0`, so `mu_bar < 7/4` for every family.
* **H2 and H3.** These transfer through `q = xi^a` and `alpha = eps xi^(a-s)`. The listed families have `s <= a`, so `alpha` is a polynomial in `xi`. This is a sketch: only the `v_2` laws are checked.

## 4. Euler's function and pentagonal words (Theorem E)

**Identity.**

`(rho;rho)_inf = sum_(k>=0) (-1)^k rho^(k(k+1)/2)/(rho;rho)_k = F_q(-1; 1)`, with `q = 1/rho`.

This is KRVZ's `q`-exponential case, `lambda = 1` and `alpha = -1`. KRVZ's hypotheses hold: `lambda` is not in `q^(Z>0)`, and `alpha` is not in `-lambda q^(Z>0)`. The product, Euler's `q`-series and the pentagonal sum agree modulo `2^4000` (`.out` E0).

**Theorem E.** `(rho;rho)_inf` is irrational in `Q_2` if `mu_bar < 3/2` (E1, elementary), `< 2` (E2) or `< 1/(1/2 - C_cyc) = 2.23719` (E3).

*Proof.*
1. **Tails.** The tails `v_n = (-1)^n sum_(m>=1) (-1)^m rho^(mn + m(m+1)/2) prod_(j=n+1)^(n+m)(1 - rho^j)^(-1)` expand as an exponential sum `(-1)^n sum_(f>=1) C_f (rho^f)^n`.
2. **Coefficients.** `C_f = sum_(m<=f) (-1)^m rho^(m(m+1)/2) h_(f-m)(rho, ..., rho^m)`, where `h` is the complete homogeneous symmetric polynomial. It has exact valuation `L f`, attained only at `m = 1`.
3. **Valuation.** Cauchy–Binet with the unique minimum at `{1, ..., n}` gives `v_2(V_n) = L n(n+1)(2n+1)/6`.
4. **Integrality.** `V_n` is a polynomial in `q` of degree `<= D_n = n(n-1)(4n+1)/6` ([KRVZ] (2.6)). So `N_n = 2^(L D_n) b^n V_n` is an integer, and `v_2(N_n) = L n^3` exactly.
5. **Archimedean bound.** `|N_n| <= n! (|a| + 2|b|)^n 2^(n(n-1)) M_0^(D_n)`.
6. **Conclusion.** Comparing `L n^3` with `(2/3) n^3 log2 M_0` gives E1.
7. **Refinements.** E2 adds `M_0^(e_0(n))` with `e_0 = n(n-1)(n-2)/6` ([KRVZ] Prop. 1, `lambda != 0`). E3 adds `prod Phi_l^hom(M_0, 2^L)^(e_l(n))` (Prop. 4 with `b(x) = x - 1`). ∎

These are the `q`-exponential analogues of Stihl's `gamma > 3`, Choulet's `gamma > 2` and KRVZ's `gamma > 1.80828`.

**Checks** (`.out` E1–E4).
* The valuation law holds for `n <= 12` on 4 maps.
* `N_n` is an integer and the divisibilities hold.
* At `rho = 2^10/3^9` the leading `n^3` coefficients are +0.49 (E1), +2.87 (E2) and +3.62 (E3). The elementary level excludes heights `2^h` at `n` = 10, 137 and 1432 for `h` = 100, `10^4` and `10^6`.
* At `rho = 2/3` (`mu_bar = 1.585`) E2 is needed.
* At `rho = 2/5` (`mu_bar = 2.32`) the case is open.

**Pentagonal swap words.**
* **Setup.** Put `B'` at the generalized pentagonal numbers with even `k`, `B''` at those with odd `k`, and `B` elsewhere, with `R_(B'') = 2 R_B - R_(B')`. The block identity then gives Bernstein number `= -(1/M_0)[R_B/(1-rho) + (R_B' - R_B)(rho;rho)_inf]`.
* **Existence.** Such triples exist for (L, A) = (9, 6), (10, 7), (11, 7), (11, 8), (12, 8), ... (`.out` E5). An example is `B' = 111110100`, `B = 111011100`, `B'' = 111101001`.
* **Coverage.** All of these are covered by E1. Under `3x+r` every such word is covered by E2, since `mu_bar <= log2 3 < 2`.
* **The case (10, 9).** It admits no triple, so no swap word has exactly `rho = 2^10/3^9`. The number `(rho;rho)_inf` itself is still irrational by E1.

**Unsigned bilateral pentagonal sum.**
* **Method.** The two families give tails with frequencies `rho^(3m)` and `rho^(3m+1)` and coefficients `rho^(P_1(m))` and `rho^(P_2(m))`, where `P_1(m) = m(3m-1)/2` and `P_2(m) = m(3m+1)/2`.
* **Result.** `v_2(T_n) = L[sum g_i + 2 sum (n-i) e_i] ~ (5/8) L n^3`, verified. Against the clearing `~ 2 n^3 log2 M_0` this gives **irrationality for `mu_bar < 21/16`** (elementary; covers `3x+1` words with `A/L < 0.828`).
* **Proof of the clearing.**
  * Put `S_s = sum_(|k|<=s)`. Then `b t_s = +-G_s/(2^(L P_1(s)) M_0^s)` with `G_s` an integer, `|G_s| <= C M_0^(P_2(s))`.
  * Scale row `i` by `2^(L R_i) M_0^i` and column `j` by `2^(L R_j) M_0^j`, with `R_i = 3i^2`. Then `R_i + R_j - P_1(i+j) = (3(i-j)^2 + i + j)/2` is a nonnegative integer, and every entry is an integer bounded by `C M_0^(R_i + R_j + i + j)`.
  * Hence `N_n = 2^(2L sum R_i) M_0^(n(n-1)) b^n T_n` is a nonzero integer with `v_2(N_n) = 6L sum i^2 + v_2(T_n)` and `|N_n| <= n! C^n M_0^(n(n-1) + 6 sum i^2)`.
  * Leading terms: `(2 + 5/8) L` against `2 log2 M_0`. ∎
* **Why it is weaker.** The unsigned sum is a product of three `q`-Pochhammer symbols, with no single `q`-exponential representation. Even for the signed sum, the `q`-exponential route (3/2) beats the two-family one (21/16).

## 5. The owner's principle: calculable defects as levers (`.out` D1–D6, K7)

* **(D1) Cube partial sums.**
  * `gcd(A_K, 3^(9K^3)) = 1`, because `A_K = 2^(10K^3) mod 3`.
  * The exact height is `3^(9K^3)(1.052...)`. The triangle bound `(K+1) 3^(9K^3)` loses only 0.9 to 3.1 bits for `K <= 12`.
  * The ratio `v_2/log2 H` tends to `1/mu_bar = 0.7010` with **zero calculable defect**.
* **(D2) Periodic approximants of Y3.**
  * The best approximants per period end have `gcd(P,Q)` up to `2^282`.
  * **All** of this comes from non-canonical representations: a non-primitive period (factor `(2^|v| - M_v)/(2^|w| - M_w)`) or a non-minimal preperiod (factor `M_z`).
  * After canonicalization the residual gcd is at most 28.6 bits, usually 0, and the gains stay bounded. These approximants are the partial sums in disguise.
* **(D3) Square forms and the exponent ladder.**
  * Zudilin's forms carry a calculable cyclotomic gcd, `18659^(n/2) * (...)`, with `18659 = 3^9 - 2^10`. It is about 10% of the height in range (399 bits at `n = 12`).
  * The height exponent per unit of 2-adic valuation, at `mu_bar = 1.4265`:

    | construction | exponent |
    |---|---|
    | periodic approximants | 1.4265 |
    | Zudilin, `m = 0` | 0.951 |
    | Zudilin, `m = n/phi` | 0.882 |
    | Hankel H1 | 0.815 |
    | H2 | 0.560 |
    | H3 | 0.496 |

  * Each step exploits a calculable defect: Padé cancellation, Vandermonde rank, then forced odd divisors.
* **(D4, D6) The exact forced divisor, squares.** For `n <= 14`, the gcd over `(a, b)` of `det F(xi)`, as a polynomial in `xi`, is `xi^(e(n))` times cyclotomic `Phi_d(xi)` powers only.
  * The tied (theta-null) case has **more** forced factors than KRVZ prove: `Phi_1(xi)`, `Phi_3(xi)`, `Phi_4(xi)`, ... with larger exponents. At `n = 9` there are 11 extra copies of `Phi_1(xi)`, i.e. `18659^11`.
  * The excess is `n^2`-sized (`excess/n^2` goes from 0.63 to 1.15 for `n = 4..14`, while `excess/n^3` decreases). The proven part is `n^3`-sized. So there is **no extra cubic lever in range**, and the theta-null tie does not move the threshold numerically.
* **(K7, D5, D6) Cubes.**
  * `v_2(T_n^cube) = L(3n^3 - 3n^2 + n)` exactly. **Proof:** the minimal term is the anti-diagonal permutation with every entry at level `m = 1`, unique by strict convexity of `s -> 3s^2 + 3s + 1`, and higher levels add at least `7L`.
  * The clearing is quartic, `K3 = 8 sum i^3 ~ 2n^4`.
  * The forced divisor is `xi^(3n^3 - 18n^2 + 34n - 19)` (law checked for `3 <= n <= 12`) times small cyclotomic factors, cubic in total.
  * Closing the 0.299 gap would need a forced divisor of about `M_0^(0.6 n^4)`. The `11`- and `1089`-parts of every forced divisor are trivial (`v_11 = 0` throughout).
  * **Verdict: no structured cancellation, mod 1089, mod 11 or otherwise, closes the cube gap for determinants of cube tails.**

## 6. "33^2 = 1089": the test (`.out` N1–N4)

* **Real facts.**
  * `1089 = 3^2 11^2`, `33 = 2^5 + 1`, `3^5 = 2*11^2 + 1` (so `3^10 = 1 mod 121`: 11 is a base-3 Wieferich prime), `3^9 = 18*1089 + 81`, `2^10 = 1089 - 65`, `3^7 - 2^11 = 139 = 2^7 + 11`.
  * `rho = 47 mod 121`, of order 55. Separately, `3^9 - 2^10 = 47*397`; the shared 47 is a coincidence.
  * `rho^5 = 2^50 mod 121` (because `3^45 = 1 mod 121`).
* **Tests.**
  * The terms `rho^(k^3) mod 121` are periodic with period 55. The partial sums drift by 0 per period, but that is generic: cubing permutes `Z/55`, so a full period sums to `(rho^55 - 1)/(rho - 1) = 0 mod 121`. The series does not converge 11-adically, and nothing 11-adic constrains a putative rational `a/b`.
  * **Real value `X_R`.** 2700 reliable partial quotients (max 3929). 41 of them are `>= 100`, against 38.7 expected. Lévy estimate 1.1742 against 1.1866. PSLQ finds no relation of degree `<= 6` with coefficients `<= 10^40`.
  * **2-adic `X`.** For `d = 1089^j, 33^j, 11^j, 121^j, 3^j` (`j < 60`), `X != c/d` for all `|c| < 2^99988`.
  * No forced divisor in section 5 has an 11-part.
  * The HYP-9131 proof (H1) needs only `mu_bar < 7/4`, so 33 enters only through `mu_bar`.
* **Label: NUMEROLOGY.** The facts are real, but no mechanism links them to the 2-adic irrationality questions.

## 7. HYP-9130 in restricted form (`.out` Z1–Z4)

**Setting.** `Lambda_(E,W) = {A in Z^(E+1) : c_j . A = 0 for E < j < W}`, with `c_j(i) = 1_cube(j - i)`.

**Lemma F2 (PROVED).** If `c_(E+1), ..., c_W` are linearly independent over `F_2`, then some `A` in `Lambda_(E,W)` has `r_W(A) = c_W . A` odd. That is, its remainder has order exactly `W` with odd leading coefficient.

*Proof.*
1. `F_2`-independence gives a unit maximal minor over `Z_2`. So `A -> (c_j . A)_(E<j<=W)` is surjective from `Z_2^(E+1)`.
2. `Lambda_(E,W) (x) Z_2` is the `Z_2`-kernel of the first `W - E - 1` rows, and `r_W` maps it onto `Z_2`.
3. Hence `r_W(Lambda_(E,W)) = dZ` with `d` odd. ∎

**Lemma E (PROVED).** For `E = N^3 - 1` with `N >= 4`, the rows `c_j` for `N^3 <= j <= (N+1)^3 - 1` have distinct minimal columns `j - N^3`. Hence they are independent. This gives the parity half for `W <= (N+1)^3 - 1`, but only with `w -> 1`.

**FINITE-EXACT.**
* **`F_2` frontier.** `W_2(E) >= 1.4265 E` for all 446 tested `E` in `[10, 1500]`. For `E >= 100`, `W_2(E) >= 1.938 E` and `2E - W_2(E)` lies in `[-1, 9]`.
* **Berlekamp–Massey profile** of the cube indicator to `N = 40000`: the maximum deviation of `L_n` from `n/2` is 8, and the increments (partial-quotient degrees over `F_2`) have frequencies 0.499, 0.255, 0.123, 0.057, exactly those of random sequences.
* **Witnesses.** Every tested `E <= 600` and `w` in `{1.5, 1.75, 1.9}` has a witness with `v_2(lead) < 10`. The largest height is `2^19.3` at `w = 1.75` and `2^92.2` at `w = 1.9` (`E = 600`, consistent with Siegel's `E^((w-1)/(2-w))`). All margins are positive.

**Exact missing estimate for HYP-9130.** Both parts are purely combinatorial.
* **(a) Parity.** The `F_2` partial quotients of `sum x^(k^3)` should have degree `o(n)`; the data say "bounded by about 16".
* **(b) Heights.** `lambda_d(Lambda_(E, wE)) = 2^(o(E))`; the data say the growth is polynomial in `E`.

**The Hankel route does not help for cubes (K7).**

## 8. Hypothesis candidates (for `INDEX.md`; no HYP files created)

* **HYP-9131 → PROVED** (Theorem H1). Suggested new status line: "PROVED (collatz_procgen_20260923_theta_beyond_phi, Theorem H1; elementary): `sum (2^33/5^23)^(k^2)` is irrational in `Q_2`; every 5x+1 square-swap word likewise (H2)."
* **HC-TH1 (the next square instance; OPEN).** `sum (2^11/9^10)^(k^2)` is irrational in `Q_2`. This is 9x+1 with 10 ones in 11 letters, `mu_bar = 2.88175`, just past H3.
* **HC-TH2 (theta non-quadraticity; OPEN, dictionary).** KRVZ Theorem 1 (`gamma > 3.27694`) transcribes to non-quadraticity for `mu_bar < 3.27694/2.27694 = 1.4392`. That would cover `rho = 2^10/3^9` (1.42647) with a margin of 0.013, and would be partial progress on HYP-9132. The 2-adic conjugate bookkeeping (KRVZ Lemma 5) was **not** done.
* **HC-TH3 (unsigned pentagonal; OPEN).** `sum_(k in Z) rho^(k(3k-1)/2)` is irrational for `mu_bar` in `[21/16, 2)`. A two-family `q`-order or cyclotomic analysis is needed.
* **HC-TH4 (tied forced divisor; OPEN, FINITE-EXACT).** The exact forced divisor of the tied Hankel determinants is `xi^(e(n)) prod_d Phi_d(xi)^(E_d(n))`, with exponents `E_d(n)` exceeding KRVZ's by `O(n^2 log n)` in total. A proof would give the sharp 2-adic constant. The data say this constant equals H3's.
* **HC-TH5 (cube determinant; OPEN).** Any `n x n` determinant of linear forms in `(1, X)`, for `X = sum rho^(k^3)`, whose forced odd divisor grows quartically (like `M_0^(c n^4)`) with `c > (1 - 1/mu_bar) * 2` would prove HYP-9127. Hankel matrices of tails achieve only cubic growth (K7).
* **HC-TH6 (the `F_2` shadow of HYP-9130; OPEN).** The partial quotients of `sum x^(k^3)` in `F_2((x))` have bounded degree. The data support it to degree 20000.

## 9. Reproduction

```bash
cd <worktree>
bash 04-computation/experiments/collatz_procgen_20260923_theta_run.sh          # full run, about 30 s, peak RSS about 100 MB
#   writes 05-knowledge/results/collatz_procgen_20260923_theta.out
bash 04-computation/experiments/collatz_procgen_20260923_theta_run.sh --quick  # smoke run, writes scratch/procgen_theta/theta_quick.out
```

Requirements: python3 with gmpy2, python-flint, sympy and mpmath.

**Sections of the `.out`:**
* **Hankel.**
  * K0: thresholds and the dictionary.
  * K1: `v_2` laws.
  * K2: integer identity and archimedean bound.
  * K3: divisibility.
  * K4: margins.
  * K5: coverage.
  * K6: quadratic families.
  * K7: cubes.
* **Euler.** E0–E6.
* **Defects.** D1–D6.
* **HYP-9130.** Z1–Z4.
* **1089.** N1–N4.

**Scratch exploration:** `scratch/procgen_theta/{hankel_check,xorder,krvz_identity,cube_xorder*,content_poly*}.py`.

## Sources

`[P]` means the primary text was read, `[R]` that a named secondary source was used.

* **[KRVZ]** C. Krattenthaler, I. Rochev, K. Väänänen, W. Zudilin, *On the non-quadraticity of values of the q-exponential function and related q-series*, Acta Arith. 136 (2009) 243–269, arXiv 0812.2921v1. [P: Introduction (the history `gamma > (3+sqrt5)/2` → Bézivin `28/15` → Choulet `28/17`), Theorems 1 and 2, section 2 (Bézivin's construction, (2.1)–(2.6)), Lemma 1, Propositions 1 and 2 with proofs, Proposition 4 statement and (4.1), Lemma 3, (4.17), Proposition 5, references.]
* J.-P. Bézivin, *Sur les propriétés arithmétiques d'une fonction entière*, Math. Nachr. 190 (1998) 31–42. [R via KRVZ]
* R. Choulet, *Des résultats d'irrationalité pour deux fonctions particulières*, Collect. Math. 52 (2001) 1–20. [R via KRVZ]
* Th. Stihl, Arch. Math. 41 (1983) 531–537. [R via KRVZ]
* L. Tschakaloff, Math. Ann. 80 (1919), 84 (1921); P. Bundschuh, Invent. Math. 9 (1970) and later papers. [R via KRVZ]
* W. Zudilin, arXiv math/0506086. [P, previous note]
* The [cube-theta note](collatz_procgen_20260923_cube_theta.md) (Proposition F, Corollary F3, Lemma P recap) and the [hard-class note](collatz_procgen_20260922_hard_class.md) (block identity, Theorem Y).
* Hypotheses HYP-9127, HYP-9130, HYP-9131 and HYP-9132 in `05-knowledge/hypotheses/`.
