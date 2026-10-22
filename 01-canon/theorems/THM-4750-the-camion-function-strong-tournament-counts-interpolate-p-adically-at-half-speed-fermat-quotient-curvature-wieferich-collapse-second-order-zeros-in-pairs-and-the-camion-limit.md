---
id: THM-4750
title: "The Camion function. Let p be an odd prime, eps = (2/p), s(n) the strong labelled tournaments (A054946), D = s - c with c(n) the circulant count, and A_m(l) = eps^(-l) D(m + l(p-1)). (1) Half-speed higher Kummer congruences: for k >= 1 and m >= 2k-1 the Mahler coefficients a_i(m) = Delta^i A_m(0) satisfy a_i(m) = 0 mod p^min(k, ceil(i/2)); in particular sum_j C(2k-1,j)(-eps)^(2k-1-j) D(m+j(p-1)) = 0 mod p^k. The exponent ceil(i/2), half of the Bernoulli/Kummer rate, comes from the edge count C(n,2) having degree 2. (2) At a Wieferich prime the speed is full: a_i(m) = 0 mod p^min(k,i) for m >= k. (3) Fermat-quotient curvature: for m >= 3, a_2(m)/p = -q_2(p) X(m-2) mod p, where q_2 is the Fermat quotient of 2 and X has EGF T''/T^2 = S'' + 2 S'^2 T; so mod p^2 every class sequence is a quadratic polynomial in l whose curvature vanishes identically exactly at Wieferich primes (for X(m-2) not 0 mod p). (4) Second-order zeros come in pairs: in a zero class the l mod p with p^2 | D(m + l(p-1)) are the roots of an explicit quadratic over F_p; for p >= 5, Weierstrass preparation gives exactly two p-adic zeros (split in Z_p or an inert conjugate pair) unless q_2(p)X(m-2) = 0 mod p; the triangle class m = 3 is always split. (5) The Camion limit Lambda_p = lim eps^a D(p^a) exists, eps^a D(p^a) = Lambda_p mod p^a, Lambda_p is the value at n = 0 of the interpolation along n = 1 mod (p-1), and the Camion tower height is h(p) = v_p(Lambda_p); Lambda_p mod p^K is computable from s(n), n <= (2K-1)p. (6) The prime 2 (Moon: |Aut| odd): f(n) = s(n)/n! = sum over strong classes of 1/|Aut| lies in Z_2, the tournament EGF is Phi_3(x) = 1 + x + x^2 mod 2, so the number of strong tournaments on n vertices up to isomorphism is odd iff n is not 2 mod 3 (and the number of all tournaments is even for n >= 3)."
status: >
  PROVED: statements 1, 2, 3, 4 (structure; the p-adic zero count needs p >= 5), 5, 6(a), 7 (r-uniform speed 1/r).
  Statements 1-3 are proved from the Session25 filtered-Hurwitz representation (credited) by a Leibniz estimate of
  Mahler coefficients; 4 by reduction mod p^2 and Weierstrass preparation; 5 from 1 (Lipschitz continuity) and the
  THM-4746 (9) freezing lemma; 6(a) from Moon's odd-order theorem and T = Phi_3 mod 2.
  FINITE-EXACT:
  - sharpness of ceil(i/2) for p in {11, 13, 17, 1091, 1097, 3517} at K = 3, 4, attained by all but O(1) classes;
  - extra divisibility exactly where p | (2k-1)!! (p = 5, 7);
  - full speed at the Wieferich primes 1093 (K = 4) and 3511 (K = 3);
  - the curvature formula (3) in every residue class of every prime p <= 10^4 (0 failures);
  - the quadratic law mod p^2 directly for p <= 59, all classes, l <= p + 1 (17158 checks);
  - second-order root statistics for 1226 primes 7 <= p <= 10^4 (4906 zero classes): classes 1, 2 and the sporadic classes split about 1:1 into 0 and 2 roots; class 3 always has 2; at 1093 and 3511 every zero class has exactly one root;
  - the Camion-Wieferich primes p <= 10^5 are exactly 3, 13, 22409 (two independent codes), and h(22409) = 2 by statement 5;
  - the 2-adic quasi-period f(n + 3*2^j) = f(n) mod 2^(j+2) (n <= 600; failures only at n <= 3 for j = 5..7);
  - r = 3 root statistics of connected 3-uniform hypergraphs (60 primes <= 2200): random-cubic proportions, forced root in class 2;
  - the simple 2-adic zero n* of f on n = 2 mod 3, v_2(f(n)) = 2 + v_2(n - n*), with n* = 16221783116972227762 mod 2^64 (n* = 3250 mod 4096), and not a rational with numerator and denominator below about 2^31.
  HEURISTIC: the reading of n* as a ratio of logarithms in the norm-one torus of Q_2(sqrt(-3)); the analogy of
  Lambda_p with the Fermat-side log_p<2> and with Kubota-Leopoldt values.
  PRE-AUDIT (mac-mini-2026-10-21-fanopsp, 2026-10-21).
  NOT CHECKED in the literature: higher Kummer congruences for these labelled counts, the curvature formula, and
  the 2-adic statements. Specker-Blatter, Mani-Stones and the Kummer-congruence literature are the nearest prior art.
session: mac-mini-2026-10-21-fanopsp (owner prompt: hunt h(p) >= 3, HYP-9357 at 27 and 31, repunit rarity, wide-net loops)
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md (section 13)
scripts:
  - 04-computation/experiments/camion_wieferich_20261021/higher_kummer.py (+ higher_kummer.out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_speed.c (+ camion_speed.out), camion_speed_base.c (+ camion_speed_base.out)
  - 04-computation/experiments/camion_wieferich_20261021/curvature_meaning.py (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_second_order.c (+ second_order_*.out)
  - 04-computation/experiments/camion_wieferich_20261021/tower_height.c (+ tower_22409_K3.out, tower_height_direct_checks.out)
  - 04-computation/experiments/camion_wieferich_20261021/cw_fast.c (+ cw_*.out, cw_fast_validation_5_10000.out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_irregular.c (+ camion_irregular_10000.out)
  - 04-computation/experiments/camion_wieferich_20261021/two_adic_face.py, two_adic_zero.py, two_adic_zero_fast.py (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/hyper3_second_order.c (+ hyper3_second_order_5_1400.out, hyper3_second_order_1401_2200.out)
depends_on:
  - 05-knowledge/results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md (the filtered algebra F and S = e - T^(-1) in F; PROVED, independently audited)
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md (statement 9, towers)
  - J. W. Moon, Topics on Tournaments (1968): the automorphism group of a tournament has odd order (CITED)
related:
  - 05-knowledge/hypotheses/HYP-9361-camion-function-open-laws-split-inert-equidistribution-tower-heights-uniformity-k-speed.md
  - 05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-wieferich-type-primes-sporadic-zero-classes-prime-squares.md
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md (Wieferich primes in circulant tournaments: a different appearance of q_2)
---

# THM-4750 — The Camion function

**One line.** Along each residue class mod `p − 1`, the strong tournament count is a `p`-adic analytic function of the class index. It converges at half the Kummer rate, because the edge count `C(n,2)` is quadratic. The half-speed part is the Fermat quotient of 2 times a tournament statistic. So the Wieferich primes 1093 and 3511 are exactly where the curvature disappears mod `p²`, and there every Camion zero class has one second-order root instead of 0 or 2.

## Notation

* `p` is an odd prime and `ε = (2/p)`.
* `s(n)` counts strong labelled tournaments (`s(1) = 1`).
* `c(n)` is the circulant count: `2^((n−1)/2)` for odd `n`, 0 for even `n`.
* `D = s − c`.
* `T(n) = 2^C(n,2)`, with EGFs `T(x)` and `S(x) = 1 − 1/T(x)`.
* For `m ≥ 1` the **class sequence** is `A_m(l) = ε^(−l)·D(m + l(p−1))`, `l ≥ 0`, with **Mahler coefficients** `a_i(m) = Δ^i A_m(0) = Σ_j (−1)^(i−j) C(i,j) A_m(j)`. Then `A_m(l) = Σ_i a_i(m) C(l,i)`.
* `q_2(p) = (2^(p−1) − 1)/p` is the Fermat quotient; `p` is **Wieferich** if `p | q_2(p)`.
* `X(n) = n!·[x^n] T''(x)/T(x)²`, with values `X = Δ²T ∗ T^(∗−2)` in the Hurwitz ring of binomial convolution.
  * Since `T = 1/(1−S)`, also `T''/T² = S'' + 2S'²T`.
  * Equivalently `C(n,2)·X(n−2) = 2·(d/da) s_a(n)|_(a=2)`, where `s_a` is the strong count built from `a^C(n,2)`: the derivative of the strong count in the base.
  * **Combinatorially (PROVED), `X(m−2)` is the number of labelled tournaments on `[m]` that become strongly connected when the arc at one fixed pair `{u,v}` is made two-way.**
    * *Proof.* Make `{u,v}` two-way and orient the other pairs. The condensation of such a digraph is transitive, and `u, v` lie in one strong component. So the `T(m)/2` such digraphs factor as tournament ∗ (strong block containing `u, v`) ∗ tournament. Hence the strong ones number `(Δ²T/2) ∗ T^(∗−2) = X/2`. Each is counted twice, once for each orientation of `uv` in the tournament.
    * Values: `2, 4, 36, 680, 25152, 1784416` for `m = 2..7`. Brute force agrees for `m ≤ 7` (`curvature_meaning.out`).

**The representation used** (Session25, PROVED and audited). Work over `R_0 = Z/p^k` and its unramified quadratic extension `R`.

* `θ` is a root of `θ^(p−1) = ε`, `G = θR_0`, and `ω` is the Teichmüller lift of 2 with `2 = ω(1 + pu)`. Then `u ≡ −q_2(p) (mod p)`, because `2^(p−1) = (1+pu)^(p−1) ≡ 1 − pu (mod p²)`.
* `F` is the `R`-module generated by `p^j δ_d ∗ E_β` with `0 ≤ j < k`, `0 ≤ d ≤ 2j` and `β ∈ G`. Here `E_β = (β^n)_n` and `(δ_d ∗ E_β)(n) = C(n,d) β^(n−d)`.
* Session25 proves that `F` is a ring and that `S = e − T^(−1) ∈ F`.
* The grading comes from the Newton expansion `(1+pu)^C(n,2) = Σ_j (pu)^j C(C(n,2), j)`: the `j`-th term has degree `2j` in `n`.

## 1. Half-speed higher Kummer congruences (PROVED)

**Statement 1.** For every `k ≥ 1`, every `m ≥ 2k − 1` and every `i ≥ 0`,

    a_i(m) ≡ 0   (mod p^min(k, ⌈i/2⌉)).

So for `m ≥ 2k − 1` the class sequence is, mod `p^k`, a polynomial in `l` of degree at most `2k − 2`, with `v_p` of the `C(l,i)` coefficient at least `⌈i/2⌉`. The case `i = 2k − 1` is the **higher Kummer congruence**

    Σ_(j=0)^(2k−1) C(2k−1, j) (−ε)^(2k−1−j) D(m + j(p−1)) ≡ 0   (mod p^k)   (m ≥ 2k − 1),

and the same holds with `D` replaced by `s`.

**Proof.** By linearity it suffices to treat one generator `g = p^j δ_d ∗ E_β` (`j < k`, `d ≤ 2j`) together with `c`. For `n ≥ m ≥ 2k − 1 > d`, `g(n) = p^j C(n,d) β^(n−d)`.

* **Nonunit `β`** (`β ∈ pR`, including 0). Then `v(g(n)) ≥ j + n − d ≥ n − j ≥ (2k−1) − (k−1) = k`, so `g` vanishes mod `p^k` on the whole class.
* **Unit `β = θb`.** Then `β^(p−1) = ε(1 + py)` with `y ∈ R_0`, and

      ε^(−l) g(m + l(p−1)) = p^j β^(m−d) · P(l) · Q(l),

  where `P(l) = C(m + l(p−1), d)` and `Q(l) = (1+py)^l`.
  * `P` is integer-valued of degree `d`, so `Δ^a P(0) ∈ Z` and `Δ^a P(0) = 0` for `a > d`.
  * `Δ^c Q(l) = (py)^c Q(l)`.
  * The discrete Leibniz rule `Δ^i(PQ)(0) = Σ_a C(i,a) Δ^a P(0) Δ^(i−a)Q(a)` gives `v(Δ^i(gε^(−l))(0)) ≥ j + max(0, i − d)`.
  * With `d ≤ 2j` this is at least `max(j, i − j) ≥ ⌈i/2⌉`.
* **Circulants.** For odd `m`, `ε^(−l) c(m + l(p−1)) = 2^((m−1)/2) (1 + pz)^l` with `2^((p−1)/2) = ε(1+pz)` (Euler). Its Mahler coefficients `2^((m−1)/2)(pz)^i` have valuation at least `i`. For even `m` the sequence is zero.

Since `S ∈ F`, the bound holds for `s`, and hence for `D`. The congruence is an integer congruence because `R_0 → R` is faithful. ∎

**What the exponent means.** For a sequence such as `n ↦ ε^(−l) a^n` (a single exponential, linear exponent) the Mahler coefficients have valuation at least `i`. That is the rate of the classical Kummer congruences for Bernoulli numbers, which underlies Kubota–Leopoldt interpolation. Here each power of `pu` arrives with two powers of `n`, one edge being two endpoints, so the rate halves. Hypergraphs built from `2^C(n,r)` should run at rate `1/r` (HYP-9361 (d)).

**The same proof applies to every sequence in Session25's closure of `F`.** Examples:

* strong semicomplete digraphs, with `ε = (3/p)` and `p ≠ 3`;
* connected labelled graphs `b(n)`, from `m ≥ 2k` (one shift: `Δb ∈ F`).

**FINITE-EXACT** (`higher_kummer.out`): the order-`(2k−1)` congruence for 13 pairs `(p,k)` with `k ≤ 5`. Failures occur only below `m = 2k − 1`.

## 2. Sharpness, small primes, and full speed at Wieferich primes

**Statement 2 (PROVED).** If `p` is Wieferich, then `a_i(m) ≡ 0 (mod p^min(k,i))` for every `m ≥ k`: full speed.

*Proof.* Now `pu = p²u'`, so the grade-`j` Newton term carries `p^(2j)` at degree `≤ 2j`.

* Let `F'` be generated by `p^J δ_d ∗ E_β` with `d ≤ J`. It is a ring (exponents and degrees add), and `T ∈ F'`.
* The inversion argument of Session25 goes through unchanged, so `S ∈ F'`.
* For a generator, the Leibniz estimate gives valuation at least `J + max(0, i − d) ≥ i`.
* Nonunit nodes contribute valuation at least `J + n − d ≥ n ≥ k`.
* The circulant part is full speed in any case. ∎

**General base (PROVED by the same argument).** Session25's representation holds for `T_a(n) = a^C(n,2)` with any `a` prime to `p`, with multiplier `(a/p)`. Statements 1–3 then hold with `q_2` replaced by the Fermat quotient `q_a(p) = (a^(p−1) − 1)/p` and `X` by `X_a = T_a''/T_a²`. So **full speed occurs exactly at the base-`a` Wieferich primes**.

* For strong semicomplete digraphs (`a = 3`, three states per pair) these are `p = 11` and `1006003`.
* `camion_speed_base.out`:
  * base 3, `p = 11` (`3^10 ≡ 1 mod 121`): min `v_p(a_i) = 0 1 2 3 4 …`, full speed;
  * base 3, `p = 13, 17`: half speed, attained;
  * base 7, `p = 5` (`7^4 ≡ 1 mod 25`): full speed;
  * base 5, `p = 7, 11`: half speed.

**FINITE-EXACT** (`camion_speed.out`). The table gives the minimum of `v_p(a_i(m))` over all `p − 1` classes, with `m ≥ 2K − 1`.

| `p` | `K` | `min v_p(a_i)`, `i = 0, 1, …` | bound `⌈i/2⌉` |
|---|---|---|---|
| 11, 13, 17 | 4 | `0 1 1 2 2 3 3 4` | equal |
| 1091, 1097, 3517 | 3 | `0 1 1 2 2 3` | equal |
| **1093** (Wieferich) | 4 | `0 1 2 3 4 4 4 4` (cap 4) | full speed |
| **3511** (Wieferich) | 3 | `0 1 2 3 3 3` (cap 3) | full speed |
| 5 | 6 | `0 1 1 2 2 4 4 5 5 6 6 6` | +1 from `i = 5` |
| 7 | 5 | `0 1 1 2 2 3 3 5 5 5` | +1 from `i = 7` |

* For non-Wieferich `p > 2k − 1`, the half-speed bound is attained by all but at most 6 classes in every tested case. At `p = 1091`, 1089 of the 1090 classes attain `v(a_2) = 1`; the exception is a class with `X(m−2) ≡ 0`.
* For small primes the extra divisibility appears exactly at `i ∈ {2k−1, 2k}` with `p | (2k−1)!!`. The top Newton coefficient of `C(C(n,2),k)` at degree `2k` is `(2k)!/(2^k k!) = (2k−1)!!`.
* This is the mechanism conjectured heuristically in HYP-9358 (e2); the full threshold statement there remains open.

## 3. The Fermat-quotient curvature (PROVED)

**Statement 3.** For every `m ≥ 3`,

    a_2(m)/p ≡ −q_2(p) · X(m − 2)   (mod p),

and mod `p²` the class sequence is exactly the quadratic

    ε^(−l) D(m + l(p−1)) ≡ a_0(m) + a_1(m)·l + a_2(m)·C(l,2)   (mod p²),

with `p | a_1(m)` and `p | a_2(m)`.

**Proof.** Take `k = 2` in §1. Only grades `j ≤ 1` survive mod `p²`.

* The grade-0 part is in the group algebra. Its second differences are `β^m (py)² ≡ 0 (mod p²)`.
* The grade-1 part of `T^(−1)` is `−A^(−1)U_1A^(−1)`, where `A` is the Teichmüller part and `U_1(n) = pu·C(n,2)·t_0(n) = pu·(δ_2 ∗ Δ²t_0)(n)`.
* So the grade-1 part of `S` is `pu·δ_2 ∗ X̃`, with `X̃ = Δ²t_0 ∗ t_0^(∗−2) ≡ X (mod p)`.
* On the class, `ε^(−l) X̃(m − 2 + l(p−1)) ≡ X(m−2) (mod p)` for `m − 2 ≥ 1`.
* `Δ_l² C(m + l(p−1), 2) = (p−1)² ≡ 1`.
* Hence `a_2 ≡ pu·X(m−2) ≡ −p q_2(p) X(m−2) (mod p²)`.
* The circulant part has second difference `2^((m−1)/2)(pz)² ≡ 0`.
* The vanishing of `a_i` mod `p²` for `i ≥ 3` is §1 with `k = 2`. ∎

**FINITE-EXACT.**
* The formula holds in every residue class `m ∈ [3, p+1]` of every prime `p ≤ 10^4`, with 0 failures (`second_order_*.out`).
* The quadratic law was checked directly against `D(m + l(p−1)) mod p²` for `l ≤ p + 1`, all classes, `p ≤ 59` (17158 checks, `camion_second_order.c` with `VMAX = 60`).

**Reading.** Mod `p²` the strong count is "Teichmüller part + Fermat drift × base-derivative". Here the base-derivative is the count of tournaments made strong by one two-way pair. The Teichmüller part `2 ↦ ω(2)` makes `T` a pure exponential sum, linear in `l`. The drift `⟨2⟩ = 1 + pu` enters through `⟨2⟩^C(n,2)`, and its first-order effect is the base-derivative `d s_a/da`. The curvature in `l` is the quadratic edge count. At Wieferich primes `2` is Teichmüller to second order, so the curvature vanishes.

## 4. Second-order zeros come in pairs (PROVED structure; FINITE-EXACT statistics)

Call `m ≥ 3` a **zero class** if `p | D(m)` (equivalently `p | D` on the whole class). This includes the universal classes `1, 2, 3` (THM-4746). For a zero class put `c_i = a_i(m)/p mod p` and

    g_m(l) = c_0 + c_1 l + c_2 C(l,2) ∈ F_p[l],   c_2 = −q_2(p) X(m−2).

**Statement 4.**

* **(a)** For `l ≥ 0`, `p² | D(m + l(p−1))` iff `g_m(l mod p) = 0`.
* **(b)** If `q_2(p)X(m−2) ≢ 0`, then `g_m` is a quadratic. It has `1 + χ(Δ)` roots in `F_p`, where `Δ = (c_1 − c_2/2)² − 2c_0c_2` and `χ` is the Legendre symbol (`χ(0) = 0`).
* **(c)** If `p` is Wieferich, every `g_m` is linear, with exactly one root unless `c_1 = 0`.
* **(d)** The triangle class `m = 3` always has the root `l = 0` (`D(3) = 0`). Under (b) it therefore has two roots, unless the root is double.
* **(e)** For `p ≥ 5`, let `F_m : Z_p → Z_p` be the **Camion function** of the class: the unique continuous function with `F_m(l) ≡ A_m(l) (mod p^k)` whenever `m + l(p−1) ≥ 2k − 1`.
  * `F_m` exists and is unique. By §1, `A_m` mod `p^k` is a polynomial on the integers `l` with `m + l(p−1) ≥ 2k − 1`, and those integers are dense in `Z_p`.
  * Its Mahler coefficients satisfy `v(α_i) ≥ ⌈i/2⌉`. Re-expanding the §1 polynomials at the origin keeps the bound.
  * Hence `v(α_i) − v(i!) ≥ (i+1)/4 → ∞`, so `F_m` is a restricted power series in `l`. In a zero class all its coefficients are divisible by `p`, and `F_m/p ∈ Z_p⟨l⟩` has reduction `g_m`.
  * Weierstrass preparation gives exactly `deg g_m` zeros, with multiplicity, in the closed unit disc of `C_p`.
  * Under (b) these are two zeros, either both in `Z_p` (`χ(Δ) = 1`) or a Frobenius-conjugate pair generating the unramified quadratic extension (`χ(Δ) = −1`).
  * Each simple root lifts to a unique `ℓ ∈ Z_p` with `v_p(D(m + l(p−1))) = 1 + v_p(l − ℓ)`. This holds for every integer `l` near `ℓ` for which the right side is below `(m + l(p−1) + 1)/2`, the range of §1.

**Proof.**
* (a)–(d) follow from Statement 3, since `D(m + l(p−1)) ≡ ε^l p·g_m(l) (mod p²)`.
* (e): the power-series coefficients of `F_m` are `b_j = Σ_(i≥j) α_i s(i,j)/i!`, with `s(i,j)` the Stirling numbers of the first kind. So `v(b_j) ≥ min_(i≥j)(⌈i/2⌉ − v(i!)) ≥ (j+1)/4`, using `v(i!) ≤ (i−1)/4` for `p ≥ 5`. In a zero class `v(α_0) ≥ 1`. For `3 ≤ i`, `⌈i/2⌉ − v(i!) ≥ 2`. So `F_m/p ∈ Z_p⟨l⟩` and its reduction is `c_0 + c_1 l + c_2 C(l,2) = g_m`. Weierstrass preparation in `Z_p⟨l⟩` gives `F_m/p = P·U`, with `P` a polynomial of degree `deg g_m` reducing to `g_m` up to a unit and `U` a unit (so `|U(l)| = 1` on `Z_p`). Hensel lifts the simple roots. The valuation formula holds for the integer sequence wherever it agrees with `F_m`. ∎

So, generically, the second-order Camion divisor of every zero class has degree 2, the degree of the edge count, and it is split or inert. At Wieferich primes the degree drops to 1, as in the Bernoulli/Kubota–Leopoldt case, where `λ = 1` is generic.

**FINITE-EXACT** (`second_order_*.out`; 1226 primes `7 ≤ p ≤ 10^4`, 4906 zero classes).

| class | 0 roots | 1 root | 2 roots |
|---|---|---|---|
| 1 | 633 | 2 (both Wieferich) | 591 |
| 2 | 612 | 5 (2 Wieferich) | 609 |
| 3 | 0 | 2 (both Wieferich) | 1224 |
| sporadic | 618 | 5 (2 Wieferich) | 605 |

* **The two Wieferich primes.** At 1093 the zero classes are 3, 468, 1, 2 and each has exactly one root; at 3511 they are 3, 15, 1, 2, again one root each. Every other prime has root counts in {0, 2} apart from 6 classes.
* **Linear and double classes.** There are 13 linear classes: the 8 Wieferich ones and 5 with `X(m−2) ≡ 0`, where about 5.5 are expected. There is 1 double root, where about 5.5 are expected. That deficit is significant at about 3%, which is not significant after multiple comparisons, and is recorded as an observation only.
* **`λ_C(p)`.** The total number of second-order roots, summed over zero classes, is even for 1222 of the 1226 primes. The 4 odd ones, `p = 7, 97, 389, 3989`, each have exactly one non-Wieferich linear or double class. At 1093 and 3511, `λ_C = 4`: four linear classes with one root each.
* **Split versus inert.** In classes 1, 2 and the sporadic classes the split/inert ratio is 1:1 within sampling error (HYP-9361 (a)).

## 5. The Camion limit and the tower height (PROVED)

**Statement 5.** Let `F_1` be the continuous interpolation of `l ↦ ε^(−l) D(1 + l(p−1))`. It agrees with the sequence mod `p^k` whenever `1 + l(p−1) ≥ 2k − 1`, and it is Lipschitz-1 by §1. Then:

* **(a)** `Λ_p := lim_(a→∞) ε^a D(p^a)` exists in `Z_p`, and `ε^a D(p^a) ≡ Λ_p (mod p^a)` for every `a ≥ 1`.
* **(b)** `Λ_p = F_1(−1/(p−1))`, the value of the class-1 Camion function at the `p`-adic point corresponding to `n = 0`. Indeed `p^a = 1 + l_a(p−1)` with `l_a = 1 + p + … + p^(a−1) → −1/(p−1)` in `Z_p`, and `ε^(l_a) = ε^a`.
* **(c)** The Camion tower height of THM-4746 (9) is `h(p) = v_p(Λ_p)`, with `h = ∞` iff `Λ_p = 0`.
* **(d)** For `p ≥ 5`, with base `m = p`: `h(p) ≥ 2` iff `g_p(0) = 0` (Camion–Wieferich). If that root is simple, then `h(p) = 1 + v_p(ℓ + p/(p−1))` for the corresponding `p`-adic zero `ℓ`. So `h(p) ≥ 3` iff `ℓ ≡ p (mod p²)`, a further condition of probability `1/p`.
* **(e)** `Λ_p mod p^K` is determined by `s(n) mod p^K` for `n ≤ p + (2K−2)(p−1)`: evaluate the degree-`(2K−2)` interpolation at `l = −p/(p−1)`. The cost is `O(K²p²)` instead of `O(p^(2K))`.

**Proof.**
* (a), (b): `|l_a − l_∞| = p^(−a)`, Lipschitz continuity, and `F_1(l_a) ≡ ε^a D(p^a) (mod p^a)`, which holds because `p^a ≥ 2a − 1`.
* (c): if `v(Λ_p) < a`, then `v(D(p^a)) = v(Λ_p)`. This is the frozen valuation of THM-4746 (9).
* (d): Statement 4 (e) at the point `l_∞ = −p/(p−1) ≡ p (mod p²)`.
* (e): §1. ∎

**The Fermat parallel, made exact.** For the base-2 Fermat test the analogous limit is `lim 2^(p^a − 1) = ⟨2⟩^(−1)`. The height of `p^a` as a Fermat pseudoprime is `v_p(⟨2⟩ − 1) = 1 + v_p(q_2(p))`, i.e. `v_p(log_p ⟨2⟩)`. So `Λ_p` plays for strong tournaments the role that `log_p⟨2⟩` plays for Fermat. The Camion–Wieferich primes `3, 13, 22409` are the analogues of `1093, 3511`. Through §3 the two lists interact: the Wieferich primes are where the Camion curvature vanishes, not where `Λ_p` is small. Neither 1093 nor 3511 is Camion–Wieferich.

**FINITE-EXACT.**
* The Camion–Wieferich primes `p ≤ 10^5` are exactly **3, 13, 22409**.
  * `cw_fast.c` uses `s(p) ≡ ε + p·u_p (mod p²)` with `u_p = [x^p] 1/T_(<p)(x)`. It reproduces `{13}` on `[5, 10^4]`.
  * 22409 was confirmed independently by `tower_height.c`, which computes `D(p) mod p³` directly: `v_p(D(22409)) = 2`.
* `h(22409) = 2`: `v_p(D(p^a)) = 2` for `a = 1..5` and `v_p(Λ_p) = 2` (`tower_22409_K3.out`).
* The interpolation tool matches direct `D(p²) mod p³` for `p = 7, 11, 13, 17, 19, 23` (`tower_height_direct_checks.out`).
* With `h(3) = h(13) = 2` (THM-4746 (9)), every prime `p ≤ 10^5` has `h(p) ≤ 2`.
* The heuristic expected number of primes with `h ≥ 3` is `Σ_(p≥5) 1/p² ≈ 0.09` (HYP-9358 (a′)).

## 6. The prime 2: Moon parity and the third cyclotomic polynomial

Every tournament automorphism group has odd order (Moon). Hence `f(n) := s(n)/n! = Σ_(strong classes) 1/|Aut|` lies in `Z_(2)`. Moreover `f(n) ≡ #(strong classes on n vertices) (mod 2)`, and similarly `2^C(n,2)/n!` is the count of all tournament classes mod 2.

**Statement 6.**

* **(a) PROVED.** The number of strong tournaments on `n` vertices up to isomorphism is odd iff `n ≢ 2 (mod 3)` (`n ≥ 1`). The number of all tournaments up to isomorphism is even for every `n ≥ 3`.
  * *Proof.* `v_2(2^C(n,2)/n!) = C(n,2) − v_2(n!) ≥ (n−1)(n−2)/2 ≥ 1` for `n ≥ 3`, so `T(x) ≡ 1 + x + x² = Φ_3(x) (mod 2)` in `Z_(2)[[x]]`.
  * Then `S = 1 − 1/T ≡ 1 − (1+x)/(1+x³) (mod 2)`, whose coefficients are 1 exactly at `n ≢ 2 (mod 3)`, `n ≥ 1`. ∎
  * The data of A051337 and A000568 agree to `n = 18`. OEIS records no parity comment on either.
* **(b) FINITE-EXACT, with a structural reason.** `f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))` holds for all `1 ≤ n ≤ 600` when `j ≤ 4`. For `j = 5..7` it fails only at `n ≤ 3`. It holds mod `2^(j+3)` on the class `n ≡ 1 (mod 3)`.
  * *Reason.* Mod `2^K`, `1/T` is the reciprocal of a polynomial `Q` of degree about `√(2K)` with `Q ≡ Φ_3 (mod 2)`.
  * Two roots `ρ, σ(ρ)` of `Q` are units near the primitive cube roots of unity in `W(F_4) = Z_2[ζ_3]`. The others have `|ρ| > 1`.
  * So `f(n) = −Tr_(Q_2(ζ_3)/Q_2)(c·ρ^(−n))` plus terms that vanish mod `2^K` for large `n`.
  * At `p = 2` the role of `F_p^×` (order `p − 1`) is played by `F_4^×` (order 3), because the reduction of the tournament EGF is the irreducible `Φ_3`.
* **(c) FINITE-EXACT.** On the zero class `n ≡ 2 (mod 3)`, `v_2(f(n)) = 2 + v_2(n − n*)` for a 2-adic integer
      n* ≡ 16221783116972227762  (mod 2^64),   n* ≡ 3250 (mod 4096).
  * This holds exactly for all 499 values `5 ≤ n ≤ 1500` (`two_adic_zero_1500.out`). The 64 bits come from the companion-matrix evaluation at large `n` (`two_adic_zero_fast.out`).
  * Rational reconstruction at `2^64` finds no fraction with numerator and denominator below about `2^31`.
  * The universal class "2" of odd `p` (no strong K2) is thus the only 2-adic zero class, and it carries a single simple zero at a "virtual order" `n*`. The exact zero `s(2) = 0` is below the threshold, as `D(1) = 0` is for class 1 at odd `p`.
  * HEURISTIC: by (b), `n*` is a ratio of 2-adic logarithms of elements of the norm-one torus of `Q_2(√−3)`.

## 6b. Degree = uniformity: r-uniform hypergraphs (PROVED; statistics FINITE-EXACT)

Let `r ≥ 2` and `T_r(n) = 2^C(n,r)`. Let `b_r` be the connected labelled `r`-uniform hypergraphs, given by `Δb_r = ΔT_r ∗ T_r^(−1)`.

**Statement 7.** Let `p` be a prime with `gcd(ord_p(2), r!) = 1`.

* **(a)** `ω^C(n,r)` has period dividing `p − 1`. By Vandermonde, `C(n+P, r) − C(n, r) = Σ_(i≥1) C(P,i) C(n, r−i)`, and `ord_p(2) | C(P,i)` when `ord_p(2) | P` and `gcd(ord_p(2), i!) = 1`.
* **(b)** Session25's argument therefore works over `Z/p^k` itself, with all nodes in `Z/p^k` and multiplier 1. The filtration now has degree `d ≤ rj` at grade `j`.
* **(c)** For `m ≥ r(k−1) + 2`, the class sequence `l ↦ b_r(m + l(p−1))` has Mahler coefficients `a_i ≡ 0 (mod p^min(k, ⌈i/r⌉))`. So the Mahler speed is `1/r`.
* **(d)** Mod `p²` this is a polynomial of degree `≤ r` in `l`. Its coefficients of `C(l,2), …, C(l,r)` are `p·u·(…)`, with `u ≡ −q_2(p)`. So they vanish mod `p²` at Wieferich primes.

*Proof.* Repeat §1 with `d ≤ rj`.

* Unit nodes give valuation at least `j + max(0, i − d) ≥ max(j, i − (r−1)j) ≥ ⌈i/r⌉`.
* Nonunit nodes vanish once `n − (r−1)j ≥ k`, i.e. for `n ≥ r(k−1) + 1`. One more index is needed for the shift in `Δb_r`.
* Mod `p²` only grades `j ≤ 1` survive. The grade-0 part is linear in `l`, and the grade-1 part is `pu·C(n,r)` times grade-0 terms. ∎

**Root statistics** (FINITE-EXACT, `hyper3_second_order_*.out`): `r = 3`, all 60 primes `p ≤ 2200` with `gcd(ord_p(2), 6) = 1`, 60444 classes.

* The cubic law holds in every class: `p | a_1, a_2, a_3` and `p² | a_4`. For `p ≤ 127` it was also checked directly against `b_3(m + l(p−1)) mod p²` for `l ≤ p + 1` (33166 checks).
* 58 classes are not genuinely cubic (`p² | a_3`), about the `1/p` rate.
* Zero classes have the second-order root counts below. A random cubic over `F_p` has 0, 1, 3 roots with probabilities `1/3, 1/2, 1/6`.

| zero class | 0 roots | 1 | 2 | 3 | random cubic (0, 1, 3) |
|---|---|---|---|---|---|
| 1 (universal: `b_3(p) ≡ 2^((p−1)(p−2)/6) − 1 ≡ 0`) | 22 | 25 | 1 | 12 | 20, 30, 10 |
| 2 (exact zero `b_3(2) = 0`: forced root) | 0 | 31 | 0 | 29 | forced root + random quadratic: 30, —, 30 |
| sporadic | 21 | 21 | 0 | 14 | 18.7, 28, 9.3 |

For tournaments the root counts are 0 or 2 (quadratic). For 3-uniform hypergraphs they are 0, 1 or 3, in random-cubic proportions. This confirms the "degree = uniformity" principle: the second-order Camion divisor has degree `r`, the degree of the edge count `C(n,r)`.

## 7. What this does and does not say

* **Says.**
  * Camion pseudoprimality is the shadow of a `p`-adic analytic object. That object is the class functions `F_m` on `Z_p`, at half the Kummer rate, twisted by `(2/p)`.
  * Its second-order structure is governed by the Fermat quotient through the quadratic edge count.
  * The tower heights `h(p)` are valuations of one `p`-adic number, `Λ_p`, the class-1 branch evaluated at `n = 0`.
  * The prime 2 behaves differently: Moon's odd-order theorem makes the groupoid counts `Σ 1/|Aut|` the natural 2-adic object, and `Φ_3` replaces `F_p^×`.
* **Does not say.**
  * Nothing about the infinitude or density of Camion–Wieferich primes (HYP-9358 (a)).
  * Nothing about the 1:1 split/inert law beyond data (HYP-9361 (a)).
  * Nothing about the repunit rarity in the `H`-spectrum (opus S51's lane, THM-4748/4749).
  * The "virtual order" `n*` is not identified.
