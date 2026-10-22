---
id: THM-4750
title: "The Camion function. Let p be an odd prime, eps = (2/p), s(n) the strong labelled tournaments (A054946), D = s - c with c(n) the circulant count, and A_m(l) = eps^(-l) D(m + l(p-1)). (1) Half-speed higher Kummer congruences: for k >= 1 and m >= 2k-1 the Mahler coefficients a_i(m) = Delta^i A_m(0) satisfy a_i(m) = 0 mod p^min(k, ceil(i/2)); in particular sum_j C(2k-1,j)(-eps)^(2k-1-j) D(m+j(p-1)) = 0 mod p^k. The exponent ceil(i/2), half of the Bernoulli/Kummer rate, comes from the edge count C(n,2) having degree 2. (2) Speed is full exactly at Wieferich primes: then a_i(m) = 0 mod p^min(k,i) for m >= k (for base a: exactly at base-a Wieferich primes, e.g. strong semicomplete digraphs at p = 11). (3) Fermat-quotient curvature: for m >= 3, a_2(m)/p = -q_2(p) X(m-2) mod p, where q_2 is the Fermat quotient of 2 and X(m-2) counts the tournaments on [m] made strong by making one fixed pair two-way (EGF T''/T^2); so mod p^2 every class sequence is a quadratic in l. (4) Second-order zeros: in a zero class the l mod p with p^2 | D(m + l(p-1)) are the roots of an explicit quadratic over F_p (linear at Wieferich primes); for p >= 5 the class function has exactly two p-adic zeros, split, inert, or (for a double root mod p) congruent and possibly ramified; the triangle class is split whenever p is not Wieferich and its root l = 0 is simple (ramified at p = 5). (5) The Camion limit Lambda_p = lim eps^a D(p^a) exists, eps^a D(p^a) = Lambda_p mod p^(a+1), Lambda_p = sigma_p - <2>^(-1/2) where sigma_p = lim eps^a s(p^a) depends only on 2 mod p, and the tower height is h(p) = v_p(Lambda_p) = v_p(2 - kappa(2)) with the Camion lift kappa(2) = omega(2) sigma_p^(-2), the unique base above 2 mod p for which every p^k is a Camion pseudoprime (Teichmueller-lift analogue). (6) The prime 2 (Moon: |Aut| odd): the number of strong tournaments on n vertices up to isomorphism is odd iff n is not 2 mod 3 (T = Phi_3 mod 2), the number of all tournaments is even for n >= 3, and the groupoid count sum 1/|Aut| has a simple 2-adic zero n* on n = 2 mod 3, a ratio of logarithms in the norm-one torus of Q_2(sqrt(-3)). (7) r-uniform hypergraphs run at speed 1/r with degree-r class polynomials mod p^2 (gcd(ord_p 2, r!) = 1)."
status: >
  AUDITED 2026-10-21 (independent adversarial audit, 05-knowledge/results/camion_function_audit_20261021.md; MISTAKE-654).
  The audit found the core proofs sound and reproduced every FINITE-EXACT number with independent code. It made two
  MAJOR corrections, both applied here:
  - M1: a double root mod p can be a ramified pair of zeros, so "split or inert" and "the triangle class is always split"
    were false at p = 5 (triangle) and p = 7 (class 2);
  - M2: the rational-reconstruction evidence for n* contradicted the session's own output; it is replaced by the audit's
    2^140 computation.
  Sixteen minor fixes were also applied, together with the audit's upgrades (the congruence mod p^(a+1), the cheaper
  Lambda_p formula, the exact double-factorial excess, and the proof sketch for 6(b), 6(c)).
  PROVED:
  - statements 1, 2 (with the converse for "exactly"), 3, 4(a)-(e), 5(a)-(e), 6(a), 7(a)-(d);
  - 6(b) (with its range n >= n_0(j), linear in j) and 6(c) in outline (entire-function argument, first audit finding 23,
    second audit finding 13).
  SECOND AUDIT PASS (2026-10-21, 05-knowledge/results/camion_function_audit_pass2_20261021.md; MISTAKE-655) on the
  post-audit continuation:
  - §5b (Teichmueller rigidity, Lambda_p = sigma_p - <2>^(-1/2), the Camion lift, the decomposition) is PROVED and
    confirmed numerically (direct rigidity to mod 3^10; tower height = v_p(a - kappa(a)) for all 794 integer bases
    below 5^4 and 7^3);
  - the r = 4 data, the base-3 list and the connected list (extended to 10^4) are confirmed;
  - one MAJOR fix: 6(b) needs its range (P2-M1);
  - the connected-graph list is not a tower-height list (drifted tower).
  FINITE-EXACT:
  - sharpness of ceil(i/2) for every tested non-Wieferich p >= 7;
  - extra divisibility exactly v_p((2k-1)!!) at p = 3, 5, 7;
  - full speed from m >= k at 1093, 3511 (base 2), at 11 (base 3) and at 5 (base 7);
  - the curvature formula in every class of every prime p <= 10^4;
  - the quadratic law directly for p <= 59 (17158 checks);
  - the root table over 1226 primes 7 <= p <= 10^4 (4906 zero classes);
  - the Camion-Wieferich primes p <= 10^5 are exactly 3, 13, 22409 (cw_fast.c and the audit's independent FLINT search
    b6), with h(22409) = 2 (two codes);
  - the 2-adic quasi-period;
  - n* mod 2^140, not a rational of height <= 2^69.5, equal to the log-ratio closed form to 90 bits;
  - the r = 3 and r = 4 root statistics (random-cubic and random-quartic proportions).
  LITERATURE: the audit found no prior statement of the higher Kummer congruences, the curvature formula, the
  Wieferich collapse, the Camion limit or lift, the strong-class parity law, or the 2-adic zero. Nearest prior art:
  Mani-Stones 2016 (first-order quasi-periods for connected graphs), Specker-Blatter and Fischer-Makowsky (periodicity),
  O'Desky-Richman 2023 (p-adic continuity of derangement-type sequences). The parity law is a two-line consequence of
  Moon's theorem and may be folklore.
session: mac-mini-2026-10-21-fanopsp (owner prompt: hunt h(p) >= 3, HYP-9357 at 27 and 31, repunit rarity, wide-net loops)
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md (section 13)
scripts:
  - 04-computation/experiments/camion_wieferich_20261021/higher_kummer.py (+ higher_kummer.out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_speed.c (+ camion_speed.out), camion_speed_base.c (+ camion_speed_base.out)
  - 04-computation/experiments/camion_wieferich_20261021/curvature_meaning.py (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_second_order.c (+ second_order_*.out, second_order_direct_3_59.out)
  - 04-computation/experiments/camion_wieferich_20261021/tower_height.c (+ tower_22409_K3.out, tower_height_direct_checks.out)
  - 04-computation/experiments/camion_wieferich_20261021/cw_fast.c (+ cw_*.out, cw_fast_validation_5_10000.out); cw_flint_np.py (+ cw_flint_*.out, search to 10^6)
  - 04-computation/experiments/camion_wieferich_20261021/camion_lift.py, camion_lift_check.py, camion_lift_map.py, camion_residue.c, camion_residue_fish.py, camion_residue_symmetry.py (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/camion_irregular.c (+ camion_irregular_10000.out)
  - 04-computation/experiments/camion_wieferich_20261021/two_adic_face.py, two_adic_zero.py, two_adic_zero_fast.py, two_adic_regulator.py (+ .out)
  - 04-computation/experiments/camion_wieferich_20261021/hyper3_second_order.c, hyper_r_second_order.c, hyper_r_tabulate.py, hyper3_class2_check.py (+ hyper3_second_order_*.out, hyper4_second_order_5_2200.out, hyper3_class2_check.out)
  - 04-computation/experiments/camion_wieferich_20261021/cw_base.c, cw_connected.c (+ cw_base_3000.out, cw_base3_*.out, cw_connected_4000.out)
  - 04-computation/experiments/camion_wieferich_20261021/audit/ (independent audit b1-b11)
depends_on:
  - 05-knowledge/results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md (the filtered algebra F and S = e - T^(-1) in F; PROVED, independently audited)
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md (statement 9, towers)
  - J. W. Moon, Tournaments with a given automorphism group, Canad. J. Math. 16 (1964); Topics on Tournaments (1968): the automorphism group of a tournament has odd order (CITED)
related:
  - 05-knowledge/hypotheses/HYP-9361-camion-function-open-laws-split-inert-equidistribution-tower-heights-uniformity-k-speed.md
  - 05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-wieferich-type-primes-sporadic-zero-classes-prime-squares.md
  - 05-knowledge/results/camion_function_audit_20261021.md (independent audit)
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md (Wieferich primes in circulant tournaments: a different appearance of q_2)
---

# THM-4750 — The Camion function

**One line.**
* Along each residue class mod `p − 1`, the strong tournament count is a `p`-adic function of the class index. It is a restricted power series for `p ≥ 5`.
* It converges at half the Kummer rate, because the edge count `C(n,2)` is quadratic, and its curvature is the Fermat quotient of 2 times a tournament statistic.
* So the Wieferich primes 1093 and 3511 are exactly where the curvature vanishes mod `p²`. There every Camion zero class has one second-order root instead of 0 or 2.
* The tower heights of THM-4746 measure the distance from 2 to a canonical "Camion lift", as the Wieferich order measures the distance to the Teichmüller lift.

## Notation

* `p` is an odd prime and `ε = (2/p)`.
* `s(n)` counts strong labelled tournaments (`s(1) = 1`).
* `c(n)` is the circulant count: `2^((n−1)/2)` for odd `n`, 0 for even `n`.
* `D = s − c`.
* `T(n) = 2^C(n,2)`, with EGFs `T(x)` and `S(x) = 1 − 1/T(x)`.
* For `m ≥ 1` the **class sequence** is `A_m(l) = ε^(−l)·D(m + l(p−1))`, `l ≥ 0`, with **Mahler coefficients** `a_i(m) = Δ^i A_m(0) = Σ_j (−1)^(i−j) C(i,j) A_m(j)`. Then `A_m(l) = Σ_i a_i(m) C(l,i)`.
* `q_a(p) = (a^(p−1) − 1)/p` is the Fermat quotient. `p` is **Wieferich** (to base 2) if `p | q_2(p)`.
* `X(n) = n!·[x^n] T''(x)/T(x)²`, with values `X = Δ²T ∗ T^(∗−2)` in the Hurwitz ring of binomial convolution.
  * Since `T = 1/(1−S)`, also `T''/T² = S'' + 2S'²T`.
  * Equivalently `C(n,2)·X(n−2) = 2·(d/da) s_a(n)|_(a=2)`, where `s_a` is the strong count built from `a^C(n,2)`: the derivative of the strong count in the base.
  * **Combinatorially (PROVED), `X(m−2)` is the number of labelled tournaments on `[m]` that become strongly connected when the arc at one fixed pair `{u,v}` is made two-way.**
    * *Proof.* Make `{u,v}` two-way and orient the other pairs. The condensation of such a digraph is transitive, and `u, v` lie in one strong component. So the `T(m)/2` such digraphs factor as tournament ∗ (strong block containing `u, v`) ∗ tournament, and the strong ones number `(Δ²T/2) ∗ T^(∗−2) = X/2`. Each is counted twice, once for each orientation of `uv` in the tournament.
    * Values: `2, 4, 36, 680, 25152, 1784416, 244023264` for `m = 2..8`.
    * Brute force agrees for `m ≤ 7` (`curvature_meaning.out`; audit b3 for `m ≤ 6`).
    * In particular `X(1) = 4 ≢ 0 (mod p)`, so the triangle class is never linear at a non-Wieferich prime.

**The representation used** (Session25, PROVED and audited). Work over `R_0 = Z/p^k` and its unramified quadratic extension `R`.

* `θ` is a root of `θ^(p−1) = ε`, `G = θR_0`, and `ω` is the Teichmüller lift of 2 with `2 = ω(1 + pu)`. Then `u ≡ −q_2(p) (mod p)`, because `2^(p−1) = (1+pu)^(p−1) ≡ 1 − pu (mod p²)`.
* `F` is the `R`-module generated by `p^j δ_d ∗ E_β` with `0 ≤ j < k`, `0 ≤ d ≤ 2j` and `β ∈ G`. Here `E_β = (β^n)_n` and `(δ_d ∗ E_β)(n) = C(n,d) β^(n−d)`.
* Session25 proves that `F` is a ring and that `S = e − T^(−1) ∈ F`.
* The grading comes from the Newton expansion `(1+pu)^C(n,2) = Σ_j (pu)^j C(C(n,2), j)`: the `j`-th term has degree `2j` in `n`.

## 1. Half-speed higher Kummer congruences (PROVED)

**Statement 1.** For every `k ≥ 1`, every `m ≥ 2k − 1` and every `i ≥ 0`,

    a_i(m) ≡ 0   (mod p^min(k, ⌈i/2⌉)).

* Equivalently, `v(a_i(m)) ≥ min(⌊(m+1)/2⌋, ⌈i/2⌉)` for every `m ≥ 1`.
* So for `m ≥ 2k − 1` the class sequence is, mod `p^k`, a polynomial in `l` of degree at most `2k − 2`, and `v_p` of its `C(l,i)` coefficient is at least `⌈i/2⌉`.
* The case `i = 2k − 1` is the **higher Kummer congruence**

      Σ_(j=0)^(2k−1) C(2k−1, j) (−ε)^(2k−1−j) D(m + j(p−1)) ≡ 0   (mod p^k)   (m ≥ 2k − 1).

  The same holds with `D` replaced by `s`.

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

**What the exponent means.**
* For a single exponential `n ↦ a^n` the Mahler coefficients along a class have valuation at least `i`. That is the rate of the classical Kummer congruences for Bernoulli numbers, which underlies Kubota–Leopoldt interpolation.
* Here each power of `pu` arrives with two powers of `n` (one edge, two endpoints), so the rate halves.
* Hypergraphs built from `2^C(n,r)` run at rate `1/r` (§6b).

**The same proof applies to every sequence in Session25's closure of `F`.** Examples:

* strong semicomplete digraphs, with `ε = (3/p)` and `p ≠ 3`;
* connected labelled graphs `b(n)`, from `m ≥ 2k` (one shift: `Δb ∈ F`).

**FINITE-EXACT.**
* `higher_kummer.out`: 13 pairs `(p,k)` with `k ≤ 5`. Failures occur only below `m = 2k − 1`.
* Audit b1, using exact big integers through the Hurwitz inverse: 34400 coefficient checks over 28 pairs, 0 violations. The `⌊(m+1)/2⌋` form holds for every `m ≥ 1`.

## 2. Sharpness, small primes, and full speed exactly at Wieferich primes

**Statement 2 (PROVED).** `p` is Wieferich iff the speed is full. If `p` is Wieferich, then `a_i(m) ≡ 0 (mod p^min(k,i))` for every `m ≥ k`.

*Proof that Wieferich gives full speed.* Now `pu = p²u'`, so the grade-`j` Newton term carries `p^(2j)` at degree `≤ 2j`.

* Let `F'` be generated by `p^J δ_d ∗ E_β` with `d ≤ J`. It is a ring, since exponents and degrees add, and `T ∈ F'` (take `J = 2j`).
* The inversion argument of Session25 goes through unchanged: `A` is the same group-algebra unit, and `U` has `J ≥ 2`. So `S ∈ F'`.
* For a generator, the Leibniz estimate gives valuation at least `J + max(0, i − d) ≥ i`.
* Nonunit nodes contribute valuation at least `J + n − d ≥ n ≥ k`.
* The circulant part is full speed in any case. ∎

*Converse.* If `p` is not Wieferich, statement 3 gives `v(a_2(3)) = 1`, because `X(1) = 4 ≢ 0`. So the speed is not full.

**General base (PROVED by the same argument).** Session25's representation holds for `T_a(n) = a^C(n,2)` with any `a` prime to `p`, with multiplier `(a/p)`. Statements 1–3 then hold with `q_2` replaced by `q_a(p)` and `X` by `X_a = T_a''/T_a²`.

**Full speed occurs exactly at the base-`a` Wieferich primes.** The converse uses `X_a(1) = a³ − 2a` and `X_a(2) = a⁶ − 4a³ − 2a² + 6a` (audit finding 5).

* Suppose `p ∤ a` and `X_a(1) ≡ 0`. Then `a² ≡ 2`, so `X_a(2) ≡ 4 − 2a`.
* That vanishes only if `a ≡ 2`, and then `4 ≡ 2`, which is impossible for odd `p`.
* So class 3 or class 4 has `v(a_2) = 1` whenever `q_a(p) ≢ 0`.

For strong semicomplete digraphs (`a = 3`, three states per pair) the full-speed primes are `p = 11` and `1006003`.

**FINITE-EXACT** (`camion_speed.out`, `camion_speed_base.out`, audit b2). The table gives the minimum of `v_p(a_i(m))` over all `p − 1` classes, with `m ≥ 2K − 1`.

| `p` (base) | `K` | `min v_p(a_i)`, `i = 0, 1, …` | bound `⌈i/2⌉` |
|---|---|---|---|
| 11, 13, 17 | 4 | `0 1 1 2 2 3 3 4` | equal |
| 1091, 1097, 3517 | 3 | `0 1 1 2 2 3` | equal |
| **1093** (Wieferich) | 4 | `0 1 2 3 4 4 4 4` (cap 4) | full speed |
| **3511** (Wieferich) | 3 | `0 1 2 3 3 3` (cap 3) | full speed |
| **11** (base 3; `3^10 ≡ 1 mod 121`) | 4 | `0 1 2 3 4 4 4 4` | full speed |
| **5** (base 7; `7^4 ≡ 1 mod 25`) | 5 | `0 1 2 3 4 5 5 5 5 5` | full speed |
| 13, 17 (base 3); 7, 11 (base 5) | 4 | half speed, attained | equal |
| 5 | 6 | `0 1 1 2 2 4 4 5 5 6 6 6` | +1 from `i = 5` |
| 7 | 5 | `0 1 1 2 2 3 3 5 5 5` | +1 from `i = 7` |

* **Threshold `m ≥ k` is sharp.** Audit b2 found full speed from `m ≥ K` at 1093 (`K ≤ 5`), 3511 (`K ≤ 4`), base 3 at 11 and base 7 at 5, with failures exactly at `m = 1..K−1`.
* **Sharpness of half speed.** For non-Wieferich `p ≥ 7` the half-speed bound is attained (audit b1). Each tested case has at most 6 exceptional classes. At `p = 1091`, 1089 of the 1090 classes attain `v(a_2) = 1`. The exception is `m = 603`, where `X(601) ≡ 0`.
* **Small primes.** The excess over `⌈i/2⌉` at `i ∈ {2k−1, 2k}` is exactly `v_p((2k−1)!!)` wherever the cap allows (FINITE-EXACT; audit finding 28):
  * `p = 3`: +1, +1, +1, +3 at `k = 2, 3, 4, 5`;
  * `p = 5`: +1 at `k = 3, 4, 5`;
  * `p = 7`: +1 at `k = 4`.
* **Mechanism.** The top Newton coefficient of `C(C(n,2),k)` at degree `2k` is `(2k)!/(2^k k!) = (2k−1)!!`. This is the mechanism conjectured in HYP-9358 (e2); the threshold statement there remains open.

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
* Audit b3, computing `X` through the other identity `S'' + 2S'²T`: 13806 classes, `p ≤ 400`.
* The quadratic law was checked directly against `D(m + l(p−1)) mod p²` for `l ≤ p + 1` in all classes:
  * `p ≤ 59`: 17158 checks (`second_order_direct_3_59.out`);
  * `p ≤ 61`: 20928 checks (audit b3).

**Reading.** Mod `p²` the strong count is "Teichmüller part + Fermat drift × base-derivative". The base-derivative is the count of tournaments made strong by one two-way pair.
* The Teichmüller part `2 ↦ ω(2)` makes `T` a pure exponential sum, linear in `l`.
* The drift `⟨2⟩ = 1 + pu` enters through `⟨2⟩^C(n,2)`, and its first-order effect is the base-derivative.
* The curvature in `l` is the quadratic edge count. At Wieferich primes `2` is Teichmüller to second order, so the curvature vanishes.

## 4. Second-order zeros (PROVED structure; FINITE-EXACT statistics)

Call `m ≥ 3` a **zero class** if `p | D(m)` (equivalently, `p | D` on the whole class). This includes the universal classes `1, 2, 3` (THM-4746). For a zero class put `c_i = a_i(m)/p mod p` and

    g_m(l) = c_0 + c_1 l + c_2 C(l,2) ∈ F_p[l],   c_2 = −q_2(p) X(m−2).

**Statement 4.**

* **(a)** For `l ≥ 0`, `p² | D(m + l(p−1))` iff `g_m(l mod p) = 0`.
* **(b)** If `q_2(p)X(m−2) ≢ 0`, then `g_m` is a quadratic. It has `1 + χ(Δ)` roots in `F_p` (a double root when `χ(Δ) = 0`), where `Δ = (c_1 − c_2/2)² − 2c_0c_2` and `χ` is the Legendre symbol.
* **(c)** If `p` is Wieferich, every `g_m` is linear, with exactly one root unless `c_1 = 0`.
* **(d)** The triangle class `m = 3` always has the root `l = 0` (`D(3) = 0`), and it is never linear at a non-Wieferich prime (`X(1) = 4`). So for non-Wieferich `p` it has two distinct roots in `F_p` (split), unless `l = 0` is a double root.
* **(e)** For `p ≥ 5`, let `F_m : Z_p → Z_p` be the **Camion function** of the class: the unique continuous function with `F_m(l) ≡ A_m(l) (mod p^k)` whenever `m + l(p−1) ≥ 2k − 1`.
  * `F_m` exists and is unique. By §1, `A_m` mod `p^k` is a polynomial on the integers `l` with `m + l(p−1) ≥ 2k − 1`, and these integers are dense in `Z_p`.
  * Its Mahler coefficients satisfy `v(α_i) ≥ ⌈i/2⌉`; re-expanding the §1 polynomials at the origin keeps the bound.
  * Hence `F_m` is a restricted power series in `l`. In a zero class `F_m/p ∈ Z_p⟨l⟩`, with reduction `g_m`.
  * Weierstrass preparation gives exactly `deg g_m` zeros, with multiplicity, in the closed unit disc of `C_p`. When `g_m` is a genuine quadratic there are three cases:
    * **split**, `χ(Δ) = 1`: two zeros in `Z_p`, distinct mod `p`;
    * **inert**, `χ(Δ) = −1`: a Frobenius-conjugate pair generating the unramified quadratic extension;
    * **congruent**, `χ(Δ) = 0`: two zeros congruent mod `p`. They may be split, unramified or **ramified**; the valuations of `D` on the disc decide which.
  * **Ramified examples** (audit b11, exact integers):
    * `p = 5`, triangle class: `g_3 = 4l + 3C(l,2)`, double root at `0`, and `v_5(D(3 + 4l)) = 2` for all 39 valid `l ≡ 0 (mod 5)`;
    * `p = 7`, class 2 (`m = 8`): `v_7(D(8 + 6l)) = 2` for all 19 valid `l ≡ 6 (mod 7)`.
    * In both, the local factor is Eisenstein, so the zeros lie at distance `p^(−1/2)` from the root mod `p`. At `p = 5` the exact zero `D(3) = 0` lies below the level-3 threshold and forces no rational zero.
  * **Valuation law.** Each simple root of `g_m` lifts to a unique `ℓ ∈ Z_p`. For integers `l ≡ ℓ (mod p)` with `n = m + l(p−1)`,

        v_p(D(n)) = 1 + v_p(l − ℓ)   whenever   1 + v_p(l − ℓ) < ⌈n/2⌉.

    Off the root discs, `v_p(D(n)) = 1`.
  * Audit b10 confirms the law with exact integers for `n ≤ 800` at `p = 7, 11, 13`.

**Proof.**
* (a)–(d) follow from statement 3, since `D(m + l(p−1)) ≡ ε^l p·g_m(l) (mod p²)`.
* (e): the power-series coefficients of `F_m` are `b_j = Σ_(i≥j) α_i s(i,j)/i!`, with `s(i,j)` the Stirling numbers of the first kind. So `v(b_j) ≥ min_(i≥j)(⌈i/2⌉ − v(i!)) ≥ (j+1)/4`, using `v(i!) ≤ (i−1)/4` for `p ≥ 5`.
  * In a zero class `v(α_0) ≥ 1`, and `⌈i/2⌉ − v(i!) ≥ 2` for `i ≥ 3`. So the reduction of `F_m/p` is `c_0 + c_1 l + c_2 C(l,2) = g_m`.
  * Weierstrass preparation in `Z_p⟨l⟩` gives `F_m/p = P·U`. Here `P` is a polynomial of degree `deg g_m` reducing to `g_m` up to a unit, and `U` is a unit, so `|U(l)| = 1` on `Z_p`.
  * Hensel lifts the simple roots.
  * `F_m ≡ A_m (mod p^k)` needs `n ≥ 2k − 1`, i.e. agreement to precision `⌈n/2⌉`. That gives the stated range. ∎

So, generically, the second-order Camion divisor of a zero class has degree 2, the degree of the edge count, and it is split or inert. Ramified pairs occur only on the rare double roots. At Wieferich primes the degree drops to 1, as in the Bernoulli/Kubota–Leopoldt case, where `λ = 1` is generic.

**FINITE-EXACT** (`second_order_*.out`; 1226 primes `7 ≤ p ≤ 10^4`, 4906 zero classes; reproduced exactly by audit b4).

| class | 0 roots | 1 root | 2 roots |
|---|---|---|---|
| 1 | 633 | 2 (both Wieferich) | 591 |
| 2 | 612 | 5 (2 Wieferich) | 609 |
| 3 | 0 | 2 (both Wieferich) | 1224 |
| sporadic | 618 | 5 (2 Wieferich) | 605 |

* **The two Wieferich primes.** At 1093 the zero classes are 3, 468, 1, 2 and each has exactly one root; at 3511 they are 3, 15, 1, 2, again one root each. This was predicted before the test. Every other prime has root counts in {0, 2} apart from 6 classes.
* **Linear and double classes.** There are 13 linear classes: the 8 Wieferich ones and 5 with `X(m−2) ≡ 0`, at `(23: 16)`, `(23: 2)`, `(97: 24)`, `(389: 2)` and `(3989: 2275)`.
  * About 4.0 linear classes are expected, since class 3 can never be linear.
  * There is 1 double root, at `(7: 2)`, against about 5.5 expected (`P(≤ 1) ≈ 0.027`). This is an observation only, given multiple comparisons.
  * The only class-3 double root is at `p = 5`, outside the table.
* **`λ_C(p)`.** The total number of second-order roots over the zero classes is even for 1222 of the 1226 primes. The 4 odd ones, `p = 7, 97, 389, 3989`, each have exactly one non-Wieferich linear or double class. At 1093 and 3511, `λ_C = 4`.
* **Split versus inert.** In classes 1, 2 and the sporadic classes the ratio is 1:1 within sampling error.
  * Classes 1 and 2 at the same prime are independent: `χ² = 0.14`.
  * There is no dependence on `ε`: 0.483 at `ε = +1` against 0.501 at `ε = −1`.
  * See HYP-9361 (a).

## 5. The Camion limit and the tower height (PROVED)

**Statement 5.** Let `F_1` be the class-1 Camion function, the continuous interpolation of `l ↦ ε^(−l) D(1 + l(p−1))`. It agrees with the sequence mod `p^k` whenever `1 + l(p−1) ≥ 2k − 1`.

* **(a)** `Λ_p := lim_(a→∞) ε^a D(p^a)` exists in `Z_p`, and `ε^a D(p^a) ≡ Λ_p (mod p^(a+1))` for every `a ≥ 1`.
* **(b)** `Λ_p = F_1(−1/(p−1))`, the value of the class-1 Camion function at the `p`-adic point corresponding to `n = 0`. Indeed `p^a = 1 + l_a(p−1)` with `l_a = 1 + p + … + p^(a−1) → −1/(p−1)` in `Z_p`, and `ε^(l_a) = ε^a`.
* **(c)** The Camion tower height of THM-4746 (9) is `h(p) = v_p(Λ_p)`, with `h = ∞` iff `Λ_p = 0`.
* **(d)** Take `p ≥ 5` and base `m = p`, with `z = −p/(p−1) = p + p² + … ≡ p (mod p²)`.
  * `h(p) ≥ 2` iff `g_p(0) = 0`, i.e. `p` is Camion–Wieferich.
  * If that root is simple, `h(p) = 1 + v_p(z − ℓ)` for the corresponding `p`-adic zero `ℓ`. So `h(p) ≥ 3` iff `ℓ ≡ p (mod p²)`, a further condition of probability `1/p`.
* **(e)** For `p ≥ 2K − 1`, `Λ_p mod p^K` is determined by:
  * `D(p) mod p^K`;
  * `D(p + j(p−1)) mod p^(K−1)` for `1 ≤ j ≤ 2K − 4`.

  Evaluate the class polynomial at `z`. For `K = 3`: `Λ_p ≡ ε[D(p) + a_1 z + a_2 C(z,2)] (mod p³)`, with `a_1, a_2` needed only mod `p²` (audit finding 15). The cost is `O(K²p²)`, instead of the `O(p^(2K))` of computing `D(p^K)` directly.

**Proof.**
* (a), (b): `F_1` is a contraction, `|F_1(x) − F_1(y)| ≤ p^(−1)|x − y|`, because `⌈i/2⌉ ≥ 1 + ⌊log_p i⌋` for `i ≥ 1` (audit finding 14). With `|l_a − l_∞| = p^(−a)` and `F_1(l_a) ≡ ε^a D(p^a) (mod p^(a+1))`, which holds as `p^a ≥ 2a + 1`, both follow.
* (c): if `v(Λ_p) ≤ a`, then `v(D(p^a)) = v(Λ_p)`. This is the frozen valuation of THM-4746 (9).
* (d): statement 4 (e) at `z`.
* (e): `v(C(z,i)) ≥ 1` for `1 ≤ i < p`, combined with §1. ∎

### 5b. Teichmüller rigidity and the Camion lift (PROVED; post-audit continuation)

**Notation.**
* For a base `a` prime to `p`, let `s_a = e − T_a^(−1)` with `T_a(n) = a^C(n,2)`.
* Put `c_a(n) = a^((n−1)/2)` for odd `n`, `c_a(n) = 0` for even `n`, `D_a = s_a − c_a` and `ε_a = (a/p)`.
* Call `n` a **Camion pseudoprime to base `a`** if `s_a(n) ≡ c_a(n) (mod n)`. Tournaments are `a = 2`.
* Let `h_a(p) = max{k : p^k | D_a(p^k)}`.
* **`p`-adic bases.** For `a ∈ Z_p^×` the same recurrence defines `s_a`, and `s_a(n) mod p^r` depends only on `a mod p^r`, because it is an integer polynomial in `a`. "`p^k` is a Camion pseudoprime to base `a`" then means `p^k | D_a(p^k)` in `Z_p`.
* With this convention `κ(a)` below is a genuinely `p`-adic object. An integer base `≡ κ(a) (mod p^K)` has height at least `K`, and in general no integer base passes every `p^k`.

**Statement 5′.**

* **(a) Rigidity.**
  * Finite level: if `a ≡ a′ (mod p)`, then `s_a(p^k) ≡ s_(a′)(p^k) (mod p^(k+1))` for every `k ≥ 1`.
  * Consequently `σ_p(a) := lim_k ε_a^k s_a(p^k)` exists, `σ_p(a) ≡ 1 (mod p)`, and `σ_p(a)` depends only on `a mod p`.
* **(b) The limit splits.** `Λ_p(a) := lim_k ε_a^k D_a(p^k) = σ_p(a) − ⟨a⟩^(−1/2)`, where `⟨a⟩ = a/ω(a)`.
* **(c) Height = distance to a canonical lift.** Put `κ(a) := ω(a)·σ_p(a)^(−2)`, the **Camion lift** of `a mod p`. Then
      h_a(p) = v_p(Λ_p(a)) = v_p(a − κ(a)),
  in parallel with the Wieferich order `v_p(a^(p−1) − 1) = v_p(a − ω(a))`.
  * `κ(a)` is the unique `p`-adic base `≡ a (mod p)` for which every power `p^k` is a Camion pseudoprime.
  * Likewise `ω(a)` is the unique such base for which every `p^k` is a Fermat pseudoprime.
* **(d) The quotient splits.** `ε_a D_a(p)/p ≡ ρ_C(a) − q_a(p)/2 (mod p)`, with the lift-independent **Camion residue** `ρ_C(a) := (σ_p(a) − 1)/p mod p`. So `p` is Camion–Wieferich iff `q_2(p) ≡ 2ρ_C(2) (mod p)`.

**Proof.**

* **(a)** In Session25's expansion `T^(−1) = Σ_k (−A^(−1)U)^k A^(−1)`, every grade-`≥1` term is a combination of `p^j δ_d ∗ E_β` with `d ≥ 1`.
  * The Newton coefficients satisfy `a_(j,0) = C(C(0,2), j) = 0` for `j ≥ 1`.
  * `C(n,d)t_0(n) = Σ_α c_α α^d (δ_d ∗ E_α)(n)`.
  * `δ_d ∗ δ_e = C(d+e,d) δ_(d+e)`.
  * Along `n = p^k`: for unit `β`, `ε^k p^j C(p^k,d) β^(p^k−d)` has valuation at least `k − v_p(d!) → ∞`. For nonunit `β` the valuation is at least `j + p^k − d → ∞`.
  * A sharper bound holds: `v_p(C(p^k,d)) = k − v_p(d)`, and `v_p(d) ≤ j − 1` because `d ≤ 2j < p^j`. So each `d ≥ 1` term has valuation at least `k + 1` at `n = p^k`. This is the finite-level statement (second audit, finding 1).
  * So only the grade-0 part `e − A^(−1)` survives. `A = ω(a)^C(n,2)` depends only on `a mod p`, and the image `e − A^(−1) = s_(ω(a))` is unique.
  * *Existence.* Session25 §7 gives the tower relation `ε^k s_a(p^k) ≡ ε^(k−1) s_a(p^(k−1)) (mod p^k)`, so the sequence is Cauchy.
  * *`σ_p ≡ 1`.* `σ_p ≡ ε s_a(p) (mod p²)` and `s_a(p) ≡ ε (mod p)`.
* **(b)** `c_a(p^k) = ω(a)^((p^k−1)/2) ⟨a⟩^((p^k−1)/2)`. Here `ω(a)^((p^k−1)/2) = ε_a^(1+p+…+p^(k−1)) = ε_a^k`, and `(p^k−1)/2 → −1/2` in `Z_p`.
* **(c)**
  * Both `σ_p` and `⟨a⟩^(−1/2)` are `≡ 1 (mod p)`, and for `x ≡ 1` and odd `p`, `v(x − 1) = v(x² − 1)`. So `v(σ_p − ⟨a⟩^(−1/2)) = v(σ_p²⟨a⟩ − 1) = v(a − ω(a)σ_p^(−2))`.
  * Every lift `a′` of `a mod p` has the same `σ_p` and `ω`, hence the same `κ`. So `h_(a′)(p) = v(a′ − κ(a))`, which is infinite only for `a′ = κ(a)`.
  * The initial-segment and freezing lemma of THM-4746 (9) holds verbatim for base `a`.
* **(d)** `⟨a⟩^(−1/2) ≡ 1 − p·u_a/2` with `u_a ≡ −q_a(p)`. The tower relation gives `εD(p) ≡ Λ_p (mod p²)`. ∎

**FINITE-EXACT** (`camion_lift.out`, `camion_lift_check.out`, `camion_lift_map.out`, `camion_residue_2000.out`, `camion_residue_fish.out`, `camion_residue_symmetry.out`).

* **Rigidity, by interpolation.** `σ_p(a)` is identical mod `p³` for the lifts `2, 2+p, 2+2p, 2+5p, 2−p+p²` at every prime `7 ≤ p ≤ 43`. The same holds for `a = 3, 5, p−1` at `p = 7, 11, 13`.
* **Rigidity, directly** (second audit c1). `s_a(p^k) mod p^(k+1)` agrees across 7 lifts per residue: 115 level checks, up to mod `3^10`, `5^7`, `7^6`, `11^5` and `13^5`.
* **Height formula.**
  * `v_p(Λ_p) = v_p(2 − κ(2))` at every tested prime: 2 at `p = 3, 13`, 1 otherwise.
  * Direct `Λ_p` from `D(p^k)` equals `σ_p(ω(2)) − ⟨2⟩^(−1/2)` at every level for 13 primes `3 ≤ p ≤ 43` (second audit c2). This is the genuine test; the session's interpolation cross-check of `D` against `s` is nearly tautological.
* **Designed and all bases** (second audit c2, c3, c3b).
  * For the base `κ(a) + 3p^h`, `v(D(p^k)) = k + 1` for `k < h` and `= h` from `k = h` on, for every `h` tested and residues `2, 3, 5, p−1`.
  * For every integer base below `5^4` and `7^3` (794 in all), the tower height equals `v_p(a − κ(a))`. The heights are distributed `(1 − 1/p)p^(−h+1)`.
* **Every `p^k` passes for the base `κ(2)`.** For the integer base `≡ κ(2) (mod p^K)`, `v_p(D_κ(p^k)) = 2, 3, 4` at `p = 11, 13`, and `2, 3, 3` at `p = 7` (capped at `K = 3`). The second audit gives `2, 3, 4, 5, 6` at `p = 7` with `K = 6`. Bases 2 and `2 + p` give 1, except base 2 at the Camion–Wieferich prime 13, which gives 2.
* **Decomposition (d).** It holds for all 301 primes `5 ≤ p ≤ 2000`.
  * `ρ_C(2) = 0` exactly at `p = 199, 809, 929`. There the two lifts agree mod `p²`, so Camion–Wieferich ⟺ Wieferich at those primes.
  * No relation `ρ_C ≡ αq_2 + βw_p + γB_(p−3) + δ` with small rationals survives past `p = 67`.
* **Structure of `κ` and `ρ_C`.** `κ` is not multiplicative. As a function on `F_p^×` (`p ≤ 19`), `ρ_C` has no symmetry under `a ↦ 1/a` or `a ↦ −a`, and its Lagrange degree is the generic `p − 2`.
* **Explicit formula** (second audit finding 4). `ρ_C(a) ≡ ε_a·[x^p] (Σ_(k<p) a^C(k,2) x^k/k!)^(−1) (mod p)`, the quantity `u_p` of `cw_fast.c`. It visibly depends only on `a mod p`. Checked at 301 primes.

**Reading.**
* Each residue `a mod p` has two canonical `p`-adic lifts:
  * the Teichmüller lift `ω(a)`, where every `p^k` is a Fermat pseudoprime;
  * the Camion lift `κ(a)`, where every `p^k` is a Camion pseudoprime.
* The Wieferich and Camion–Wieferich primes are where the integer 2 agrees with the respective lift to second order.
* Neither 1093 nor 3511 is Camion–Wieferich. Wieferich primes make the curvature of §3 vanish; Camion–Wieferich primes make `Λ_p` small.

**The Fermat parallel.**
* For the base-2 Fermat test the analogous limit is `lim 2^(p^a − 1) = ⟨2⟩^(−1)`.
* The height of `p^a` as a Fermat pseudoprime is `v_p(⟨2⟩ − 1) = 1 + v_p(q_2(p))`.
* So `Λ_p` plays the role for strong tournaments that `log_p⟨2⟩` plays for Fermat.

**FINITE-EXACT (search).**
* **Camion–Wieferich primes `p ≤ 10^5`.** They are exactly **3, 13, 22409**.
  * `cw_fast.c` uses `s(p) ≡ ε + p·u_p (mod p²)` with `u_p = [x^p] 1/T_(<p)(x)`.
  * The audit's independent quasi-linear FLINT search (b6) covered all 9590 primes `5 ≤ p ≤ 10^5`.
* **`h(22409) = 2`.** It was found by `tower_height.c` and confirmed by audit b5 through `D(22409) mod p³`. The root `l = 0` is simple.
* **The search to `10^6`** (`cw_flint_np.py`, an `O(p log p)` method; validated by the second audit, c6): IN PROGRESS. Results will appear in `cw_flint_*.out`.
* **Every prime `p ≤ 10^5` has `h(p) ≤ 2`.** The heuristic expected number with `h ≥ 3` is `Σ_(p≥5) 1/p² ≈ 0.09`.
* **Other structures.**
  * **Strong semicomplete digraphs (base 3).** The Camion–Wieferich primes `p ≤ 10^4` are 17 and 31 (two codes). §5b applies with `a = 3`, so `h_3(p) ≥ 2` exactly at those primes. The base-3 Wieferich prime 11 is not on the list.
  * **Connected labelled graphs.** The primes `p ≤ 10^4` with `p² | E(p)`, where `E(n) = cg(n) − 2^((n−1)/2) + 1` (rotation: `cg(p) ≡ 2^((p−1)/2) − 1`), are 5 and 197.
    * This is NOT a tower-height list (second audit finding 20).
    * The target has a constant term, so the tower drifts: `E(p^r) ≡ εE(p^(r−1)) + (1 − ε) (mod p^r)`.
    * For `(2/p) = −1`, which holds at both 5 and 197, `p²` never passes. The passing powers need not form an initial segment: at `p = 5`, 5 and 125 pass and 25 fails.
    * The rigidity of `b_a(p^k)` holds numerically but is unproved.

## 6. The prime 2: Moon parity and the third cyclotomic polynomial

Every tournament automorphism group has odd order (Moon): an involution would have to reverse the arc it swaps.
* Hence `f(n) := s(n)/n! = Σ_(strong classes) 1/|Aut|` lies in `Z_(2)`, and `f(n) ≡ #(strong classes on n vertices) (mod 2)`.
* Similarly `2^C(n,2)/n!` is the number of all tournament classes, mod 2.
* Structures whose automorphisms can be involutions, such as semicomplete digraphs or hypergraphs, have no such 2-adic face.

**Statement 6.**

* **(a) PROVED.** The number of strong tournaments on `n` vertices up to isomorphism is odd iff `n ≢ 2 (mod 3)` (`n ≥ 1`). The number of all tournaments up to isomorphism is even for every `n ≥ 3`.
  * *Proof.* `v_2(2^C(n,2)/n!) = C(n,2) − v_2(n!) ≥ (n−1)(n−2)/2 ≥ 1` for `n ≥ 3`, so `T(x) ≡ 1 + x + x² = Φ_3(x) (mod 2)` in `Z_(2)[[x]]`.
  * Then `S = 1 − 1/T ≡ 1 − (1+x)/(1+x³) (mod 2)`, whose coefficients are 1 exactly at `n ≢ 2 (mod 3)`, `n ≥ 1`. ∎
  * The OEIS data of A051337 and A000568 agree up to `n = 19`. OEIS records no parity comment on either sequence.
* **(b) PROVED in outline, with its range; FINITE-EXACT.**

      f(n + 3·2^j) ≡ f(n)  (mod 2^(j+2))   whenever v_2(r_n) ≥ j + 2,

  where `r_n` is the remainder below. In particular this holds for all `n ≥ 2(j−1)/3`, using the observed bound `v_2(r_n) ≥ 3n/2 + 3`.
  * *Failures below the range.* They occur for every `j ≥ 5`, and the exception set grows linearly in `j` (second audit c8, MISTAKE-655):
    * `n ≤ 2` for `j = 5, 6`;
    * `n ≤ 3` for `j = 7`;
    * `n ≤ 4` for `j = 8..10`;
    * `n ≤ 6` for `j = 11, 12`;
    * `n ≤ 8` for `j = 14..16`;
    * `n ≤ 10` for `j = 17, 18`.
  * The exact rule "congruence ⟺ `v_2(r_n) ≥ j + 2`" is FINITE-EXACT for `j ≤ 18` and `n ≤ 60`.
  * *Class `n ≡ 1 (mod 3)`.* The congruence holds mod `2^(j+3)` for `7 ≤ n ≤ 1500 − 3·2^j` and `j ≤ 8`. Its exceptions also grow: `n = 1` for `j ≥ 4`, `n = 4` for `j ≥ 7`, `n = 7` for `j ≥ 12`, `n = 10` for `j ≥ 16`.
  * *Reason.* `T` is entire on `C_2`, since `v_2(τ_n) ≈ n²/2`. By its Newton polygon its only zeros in the closed unit disc are two conjugate units `ρ, σ(ρ) ∈ Z_2[ζ_3]` with `ρ ≡ ζ_3`; all other zeros satisfy `|x| ≥ 2^(3/2)`.
  * Hence `f(n) = Tr_(Q_2(ζ_3)/Q_2)(C ρ^(−n)) + r_n`, with `C = 1/(ρT′(ρ))`.
  * The remainder satisfies `v_2(r_n) ≥ (3/2)n − O(1)` from the next zeros on `|x| = 2^(3/2)`. Observed exactly: `v_2(r_n) = 3n/2 + 3, 9/2, 3, 7/2` by `n mod 4`, for `n ≤ 120` (second audit c8). The Newton polygon has slopes `0, 3/2, 4, 5, 8, 19/2` (2 zeros each).
  * Writing `ρ = ζw`, `v_2(w³ − 1) = 2` gives the exponent `j + 2`.
  * *Checks.* `f(n) ≡ Tr(Cρ^(−n))` mod `2^72` at `n = 10^5, …, 3·10^6 + 5` (`two_adic_regulator.out`), and mod `2^100` for `n ∈ [1500, 2200]` (audit b9).
  * At `p = 2` the role of `F_p^×` (order `p − 1`) is played by `F_4^×` (order 3), because the reduction of the tournament EGF is the irreducible `Φ_3`.
* **(c) PROVED in outline; FINITE-EXACT.** On the zero class `n ≡ 2 (mod 3)`:

      v_2(f(n)) = 2 + v_2(n − n*),   n* = log t / log r,

  with `r = σ(w)/w` and `t = −ζσ(C)/C`, both in the norm-one torus of `Q_2(√−3)`.
  * *Proof.* On this class `Tr(Cρ^(−n)) = Cζw^(−n)(1 − t r^(−n))`, with `v(C) = 0` and `v(r − 1) = 2`.
    * `t ∈ r^(Z_2)`: the group `N^1 ∩ (1 + 4O)` is torsion-free of rank 1, `r` generates it (`v(r − 1) = 2` is minimal), and `v(t − 1) = 3`. Hence `v_2(n*) = 1`.
    * The valuation law holds whenever `2 + v_2(n − n*)` is below the decay rate of `r_n`.
  * *Values:*

        n* ≡ 16221783116972227762  (mod 2^64),   n* ≡ 3250 (mod 4096),
        n* ≡ 1090019148193376184792522869885760632192178  (mod 2^140)   (audit b7, Kitamasa evaluation).

  * The law holds exactly for all 499 values `5 ≤ n ≤ 1500` (`two_adic_zero_1500.out`, audit b7), and at 300 random huge `n` (audit b7).
  * The closed form `log t / log r` agrees with `n*` to 90 bits (audit b9). The check `r^N ≡ t (mod 2^66)` with `N ≡ n* (mod 2^64)` matches the predicted precision `64 + v(r − 1)` (`two_adic_regulator.out`).
  * **`n*` is not a small rational.** Rational reconstruction at `2^64` returns a generic fraction of maximal height, `765490658/1598160473`, which fails mod `2^80`. At `2^140` no valid reconstruction exists. So `n*` is not a rational with numerator and denominator below `2^69.5` (audit b7).
  * The universal class "2" of odd `p` (no strong K2) is the only 2-adic zero class, and it carries a single simple zero at the "virtual order" `n*`. The exact zero `s(2) = 0` lies below the analytic range.

## 6b. Degree = uniformity: r-uniform hypergraphs (PROVED; statistics FINITE-EXACT)

Let `r ≥ 2` and `T_r(n) = 2^C(n,r)`, and let `b_r` count connected labelled `r`-uniform hypergraphs: `Δb_r = ΔT_r ∗ T_r^(−1)`.

**Statement 7.** Let `p` be a prime with `gcd(ord_p(2), r!) = 1`.

* **(a)** `ω^C(n,r)` has period dividing `p − 1`. By Vandermonde, `C(n+P, r) − C(n, r) = Σ_(i≥1) C(P,i) C(n, r−i)`, and `ord_p(2) | C(P,i)` when `ord_p(2) | P` and `gcd(ord_p(2), i!) = 1`.
* **(b)** Session25's argument therefore works over `Z/p^k` itself, with all nodes in `Z/p^k` and multiplier 1. The filtration now has degree `d ≤ rj` at grade `j`.
* **(c)** For `m ≥ r(k−1) + 2`, the class sequence `l ↦ b_r(m + l(p−1))` has Mahler coefficients `a_i ≡ 0 (mod p^min(k, ⌈i/r⌉))`. So the Mahler speed is `1/r`.
* **(d)** Mod `p²` this is a polynomial of degree `≤ r` in `l`. Its coefficients of `C(l,2), …, C(l,r)` are `p·u·(…)` with `u ≡ −q_2(p)`, so they vanish mod `p²` at Wieferich primes.
* **(e)** Class 1 is a universal zero class. For example `b_3(p) ≡ 2^((p−1)(p−2)/6) − 1 ≡ 0 (mod p)`: rotations act freely except on unions of the free orbits of 3-sets, and every nonempty union is connected.

*Proof.* Repeat §1 with `d ≤ rj`.
* Unit nodes give valuation at least `j + max(0, i − d) ≥ max(j, i − (r−1)j) ≥ ⌈i/r⌉`.
* Nonunit nodes vanish once `n − (r−1)j ≥ k`, i.e. for `n ≥ r(k−1) + 1`. One more index is needed for the shift in `Δb_r`.
* Mod `p²` only grades `j ≤ 1` survive. The grade-0 part is linear in `l`, and the grade-1 part is `pu·C(n,r)` times grade-0 terms. ∎

**Root statistics** (FINITE-EXACT). The data cover the 60 primes `p ≤ 2200` with `gcd(ord_p(2), 6) = 1`, 60444 classes in each case. Root counts are compared with a random degree-`r` polynomial over `F_p`, whose root count is the number of fixed points of a uniform permutation in `S_r`.

* **`r = 3`** (`hyper3_second_order_*.out`; reproduced by audit b8).
  * The cubic law holds in every class: `p | a_1, a_2, a_3` and `p² | a_4`. For `p ≤ 127` it was also checked directly (33166 checks).
  * 58 classes are not genuinely cubic.
  * Observed thresholds are lower than proved: at `k = 2` the law fails only at the bases `m ∈ {1, 3}`, and at `k = 3` only at `m ≤ 6` (audit b8).
* **`r = 4`** (`hyper4_second_order_5_2200.out`, `hyper_r_second_order.c`).
  * The quartic law holds in all 60444 classes.
  * 50 classes are not genuinely quartic.

| `r` | zero class | observed root counts `0/1/2/3/4` | prediction |
|---|---|---|---|
| 3 | 1 | 22 / 25 / 1 / 12 | random cubic 20 / 30 / 0 / 10 |
| 3 | 2 (exact zero `b_3(2) = 0`) | 0 / 31 / 0 / 29 | forced root + random quadratic 0 / 30 / 0 / 30 |
| 3 | sporadic | 21 / 21 / 0 / 14 | 18.7 / 28 / 0 / 9.3 |
| 4 | 1 | 23 / 19 / 16 / 0 / 2 | random quartic 22.5 / 20 / 15 / 0 / 2.5 |
| 4 | 2 (`b_4(2) = 0`) | 0 / 18 / 32 / 1 / 9 | forced + random cubic 0 / 20 / 30 / 0 / 10 |
| 4 | 3 (`b_4(3) = 0`) | 0 / 21 / 26 / 0 / 13 | forced + random cubic 0 / 20 / 30 / 0 / 10 |
| 4 | sporadic | 21 / 19 / 13 / 0 / 3 | 21 / 18.7 / 14 / 0 / 2.3 |

* **The forced roots are FINITE-EXACT, not a consequence of the theorem.** Their mechanism is that the class polynomials already hold at the exact zeros `n = 2` (and `n = 3` for `r = 4`), below the proved threshold.
  * For `r = 3` the class-2 cubic vanishes at `n = 2` in 60 of 60 cases, and base `m = 2` obeys the law at all 60 primes (audit b8).
  * The extrapolation check `hyper3_class2_check.out` reproduces `b_3(n)` at `n = 2, 4, 5` but not at `n = 1, 3`.
* **The `S_r` law (HYP-9361 (d)).** For tournaments (`r = 2`) this is the split/inert 1:1 law, with the triangle forced split. For `r = 3, 4` the second-order zeros behave like the roots of a random degree-`r` polynomial over `F_p`.

## 7. What this does and does not say

* **Says.**
  * Camion pseudoprimality is the mod-`p` shadow of a `p`-adic object: the class functions `F_m` on `Z_p` (restricted power series for `p ≥ 5`). They run at half the Kummer rate and are twisted by `(2/p)`.
  * Their second-order structure is governed by the Fermat quotient through the quadratic edge count.
  * Their zero divisor has degree 2, the degree of `C(n,2)`. In general the degree is `r`, with `S_r` statistics.
  * The tower heights `h(p)` are the valuations of `Λ_p`, the class-1 branch at `n = 0`, equivalently the distance from 2 to its Camion lift.
  * **Reading.** The universal residues of THM-4746 can be seen as having `p`-adic avatars. Each sub-claim is proved; the word "avatar" and the `K2` identification across `p = 2` and odd `p` (through `F_4^×` in place of `F_p^×`) are interpretive.
    * K1 (class 1) carries the Camion limit `Λ_p` and the tower heights.
    * The absent strong K2 (class 2) is the only 2-adic zero class, with its simple zero `n*`.
    * The triangle (class 3) is the one exact zero with `n ≥ 3`. For `p ≥ 5` it forces a split second-order zero at every non-Wieferich prime where the root is simple.
  * The prime 2 behaves differently. Moon's odd-order theorem makes the groupoid counts `Σ 1/|Aut|` the natural 2-adic object, and `Φ_3` replaces `F_p^×`.
* **Does not say.**
  * Nothing about the infinitude or density of Camion–Wieferich primes (HYP-9358 (a)).
  * Nothing about the 1:1 split/inert law, or the `S_r` statistics, beyond data (HYP-9361 (a), (d)).
  * Nothing about the repunit rarity in the `H`-spectrum (opus S51's lane, THM-4748/4749).
  * The "virtual order" `n*` has a closed form, but no combinatorial meaning.
