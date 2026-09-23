# The Collatz trunk `(4^i-1)/3` and the Riemann hypothesis: one exact Iwasawa coordinate, four different "1/2"s, no bridge

**Status. RH is OPEN and nothing below bears on it.**

* **PROVED** (hand proofs below; each also checked by machine, cited results marked):
  1. `T(i) = (4^i-1)/3` is a bijective isometry of `Z_3`. Its fixed points are exactly `{0, 1, -1/2}` (by Strassmann's theorem). It takes the values `T(1/2) = -1`, `T(n/2) = ((-2)^n-1)/3` and `T(i*) = 1/2`, where `i* = log(5/2)/log 4` is irrational. The minus-sheet analogues are §1.3.
  2. In logarithmic charts, every reverse move of the backward E-game is a translation. The foundry's kappa formula is the trunk isometry, since `kappa(3m+1) = 2 i(m)`. The nonlinearity sits entirely in a 3-fold expanding re-chart map and in the archimedean price.
  3. In the trunk chart, `s -> 1-s` is the Möbius involution `M(m) = (1-m)/(3m+1)`. Its fixed point in `Z_3` is `-1`, and it does not preserve the hostile set `{1, 1/2}`.
  4. Using the cited Kubota–Leopoldt construction: `zeta_3(1-s) = G(3T(s)) / (3T(s))`, with `G` a **unit** of `Z_3[[X]]`. Hence `zeta_3` has **no zeros** on its weight disc, and `|(s-1) zeta_3(s)|_3 = 3` on `Z_3`.
  5. The trunk Dirichlet series `sum T_i^-s` has no zeros in `Re s > 0.459366`, so none on the critical line.
  6. Euler-factor neutrality: every trunk identity with `zeta` multiplies it by factors whose zeros lie on `Re s = 0` or `Re s = 1`.
  7. The trunk is an Iwasawa coordinate for every non-Wieferich `qx+1` map as well, so the identity is blind to DRIFT and to SHEET.
* **AUDIT (orchestrator, 2026-09-23).** Theorem 1.2's Strassmann count was re-derived: `a_1 = lambda/3 - 1 = 3/2 + O(27)`, so `v_3 = 1`; `v_3(a_2) = v_3(a_3) = 1`; `v_3(a_n) >= (n-1)/2` for `n >= 4`, so there are at most 3 zeros, and `0, 1, -1/2` are zeros. The Jacobsthal identification `J_2n = T(n)`, `J_(2n+1) = S(n)` was also re-checked. No gap found.
* **FINITE-EXACT.** Exact 3-adic and rational checks (sections P1–P12 of the `.out`), and exact counts `N_k`. Items marked (fp) are floating-point numerics with a validated evaluator, not interval arithmetic:
  * the Davenport–Heilbronn and Epstein controls (functional equations; off-line zeros counted by the argument principle);
  * the zeros of the trunk series;
  * the dynamical zeta of the canonical escape `Psi`;
  * the data tests T1–T3, over `m <= 10^7` (T1) and `m <= 10^6` (T2, T3).
* **CITED:** the primary sources listed in §8, each marked read (abstract, or the named sections) or UNVERIFIED. The Kubota–Leopoldt statements come from Rodrigues Jacinto–Williams (sections read).
* **OPEN:** RH; E-SCC Q1 and Q2; the hypothesis candidates of §6.
* **REFUTED:**
  * "Under the trunk/Iwasawa map the hostile point `1/2` sits at the critical point `1/2`." In fact `m = 1/2` corresponds to `s = i*`, and `s = 1/2` corresponds to `m = -1`.
  * "The trunk series has zeros on `Re s = 1/2`": it is provably zero-free there.
  * "E-SCC or Collatz Dirichlet series satisfy a degree-1 or degree-2 functional equation" (fp; conductor `<= 300`, resp. level `<= 150` and weight `<= 12`).
  * "The E-game's own zeta satisfies an Ihara-type Riemann hypothesis" (fp, two weights). No one had claimed it; it was the natural test.

Session `collatz-procgen-20260922` (mac-mini), trunk/RH lane, 2026-09-23.

* **Scripts:** `04-computation/experiments/collatz_procgen_20260923_trunk_rh_{padic,dirichlet,data}.py`.
* **Output:** [collatz_procgen_20260923_trunk_rh.out](collatz_procgen_20260923_trunk_rh.out). Its sections are P1–P12 (3-adic), D0–D5 (archimedean), and C0 plus T1–T3 (data).
* **Inputs:**
  * [synthesis](collatz_procgen_20260922_synthesis.md);
  * [choice ladder](collatz_procgen_20260922_choice_ladder.md) §5 (the trunk as the exit set of `1/2`);
  * [Q2 endgame](collatz_procgen_20260922_q2_endgame.md): Lemma 1.1, `Psi`, §3.4;
  * [transversality foundry](collatz_procgen_20260922_transversality_foundry.md) §3 (the kappa formula);
  * MISTAKE-497, for the character-branch typing of `p`-adic zeta values.

## 0. The answer in brief

The trunk is related to zeta functions in two exact ways. Neither one reaches the Riemann hypothesis.

1. **3-adically (REAL).** `u = 4 = 1+3` is the standard topological generator of `1+3Z_3`, and `X = u^s - 1` is Iwasawa's coordinate on the disc of weights `x -> <x>^s`. So the trunk is literally `X/3`. Moreover `3T(s)` is the pole factor of the Kubota–Leopoldt function: `zeta_3(1-s) = G(3T(s))/(3T(s))`. Here `G` is a unit power series, so `zeta_3` has no zeros anywhere, and there is no 3-adic "critical line" to talk about.
2. **Archimedeanly (REAL).** `1 - 4^i = -3T(i)` is the factor at the prime 2 that turns `zeta(1-2i)` into `eta(1-2i)` and `B_2i` into the Genocchi number `G_2i = -6T(i)B_2i`. Its complex zeros lie on `Re s = 1`.

The E-game's hostile point `1/2` and RH's `1/2` are different objects. The first is a point of the 3-adic phase space `Z_3`, namely the unit with `2x = 1`. The second is the real part of a spectral parameter. The only exact bridge, the trunk chart, does not match them:

* it sends `m = 1/2` to `s = i* = log(5/2)/log 4`, which is irrational;
* it sends `s = 1/2` to `m = -1`. This is because `4^(1/2) = -2` in `1+3Z_3`, and `-1` happens to be the hostile point of the *forward* game and of the minus sheet.

No candidate series behaves like an L-function on the critical line:

* **The trunk series** is provably zero-free for `Re s > 0.4594`.
* **The E-SCC descent statistics and Collatz stopping times** give Dirichlet series with no functional equation: best residual `>= 9*10^-3`, against `<= 5*10^-11` for genuine L-functions.
* **The fluctuation exponents** of those statistics are halved exceptional dimensions, about `0.25` for `Psi` descent and between `0.40` and `0.50` for stopping times. The Möbius function and shuffled data give `0.50`.

The one place where "zeros on a line" and the hostile point `1/2` genuinely meet is the E-game's own dynamical zeta function. That is question R in §5.

## 1. The trunk as a 3-adic object (PROVED)

**Setting.** Put `lambda = log_3 4`. Then `v_3(lambda) = 1` and `v_3(lambda - 3) = 2`. For `s` in `Z_3`, define `4^s = exp(s lambda)`, an element of `1+3Z_3`. For `p` odd, `log` and `exp` are mutually inverse isometries between `1+pZ_p` and `pZ_p`.

**Theorem 1.1 (isometry).** `T(s) = (4^s-1)/3` is a bijection `Z_3 -> Z_3` with `v_3(T(s)-T(s')) = v_3(s-s')`. Its inverse is the **trunk chart** `i(m) = log(3m+1)/log 4`.

*Proof.* `T(s) - T(s') = 4^(s')(4^(s-s') - 1)/3`, and `v_3(4^x - 1) = v_3(x lambda) = 1 + v_3(x)`. The map `s -> s lambda` sends `Z_3` onto `3Z_3`, and `exp` sends `3Z_3` onto `1+3Z_3`. ∎

*Checks (P2).* 22,800 pairs, including `1/2, -1/2, 1/4, 1/5, 5/7`; round trips `i(T(s)) = s` on 310 samples; agreement with the integers `(4^i-1)/3` for `i < 200`.

In particular `v_3(T(i)) = v_3(i)`, so `T(i)` is a unit iff `3` does not divide `i`. This is exactly the legality rule for the exits of `1/2` (choice ladder §5): the move `k = 2i+1` is legal iff `i` is not `0 mod 3`.

**Theorem 1.2 (fixed points and special values).**
* `Fix(T) = {0, 1, -1/2}` exactly.
* `T(n/2) = ((-2)^n - 1)/3` for every integer `n`; in particular `T(1/2) = -1` and `T(-1) = -1/4`.
* `T(i*) = 1/2` for `i* = log(5/2)/log 4`. This number is irrational and satisfies `v_3(i* - 1/2) = 1`.

*Proof.*
1. Write `f(i) = T(i) - i = sum_n a_n i^n`, with `a_1 = lambda/3 - 1` and `a_n = lambda^n/(3 n!)` for `n >= 2`.
2. Then `v_3(a_1) = v_3(a_2) = v_3(a_3) = 1`. For `n >= 4`, `v_3(a_n) = n - 1 - v_3(n!) >= (n-1)/2 > 1` by Legendre's formula.
3. By Strassmann's theorem `f` has at most 3 zeros in `Z_3`.
4. The points `0`, `1` and `-1/2` are zeros. For `-1/2`, use `4^(-1/2) = (-2)^(-1)`, since `-2` is the square root of `4` lying in `1+3Z_3`. The same fact gives the formula for `T(n/2)`.
5. If `i* = a/b` with `b != 0`, then `(5/2)^b = 4^a`, which is impossible.
6. `v_3(i* - 1/2) = v_3(T(i*) - T(1/2)) = v_3(3/2) = 1`. ∎

*Checks (P3).*
* The count `#{i mod 3^r : T(i) = i mod 3^r}` is `3, 6, 12, 21, 21, ..., 21` for `r = 1..9`.
* `21 = 3^1 + 3^2 + 3^2`. This matches `v_3(T'(z) - 1) = 1, 2, 2` at the three fixed points. (The lead listed only `1` and `-1/2`; `0` is the third fixed point.)

| `s` | `T(s)` | role |
|---|---|---|
| `0` | `0` | the pole of `zeta_3` (§2); not a unit |
| `1` | `1` | Q2 hostile point `1` (fixed) |
| `-1/2` | `-1/2` | minus-sheet hostile point (fixed) |
| `1/2` | `-1` | hostile point of the forward game Q1 and of the minus sheet |
| `n + 1/2` | `-(2^(2n+1)+1)/3` | the exits of `-1` (the negated minus trunk; P12) |
| `i*`, `= 101220102010212022112001200102` mod `3^30` | `1/2` | Q2 hostile point `1/2` |

**1.3 The minus sheet (PROVED).**
* Negation conjugates `E` to `E_-`, whose arrows are `n -> n/2` and `n -> 3n-1`. The minus-sheet hostile points are `-1` and `-1/2`. These are the classes `26, 13 mod 27` excluded by the Lean mirror lemma.
* Their exit set is `-T(n) = (1-4^n)/3`, with chart `m -> i(-m)`. The conjugate chart `T~(i) = -T(-i)` fixes `{0, -1, 1/2}` (P4).
* The proper minus trunk is the set of positive integers `y` with `3y-1` a power of 2, namely the `E_-` reverse images of the root 1. It is `S(n) = (2^(2n+1)+1)/3 = 2T(n)+1`, which gives `1, 3, 11, 43, 171, ...`.
  * `S` is an isometry with exactly one fixed point. Its Strassmann valuations are `0, 0, 1, 1, 2, ...`, and the fixed point satisfies `i_S = 2 mod 3`.
  * `S(-1) = 1/2`, `S(1/2) = -1` and `S(-1/2) = 0`.
* **The two trunks are one object.** They interleave into the Jacobsthal numbers `J_t = (2^t - (-1)^t)/3`: `J_2n = T(n)` and `J_(2n+1) = S(n)`. Also `((-2)^t - 1)/3 = (-1)^t J_t`. So both trunks are the lattice points of the single chart `t(m) = log(3m+1)/log(-2) = 2 i(m)`, the chart of the kappa formula.
* Since `4^(1/2) = -2`, `T(n + 1/2) = -S(n)`: the 3-adic interpolation of the plus trunk passes through the negated minus trunk at the half-integers.

**1.4 The backward E-game in logarithmic charts (PROVED; checked on 5000 moves, P5).** For a 3-adic unit `x`:
* let `eps(x)` be `0` if `x = 1 mod 3` and `1` if `x = 2 mod 3`;
* let `h(x) = 2^(-eps(x))`, which is `1` or `1/2`;
* the **departure chart** is `delta(x) = log(2^eps(x) x)/log 4`;
* the **arrival chart** is `i(y) = log(3y+1)/log 4`.

**Proposition 1.4.**
* **(a) Charts of the hostile points.** `delta(1) = delta(1/2) = 0`, `delta(-1) = 1/2` and `delta(-1/2) = -1/2`. The precision of Lemma 1.1 of the endgame lane is `k(x) = v_3(x - h(x)) = 1 + v_3(delta(x))`.
* **(b) Translation law.** The legal moves at `x` are `k = 2n + eps(x)` with `n >= 0`. Such a move lands at `y = (2^k x - 1)/3` with `i(y) = n + delta(x)`, because `3y+1 = 4^n 2^eps x`. The landing is legal iff `3` does not divide `n + delta(x)`.
* **(c) Landing law.** For every `h` in `Z_3`, `v_3(y - h) = v_3(n + delta(x) - i(h))`. The targets are `i(1) = 1` and `i(1/2) = i*`.
* **(d) The kappa formula is (c).** The formula `v_3(2^K - r) = [K = e mod 2](1 + v_3(K - kappa(r)))` is (c) read in the doubled chart, since `kappa(3m+1) = 2 i(m)` (checked on 16,000 and 200 cases).
* **(e) Re-chart.** `delta(y) = Theta(i(y))` with `Theta = delta o T`. Here `Theta' = 4^i/(4^i-1)`, so `Theta` multiplies 3-adic distances by 3. Its zeros are exactly `i = 1` and `i = i*`.

The game therefore reads `i_(t+1) = n_t + Theta(i_t)`. The shift `n_t >= 0` is free, the move is legal iff `3` does not divide `i_(t+1)`, and it costs the real price `2^(2n_t + eps_t)/3`.

**Answer to the coordinator's question.** Yes, the kappa formula is linear in the trunk chart: every move is a translation, and every landing precision is a 3-adic distance between chart values. But this is only the logarithm of multiplication by powers of 2. The "+1" of `3n+1` is moved into the expanding map `Theta`, and hostility is decided by the archimedean price, which the chart cannot see. The linearization relocates the endgame problem of the Q2 lane; it neither removes nor weakens it.

**1.5 The involution `s -> 1-s` (PROVED; P6).**
* `T(1-s) = M(T(s))` with `M(m) = (1-m)/(3m+1)`, an involution.
* In the variable `X = 3T`, `M` is `X -> u(1+X)^(-1) - 1` with `u = 4`, the Iwasawa-algebra involution `gamma -> u gamma^(-1)`.
* The fixed points of `M` are `-1 = T(1/2)` and `1/3`, and `1/3` is not in `Z_3`.
* `M(1) = 0` (not a unit), `M(1/2) = 1/5` (precision 1, a generic point), `M(0) = 1` and `M(-1/2) = -3`.

So the hostile set `{1, 1/2}` is not `M`-invariant. Nothing in the E-game corresponds to the symmetry of the functional equation.

## 2. The Iwasawa / Kubota–Leopoldt connection

**2.1 What is cited.** [RJW] is J. Rodrigues Jacinto and C. Williams, *An introduction to p-adic L-functions*, Ess. Number Th. 4 (2025) 101–216, arXiv:2309.15692. I read the statements used:
* **Thm 4.1:** there is a unique pseudo-measure `zeta_p` on `Z_p^x` with `int x^k zeta_p = (1-p^(k-1)) zeta(1-k)` for all `k > 0`.
* **Prop 4.4:** `mu_a` has Mahler transform `F_a(T) = 1/T - a/((1+T)^a - 1)`, which lies in `Z_p[[T]]`.
* **Prop 4.8:** `int_(Z_p^x) x^k mu_a = (-1)^k (1-p^k)(1-a^(k+1)) zeta(-k)`.
* **Lemma 3.36(i):** a measure on `Z_p^x` all of whose moments `int x^k` (`k > 0`) vanish is zero.
* **Def 4.10 and Prop 4.11:** `zeta_p = x^(-1) Res mu_a / ([a] - [1])`.
* **Remark 3.47:** evaluation at a topological generator of `1+pZ_p` identifies `Hom_cts(1+pZ_p, C_p^x)` with the open unit disc.
* **Thm 5.1:** `int chi(x) x^k zeta_p = L(chi, 1-k)` for primitive `chi` of conductor `p^n`, `n >= 1`.
* **Def 5.16, Thm 5.17, Remark 5.21:** the branches `zeta_(p,i)(s) = int omega(x)^i <x>^(1-s) zeta_p = L_p(omega^i, s)` satisfy `zeta_(p,i)(1-k) = (1 - omega^(i-k)(p) p^(k-1)) L(omega^(i-k), 1-k)`, and `zeta_(p,i)` vanishes identically for odd `i`.
* **Thm 7.1(ii):** `zeta_(p,p-1)` has a simple pole at `s = 1` with residue `1 - 1/p`.
* **§11:** Iwasawa's theorem on the zeros of `zeta_p`. Used for context only.

**Normalization.** Set `zeta_3(1-s) = int_(Z_3^x) <x>^s zeta_p` for `s` in `Z_3 \ {0}`, the trivial-branch disc. This is RJW's branch `zeta_(3,2)(1-s) = L_3(omega^2, 1-s) = L_3(1, 1-s)`, since `omega^2 = 1` when `p = 3` (Def 5.16, Remark 5.21). Following MISTAKE-497, the branch is typed explicitly. `zeta_3(1-k)` equals `(1-3^(k-1)) zeta(1-k)` for even `k`. For odd `k` it equals `L(omega, 1-k)`, where `omega` is the character mod 3 (by Thm 5.1). For example `zeta_3(-1) = 1/6`, `zeta_3(0) = 1/3` and `zeta_3(-2) = -2/9`.

**Proposition 2.1 (the trunk is the pole factor; PROVED from the cited results).** For `s` in `Z_3 \ {0}`,

    zeta_3(1-s) = G(X)/X,   X = 4^s - 1 = 3 T(s),   G(X) = int_(Z_3^x) (1+X)^(l(x)) x^(-1) mu_4(x),

where `l(x) = log<x>/log 4` lies in `Z_3`, and `G` has coefficients in `Z_3`.

*Proof.*
1. By Prop 4.8 and `zeta(1-k) = -B_k/k`, the measure `nu = x^(-1) Res mu_4` has moments `int x^k nu = (4^k-1)(1-3^(k-1)) zeta(1-k)` for even `k`, and `0` for odd `k`.
2. These are the moments of `([4]-[1]) zeta_p`, so the two measures are equal by Lemma 3.36(i).
3. Integrate the character `<x>^s`. For `s != 0`, `[4]-[1]` contributes `4^s - 1`, which is nonzero.
4. Finally `<x>^s = 4^(s l(x)) = sum_n binom(l(x), n) X^n`, and each `binom(l(x), n)` lies in `Z_3`. ∎

*Checks (P7, P8).*
* The script realizes `mu_4` as Mazur's regularized Bernoulli measure `E_(1,4)`. Its Riemann sums reproduce the moments to 3-adic precision `>= 9` for `k <= 12`.
* The coefficients `c_0..c_9` were computed mod `3^(10 - v_3(n!))`.
* `G(4^k-1)/(4^k-1)` reproduces `zeta_3(1-k)` mod `3^6` for `k = 1..12`, both parities.

**Theorem 2.2 (`zeta_3` has no zeros; PROVED).** `G(0) = -(2/3) log_3 4` is a 3-adic unit (`= 1 mod 3`). The value follows from RJW Thm 7.1(ii), since `(4^s - 1) ~ s log 4` and the residue at `s = 1` is `1 - 1/3`. The unit property is also proved directly below. Hence:
* `G` is a unit of `Z_3[[X]]`, and its `mu`- and `lambda`-invariants are `0`;
* `G(X) != 0` on the whole open unit disc of `C_3`;
* the trivial branch of `zeta_3` has no zeros, and `|(s-1) zeta_3(s)|_3 = 3` for all `s` in `Z_3`.

*Proof.*
1. For even `k >= 2` we have `(3-1) | k`, so von Staudt–Clausen gives `3B_k = -1 mod 3`.
2. Hence `v_3(zeta_3(1-k)) = -1 - v_3(k)`, while `v_3(4^k-1) = 1 + v_3(k)`.
3. So `G(4^k-1)` is a unit for every even `k`.
4. Since `G(X) = G(0) mod 3` for `X` in `3Z_3`, `G(0)` is a unit. ∎

*Checks (P8, P9).*
* `c_0 = 33664 mod 3^10`, which equals `-(2/3) log 4` there, as Thm 7.1(ii) predicts.
* `v_3((1-3^(k-1))B_k) = -1` for all even `k <= 300`.
* `v_3(B_(k,omega)) = -1` for all odd `k <= 119`.

*Remarks.*
* **Only one branch.** For `p = 3` the trivial branch is the only even branch, and the odd branch `zeta_(3,1)` vanishes identically (RJW Thm 5.17). So the Kubota–Leopoldt `zeta_3` has no zeros anywhere, and "3 is regular" is automatic rather than an input. The analogous statement for every odd `p` (the trivial branch is zero-free, and zeros live only on the branches `omega^(2j)` of irregular primes) is standard background; I did not read it here, so it is UNVERIFIED and not used.
* **Values at the two points of interest.**
  * At `s = 1/2`: `X = 3T(1/2) = -3`, so `zeta_3(1/2) = G(-3)/(-3)`. Here `G(-3) = 1789 mod 3^7` is a unit, so `|zeta_3(1/2)|_3 = 3`.
  * At the E-game point `s = i*`: `X = 3/2`, so `zeta_3(1-i*) = (2/3) G(3/2)`, with `G(3/2) = 826 mod 3^7`. It is nonzero.

**2.3 What the connection says.** Everything in this list is REAL.
* The trunk is Iwasawa's normalized coordinate `X/3` for `u = 1+3`, on the closed disc `|X| <= 1/3` of the weights `x -> <x>^s` (RJW Remark 3.47).
* The coincidence behind it is `3*1 + 1 = 4 = 1 + 3`: the Collatz map sends the root `1` to the canonical generator.
* `3T(s)` is the pole factor of `zeta_3`, and it is also Mazur's regularizer `c^s - 1` for `c = 4`. The product `3T(s) zeta_3(1-s)` is a unit power series in `3T(s)`.

**2.4 What the connection does not say.**
* **Nothing about the zeros of `zeta(s)`.** The connection is 3-adic, `zeta_3` is zero-free, and `C_3` has no "real part", hence no critical line.
* **Nothing Collatz-specific (P11).** For every prime `q` with `2^(ord_q 2)` not `1 mod q^2`, the `qx+1` trunk `{(u_q^i - 1)/q}`, with `u_q = 2^(ord_q 2)`, is a bijective isometry of `Z_q` and Iwasawa's coordinate for the generator `u_q`.
  * For the Mersenne primes `q = 3, 7, 31, 127`, `u_q = 1+q` exactly. Among these is `7x+1`, the repository's DRIFT control, which is expected to have divergent orbits.
  * The minus sheet uses the negated trunk.
  * Only Wieferich primes (1093, 3511) break the identity; for them the trunk lands in `qZ_q`.
  * So the identity is DRIFT-blind and SHEET-blind in the sense of the foundry.
* **Nothing about the point 1/2 of the E-game.** It sits at `X = 3/2`, where `G` is a unit, which is an unremarkable point.
* **The choice of generator is a convention.** Any `u' = 4^a` with `a` in `Z_3^x` gives an equivalent coordinate, and `-2 = 4^(1/2)` gives the Jacobsthal/kappa chart.

**2.5 The archimedean shadow (PROVED; exact checks P10, D2).**
* `eta(s) = (1-2^(1-s)) zeta(s)`. At `s = 1-2i` this gives `eta(1-2i) = (1-4^i) zeta(1-2i) = -3T(i) zeta(1-2i)`.
* The Genocchi numbers satisfy `G_2i = 2(1-4^i)B_2i = -6T(i)B_2i`; checked for `i <= 30` from `2t/(e^t+1)`.
* On `Z_3` the `c = 2` Mazur measure has moments `(1-2^k)B_k/k = G_k/(2k)`.

These are values at negative integers, where `zeta(1-2i) != 0`. The factor `1 - 2^(1-s)` vanishes exactly at `s = 1 + 2 pi i k/log 2`, on `Re s = 1`. That is also where the complex trunk `T(z) = (4^z-1)/3` vanishes under `s = 1-2z`: numerically `|eta| ~ 10^-30` while `|zeta| = 1.35..2.34` for `k = 1..5`.

So the trunk is the Euler-type factor at the prime 2, sampled on a lattice. Any vertical line can be moved to `Re s = 1/2` by an affine change of variable, so "its zeros lie on a line" has no content.

**Riesz-type criteria** (RH in terms of `1/zeta(2k)`; UNVERIFIED, not read) see this factor only as `1/(1-2^(1-2k))`, which tends to 1.

## 3. An RH foundry: approaches, controls, and where each trunk bridge fails

**3.1 Controls.** A valid proof mechanism must not apply to the controls on which RH fails while the listed structure is present: Davenport–Heilbronn, Epstein, Beurling systems, and the de Bruijn–Newman deformation `H_t` for `t < 0`. It should also be compatible with the proved function-field case.

| control | what it has | RH there | an argument that works here would prove something false unless it uses... | source |
|---|---|---|---|---|
| **Davenport–Heilbronn** `f = ((1-i kappa)/2) L(s,chi) + ((1+i kappa)/2) L(s,chi-bar)`, `chi mod 5`, `chi(2) = i` | the functional equation `Lambda(s) = Lambda(1-s)` of degree 1 and conductor 5; a Dirichlet series with bounded coefficients; **no Euler product** | **false**. (fp) Exactly 4 zeros in `0.52 < Re s < 1.6`, `t <= 200`: `0.808517+85.699348i`, `0.650830+114.163343i`, `0.574356+166.479306i`, `0.724258+176.702461i`, plus their mirrors. `kappa = 0.284079043840412` from the root number, equal to the closed form; functional equation checked to `10^-16` | multiplicativity | D3; Balanzario–Sánchez-Ortiz 2007 (abstract read) |
| **Epstein** `sum' (m^2+5n^2)^-s`, `h(-20) = 2` | functional equation of degree 2; nonnegative coefficients; a modular theta function; no Euler product | **false**. (fp) 60 zeros in `0.52 < Re s < 1.6`, `t <= 300`; the first is `0.932970+15.668250i`. The decomposition `zeta L(chi_-20) + L(chi_-4) L(chi_5)` agrees with the lattice sum to `10^-7`; functional equation checked to `10^-20` | Euler product (positivity plus modularity is not enough) | D4; Lamzouri 2021; Rezvyakova 2024 (positive proportion on the line) (abstracts read) |
| **Beurling generalized primes** | an Euler product by construction; integers with `N(x) = ax + O(x^beta)` | **false in general**: the de la Vallée Poussin error term is optimal (DMV 2006) | the additive structure of `Z` | BDR 2025, introduction (read); DMV not read |
| **function fields** | Frobenius, cohomology, the index theorem | **PROVED** (Hasse `g = 1`; Weil; Deligne) | – (a positive control: a mechanism should have an analogue here) | Bombieri, Clay description (read) |
| **de Bruijn–Newman** | the heat flow `H_t` | RH iff `Lambda <= 0`. Rodgers–Tao: `Lambda >= 0`. Polymath: `Lambda <= 0.22` | no slack: any argument stable under `H_t` for small `t < 0` fails | abstracts read |
| **Nyman–Beurling / Báez-Duarte** | the `L^2` closure of the dilations `rho(1/(ax))`, `a` in `N` | equivalent to RH | – (a reformulation; it uses every `a` in `N`) | Báez-Duarte (abstract read) |
| **Li / Weil positivity** | positivity of `lambda_n`, respectively of the explicit-formula functional | equivalent to RH. For any multiset symmetric under `rho -> 1-conj(rho)`, Li positivity is equivalent to all points lying on the line, so the criterion is arithmetic-blind | the function-field proof of positivity uses the index theorem | Lagarias 2007 (abstract, §1); Bombieri |
| **random matrices** | Montgomery pair correlation, Odlyzko, Katz–Sarnak | evidence; an input to Rodgers–Tao | a spectral model must reproduce GUE statistics | Bombieri |
| numerics | – | RH verified to height `3*10^12` | – | Platt–Trudgian 2021 (abstract) |
| Goss zeta (characteristic `p`) | – | an "RH" (zeros in `F_q((1/T))`) is PROVED | shows that non-archimedean RH analogues exist but have a different content | Sheats 1998; Kramer-Miller–Upton (abstracts) |

**In-house controls.** These are the cheap tests every trunk bridge must pass.
* **FIN (Euler-factor neutrality; PROVED).** Let `Phi` be a finite product of factors `(1 - c p^-s)^(+-1)` with `|c| = 1`, and `(1 - c^(1-s))^(+-1)` with `c > 1`. Then `Phi` has zeros and poles only on `Re s = 0` and `Re s = 1`. So `Phi zeta` and `zeta` have the same zeros in `0 < Re s < 1`: any statement whose only arithmetic input is finitely many primes, such as Collatz's 2 and 3, is RH-neutral.
* **TRUNC.** Truncated Dirichlet series, including the partial sums of `zeta` itself, do not have their zeros on the line. Partial sums have zeros in a strip `alpha < sigma < 1.73` with `alpha > -X` (Gonek–Ledoan, Thm 1), and zeros with `Re s > 1` for every `N` except 1–18, 20, 21, 28 (Platt–Trudgian 2016). See also T3.
* **SHUF.** Square-root cancellation (`H = 1/2`) is produced by every weakly dependent sequence. See T1.
* **GENERIC.** Zeros on a prescribed line are non-generic unless a symmetry makes a normalized function real there, as the functional equation does for Hardy's `Z`. So T2's negative result already disposes of zero alignment for the data series.
* **DRIFT and SHEET.** The foundry's barriers, applied through P11.

**3.2 Where each proposed bridge fails.**

| # | bridge | type | fails |
|---|---|---|---|
| B1 | trunk = Iwasawa coordinate and pole factor of `zeta_3` (§2) | REAL identity | FIN (it removes the Euler factor at 3 and applies the `c = 4` regularizer); `zeta_3` is zero-free; DRIFT and SHEET (P11) |
| B2 | trunk = the `eta` / Genocchi factor | REAL identity | FIN (the zeros lie on `Re s = 1`) |
| B3 | E-game hostile point `1/2` = critical line `1/2` | NUMEROLOGY | type mismatch; the trunk chart sends `1/2` to `i*` and `s = 1/2` to `-1`; `M` does not preserve `{1, 1/2}` |
| B4 | `s <-> 1-s` becomes `M` on the trunk chart | REAL identity, empty | `zeta_3` has no zeros for the symmetry to constrain; the E-game is not `M`-symmetric |
| B5 | `D_T(s) = sum_(i>=1) T_i^-s` | object | zero-free for `Re s > 0.4594` (PROVED); poles on every line `Re s = -j`; no functional equation (T2) |
| B6 | Dirichlet series of the E-SCC data (`d`, `L1`, `H`) and of Collatz stopping times | object | GENERIC together with T2 (no functional equation); TRUNC (T3); SHUF (the exponents are halved dimensions, T1) |
| B7 | Collatz dilations `x -> 2^k x/3` compared with Nyman–Beurling dilations | ANALOGY | FIN. The two share only the Mellin-diagonal part, dilation by `2^k 3^-j`, which involves only the primes 2 and 3; the Collatz translation `-1/3` is not a Mellin multiplier. A Nyman–Beurling statement restricted to `<2,3>` would imply RH, since its closure lies inside Báez-Duarte's, so it is at least as hard as RH |
| B8 | square-root cancellation (Littlewood: RH iff `M(x) = O(x^(1/2+eps))`, UNVERIFIED primary; the `pi(x)` form is stated in BDR) | ANALOGY | SHUF; and the E-SCC data do not even show it (`H` about `0.25`) |
| B9 | the dynamical zeta of `Psi` (§4.3) | REAL object | it is not `zeta(s)` (FIN); no Ihara-type alignment was found |
| B10 | "F2[x] Collatz is provable and function-field RH is proved" | ANALOGY | the mechanisms differ (no carries, versus Frobenius and the index theorem); the two share only the absence of an archimedean place |

## 4. Numerical tests of the candidate bridges (exact scope)

**4.1 The trunk series (D1).**
* Put `sum_(i>=2) T_i^(-sigma_0) = 1`. This gives `sigma_0 = 0.459365562137` for all `i >= 1`, and `0.353857453736` for the legal exits (`3` not dividing `i`); the tail bound is `< 4*10^-83`.
* For `Re s > sigma_0`, `|D_T(s)| >= 1 - sum_(i>=2) T_i^(-Re s) > 0`. The margin at `Re s = 1/2` is `0.1178`.
* **Continuation.** `D_T(s) = 3^s sum_(j>=0) (s)_j/j! (4^(s+j) - 1)^(-1)` agrees with the direct sum to `10^-13`. The poles are at `s = -j + i pi k/log 2`. The residue at `i pi/log 2` is `3^s/log 4`, checked.
* **Zeros** (fp). In the box `[-3.9, 1.2] x [0.15, 30]`, 30 zeros were located (grid cells of winding `+1`), each polished to `|D_T| < 10^-20` at 30 digits, with real parts in `[-3.79, 0.4485]`. The 24 cells of winding `-1` are exactly the predicted poles. No zero comes within 0.05 of `Re s = 1/2`, as the proof requires.

**4.2 Data (C0) and the tests T1–T3.** The data are:
* the canonical-escape statistics `d(m)` (`Psi` steps until the value falls below `m`), `L1(m)` (expensive `1/2`-transfers) and `H(m)` (steps taken from precision `k >= 3`), for all `m` prime to 3;
* the Collatz shortcut stopping time `sigma(m)` and total stopping time `tau(m)`.

Summary facts:
* All 666,666 units `m <= 10^6` descend.
* The mean of `d` is `1.0714`. On `m = 14 mod 27` it is `2.2846`, against the endgame lane's `2.2848` to `10^11`; this is an independent code path.
* The maximum of `d` is 9 at `m = 221576` (below `10^6`) and 11 at `m = 8751065` (below `10^7`).
* The mean of `sigma` is `3.4849`, with maximum 176 at `m = 626331`. The maximum of `tau` is 329 at `m = 837799`.
* Sanity checks: `M(10^6) = 212` and `L(10^6) = -530`.

**T1: fluctuation exponents**, `m <= 10^7`. Two estimators were used. The first is the path exponent `alpha`, from `max|S|` and from `rms S`, on `[10^4, 10^7]`; this is the literal abscissa of convergence. The second is the Hurst exponent `H` from variances of block sums, with block lengths `2^4..2^18`. The null distribution comes from 20 shuffles.

| sequence | `alpha` (max, rms) | `H` | `H` on `2^4..2^11` / `2^11..2^18` | shuffle null: `alpha`(max), `H` |
|---|---|---|---|---|
| `Psi` steps `d` | 0.256, 0.290 | **0.255** | 0.208 / 0.312 | 0.467 ± 0.064, **0.496 ± 0.006** |
| `L1` links | 0.261, 0.310 | 0.219 | 0.222 / 0.251 | 0.478 ± 0.106, 0.498 ± 0.008 |
| hostile landings `H` | 0.250, 0.279 | 0.253 | 0.311 / 0.253 | 0.443 ± 0.093, 0.498 ± 0.009 |
| Collatz `sigma` | 0.602, 0.640 | **0.446** | 0.400 / 0.500 | 0.467 ± 0.123, **0.498 ± 0.008** |
| Möbius `mu` | 0.477, 0.491 | 0.496 | | |
| Liouville `lambda` | 0.489, 0.501 | 0.496 | | |
| random signs | 0.565, 0.601 | 0.488 | | |
| `1_(m=1 mod 27) - 1/27` | 0.000, 0.000 | -0.012 | | |
| `v_3(m) - 1/2` | 0.082, 0.080 | 0.030 | | |

Reading the table:
* **The estimators.** Single-path `alpha` estimates are too noisy to use: the null standard deviation is 0.06–0.17. `H` is sharp, with standard deviation `<= 0.009`.
* **E-SCC.** The E-SCC statistics give `H` about `0.22–0.26` overall (`0.21–0.31` by scale), 27 to 40 standard deviations below the shuffle null.
* **Collatz.** Stopping times give `H = 0.446`, 6.5 standard deviations below the null. The value rises with scale, from `0.400` to `0.500`.
* **Only Möbius, Liouville, random signs and shuffles give 1/2.**
* **Explanation (candidate C1, §6).** Model a statistic determined by digits as `sum_k 1_(A_k)(m mod p^k)`, with the exceptional residue classes placed at random. That model gives `Var(block of length L) ~ L^delta`, hence `H = delta/2`, where `delta` is the local growth rate of the exceptional counts.
  * **Collatz.** The exact counts `N_k` of classes mod `2^k` with no multiplicative descent reproduce the choice ladder's `93,222` (`k = 22`) and `1,037,374` (`k = 26`). Their local rates predict `H = 0.387` (low scales) and `0.429` (high scales), with limit `h(log_3 2)/2 = 0.47498`. The measured values are `0.400` and `0.500`: in agreement at low scales, and at high scales within the noise of 38–2441 blocks (which also cannot exclude `1/2`).
  * **`Psi`.** The endgame lane's alive counts 24, 274, 3050 and 31438 at `n = 7, 11, 15, 19` give local rates `0.531–0.554`. These predict `H` about `0.27`, with limit `<= 0.6415/2 = 0.32075`; the measured values are `0.208` and `0.312`.

**T2: functional-equation detector.** The test asks whether `Lambda(s) = eps Lambda(1-s)` holds for some gamma factor:
* degree 1: `(q/pi)^((s+kappa)/2) Gamma((s+kappa)/2)`, with `q <= 300` and `kappa` in `{0, 1}`;
* degree 2: `(sqrt(q)/2pi)^s Gamma(s + (w-1)/2)`, with `q <= 150` and `w <= 12`.

Method:
* This is tested through the equivalent theta relations `theta(1/t) = eps t^(kappa+1/2) theta(t)`, respectively `g(1/y) = eps y^w g(y)`, fitted by least squares on 81 points of `[1/12, 12]`. Polar terms are allowed: `{1, t^(1/2), t, t^(3/2)}`, respectively `{1, y^w}`.
* The residual is normalized by the part of `theta(1/t)` orthogonal to the polar terms. Fits where that part vanishes are "degenerate" and skipped. This matters: for `a_n = 1` and large `q` the theta function is pure Poisson mass. A first version without this normalization accepted `zeta` at `q = 136`.

Results (best residual `(q, kappa or w, eps)`):

| sequence | degree 1 | degree 2 |
|---|---|---|
| control `zeta` | 6.9e-14 (1, 0, +1.0000) | – |
| control `zeta` with `a_2` perturbed by `10^-6` | 1.5e-6 (1, 0, +1.0000) | – |
| control `L(chi_-4)` | 1.2e-15 (4, 1, +1) | – |
| control `L(chi_5)` | 8.4e-15 (5, 0, +1) | – |
| control random signs | 7.2e-2 | – |
| control `Delta` | – | 4.0e-11 (1, 12, +1.0000) |
| control `Delta` with `tau(2)` perturbed by `10^-6` | – | 3.8e-4 |
| control `E_4` | – | 5.4e-12 (1, 4, +1) |
| trunk indicator | 1.3e-1 | 9.1e-3 (89, 12, +0.11) |
| `Psi` steps `d` | 1.8e-1 | 2.9e-2 |
| `L1` | 8.2e-2 | 5.9e-1 |
| `H` | 4.1e-2 | 7.4e-1 |
| Collatz `sigma` | 2.8e-1 | 4.6e-1 |
| Collatz `tau` | 1.9e-1 | 4.9e-1 |

Every genuine L-function is recovered with the correct conductor and sign `eps = +1`. No data sequence comes within six orders of magnitude of the real ones: the best, the trunk indicator in degree 2 at `9.1*10^-3`, does not even have `eps` near `+-1`.

**T3: zeros of `sum_(n<=1000) a_n n^-s` in `[-1.5, 2] x [2, 40]`** (fp; cell windings, then Newton polishing).

| sequence | zeros | range of `Re s` | fraction with `abs(Re s - 1/2) < 0.05` | median `abs(Re s - 1/2)` | distances to the 6 zeros of `zeta` with `t < 40` |
|---|---|---|---|---|---|
| partial sums of `zeta` | 42 | [0.30, 0.93] | 0.24 | 0.125 | 0.24–0.34 |
| Möbius | 39 | [0.08, 0.54] | 0.51 | 0.042 | 0.81–0.88 |
| random signs | 42 | [-0.88, 1.27] | 0.14 | 0.237 | 0.08–0.71 |
| trunk (5 terms) | 35 | [-0.29, 0.42] | 0 | 0.549 | 0.47–0.68 |
| `Psi` steps `d` | 37 | [0.21, 0.99] | 0.16 | 0.144 | 0.04–0.35 |
| `d`, 10 mean-matched shuffles | | | 0.17 ± 0.06 | 0.154 ± 0.021 | means 0.16–0.35; minima 0.03–0.27 |
| Collatz `sigma` | 35 | [-0.36, 1.52] | 0.03 | 0.414 | 0.35–0.88 |
| `sigma`, 10 shuffles | | | 0.09 ± 0.05 | 0.27 ± 0.10 | means 0.37–0.48 |

* A sequence with mean `c` contains `c zeta_N(s)`, multiplied by `(1-3^-s)` when it is supported on units, so its zeros are perturbed partial-sum zeros.
* The close approaches of the `d`-polynomial to two zeros of `zeta` (0.04 and 0.08) are reproduced by shuffles, whose minima are 0.03–0.04. Nothing here is Collatz-specific.
* The Möbius polynomial's zeros cluster near `Re s = 1/2`, the line of square-root cancellation, while staying about `0.8` away from the zeros of `zeta`. That is the SHUF "1/2", not RH.

**4.3 The E-game's own zeta (D5).** `Psi` has full branches `(h,k)`, each mapping its shell bijectively onto `Z_3^x` while expanding by `3^k`, with price `rho_h(k)`, and `rho_(1/2)(k) = 2 rho_1(k)`. Its price-weighted Ruelle (Artin–Mazur) zeta function is therefore

    zeta_Psi(s; theta) = 1/(1 - W(3^-s; theta)),    W(x; theta) = (1 + 2^theta) sum_(k>=1) rho_1(k)^theta x^k.

* **Where the hostile point `1/2` enters.** It enters as the factor `(1 + 2^theta)`.
* **The coefficients.** For `k >= 3`, `rho_1(k) = 2^(-{k log2 3})`, so the coefficients form a Hecke–Mahler-type series in the loop clock `floor(k log2 3)`.
* **Periodic points.** They are rational E-cycles. The two `k = 1` branches fix `-1/2` and `-1`.
* **Consistency (fp).** `min_theta s(theta) = 0.64150` at `theta* = 2.6735`, which is the endgame lane's `dim_H Bad_Psi <= 0.64150`. The Lundberg exponent is `theta_L = 8.6433`, against the lane's `8.6434`.
* **Resonances** (poles with `|x| < 0.8`, stable under truncation at `K = 150` and `K = 300`, polished with `K = 4000`):
  * at `theta*`: the leading pole `s = 0.64150` and one pair at `Re s = 0.2355`, where an Ihara-type line would be at `s_0/2 = 0.32075`;
  * at `theta_L`: `Re s = 1, 0.9065, 0.7508, 0.5276, 0.4239`, where the Ihara-type line would be at `0.5`;
  * at `theta = 1`: only the leading pole.

So the E-game's zeta satisfies no Ihara-type Riemann hypothesis at these weights.

## 5. Verdict

**REAL.**
* **The trunk as a 3-adic object.**
  * `T` is the 3-adic interpolation of the trunk and a bijective isometry of `Z_3`.
  * Its fixed points are exactly `{0, 1, -1/2}`. `T(n/2) = ((-2)^n-1)/3` places the minus trunk at the half-integers.
  * In the logarithmic charts the backward E-game is "translate, then apply the expanding re-chart `Theta`"; the kappa formula is the isometry.
* **The trunk as Iwasawa's coordinate.**
  * The trunk is Iwasawa's coordinate for `u = 4 = 1+3`, and `3T(s)` is exactly the pole factor of `zeta_3`, with `3T zeta_3` a unit power series in `3T`. Consequently `zeta_3` has no zeros.
  * `s -> 1-s` acts as the Iwasawa involution `M`.
  * The same coordinate exists for every non-Wieferich `qx+1` map.
* **The trunk archimedeanly.** It is the Euler-type factor at 2 in `eta` and in the Genocchi numbers, with zeros on `Re s = 1`.
* **The E-game's own zeta.** It is a real object in which the price of the hostile point `1/2` appears. Its leading pole reproduces the dimension bound 0.6415.

**ANALOGY.** "1/2" means four different things here:
* **RH's `1/2`:** the axis of `xi(s) = xi(1-s)`, equivalently the square-root error term in `pi(x)` and `M(x)`.
* **The CLT `1/2`:** the square-root cancellation of weakly dependent sequences, which Möbius and shuffles show and E-SCC data do not.
* **The E-game's `1/2`:** the 3-adic unit with `2x = 1`. It is hostile because escape prices exceed 1, an archimedean property of the costs.
* **The 3-adic `s = 1/2`:** the point `X = -3` of the weight disc. The trunk chart sends it to `-1`, because `log(-2)/log 4 = 1/2` is the exponent relating the two generators `-2` and `4` of `1+3Z_3`.

The first two are genuinely related (Littlewood). The third is related to none of them, and the fourth links only to `-1`. Two further analogies connect no mechanisms:
* the Nyman–Beurling dilations against the Collatz dilations (only the primes 2 and 3 are shared, and the Collatz translation is not Mellin-diagonal);
* "function-field Collatz and function-field RH are both provable".

**NUMEROLOGY.**
* `T(1/2) = -1` "is a hostile point": true, but it is a hostile point of other games (the forward game and the minus sheet), not of the Q2 game whose hostile point is `1/2`. It comes from `sqrt 4 = -2` in `Z_3`.
* `i* = 1/2 mod 3`: forced by the isometry and `1/2 = -1 mod 3`.
* The trunk's complex zeros "lie on a line": every single exponential does that.
* A Collatz-statistics exponent "near 1/2" is `h(log_3 2)/2 = 0.475` in the limit, and `0.40–0.50` measured.
* The Dirichlet-polynomial zeros "near" zeros of `zeta`: reproduced by shuffles.

**One research question that would make the analogy precise (question R).** Consider the canonical escape's dynamical zeta `zeta_Psi(s; theta) = 1/(1 - W(3^-s; theta))`, in which the hostile point `1/2` contributes the factor `1 + 2^theta` and the coefficients run along the loop clock `floor(k log2 3)`. Locate its poles (the resonances) in the strip `s(theta) - 1 < Re s < s(theta)`.
* Is there a spectral gap that is uniform in `theta`?
* Is the gap governed by the `1/2`-branch?
* What natural boundary does `W` have? For irrational slopes, Hecke–Mahler-type series have the unit circle as natural boundary; this is UNVERIFIED background.

This is the one well-posed setting in which "zeros on a line" and the hostile point `1/2` live in the same function. It concerns the E-game's zeta, not Riemann's. By standard renewal theory the resonances govern the finite-level fluctuations of the `Psi`-alive counts, and hence the approach of the counts `24, 274, 3050, 31438, ...` to the dimension 0.6415. The numerics of §4.3 already exclude the Ihara-type "RH" at two natural weights. The open part is the gap.

## 6. Hypothesis candidates (not filed as HYP files)

* **C1 (a dimension law for fluctuation exponents).** Suppose a statistic `a(m)` of a descent process is determined by the `p`-adic digits it consumes, and its exceptional set has dimension `delta` (normalized to base `p`). Then:
  * the summatory fluctuation `sum_(m<=x) (a(m) - mean)` grows like `x^(delta/2 + o(1))`, so the fluctuation Dirichlet series has abscissa of convergence `delta/2`;
  * the variance of block sums scales like `L^(delta+o(1))`.

  For Collatz stopping times this predicts `h(log_3 2)/2 = 0.47498`. For `Psi` descent it predicts `dim Bad_Psi / 2 <= 0.32075`.
  * **Support:** T1, where both the scale dependence and the exact `N_k` rates match at low scales; and the random-placement heuristic.
  * **Against, or open:** at high scales `sigma` gives `H = 0.500 ± (noise)`, not separated from `1/2`. The path exponents are too noisy to test the abscissa directly.
  * **Test:** block statistics of `sigma` to `2^34` in C.
* **C2 (resonances of `Psi`).** Question R above, in the form of a conjecture: every non-leading pole of `zeta_Psi(.; theta)` with `|3^-s| < 1` satisfies `Re s <= s(theta) - eta` for some `eta > 0` uniform in `theta` in compact subsets of `(0, infinity)`. Numerics: at `theta*` the gap is `0.64150 - 0.23547 = 0.406`, and at `theta_L` it is `1 - 0.90647 = 0.094`.

## 7. Reproduction

```bash
python3 04-computation/experiments/collatz_procgen_20260923_trunk_rh_padic.py      # P1-P12, about 7 s, < 50 MB
python3 04-computation/experiments/collatz_procgen_20260923_trunk_rh_dirichlet.py  # D0-D5, about 3 min, < 250 MB
python3 04-computation/experiments/collatz_procgen_20260923_trunk_rh_data.py       # C0, T1-T3, about 2.5 min, peak RSS 0.75 GB
```

The `.out` concatenates the three outputs, with fixed seeds. Run the three scripts one at a time.
* **The Euler–Maclaurin Hurwitz evaluator** is validated against mpmath to `3.6*10^-13` (D0).
* **The 3-adic computations** are exact modulo `3^56` or better.
* **The vectorized data code** is checked against the loop code on `m <= 10^6`.

Downloaded literature stays in `scratch/procgen_trunkrh/lit/`, which is git-ignored by the local `.gitignore`. No bot check or paywall was bypassed: the AMS pages and the Notices article sit behind one and were not read.

## 8. Sources

**Read.** Status: [P] means the stated sections of the primary source were read; [A] means the abstract only, via the arXiv API or Crossref metadata.
* [P] J. Rodrigues Jacinto, C. Williams, *An introduction to p-adic L-functions*, Ess. Number Th. 4 (2025) 101–216, arXiv:2309.15692. Sections read: Thm 2.13, 4.1, 5.1, 5.17, 7.1; Props 4.4–4.11; Lemma 3.36; Remarks 3.47, 5.21; Def 5.16; §11.1–11.3.
* [P] E. Bombieri, *Problems of the Millennium: the Riemann Hypothesis*, Clay Mathematics Institute official problem description, pp. 4–10. Covers the function-field case (Hasse, Weil, Deligne), Weil's explicit formula and positivity, the index theorem, Montgomery/Odlyzko/Rudnick–Sarnak/Katz–Sarnak, the result that more than 40% of zeros lie on the line (Selberg, Levinson, Conrey), and density theorems.
* [P] F. Broucke, G. Debruyne, Sz. Révész, *Some examples of well-behaved Beurling number systems*, Trans. AMS 378 (2025) 477–501, arXiv:2309.01567. Introduction pp. 1–3, including the statement of Diamond–Montgomery–Vorhauer's optimality result.
* [P] J. C. Lagarias, *Li coefficients for automorphic L-functions*, Ann. Inst. Fourier 57 (2007) 1689–1740, arXiv:math/0404394. Abstract and §1: Li's criterion, and positivity for symmetric multisets.
* [P] S. M. Gonek, A. H. Ledoan, *Zeros of partial sums of the Riemann zeta-function*, arXiv:0807.0019. Abstract, Thms 1–2.
* [A] B. Rodgers, T. Tao, *The de Bruijn–Newman constant is non-negative*, arXiv:1801.05914 (Forum Math. Pi 2020).
* [A] D. H. J. Polymath, arXiv:1904.12438 (`Lambda <= 0.22`).
* [A] L. Báez-Duarte, *A strengthening of the Nyman–Beurling criterion for the Riemann Hypothesis*, arXiv:math/0202141.
* [A] D. J. Platt, T. S. Trudgian, *The Riemann hypothesis is true up to 3·10^12*, arXiv:2004.09765.
* [A] D. J. Platt, T. S. Trudgian, *Zeroes of partial sums of the zeta-function*, LMS J. Comput. Math. 19 (2016), arXiv:1507.01340.
* [A] Y. Lamzouri, *Zeros of the Epstein zeta function to the right of the critical line*, Math. Proc. Camb. Phil. Soc. 171 (2021), arXiv:1907.06387.
* [A] I. S. Rezvyakova, arXiv:2411.18492 (a positive proportion of Epstein zeros on the line).
* [A] E. P. Balanzario, J. Sánchez-Ortiz, *Zeros of the Davenport–Heilbronn counterexample*, Math. Comp. 76 (2007) 2045–2049, doi:10.1090/S0025-5718-07-01999-0.
* [A] J. T. Sheats, arXiv:math/9801158; J. Kramer-Miller, J. Upton, arXiv:2312.01264 (RH for Goss zeta functions).
* [A] L. Zhao, arXiv:2201.08870 (context for the branch typing of MISTAKE-497).

**UNVERIFIED (not read).**
* H. G. Diamond, H. L. Montgomery, U. Vorhauer, *Beurling primes with large oscillation*, Math. Ann. 334 (2006): known only through BDR.
* H. Davenport, H. Heilbronn, *On the zeros of certain Dirichlet series* I and II, J. London Math. Soc. 11 (1936): titles only.
* A. A. Karatsuba (1991), on Davenport–Heilbronn zeros on the critical line: title only.
* J. B. Conrey, *The Riemann Hypothesis*, Notices AMS 2003: behind a bot check.
* L. C. Washington, *Introduction to Cyclotomic Fields*: the moment formula for Mazur's `E_(1,c)` used by the script. It is not needed: the proof uses RJW's `mu_4`, and P7 checks numerically that the two measures have the same moments on `Z_3^x`.
* The Riesz criterion, Littlewood's Möbius criterion (the primary sources), and the Hecke–Mahler natural boundary.
* Strassmann's theorem and von Staudt–Clausen are standard; I did not read a specific textbook page.

## 9. HYP-9133 is PROVED (orchestrator, 2026-09-23)

**Theorem (resonance gap for the canonical escape).** Let
`zeta_Psi(s; theta) = 1/(1 - W(3^-s; theta))` with
`W(x; theta) = (1 + 2^theta) sum_(k>=1) rho_1(k)^theta x^k`, where `rho_1(k) = 2^(K*(k-1))/3^k`.
For every `theta > 0` the leading pole `s(theta)` (`x_0 = 3^(-s(theta))`) is simple, and it is the only pole on
`|x| <= x_0`. Every other pole with `|3^-s| < 1` satisfies `Re s <= s(theta) - eta`. The gap `eta > 0` is uniform for
`theta` in compact subsets of `(0, infinity)`. Since every digit count `k` is an integer, `zeta_Psi` is
`2 pi i/log 3`-periodic in `s`. "Other" poles are therefore meant modulo this period, i.e. in the variable
`x = 3^-s`. The copies `s(theta) + 2 pi i m/log 3` are the same pole `x_0`.

*Proof.*
1. **Coefficients.** `c_k = (1+2^theta) rho_1(k)^theta` satisfies `(1+2^theta) 3^(-theta) <= c_k < 1 + 2^theta`,
   because `1/3 <= rho_1(k) < 1` (the loops-lane values: `rho_1(1) = 1/3`, `rho_1(2) = 4/9`, and
   `rho_1(k) = c(k-1)/3` in `(1/2, 1)` for `k >= 3`). So `W` is analytic on `|x| < 1`, with radius exactly 1,
   and `W(r) -> infinity` as `r -> 1^-`.
2. **The leading zero.** `W` is strictly increasing on `[0,1)` from `0`. So `1 - W` has a unique zero
   `x_0 in (0,1)`, and it is simple, since `W'(x_0) > 0`.
3. **Pringsheim with aperiodicity.** If `|x| <= x_0`, then `|W(x)| <= W(|x|) <= 1`. Equality throughout
   forces `|x| = x_0` and `x^k > 0` for every `k` with `c_k > 0`. Since `c_1 > 0`, this forces `x = x_0`.
   So `x_0` is the only zero of `1 - W` in the closed disc `|x| <= x_0`.
4. **The gap for fixed theta.** Zeros of `1 - W` are isolated in `|x| < 1`. So for any `r_1 in (x_0, 1)` only
   finitely many lie in `|x| <= r_1`, and all of them other than `x_0` have modulus `> x_0`. Take
   `eta = log_3(min(r_1, min |x_j|)/x_0) > 0`.
5. **Uniformity.** Suppose `theta_n` lie in a compact set `K`, `x_n != x_0(theta_n)` are zeros, and
   `|x_n|/x_0(theta_n) -> 1`. Pass to a subsequence with `theta_n -> theta_inf` and
   `x_n -> x_inf`, `|x_inf| = x_0(theta_inf)`. Here `x_0` is continuous, and `W(.; theta_n) -> W(.; theta_inf)`
   locally uniformly on `|x| < 1`, because the coefficients are continuous and uniformly bounded on `K`.
   By Hurwitz, `x_inf` is a zero, so by step 3 `x_inf = x_0(theta_inf)`. That zero is simple, so near it
   `1 - W(.; theta_n)` has exactly one zero for large `n`. That zero is `x_0(theta_n)`, which contradicts
   `x_n != x_0(theta_n)`. ∎

**Numerics** (`04-computation/experiments/collatz_procgen_20260923_hyp9133_gap.py`). The roots of the
degree-300 and degree-600 truncations agree to 6 digits inside `|x| <= 0.95`.

| theta | s(theta) | next resonance, Re s | gap |
|---|---|---|---|
| 1 | 0.74464 | 0.06148 | 0.683 |
| 2 | 0.65410 | 0.13911 | 0.515 |
| 4 | 0.67348 | 0.40793 | 0.266 |
| 8.6434 (theta_L) | 1.00000 | 0.90647 | 0.094 |

The last row matches §4.3 above. At `theta = 0.5` no other zero lies in `|x| < 0.95`.

**What this does and does not say.** The Ihara-type "Riemann hypothesis" (all non-leading resonances on one line)
is still REFUTED numerically (§4.3). The gap statement holds for the elementary reason above: positive
coefficients with `c_1 > 0`. So the one precise form of the "hostile 1/2 versus critical line" analogy is a
spectral gap, and it holds. It is a renewal-theoretic fact, not an RH-type fact.
