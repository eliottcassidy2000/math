---
id: THM-4555
title: "Uniform switches of the Collatz rewrite compiler are collisions at -1: if two reduced words satisfy f_u(-1) = f_u'(-1) (f_c(x) = (3x+1)/2^c) and |u'| = |u| + D, every odd n whose word starts 1^r u (r >= D) merges at equal length with m = (n+1)/2^D - 1 (delete D trailing ones); the reset switch is the root collision (a) ~ (2, a-2); for reset-2 sources the root collision points upward, so their downward equal-length uniform switches come only from sporadic collisions (377 of 2500 reset-2 sources below 2*10^4)"
status: "PROVED (i)-(v); FINITE-EXACT (vi); INDEPENDENTLY AUDITED (2026-10-06; corrections in MISTAKE-570; record in section 9 of the results note)"
session: mac-mini-2026-10-06-mod1819
source: 05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md (section 6b)
scripts:
  - 04-computation/experiments/mod1819_20261006_trailing_ones_switch.py (+ .out, ALL CHECKS PASSED)
related:
  - 05-knowledge/results/checked_switch_phase19_20261004.md (the compiler, the reset switch (3), the debt states (4)-(5), the 239 residual seeds)
  - 01-canon/theorems/THM-4554-backward-sieve-keeps-a-positive-proportion-moran-duality.md (the 3-adic side; Mersenne sampling)
---

# THM-4555 — uniform switches are collisions at −1

**Setting.**
* `U(x) = oddpart(3x+1)` on odd `x`. A word `w = (a_1, ..., a_j)` lists exponents, with total `A_w = Σ a_i`.
* Its sources are the odd `n` with `U^j(n) = z` along `w`, and `n = I_w(z) = (2^(A_w) z − B_w)/3^j` (`checked_switch_phase19`, (1)).
* On `Q` put `f_c(x) = (3x+1)/2^c`, and let `f_u` be the composition along `u`. A word is **reduced** if its first letter is `>= 2`.
* Every odd `n > 1` is `n = 2^(r+1) t − 1` with `t` odd. Its word is `1^r u` with `u` reduced, where `r` counts the leading exponents 1.
* Two reduced words **collide** if `f_u(−1) = f_u'(−1)`.
* A switch is **uniform** if it pairs the family `1^r u` with a family `1^(r+s) u'` (`u`, `u'` reduced, `s` fixed) for infinitely many `r`.

**(i) Endpoint progression (PROVED).** The endpoints of the word `1^r u` are exactly the positive odd `z` with `z ≡ f_u(−1) (mod 3^(r+|u|))`. For each such `z`, `I_w(z)` is a positive odd integer whose word starts with exactly `1^r u`.

**(ii) Denominators (PROVED).** For reduced `u`, `f_u(−1) = N/2^(A_u − 1)` with `N` odd. Hence colliding reduced words have equal totals.

**(iii) Uniform families are collisions (PROVED).** Fix `u`, `u'` and a shift `s`. The families `1^r u` and `1^(r+s) u'` share endpoints for infinitely many `r` iff `u` and `u'` collide.

**(iv) The switch (PROVED).** Suppose `u` and `u'` collide and `|u'| = |u| + D` with `D >= 1`. Then every odd `n` whose word starts with `1^r u`, `r >= D`, satisfies

    U^j(n) = U^j(m),     j = r + |u|,     m = (n + 1)/2^D − 1 < n,

and the word of `m` starts with `1^(r−D) u'`. In binary, `m` is `n` with `D` trailing ones deleted.

* The **reset switch** of `checked_switch_phase19` (its (3), `m = (n−1)/2`) is the root collision `(a) ~ (2, a−2)`, `a >= 3`, with `D = 1`.
* Within any family the multiples of 3 form one class of the endpoint progression mod `3^(j+1)`, and (iv) applies to them like any other source.
* **No uniform switch is special to multiples of 3.**
  * By (iii), a uniform switch firing on multiples of 3 for infinitely many `r` is a collision.
  * By (v), a longer partner fires only on `n ≡ −1 (mod 3)`.
  * A partner of length `k <= j` gives `m + 1 = 2^s 3^(j−k)(n+1)` on the whole family.
  * This answers D-α negatively in the uniform sense: multiples of 3 receive exactly the uniform switches of their word family (for example 10 of the 153 seeds).

**(v) Longer partners and reset 2 (PROVED).**
* If the partner word `1^(r+s) u'` has length `k > j` (i.e. `s > −D`), the join fires only on `n ≡ −1 (mod 3^(k−j))`. Such `n` already have the smaller predecessor `(2n−1)/3`.
* For a reset-2 source (word `1^r (2, c, ...)`) the root collision `(2, c, ·) ~ (c+2, ·)` has `D = −1`: it points upward (it is the reset switch read backwards). With the partner's run lengthened by 3 (partner `1^(r+3)(c+2, ·)`, `m = 8(n+1)/9 − 1`) it fires downward exactly on `n ≡ 8 (mod 9)`.

**(vi) Census (FINITE-EXACT; merges are first equal-length meetings on orbits stopped at 1).**
* **Collisions.** Among reduced words of length `<= 6` with letters `<= 8`, 28014 values are hit more than once. 774 of them are not explained by the root collision alone (their words fall into at least two classes of `(2, c, x) ~ (c+2, x)`; extensions `ux ~ u'x` count separately).
  * The smallest sporadic one (by total and by `|N|`, over all reduced words) is `(8, c) ~ (4, 1, 1, c+2)`, an identity for every `c >= 1` (PROVED: both sides equal `125/2^(7+c)`). Both give `125/2^(7+c)`, because `125 = 2^7 − 3 = 2^5 + 3·2^4 + 9·2^3 − 27`.
  * Composed with the root collision, it gives `(2, 6, c) ~ (4, 1, 1, c+2)`. So reset-2 sources with word `1^r (2, 6, c)` switch to `(n−1)/2` at depth `r + 3`.
* **Reset-2 sources `n < 2·10^4`.** 377 of 2500 (15.1%) have a uniform trailing-ones switch, with merge depth beyond the run between 3 and 40. With orbits stopped at 1, every equal-length merge of such an `n` with some `(n+1)/2^D − 1` is a collision; there are no non-collision merges. (Under `U(1) = 1` every partner also meets `n` at 1, e.g. `n = 7`, `m = 3` at step 5, which is not a collision; the collision counts 377, 16 and 37 are unchanged.)
* **The 239 residual seeds** of `checked_switch_phase19`: 16 have a uniform switch (10 of the 153 multiples of 3).
* **Odd Mersenne numbers** `2^a − 1` (all reset-2 sources, `t = 1`):
  * 37 of the 60 odd `a` in `[3, 121]` switch to some `2^(a−D) − 1`.
  * The least such `D` is always odd. So the partner is an even Mersenne number, which the reset switch takes one more step down.
  * All of these switches are collisions.

## Proofs

(i) `U_w(x) = (3^j x + C_w)/2^(A_w)` is affine with slope `3^j/2^(A_w)`, and as a map of `Q` it is `f_w`. Since `f_1(−1) = −1`, `f_(1^r u)(−1) = f_u(−1) =: Q`. For a source `n` of `w`,

    z − Q = U_w(n) − U_w(−1) = 3^j (n + 1)/2^(A_w),     so     n + 1 = 2^(A_w) (z − Q)/3^j.        (*)

Here `n + 1` is an integer and `Q` has a power-of-2 denominator, so `z ≡ Q (mod 3^j)`. Conversely, for positive odd `z ≡ Q`, `I_w(z)` is an integer. It is positive because each backward step `x -> (2^a x − 1)/3` keeps positive integers positive. And the compiler note's backward-integrality argument shows that its word is `w`. ∎

(ii) `f_(c_1)(−1) = −1/2^(c_1 − 1)` with `c_1 − 1 >= 1`. For `A >= 1` and `N` odd, `f_c(N/2^A) = (3N + 2^A)/2^(A+c)` with `3N + 2^A` odd. By induction the denominator exponent is `A_u − 1` and the numerator is odd. Equal values therefore have equal denominators. ∎

(iii) By (i), the two progressions are `z ≡ f_u(−1) (mod 3^(r+|u|))` and `z ≡ f_u'(−1) (mod 3^(r+s+|u'|))`. They intersect iff the two values agree modulo the smaller modulus. This happens for infinitely many `r` iff they are equal in `Z_3`, which for elements of `Z[1/2]` means equal in `Q`. By (i), every positive odd `z` in the intersection is an endpoint of both families. ∎

(iv)
1. Let `w = 1^r u` and `v = 1^(r−D) u'`. Both have length `j`, because `|u'| = |u| + D`. They have the same `Q` (collision). By (ii), `A_v = A_w − D`.
2. For a source `n` of `w`, `z = U^j(n)` satisfies `v`'s progression. So `m = I_v(z)` is an integer whose word is `v`.
3. By (*) for both words, `m + 1 = 2^(A_v)(z − Q)/3^j = (n + 1)/2^D`.
4. `n + 1 = 2^(r+1) t`, so `m + 1 = 2^(r+1−D) t >= 2`. Hence `m` is a positive odd integer, and `m < n`.
5. For the reset switch, `f_a(−1) = −2^(1−a) = f_(2,a−2)(−1)`. Its multiples of 3: `n ≡ 0 (mod 3)` iff `n + 1 ≡ 1`, iff `z ≡ Q + 3^j 2^(−A_w) (mod 3^(j+1))` by (*). ∎

(v)
* By (*) for the longer word `v`, `z ≡ Q (mod 3^k)`. So `3^(k−j)` divides `z − Q`, and hence divides `n + 1`.
* For the reset collision, `f_(2,c)(−1) = −2^(−1−c) = f_(c+2)(−1)`, where `(2, c)` is the longer word.
* Pairing `1^r (2, c)` with `1^(r+s) (c+2)` gives `m/n -> 3 (2/3)^s`, which is below 1 only for `s >= 3`. The partner length is then `k = j + s − 1`, so the join needs `n ≡ −1 (mod 3^(s−1))`. For `s = 3` this is `n ≡ 8 (mod 9)`. ∎

(vi) Script sections 2 and 5. The census of equal-length merges searches `D = 1..r` on orbits stopped at 1, depth `<= 300` (`<= 4000` for the Mersenne numbers). No orbit exceeds 102 (resp. 593) odd steps, so the caps never bind. ∎

## Remarks

* **What the reset-2 "debt" is.**
  * A reset-2 source's own root collision points upward. A uniform equal-length (trailing-ones) rewrite for it therefore needs a *sporadic* collision. That is a Diophantine coincidence `−3^(p−1) + Σ 3^(p−1−i) 2^(A_i) = −3^(p'−1) + Σ 3^(p'−1−i) 2^(A'_i)` with strictly increasing exponents `1 <= A_1 < …` below a common final exponent `E = A_u − 1`, other than the root identity `−1 = −3 + 2`. (The root collision itself gives only the longer-partner joins on `n ≡ 8 (mod 9)` of (v).)
  * In the census these give uniform switches to 15% of reset-2 sources. For the other 85%, no trailing-ones partner merges at equal length before the orbits reach 1 (all within 102 odd steps). Their rewrites must either depend on `r` or use other partners, which is the debt-state route of `checked_switch_phase19` (4)–(5).
* **2 <-> 3 mirror for Mersenne numbers.**
  * For `K >= 3`, up to total `K`, the reduced word of `2^a − 1` depends only on `a mod 2^(K−2)`. This is the 2-adic clock of 3, `ord_(2^K)(3) = 2^(K−2)`, read from `U^(a−1)(2^a − 1) = 2·3^(a−1) − 1`.
  * Backward-minimality of `2^a − 1` to depth `k` (odd `a`) instead depends only on `a mod 2·3^(k−1)` (THM-4554(v)), the 3-adic clock of 2. The forward switches and the backward sieve of the same family run on the two mirror clocks.
* **Prior art.** The root case is Theorem 2.1 of M. M. Ahmed, *The intricate labyrinth of Collatz sequences*, arXiv:1602.01617 (2016). That theorem covers the reset switch (3) and its upward reading for reset 2: `N` and `2N + 1` merge when `N`'s first reset is 2, and `4N + 3` merges with `2N + 1` otherwise. We know of no source for the general collision switch (iv).
* **Scope.** These are statements about rewrite rules. Each uniform switch transports a supplied certificate for the smaller `m`; nothing here proves a certificate exists. Collatz is OPEN.
