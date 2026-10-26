# Seven and twenty-one as the threefold iterates of the append maps 2x+1 and 4x+1: the reset switch is the sibling identity at the end of the run (Ahmed 2016), so trailing-ones switch lags pair up as {D, D+1} for every source and the least lag of every reset-2 source is odd; 2x+1 and 4x+1 generate F_21 = Aut(P_7) only at 7; finite-depth CRT independence of the Mersenne mirror clocks; collision statistics at −1; certified Mersenne switching to template total 27; the openai/math collection read for compression

2026-10-06, session opus-2026-10-06-S18 (worktree `codex/session-mersenne-f21-20261006`).

Owner's prompt (verbatim): "think of how ideas ananlgous to {7,21} relate deeply with the idea that every odd-exponent Mersenne switch has an odd shift distance D, and 37 of the 60 odd exponents in [3,121] produce collision. also explore deeply and intricately this repository of work and look to blend its ideas with ours in useful ways and extend concepts meaningfully toward proofs https://github.com/openai/math thinking in terms of information compression creatively especially https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026 explore around there very open mindenly, pursuing any possible connection between themes or similar looking equations, even if the topics seem like they could not be related". Three preprints from that repository were attached: *Weil classes and Hodge classes on abelian powers*, *Unrestricted pro-modularity at the prime two*, and *Weak mixing of triangular billiards with an irrational angle*.

**Concurrent work on the same prompt.** The mac-mini session `mac-mini-2026-10-06-mod1819` answered the same prompt in parallel and pushed during this session:

* [THM-4556](../../01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md): odd least shift on the Mersenne line, 2-adic periodicity, certified density 0.1199;
* [THM-4557](../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md): `Aut(T_k) = F_21`, known since Hanaki 2020 (MISTAKE-571);
* [HYP-9213](../hypotheses/HYP-9213-mersenne-collatz-trajectories-coalesce-plateaus-of-odd-step-time.md) and [HYP-9214](../hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md);
* the results note [seven_twentyone_mersenne_openai_math_20261006.md](seven_twentyone_mersenne_openai_math_20261006.md).

Their note reads the multiplication paper in depth: the negacyclic 3-adic clock `19·27 = 2^9 + 1`, Rader index maps, double-base representations, and Alon's orthogonal representation. It also cites this note's additions. This note cites theirs, reproduces the numbers it uses (script §A, §E, §H), and does not repeat their readings.

**Status.**

* **Restated, with credit.** Lemma 1, the sibling form of the reset switch, is Ahmed 2016, Theorem 2.1. Its general form is `U^r(2m+1) = 2^(e(m)) U^r(m) + 1`, and `checked_switch_phase19` (3)–(4) also contains it.
* **PROVED (elementary):**
  * Theorem 2, the parity law for trailing-ones lags. Part 1 is Ahmed's theorem iterated. Parts 2–3, the closure of the lag set and the odd least lag for every reset-2 source, are *new*. THM-4556 (iii) is the Mersenne case `t = 1`.
  * Proposition 3 (*new*): the two append maps generate `F_21 = Aut(P_7)`, and only at 7. That `F_21` is doubling plus translation is THM-4556's remark.
  * Proposition 5 (*new*, elementary): CRT independence of the forward and backward Mersenne clocks at every finite depth.
  * Proposition 6 (*new*): `μ_2(S) = 1` implies HYP-9213.
* **FINITE-EXACT:** every number so typed is printed by [the script](../../04-computation/experiments/mersenne_switch_parity_f21_compression_20261006.py) (output [`.out`](mersenne_switch_parity_f21_compression_20261006.out), ending `ALL CHECKS PASSED`) or by the [numba script](../../04-computation/experiments/mersenne_certified_density_numba_20261006.py) (output [`.out`](mersenne_certified_density_numba_20261006.out)).
* **NUMERICAL:** the fits.
* **HEURISTIC:** the coalescing-walk reading (HYP-9214's model).
* **ANALOGY / DICTIONARY:** all bridges to the openai/math preprints (section 7). These are unrefereed, model-written and partly unformalized; their claims are reported as claims.
* **AUDITED:** an independent adversarial audit (2026-10-06) recomputed everything with its own code. The mathematics holds; its corrections are applied below (section 10; MISTAKE-573).
* **Collatz: OPEN.**

## 0. Answers in brief

1. **Where the prompt's facts come from.**
   * "Every odd-exponent Mersenne switch has an odd shift distance `D`" and "37 of the 60 odd exponents in `[3,121]`" are THM-4555 (vi).
   * A *switch* there is an equal-time merge of `2^a − 1` with `2^(a−D) − 1`, i.e. with `2^a − 1` minus `D` trailing binary ones.
   * `{7, 21}` is the Paley tournament `P_7` and its automorphism group `F_21 = C_7 ⋊ C_3`: THM-4553's ladder, and the tower `T_k` of THM-4557.
2. **The mechanism behind "odd `D`" is one identity, and its maps are the ones whose threefold iterates give 7 and 21.**
   * Put `A_1(x) = 2x + 1` ("append a binary 1") and `R(x) = 4x + 1` ("append 01"). Then `A_1^3(0) = 7 = 111_2` and `R^3(0) = 21 = 10101_2`.
   * Collatz erases `R` in one step: `U(4x + 1) = U(x)`.
   * **Ahmed's theorem (Lemma 1).** If `n = 2m + 1` has run `r >= 1`, then `U^r(n) = 2^(e(m)) U^r(m) + 1`, where `e(m)` is `m`'s reset exponent.
     * When `e(m) = 2` this is `R(U^r(m))`, erased one step later. That is the reset switch `n ⇒ (n − 1)/2`, which happens iff `n`'s reset is at least 3.
     * When `e(m) >= 3` it creates the debt state of `checked_switch_phase19` (4).
   * **Theorem 2 (parity law).** Let `n = 2^(r+1) t − 1` with run `r >= 2`, and add the lag 0 when `n`'s reset is at least 3. Then the lags `D <= r − 1` at which `n` meets `(n + 1)/2^D − 1` at equal time come in pairs `{D, D+1}` with `D` "good". The only exception is an unpaired good top lag `r − 1`.
     * Good lags have the parity of `r + 1` if `t ≡ 1 (mod 4)`, and of `r` if `t ≡ 3`.
     * So the least lag is odd for every reset-2 source; it is 1 for every reset-`>= 3` source.
     * Mersenne numbers are the case `t = 1` (THM-4556 (iii)).
     * Checked on all 2500 sources below `2·10^4`, and in THM-4555's wider convention (all 377 reset-2 sources with a partner).
   * The parity law does not depend on 7. The link to `{7, 21}` is only that the same two maps appear: DICTIONARY.
3. **Why 7 and 21 in particular (Proposition 3; EXACT).**
   * Modulo 7, `A_1` and `R` generate exactly `F_21 = Aut(P_7)`.
   * Modulo a Mersenne prime `p = 2^k − 1` they generate `Z/p ⋊ <2>`. This is all of `Aut(P_p)` only when `<2> = QR_p`, i.e. `k = 3`, the equation of the Paley bridge note's Proposition 3.
   * `|F_21| = 21 = R^3(0)` amounts to `2^k + 1 = 3k`, i.e. `k ∈ {1, 3}`. At `k = 3` this is `2^3 + 1 = 3^2`, the Catalan–Mihăilescu solution that mac-mini's note finds behind 63.
4. **"37 of 60" is a compression statement.**
   * Mersenne switches are coincidences of the odd-step time `σ`.
     * For odd `a <= 2000`, every least-lag merge happens at a value `>= 911` (minimum at `a = 107`).
     * Even exponents merge with `a − 1` by the reset switch (at 5, 91, 205 for `a = 4, 6, 8`).
     * So no class is held together by simultaneous arrival at 1.
   * `σ(2^a − 1)` takes 68 values for `a = 2..2000` and 86 for `a = 2..4000` (HYP-9213: 104 for `a = 1..6000`).
   * Finite 2-adic templates certify 0.1199 of odd exponents at template total 20 (THM-4556 (v)). A parallel kernel reaches 0.1556 at total 27 (FINITE-EXACT). The template total at an actual least-lag merge has median 247.
   * All 932 least-lag switches of odd `a <= 2000` are collisions in THM-4555's sense (FINITE-EXACT). The first template, `a ≡ 95 (mod 128)`, is THM-4555's sporadic collision `(8,1) ~ (4,1,1,3)`, of value `125/256`.
5. **Mirror clocks (Proposition 5).**
   * The forward template switch at bounded total is periodic in `a` modulo `2^(K−2)` (THM-4556 (iv)).
   * Backward-minimality to depth `k` is periodic modulo `2·3^(k−1)` (THM-4554 (v)).
   * So the two events are independent among odd exponents at every finite depth.
   * The census is consistent: `P(both) = 0.367, 0.495, 0.550` against products `0.360, 0.487, 0.542`. This is an empirical observation, not a test of the finite-depth statement.
6. **Collisions at −1, counted (FINITE-EXACT; the "hash" reading is ANALOGY).**
   * The map `u -> f_u(−1)` is exactly 2-to-1 up to total 8 (the root collision).
   * Sporadic collisions take the image to `0.449·M` at total 18.
   * Their per-unit growth ratio is still falling: 2.36 to 2.14 between totals 13 and 18.
7. **openai/math (section 7).**
   * The multiplication paper's power saving is a rank inequality `s < Wm` for a finite exchange network built from three-element subsets of `[100]`, applied recursively.
   * Its ingredients have typed counterparts in our compiler. The Mersenne plateaus (HYP-9213, OPEN) are where our data show the same shape of saving.
   * The scout survey adds:
     * the circulant-Hadamard/Barker link (Barker 7 = `{1,2,4} + 2`);
     * the Kraft/Gibbs identity of the 9/4 paper;
     * the prefix-instruction universality class.

## 1. The append maps and Ahmed's identity (restated, with credit; FINITE-EXACT checks)

**Notation.**

* `U(x) = oddpart(3x+1)` on odd `x`.
* Every odd `n > 1` is `n = 2^(r+1) t − 1` with `t` odd. Its first `r` exponents are 1 (the run), and its next exponent `e(n)` is the reset.
* `m_D = (n + 1)/2^D − 1` is `n` with `D` trailing binary ones deleted, so `m_D + 1 = 2^(r+1−D) t`.
* `A_1(x) = 2x + 1` and `R(x) = 4x + 1`. Then `A_1^k(0) = 2^k − 1`, `R^k(0) = (4^k − 1)/3`, and `m_D = A_1^(−D)(n)`.

**Lemma 1 (Ahmed 2016, Theorem 2.1).** Let `n = 2m + 1` have run `r >= 1`. Then

    U^r(n) = 2^(e(m)) · U^r(m) + 1.

In particular `U^r(n) = 4 U^r(m) + 1 = R(U^r(m))` iff `e(m) = 2`, iff `e(n) >= 3`. In that case `U^(r+1)(n) = U^(r+1)(m)` (the reset switch).

*Proof.*

1. `U^r(n) = 2·3^r t − 1`.
2. `m = 2^r t − 1` has run `r − 1` and `U^(r−1)(m) = 2·3^(r−1) t − 1`. Then `3U^(r−1)(m) + 1 = 2(3^r t − 1)`, so `e(m) = 1 + v` with `v = v_2(3^r t − 1)`, and `U^r(m) = (3^r t − 1)/2^v`. Hence `2^(e(m)) U^r(m) + 1 = 2·3^r t − 1`.
3. `e(n) = 1 + v_2(3^(r+1) t − 1) >= 3` iff `3^(r+1) t ≡ 1 (mod 4)` iff `3^r t ≡ 3 (mod 4)` iff `v = 1`.
4. Finally `U(4y + 1) = U(y)`. ∎

**Checks (script §C).** The general form holds on every odd `n < 2·10^5` with run `>= 1`. The sibling case holds exactly when `e(n) >= 3`, and the merge occurs in all 24998 such cases.

**Readings.**

* In general, "append a 1" (`n = A_1(m)`) becomes, after the run, "append `0^(e(m)−1)1`".
  * When `e(m) = 2` the appended block is `01`, which Collatz deletes in one step.
  * When `e(m) >= 3` it becomes the debt state `(1, e(m) − 2, U^r(m))` of `checked_switch_phase19` (4).
* For Mersenne numbers with even `a`, the identity reads `2·3^(a−1) − 1 = 4·(3^(a−1) − 1)/2 + 1`. This is THM-4556 (ii)'s repunit form, seen from the other side.
* **Offsets at length 3.** `A_1^3(n) − n = 7(n + 1)` and `R^3(n) − n = 21(3n + 1)`; the latter is in the mod18/19 clock-tower table. The two maps fix the 2-adic poles `−1` (the limit of `2^a − 1`) and `−1/3` (the limit of `(4^k − 1)/3`), as in THM-4556's remark.

## 2. The parity law (PROVED; FINITE-EXACT checks)

**Definition.** For odd `n` with run `r >= 2`:

* `L(n)` is the set of lags `1 <= D <= r − 1` with `m_D > 1` such that `n` and `m_D` meet at equal time (orbits stopped at 1).
* `L+(n) = L(n) ∪ {0}` if `e(n) >= 3`, and `L+(n) = L(n)` otherwise.
* `D` is **good** if `e(m_D) >= 3`, where `m_0 = n`.

**Theorem 2.**

1. `D` is good iff `(−1)^(r−D+1) ≡ t (mod 4)`. So good and bad lags alternate, and good lags have the parity of `r + 1` when `t ≡ 1 (mod 4)` and of `r` when `t ≡ 3`.
2. For `0 <= D < D + 1 <= r − 1` with `D` good, `D ∈ L+(n)` iff `D + 1 ∈ L+(n)`. Hence `L+(n)` is a union of pairs `{D, D+1}` with `D` good, except possibly an unpaired good top lag `D = r − 1`.
3. The least element of `L(n)` is good, hence odd, for every reset-2 source. For every reset-`>= 3` source, `1 ∈ L(n)`.

*Proof.*

1. **Part 1** is Lemma 1's congruence for `m_D`, whose run is `r − D` and whose `t` is the same.
2. **No meeting during the run.** For `i <= r`, `U^i(m_D) + 1 <= (3/2)^i (m_D + 1) < (3/2)^i (n + 1) = U^i(n) + 1`, because `U(x) + 1 <= (3/2)(x + 1)` with equality on runs. So meetings happen at times `>= r + 1`.
3. **The reset pairs.** Let `D` be good with `D + 1 <= r − 1`. Then `m_D` (run `r − D >= 2`, reset `>= 3`) and `m_(D+1) = (m_D − 1)/2` meet at time `r − D + 1 <= r` by Lemma 1, and agree from then on. So `n` meets `m_D` at a time `>= r + 1` iff it meets `m_(D+1)` there.
4. **Part 3.** For a reset-2 source, 0 is bad, so the good lags are odd. A bad least element `D_0` would pair with the good lag `D_0 − 1 >= 1`, contradicting minimality. For a reset-`>= 3` source, 0 is good, and the pair `{0, 1}` gives `1 ∈ L(n)`: the reset switch. ∎

**Checks (script §B).** For the 2500 odd sources `< 2·10^4` with run `>= 2`:

* the good-parity prediction has 0 failures, and the closure of part 2 has 0 violations;
* 148 sources have an unpaired good top lag;
* each of the 202 (of 1250) reset-2 sources with a partner has odd least lag.

In THM-4555's wider convention (runs `>= 1`, lags `D <= r`), all 377 of the 2500 reset-2 sources below `2·10^4` that have a partner have odd least lag.

**Consequence for the compiler.**

* Trailing-ones searches need only good lags; each hit at `D <= r − 2` gives `D + 1` free.
* A switching reset-2 source hands its debt to the smaller reset-2 source `m_(D+1)` (total deletion `D + 1`, which is even).
* The chain of debts ends at a reset-2 source with no trailing-ones partner. For Mersenne numbers this is a class root (section 4). For random sources it is HYP-9214's unresolved fraction.

## 3. 2x+1 and 4x+1 generate F_21 only at 7 (PROVED; FINITE-EXACT)

**Proposition 3.**

1. Modulo 7, `A_1` (multiplier 2) and `R` (multiplier 4) lie in `F_21 = {x -> ax + b : a ∈ QR_7} = Aut(P_7)` and generate it.
   * `A_1` has cycles `(0 1 3)(2 5 4)` and fixes `6 = −1`.
   * `R` has cycles `(0 1 5)(3 6 4)` and fixes `2 = −1/3`.
2. Modulo a Mersenne prime `p = 2^k − 1` with `k >= 3`, they generate `Z/p ⋊ <2>`, of order `kp`. This lies in `Aut(P_p) = Z/p ⋊ QR_p`, and equals it iff `k = (p − 1)/2`, i.e. `k = 3`.
3. `kp = k(2^k − 1)` equals `R^k(0) = (4^k − 1)/3` iff `2^k + 1 = 3k`, iff `k ∈ {1, 3}`. At `k = 3` this is `2^3 + 1 = 3^2`.
4. Modulo `63 = 9·7` (level 2 of the clock tower), `A_1^6 = id` and `R^3` is translation by 21.

*Proof.*

1. `<2> ⊆ QR_p` because `p ≡ 7 (mod 8)`. Since `A_1^2(x) = 4x + 3`, the map `R ∘ A_1^(−2)` is `x -> x − 2`, which generates `Z/p`.
2. `|Aut(P_p)| = p(p − 1)/2` for prime `p` (the affine maps with square multiplier).
3. `(4^k − 1)/3 = (2^k − 1)(2^k + 1)/3`. ∎

Checked for `k = 3, 5, 7` (orders 21, 155, 889 against 21, 465, 8001), and parts 2–3 for `k <= 300` (script §D).

**Typing.**

* The identities are EXACT. The link to the dynamics is that the maps of Lemma 1 are these two: DICTIONARY.
* Lemma 1 is a pointwise intertwining `U^r(A_1(m)) = R(U^r(m))`, valid where `e(m) = 2` and with `r` depending on `m`. It is not a conjugacy.
* This explains neither the sporadic collisions (Diophantine; THM-4555, remarks) nor which odd exponents switch.
* THM-4556 (iii)'s landing point `2^(2k) − 1 = 3·R^k(0)` (with `63 = 3·21`) involves the same pair of maps.

## 4. Equal-time classes of Mersenne numbers (FINITE-EXACT, NUMERICAL, cited)

**Classes.**

* Orbits that meet at equal time stay together, so the classes form an equivalence relation on which `σ` is constant.
* For odd `a <= 2000` every least-lag merge happens at a value `>= 911`, and even exponents merge with `a − 1` by the reset switch. So the classes are exactly the level sets of `σ(2^a − 1)`, and for `a >= 3` they are unions of reset pairs.
* Class counts for exponents `2..A` are 22, 30, 36, 49, 53, 62, 68, 78, 86 at `A = 100, 200, 400, 800, 1000, 1600, 2000, 3000, 4000`. THM-4556 (vi) and HYP-9213 count from `a = 1`: 23, 37, 58 at `A = 100, 400, 1200`, and 104 at 6000.
* A least-squares log-log fit over `A ∈ [100, 4000]` gives exponent 0.36 (HYP-9213: about 0.37). NUMERICAL.
* The share of odd exponents that switch, in windows of 500 consecutive exponents (250 odd each) up to 4000: 0.847, 0.944, 0.968, 0.972, 0.988, 0.972, 0.992, 0.976.

**Certified versus actual.**

* **Certified shares (FINITE-EXACT).** Templates of total `<= K` certify 0.1199 at `K = 20` (THM-4556 (v); reproduced to `K = 18`). The numba kernel gives 0.1255, 0.1308, 0.1360, 0.1411, 0.1460, 0.1508, 0.1556 for `K = 21, …, 27`.
  * The kernel lists each exponent class twice (classes of `X mod 2^(K+1)` against `a mod 2^(K−2)`). The shares are unaffected.
* **No extrapolation (NUMERICAL).** The increments decay slowly, but their local power-law exponent drifts: about 1.0 near `K = 20` and 0.45 near `K = 27` (audit). So the data cannot distinguish a limit of 1 from saturation.
* **Actual merges (FINITE-EXACT).**
  * For the 543 switching odd `a <= 1200`, the template total at the least-lag merge has median 247; only 0.110 of odd `a` have it `<= 20`.
  * All 932 least-lag switches of odd `a <= 2000` are collisions: their exponent totals differ by exactly `D`.
* **The first template** is `a ≡ 95 (mod 128)`.
  * Post-run words `(2,6,1)` and `(4,1,1,3)`, common value `125/256`, normal forms `(8,1) ~ (4,1,1,3)`: THM-4555's smallest sporadic collision.
  * Checked for `a = 95 + 128t`, `t <= 7`; the audit checked `t < 40`.

**Proposition 6 (PROVED).** Let `S ⊂ Z_2` be the open set of odd 2-adic exponents `a` for which the post-run words of `2·3^(a−1) − 1` and `2·3^(a−1−D) − 1` collide at some finite total, for some odd `D`. If `μ_2(S) = 1`, then HYP-9213 holds.

*Proof.*

1. Given `ε > 0`, choose `K` with `μ_2(S_(<=K)) > 1 − ε`.
2. By THM-4556 (iv), membership in `S_(<=K)` is periodic modulo `2^(K−2)` for `a` larger than a constant. So the odd `a <= A` in `S_(<=K)` have density tending to `μ_2(S_(<=K))`.
3. Each such `a` has an equal-time Mersenne partner below it. So the class roots among odd exponents have upper density `<= ε`.
4. Every even exponent lies in its odd predecessor's class. Hence the number of classes is `o(A)`. ∎

The converse is not claimed.

**HEURISTIC (HYP-9214's model on the Mersenne line).**

* After the run, the orbit behaves like a random walk in `log_2` with drift `log_2(3/4)`. At equal times, the walks of exponents `a` and `a − D` start about `2D` bits apart, and one-dimensional walks recur.
* But an actual merge needs an exact sibling configuration: `U(x) = U(y)` iff `x = R^k(y)` or `y = R^k(x)`. The per-visit merge probability is HYP-9214's open point.

## 5. Mirror clocks are CRT-independent at finite depth (PROVED; empirical census)

**Proposition 5.** Fix `K` and `k`. Among odd exponents `a`:

* the event "`2^a − 1` has a certified trailing-ones switch of template total `<= K`" is periodic modulo `2^(K−2)` (THM-4556 (iv));
* the event "`2^a − 1` has a smaller ancestor through an inverse word of length `<= k`" is periodic modulo `2·3^(k−1)` (THM-4554 (v)).

Their joint natural density is the product of their densities.

*Proof.* On odd `a`, the residues `a mod 2^(K−2)` and `a mod 3^(k−1)` are jointly equidistributed, and each event is a union of classes of one of them. ∎

**Census (script §F).** Backward depth 30 reproduces THM-4554 (vi)'s 35 of 59; switching is taken at unbounded depth.

| odd `a` | `P(switch)` | `P(backward-minimal)` | `P(both)` | product |
|---|---|---|---|---|
| `[3, 121]` | 0.617 | 0.583 | 0.367 | 0.360 |
| `[3, 401]` | 0.825 | 0.590 | 0.495 | 0.487 |
| `[123, 401]` | 0.914 | 0.593 | 0.550 | 0.542 |

The class roots are backward-minimal at the base rate (19 of 35).

**Reading.** At finite depth the two clocks certify independent fractions `p` and `q`. Together they leave a fraction `(1 − p)(1 − q)` of exponents with neither a certified forward switch nor a smaller ancestor. Independence says nothing about infinite depth, and does not by itself say that a proof must use size.

## 6. Collisions at −1, counted (FINITE-EXACT; reading ANALOGY)

The reduced words (first letter `>= 2`) of total `A` number `M = 2^(A−2)`, and `f_u(−1) = N/2^(A−1)` with `N` odd (THM-4555 (ii)). Census (script §G):

| `A` | `M` | values `V` | `V/M` | colliding pairs | sporadic pairs | `H_2 − log_2 M` |
|---|---|---|---|---|---|---|
| 8 | 64 | 32 | 0.500 | 32 | 0 | −1.000 |
| 12 | 1024 | 484 | 0.473 | 624 | 28 | −1.150 |
| 15 | 8192 | 3765 | 0.460 | 5424 | 332 | −1.217 |
| 18 | 65536 | 29442 | 0.449 | 46180 | 3353 | −1.269 |

**Readings.**

* **Up to total 8 the map is exactly 2-to-1.** The root collision `(2, c, ·) ~ (c+2, ·)` is a bijection between the words that begin with 2 and those that begin with `>= 3`.
* **Sporadic pairs first appear at total 9**, as `(8,1) ~ (4,1,1,3)`.
  * Their per-unit growth ratio falls from 2.36 to 2.14 between totals 13 and 18, but stays above `M`'s ratio of 2.
  * So the image keeps falling below half, and the Rényi-2 deficit grows slowly.
* **Interpretation.** One may call this a hash with structured collisions (ANALOGY). No random-hash model was tested.

## 7. openai/math (ANALOGY / DICTIONARY; claims reported as claims)

### 7.1 *Integer multiplication below n log n* (OpenAI, 2026-09-23)

**Claim.** One fixed multitape Turing machine multiplies `n`-bit integers in `O(n (lg n)^(1−κ))` time, with `κ = 2^(−182)`. This would disprove, in that model, Schönhage and Strassen's conjecture that `n log n` is optimal. We read the LaTeX source; nothing was executed.

**Mechanism.**

* **The pipeline.** Harvey–van der Hoeven's Gaussian resampling, Bluestein's chirp, and Nussbaumer's synthetic transforms over `Z[y]/(y^r + 1)`, where roots of unity are signed shifts.
* **Two tape procedures**, needed because scanning `Θ(n)` bits at `Θ(log n)` levels costs `n log n`:
  * an address-field interchange in `O(V u^τ)` through XOR combinations;
  * simultaneous butterfly layers `H_0(u,v) = ((u+v)/2, (u−v)/2)` in `O(V d^(λ'))`.
* **The recursive linear networks.**
  * Each has `W` wires, and each wire carries one scalar.
  * Each source, gate and sink carries an `m × m` frame matrix. The rank cost is `s = Σ_edges rank(M_head − M_tail)`, and the contract is `M_out(ρ(w)) − M_in(w) = I_m` with `s < Wm`.
* **The construction.**
  * It uses all three-element subsets of `[100]`, with `m = 100^3`. Over `F_2`, the identity on triples is the incidence Gram matrix (`|S ∩ T| mod 2`, rank `<= 100`) plus side wires on the intersection-1 pairs.
  * Banks are swapped by three shears.
  * A rank-`a` edge costs `a` smaller interchanges via `A = E_1 Π E_2`, each unit of `Π` being one shear-swap-shear.
  * Scratch values are restored (Buhrman–Cleve–Koucký–Loff–Speelman).
* **Ordering by CRT.** Coefficients are placed on prime cyclic axes by the triangular CRT map `b_i = a_i + μ_i Σ_(j<i) a_j P_j (mod s_i)`, executed as prefix-controlled cyclic shifts.

**Dictionary to the Collatz compiler.** See mac-mini's §2 for the negacyclic clock, Rader and double-base readings.

| multiplication paper | Collatz (this repo) | type |
|---|---|---|
| bank swap = three shears | reset switch = Lemma 1's intertwining, then `R` erased | ANALOGY |
| restored scratch | `U ∘ R = U`: appended `01` blocks vanish; runs become powers of 3 | ANALOGY |
| rank saving `s < Wm`, iterated | Mersenne classes `o(A)` (HYP-9213, OPEN); templates certify 0.156 at total 27 | ANALOGY (shape of saving) |
| triangular CRT by prefix-controlled shifts | THM-4555 (i): the 2-adic prefix `u` fixes the 3-adic offset `z ≡ f_u(−1) (mod 3^j)` | EQUATION-LOOKALIKE (triangular change of axes) |
| roots of unity as signed shifts | the trailing-ones lag `D` as the rotation `2^(−D)` of `Z/(2^a − 1)` | DICTIONARY; it cannot see the parity of `D` (section 9) |
| Gaussian resampling between lengths | the 2 ↔ 3 mirror: THM-4554 (iv), Proposition 5 | ANALOGY |

**What might transfer.** The shape: a fixed finite gadget whose rank sum is strictly below its volume, iterated, gives a power saving.

* The Collatz counterpart would be a finite set of collision templates certifying, at each scale, a family with strictly fewer independent certificates.
* For the Mersenne family, `N(2A)/N(A)` is 1.20–1.36 (script §E), below 2 (NUMERICAL).

### 7.2 The three attached preprints

* **Unrestricted pro-modularity at the prime two** (2026-10-04).
  * *Claim:* every continuous, odd, absolutely irreducible *two-dimensional* 2-adic representation of `G_Q`, unramified outside finitely many primes, occurs in the full completed Hecke algebra `T_2(N) = lim_k Z_2 ⊗ T_(<=k)(N)` at some odd tame level. This is an occurrence result, with no de Rham hypothesis.
  * *Mechanism:* fixed-determinant pseudodeformation families and local block equivalences for `GL_2(Q_2)`.
  * *Bridge:* DICTIONARY only. Weights are interpolated 2-adically by Kummer-type congruences, and the Mersenne exponent `a` is interpolated the same way.
    * Post-run words depend only on `a mod 2^(K−2)` (THM-4556 (iv)), because `a -> 3^a` is 2-adically continuous.
    * So templates are locally constant on that "weight space", and THM-4556 (v)'s classes are open sets in it.
* **Weil classes and Hodge classes on abelian powers** (2026-09-30).
  * *Claim:* the rational Hodge conjecture for every self-power of an abelian sixfold with an imaginary-quadratic action whose Hermitian form is hyperbolic of signature `(3,3)`. Also for every self-power of an abelian variety of dimension `<= 5` admitting an imaginary-quadratic action, with no signature condition.
  * *Bridge:* CONTEXT. The S17 field `Q(sqrt(−7))` acts on `E = 49a1`. For instance `E^6`, with the action through one embedding on three factors and the conjugate on the other three, has signature `(3,3)`; as a CM power it lies in the companion CM case. No Collatz content.
* **Weak mixing of triangular billiards with an irrational angle** (2026-10-05).
  * *Claim:* weak mixing for every triangle with an angle irrational relative to `π`.
  * *Mechanism:* a measurable eigenfunction `Xf = iλf` is expanded in angular Fourier modes. Energy identities telescope in `λ` to `||Xf||^2 + ||Yf||^2 <= λ^2 ||f||^2`, so `Yf = 0`, and rigidity finishes.
  * *Bridge:* ANALOGY ("no eigenfunction" vs "no hidden clock"). The Collatz notes have found only the two mirror clocks, and Proposition 5 makes them independent at finite depth. EQUATION-LOOKALIKE: the telescoping resembles the procgen three-mirrors character identity.

### 7.3 Survey of the repository (scout pass)

The scout read all 722 abstracts and family descriptions, and about 25 sources in detail; nothing was executed. **No preprint addresses Collatz, `|2^A − 3^l|`, linear forms in logarithms, or an irrationality measure of `log_2 3`.** Ranked bridges (as claims):

1. **The circulant Hadamard conjecture** (family 179).
   * *Claim:* real circulant Hadamard matrices exist only in orders 1 and 4, so the Barker lengths are 2, 3, 4, 5, 7, 11, 13.
   * *Mechanism:* Turyn's descent, and an alternating product `Δ(x)` of character values that Kronecker's theorem forces to be a root of unity.
   * *Checked here (script §I):*
     * Barker 7's minus set is `{3,4,6} = {1,2,4} + 2`, the trivial cycle's real code and the `(7,3,1)` Fano set;
     * Barker 11's plus set is `QR_11 + 8`, the Paley `(11,5,2)` set;
     * Barker 13's minus set is a `(13,4,1)` Singer set.
   * The Collatz-side finiteness is the paragraph after the Paley bridge note's Proposition 3: a single doubling orbit is a nontrivial difference set only at length 3.
   * TECHNIQUE candidate: the alternating-product rigidity for Gauss periods of codes.
2. **The 9/4 bound for ω** (family 107).
   * The shared-leg entropy inequality runs on `max_q e^(H(q)) Π A_i^(q_i) = Σ A_i` (script §I).
   * This is the same convex duality as the refuel-bill Kraft sums and THM-4554 (iv)'s generating function. EQUATION-LOOKALIKE.
3. **Finite tensor savings and exact Fourier circuits** (family 130), with 109.
   * One finite saving amplified, with CRT and uncomputation. The existence of a finite win is proved by duality with a "price" functional.
   * TECHNIQUE (shape). Compare the Bellman/Kraft dualities.
4. **Prefix instructions and incompressible flows; incompressible box transport** (family 376).
   * Turing-universal generalized shifts `x -> B^k x + c`, with a determinant-one recorder ("incompressible" means divergence-free).
   * Collatz is outside this radix-power class, since 3 is coprime to 2. Its determinant-one completion is the `×2×3` solenoid extension. CONTEXT, but precise.
5. **Artin primitive roots for every admissible base** (family 029).
   * *Claim:* at least `c_a x/(log x)^2` primes in `(x, 2x)` with `a` as a primitive root.
   * *Bridges:*
     * the repo's Artin coordinate (HYP-9174);
     * the clock tower's primes `2·3^k + 1` (mac-mini's §2 corrects that list);
     * S17 Proposition 4's denominators `m = 7p`, conditionally on the claim holding for `p ≡ 2 (mod 3)`.
6. **The irrationality exponent of π is 2** (family 017).
   * The collection's only irrationality-exponent method; it does not transfer to `log_2 3` as written.
   * Even `μ(log_2 3) = 2` would be ineffective, and useless for finite cycle exclusion.
7. **The entropy-rate dimension formula** (family 148). Collatz block maps never coincide as maps, so the repo's collisions are coincidences of one orbit point. CONTEXT.
8. **Context only:**
   * Steinitz–Bergström prefix discrepancy (097);
   * Bernoulli convolutions (153), where `λ = 2/3` remains open;
   * Catalan's constant (005), via Frobenius certificates;
   * generalized star height (134), a prefix-code compiler cousin, though the Collatz language is not regular;
   * the Thorp shuffle (238), Thompson's group F (248), dyadic Erdős similarity (084), Deligne–Drinfeld (008), mean-payoff games (104).

## 8. Toward proofs: what moved, and next steps

**Moved.**

* The odd least lag is now a theorem for every reset-2 source (Theorem 2), with exact boundary behaviour, on top of Ahmed's identity.
* The `{7, 21}` question has an exact dictionary answer (Proposition 3).
* The Mersenne census is placed precisely:
  * the certified 2-adic part reaches 0.1556 at total 27, with proved infinite families;
  * collisions occur everywhere in the census;
  * HYP-9213 follows from the 2-adic statement `μ_2(S) = 1` (Proposition 6, one direction only).

**Not moved.** Collatz; HYP-9213; HYP-9214.

**Next steps (ranked).**

1. **A local merge lemma.** Bound from below the per-visit probability that two orbits at equal times reach an exact sibling configuration (HYP-9214).
2. **Certified density past `K = 27`**, refining only unresolved classes. Is `μ_2(S) = 1`?
3. **Theorem 2 in the compiler** (`checked_switch_phase19`): good-lag search, paired lags, debt hand-off.
4. **A certificate-rank inequality** of the multiplication paper's shape, for a template-defined subfamily.

## 9. Hostile checks and numerology

* **The parity law is independent of 7.** Proposition 3 says only that the two maps of Lemma 1 generate the full Paley symmetry at `p = 7`. `21 = |F_21| = R^3(0)` is the `k = 3` solution of `2^k + 1 = 3k` (EXACT; a two-solution coincidence).
* **The shift dictionary cannot see parity.** In `Z/(2^a − 1)` with `a` odd, the lags `D` and `D + a` give the same rotation. So the parity law is the 2-adic congruence `3 ≡ −1 (mod 4)`, not a Mersenne-ring fact.
* **7 and 21 are class roots, but so are all odd `a <= 27` except 13:** a short-orbit effect. NUMEROLOGY.
* **"37 of 60" depends on the window.** At larger `a` the share is 0.95–0.99, and no class rests on simultaneous arrival at 1 (section 4).
* **The `a ≡ 95 (mod 128)` family:** `128 = 2^(K−2)` at `K = 9`, and `125 = 2^7 − 3`. These 7s come from the clock, not from Fano.
* **The multiplication paper's triples are all three-element subsets of `[100]`, not a design.** A `{7, 21}` reading of the gadget is at best EQUATION-LOOKALIKE.
* **Propositions 5 and 6 are finite-depth or one-directional.** The infinite-depth switching set is open in `Z_2`, and integer densities need not exist.

## 10. Audit

An independent adversarial audit (2026-10-06, blind subagent, own code: Python and numba, with read-only checks of the openai/math sources) confirmed:

* Lemma 1 (no failures in 50000 sources);
* Theorem 2's proof and data (clean to `2·10^6`);
* Proposition 3, including `|Aut(P_p)| = p(p−1)/2` by backtracking for `p <= 31`;
* Proposition 5;
* every census number;
* the certified densities for `K = 9..27`, by a separate kernel, and for `K <= 15` also by actual big-integer merges;
* the cross-tab at depths 30 and 60, the collision table, and the `a ≡ 95 (mod 128)` family for `t < 40`.

**Errors found and corrected:**

1. "Every least-lag merge … `>= 911`" holds for odd `a` only.
2. Theorem 2 as first stated failed at two boundaries: the lag 0 for reset-`>= 3` sources, and an unpaired good top lag `r − 1` (148 sources). It also needs `r >= 2`.
3. The windows were 500 consecutive exponents, not 500 odd ones.
4. The `A^0.35` fit was a two-point estimate; least squares gives 0.36–0.37.
5. Smaller slips:
   * the Mersenne reading of Lemma 1 holds for even `a` only;
   * the conventions for the class counts were mixed;
   * "unions of reset pairs" needs `a >= 3`.

**Retypings and relabelings:**

* Lemma 1 relabelled as Ahmed's theorem.
* The 7/21 link retyped as DICTIONARY, with "conjugation" replaced by "pointwise intertwining".
* Proposition 5 stripped of "barrier" and "only through size", the overclaims MISTAKE-569 had removed elsewhere.
* "HYP-9213 is essentially `μ_2(S) = 1`" replaced by the one-directional Proposition 6.
* The `K^(−0.6)` extrapolation withdrawn.
* "Hash" and "birthday" typed as ANALOGY, with the growth ratio noted as falling.
* "`o(A)` certificates" marked OPEN.

**Script repairs:**

* the claims not yet in the script (911, windows, fit, roots, ratios, Gibbs, collisions) were added;
* the budget fallback was fixed.

**Citation fixes:**

* pro-modularity is for two-dimensional representations, and its method does not use patching;
* the dimension-`<= 5` Hodge statement has no signature condition;
* frames sit on vertices, not wires;
* the multiplication paper would disprove an optimality conjecture, not refute a lower bound;
* the Paley bridge reference is now located;
* THM-4556's `F_21` remark and the clock-tower table are credited.

All of this is logged as MISTAKE-573.

## 11. Reproduction and references

**Reproduction.**

* `python3 04-computation/experiments/mersenne_switch_parity_f21_compression_20261006.py` (about 30 s), ending `ALL CHECKS PASSED`.
* `python3 04-computation/experiments/mersenne_certified_density_numba_20261006.py 9 27` (needs numba; under 20 s).

**Repo inputs.**

* `checked_switch_phase19_20261004.md` (3)–(5);
* THM-4553, THM-4554, THM-4555, THM-4556, THM-4557;
* HYP-9213, HYP-9214;
* `mod18_mod19_seven_sixtythree_fractal_20261006.md`;
* `seven_twentyone_mersenne_openai_math_20261006.md`;
* `collatz_paley_bridge_20261001.md`;
* `ramanujan_heegner7_trivial_cycle_20261006.md` (S17).

**Literature.**

* M. M. Ahmed, *The intricate labyrinth of Collatz sequences*, arXiv:1602.01617 (2016), Theorem 2.1.
* A. Hanaki (2020), credited in THM-4557 (MISTAKE-571).
* A. Schönhage and V. Strassen (1971).
* D. Harvey and J. van der Hoeven, *Integer multiplication in time O(n log n)*, Ann. of Math. (2021).
* H. Nussbaumer (1980).
* H. Buhrman, R. Cleve, M. Koucký, B. Loff and F. Speelman, *Computing with a full memory: catalytic space* (STOC 2014).
* R. H. Barker (1953) and R. J. Turyn, on Barker sequences.

**openai/math preprints** (https://github.com/openai/math; unrefereed, model-written), cited by directory under `preprints/`:

* *Integer-multiplication-below-n-log-n-September-23-2026*;
* *The-circulant-Hadamard-conjecture-September-23-2026*;
* *Matrix-Multiplication-Nine-Fourths-October-2-2026*;
* *Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026*;
* *Prefix-Instructions-and-Incompressible-Flows-September-27-2026* and *Incompressible-Box-Transport-and-Finite-Computation-September-27-2026*;
* *Primitive-roots-for-every-admissible-integer-base-October-4-2026*;
* *The-irrationality-exponent-of-pi-is-2-September-24-2026*;
* *The-entropy-rate-dimension-formula-for-self-similar-measures-on-the-line-September-24-2026*;
* the three attached preprints.
