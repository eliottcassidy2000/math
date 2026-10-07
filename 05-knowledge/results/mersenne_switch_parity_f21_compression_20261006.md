# Seven and twenty-one are the two binary-append maps: the run turns "append 1" into "append 01", so trailing-ones switch lags come in pairs {D, D+1} for every source (least lag odd for every reset-2 source, not only Mersenne numbers); 2x+1 and 4x+1 generate F_21 = Aut(P_7) exactly at 7; the Mersenne mirror clocks are CRT-independent; collisions at −1 are birthday compression; the openai/math multiplication paper's rank saving read against the certificate compiler

2026-10-06, session opus-2026-10-06-S18 (worktree `codex/session-mersenne-f21-20261006`).

Owner's prompt (verbatim): "think of how ideas ananlgous to {7,21} relate deeply with the idea that every odd-exponent Mersenne switch has an odd shift distance D, and 37 of the 60 odd exponents in [3,121] produce collision. also explore deeply and intricately this repository of work and look to blend its ideas with ours in useful ways and extend concepts meaningfully toward proofs https://github.com/openai/math thinking in terms of information compression creatively especially https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026 explore around there very open mindenly, pursuing any possible connection between themes or similar looking equations, even if the topics seem like they could not be related". Three preprints from that repository were attached: *Weil classes and Hodge classes on abelian powers*, *Unrestricted pro-modularity at the prime two*, and *Weak mixing of triangular billiards with an irrational angle*.

**Concurrent work on the same prompt.** The mac-mini session `mac-mini-2026-10-06-mod1819` answered the same prompt in parallel and pushed it during this session:

* [THM-4556](../../01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md): the Mersenne line is a chain of debt states. It proves the odd least shift for Mersenne numbers and 2-adic periodicity, and certifies a density of at least 0.1199.
* [THM-4557](../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md): `Aut(T_k) = F_21`.
* [HYP-9213](../hypotheses/HYP-9213-mersenne-collatz-trajectories-coalesce-plateaus-of-odd-step-time.md): Mersenne coalescence.
* [HYP-9214](../hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md): reset-2 debt resolves with probability tending to 1.

This note does not repeat those results. It cites them, reproduces their key numbers independently where it uses them (sections A, E, H of the script), and adds the items marked *new* below.

**Status.**

* **PROVED (elementary):**
  * Lemma 1 (*new*, all sources): the sibling form of the reset switch.
  * Theorem 2 (*new*, all sources): the parity law for trailing-ones lags. THM-4556 (iii) is the Mersenne case `t = 1`, proved there by the same composition.
  * Proposition 3 (*new*): the two append maps generate `F_21 = Aut(P_7)`, and only at 7.
  * Proposition 5 (*new*): CRT independence of the forward and backward Mersenne clocks at every finite depth.
* **FINITE-EXACT:** every check in [the script](../../04-computation/experiments/mersenne_switch_parity_f21_compression_20261006.py), with output [`.out`](mersenne_switch_parity_f21_compression_20261006.out) ending `ALL CHECKS PASSED`.
  * Reproduced: THM-4555 (vi) (37 of 60); THM-4556 (v) (certified shares to `K = 18`) and (vi) (class counts); THM-4554 (vi) (35 of 59 backward-minimal, depth 30).
  * New: the parity law on all 2500 odd sources `< 2·10^4` with run `>= 2`; the sibling form on all odd `n < 2·10^5`; the cross-tab of section 5; the collision census at −1 for totals `<= 18`; template totals of actual switches; and the large-`a` check of the `a ≡ 95 (mod 128)` family.
* **HEURISTIC:** the coalescing-walk reading of the plateaus (section 4), which is HYP-9214's model.
* **ANALOGY / DICTIONARY:** every bridge to the openai/math preprints (section 7). These preprints are unrefereed, model-written, and partly unformalized (their README says so). Their claims are reported as claims.
* **Collatz: OPEN.**

## 0. Answers in brief

1. **Where the prompt's facts come from.**
   * "Every odd-exponent Mersenne switch has an odd shift distance `D`" and "37 of the 60 odd exponents in `[3,121]`" are THM-4555 (vi). There a *switch* is an equal-time merge of `2^a − 1` with `2^(a−D) − 1` = `2^a − 1` with `D` trailing binary ones deleted.
   * `{7, 21}` is the Paley tournament `P_7` and its automorphism group `F_21 = C_7 ⋊ C_3`: THM-4553's ladder, and THM-4557 for the tower `T_k`.
2. **The mechanism behind "odd `D`" is one identity, and it is about 7 and 21.**
   * Put `A_1(x) = 2x + 1` ("append a binary 1"; `7 = A_1^3(0) = 111_2`) and `R(x) = 4x + 1` ("append 01"; `21 = R^3(0) = 10101_2`).
   * Collatz erases `R` in one step: `U(4x + 1) = U(x)`, the sibling identity. It transports `A_1` along runs.
   * **Lemma 1.** The run map turns `A_1` into `R`. If `n = 2m + 1` has run `r >= 1`, then `U^r(n) = 4·U^r(m) + 1` exactly when `n`'s reset exponent is at least 3. One more step erases the `R`, which is the reset switch `n ⇒ (n − 1)/2`.
   * **Theorem 2 (parity law, every source).** For any odd `n = 2^(r+1) t − 1`, the lags `D` for which `n` meets `(n + 1)/2^D − 1` at equal time come in consecutive pairs `{D, D+1}`. The first member `D` of each pair has a fixed parity: `D ≡ r + 1 (mod 2)` if `t ≡ 1 (mod 4)`, and `D ≡ r (mod 2)` if `t ≡ 3 (mod 4)`.
     * So the least lag is odd for every reset-2 source, and equals 1 for every reset-`>= 3` source.
     * Odd-exponent Mersenne numbers are reset-2 sources with `t = 1`, and their good lags are odd: THM-4556 (iii).
     * Checked on all 2500 odd sources below `2·10^4` with run `>= 2`.
3. **Why 7 and 21 in particular (Proposition 3).**
   * Modulo 7, `A_1` and `R` are automorphisms of `P_7`, and they generate all of `F_21 = Aut(P_7)`.
   * Modulo a Mersenne prime `p = 2^k − 1` they generate `Z/p ⋊ <2>`, of order `kp`. This is all of `Aut(P_p)` only when `<2> = QR_p`, i.e. `2^(k−1) − 1 = k`, i.e. `k = 3`. This is the same equation as the Paley bridge note's Proposition 3.
   * The coincidence `|F_21| = 21 = R^3(0)` is `(4^k − 1)/3 = k(2^k − 1)`, i.e. `2^k + 1 = 3k`, whose only solutions are `k = 1, 3`.
   * So `{7, 21}` is the unique level at which "append 1" and "append 01" generate the full Paley symmetry, and the reset switch is the Collatz shadow of their relation.
4. **"37 of 60" is a compression statement, and mostly a deep one.**
   * The switches are exactly the coincidences of the odd-step time `σ` on Mersenne numbers. For `a <= 2000`, every least-lag merge happens at a value `>= 911`, never only at 1.
   * `σ(2^a − 1)` takes only 68 values for `a <= 2000` (HYP-9213: 104 for `a <= 6000`).
   * The certified fraction from finite 2-adic templates is 0.1199 at template total 20 (THM-4556 (v); reproduced). A parallel kernel takes it to 0.1556 at total 27. But the template total at the actual least-lag merge has median 247.
   * So the bulk of the 0.95–0.98 empirical switching rate is late coalescence, not short collisions. HYP-9214 sees the same split for random sources.
   * The smallest template, `a ≡ 95 (mod 128)`, runs on THM-4555's smallest sporadic collision: post-run words `(2,6,1)` and `(4,1,1,3)`, both `125/256`.
5. **Mirror clocks are CRT-independent (Proposition 5).**
   * Forward switching at bounded template total is periodic in `a` modulo a power of 2 (THM-4556 (iv)). Backward-minimality to depth `k` is periodic modulo `2·3^(k−1)` (THM-4554 (v)).
   * So at every finite depth the two events are independent among odd exponents. The census agrees: `P(both) = 0.367, 0.495, 0.550` against products `0.360, 0.487, 0.542`.
6. **Collisions at −1 are hash collisions.**
   * The map `u -> f_u(−1)` on the `2^(A−2)` reduced words of total `A` is about 2-to-1: the root collision `(2, c, ·) ~ (c+2, ·)` pairs the words exactly.
   * Sporadic pairs take the image below half: `0.449·M` at `A = 18`. Their count grows by about 2.1 per unit of total, slightly faster than the word count. The Rényi-2 entropy deficit grows slowly, from 1.00 to 1.27 bits.
7. **openai/math (section 7).**
   * The multiplication paper saves a power of `log` by a rank inequality `s < Wm` for a finite exchange network built from three-element subsets, applied recursively.
   * Its ingredients each have a Collatz counterpart in our compiler, typed in a dictionary:
     * three-shear bank swaps ↔ the reset switch;
     * triangular CRT by prefix-controlled shifts ↔ THM-4555 (i)'s endpoint progression;
     * signed shifts as roots of unity ↔ trailing-ones lags as rotations in `Z/(2^a − 1)`;
     * restored scratch ↔ `U ∘ R = U`.
   * The transferable idea is the *shape* of the saving: a strict inequality in a finite gadget, iterated. The Mersenne plateaus are the one place where our data show a power saving in certificates.
   * The three attached preprints give a dictionary (2-adic weight families ↔ Mersenne exponent families), context (`Q(sqrt(−7))`), and an analogy (no hidden clock ↔ weak mixing).

## 1. The two append maps and the sibling form of the reset switch (PROVED + FINITE-EXACT)

**Notation.**

* `U(x) = oddpart(3x+1)` on odd `x`.
* Every odd `n > 1` is `n = 2^(r+1) t − 1` with `t` odd. Its first `r` exponents are 1 (the run), and its next exponent `e(n)` is the reset.
* `m_D = (n + 1)/2^D − 1` is `n` with `D` trailing binary ones deleted, so `m_D + 1 = 2^(r+1−D) t`.
* `A_1(x) = 2x + 1` and `R(x) = 4x + 1`. Then `A_1^k(0) = 2^k − 1`, `R^k(0) = (4^k − 1)/3`, and `m_D = A_1^(−D)(n)`.

**Lemma 1 (sibling form; PROVED).** Let `n = 2m + 1` have run `r >= 1`. Then

    U^r(n) = 4·U^r(m) + 1      ⟺      e(n) >= 3,

and in that case `U^(r+1)(n) = U^(r+1)(m)`.

*Proof.*

1. During its run, `U^i(n) = 3^i 2^(r+1−i) t − 1`, so `U^r(n) = 2·3^r t − 1`.
2. `m = 2^r t − 1` has run `r − 1`, so `U^(r−1)(m) = 2·3^(r−1) t − 1`. Then `3U^(r−1)(m) + 1 = 2(3^r t − 1)`, so `m`'s reset exponent is `1 + v` with `v = v_2(3^r t − 1)`, and `U^r(m) = (3^r t − 1)/2^v`.
3. Hence `4U^r(m) + 1 = 2·3^r t − 1` iff `v = 1`, i.e. iff `m`'s reset is exactly 2. (For `v >= 2` the left side is at most `3^r t`.)
4. And `e(n) = 1 + v_2(3^(r+1) t − 1)`. One checks `e(n) >= 3` iff `3^(r+1) t ≡ 1 (mod 4)` iff `3^r t ≡ 3 (mod 4)` iff `v = 1`.
5. Finally `U(4y + 1) = U(y)`. ∎

**Reading.**

* "Append a 1" (`n = A_1(m)`) becomes, after the run, "append 01" (`U^r(n) = R(U^r(m))`), and Collatz deletes a trailing `01` in one step.
* This is the reset switch of `checked_switch_phase19` (3): root case of THM-4555 (iv), Ahmed 2016, Theorem 2.1.
* For Mersenne numbers it reads `2·3^(a−1) − 1 = 4·(3^(a−1) − 1)/2 + 1`, i.e. THM-4556 (ii)'s repunit form `U(R_(a−1)) = oddpart(R_a)` seen from the other side.
* Checked: the equivalence on every odd `n < 2·10^5` with run `>= 1`, and the merge in all 24998 reset-`>= 3` cases (script §C).

**Offsets at length 3 (PROVED, trivial).** `A_1^3(n) − n = 7(n + 1)` and `R^3(n) − n = 21(3n + 1)`. The two maps fix the two 2-adic poles `−1` (the limit of `2^a − 1`) and `−1/3` (the limit of `(4^k − 1)/3`). This is THM-4556's remark "seven and twenty-one".

## 2. The parity law for every source (PROVED + FINITE-EXACT)

**Definition.** For odd `n` with run `r >= 2`, let `L(n)` be the set of lags `1 <= D <= r − 1` with `m_D > 1` such that `n` and `m_D` meet at equal time (`U^i(n) = U^i(m_D)` for some `i`, orbits stopped at 1). Call `D` **good** if `e(m_D) >= 3`.

**Theorem 2.**

1. `D` is good iff `(−1)^(r−D+1) ≡ t (mod 4)`. So good and bad lags alternate. Good lags have the parity of `r + 1` when `t ≡ 1 (mod 4)`, and of `r` when `t ≡ 3 (mod 4)`.
2. `L(n)` is a union of pairs `{D, D+1}` with `D` good (for `D + 1 <= r − 1`). In particular, for a reset-2 source `min L(n)` is good.
3. Hence the least lag is odd for every reset-2 source. For every reset-`>= 3` source the pair `{0, 1}` is the reset switch, and `1 ∈ L(n)`.

*Proof.*

1. Part 1 is Lemma 1's congruence applied to `m_D`, whose run is `r − D` and whose `t` is unchanged.
2. **No meeting during the run.** For `i <= r`, `U^i(m_D) + 1 <= (3/2)^i (m_D + 1) < (3/2)^i (n + 1) = U^i(n) + 1`. This holds because `U(x) + 1 <= (3/2)(x + 1)`. So any equal-time meeting has `i >= r + 1`.
3. **The reset pairs.** If `D` is good, then `m_D = 2m_(D+1) + 1` has run `r − D >= 1` and reset `>= 3`. By Lemma 1, `m_D` and `m_(D+1)` meet at time `r − D + 1 <= r`, and they agree from then on.
4. **Closure.** So `n` meets `m_D` at some time `i >= r + 1` iff it meets `m_(D+1)` there. This gives `D ∈ L(n) ⟺ D + 1 ∈ L(n)`.
5. **Parity of the least lag.** For a reset-2 source, `D = 0` is bad, so good lags are odd.
6. **Mersenne case.** For `n = 2^a − 1` (`t = 1`, `r = a − 1`), good lags have the parity of `a`, as in THM-4556 (iii). ∎

**Checks (script §B).** For the 2500 odd sources `< 2·10^4` with run `>= 2`:

* the good-parity prediction has 0 failures, and the pair closure 0 violations;
* each of the 202 (of 1250) reset-2 sources that has a partner in this range has odd least lag.

THM-4555 (vi)'s count, 377 of 2500 reset-2 sources, uses a wider convention (runs `>= 1`, all `D <= r`).

**Consequence for the compiler.**

* A search for trailing-ones switches needs only good lags. Each hit at `D` gives `D + 1` for free.
* Every reset-2 source that switches hands its debt to the smaller reset-2 source `m_(D+1)`, with even total deletion `D + 1 >= 2`.
* The chain of debts ends at a reset-2 source with no trailing-ones partner. For Mersenne numbers this is a root of a `σ`-class; for random sources it is HYP-9214's unresolved fraction.

## 3. Seven and twenty-one as the generators of F_21 (PROVED + FINITE-EXACT)

**Proposition 3.**

1. Modulo 7, `A_1` (multiplier 2) and `R` (multiplier 4) lie in `F_21 = {x -> ax + b : a ∈ QR_7} = Aut(P_7)` and generate it.
   * `A_1` has cycles `(0 1 3)(2 5 4)` and fixes `6 = −1`.
   * `R` has cycles `(0 1 5)(3 6 4)` and fixes `2 = −1/3`.
2. Modulo a Mersenne prime `p = 2^k − 1`, `k >= 3`, the maps `A_1` and `R` generate `G_k = Z/p ⋊ <2>`, of order `kp`. This lies in `Aut(P_p) = Z/p ⋊ QR_p` (`p ≡ 7 (mod 8)`), and equals it iff `k = (p − 1)/2`, i.e. `k = 3`.
3. `|G_k| = kp = k(2^k − 1)` equals `R^k(0) = (4^k − 1)/3` iff `2^k + 1 = 3k`, iff `k ∈ {1, 3}`.
4. Modulo `63 = 2^6 − 1 = 9·7` (the clock tower's level 2), `A_1^6 = id` and `R^3` is translation by 21.

*Proof.*

1. **The group.** Both multipliers lie in `<2>`, which is contained in `QR_p` because `p ≡ 7 (mod 8)`. Since `A_1^2(x) = 4x + 3`, the composite `R ∘ A_1^(−2)` is the translation `x -> x − 2`, which generates `Z/p` for odd prime `p`. Hence `<A_1, R> = Z/p ⋊ <2>`. `|Aut(P_p)| = p(p−1)/2` for prime `p` (affine maps with square multiplier).
2. **Part 3.** `(4^k − 1)/3 = (2^k − 1)(2^k + 1)/3`.
3. **Part 4.** Direct. ∎

Checked for `k = 3, 5, 7` (orders 21, 155, 889 against 21, 465, 8001), and for `k <= 300` in parts 2–3 (script §D).

**What this explains, and what it does not.**

* It explains why the symmetry of 7 and the trailing-ones calculus meet: the reset switch is the conjugation of `A_1` into `R` by the run (Lemma 1), and these two maps generate the Paley symmetry exactly at `p = 7`.
* It also places THM-4556 (iii)'s landing point `2^(2k) − 1 = 3·R^k(0)` ("bipolar", with `63 = 3·21`) inside the same pair of maps.
* It does not explain the sporadic collisions, which are Diophantine coincidences (THM-4555 Remarks), nor which odd exponents switch.
* Typed: EXACT identities; the link to the Collatz dynamics is the conjugation of Lemma 1; the rest is DICTIONARY.

## 4. Equal-time classes: compression of the Mersenne family (FINITE-EXACT + cited)

**Classes and their counts.**

* Orbits that meet at equal time stay together, so equal-time meeting is an equivalence relation, and `σ` is constant on its classes.
* For Mersenne numbers with `a <= 2000`, every non-root exponent is linked to a smaller one by a merge at a value `>= 911`. So the classes are exactly the level sets of `σ(2^a − 1)`, and they are unions of reset pairs.
* Class counts: 22, 36, 57, 68 for exponents `2..A`, `A = 100, 400, 1200, 2000`. THM-4556 (vi) counts 23, 37, 58 including `a = 1`; HYP-9213 gives 104 at 6000.
* In windows of 500 odd exponents, the share switching rises 0.85, 0.94, 0.97, 0.97, 0.99, 0.97, 0.99, 0.98 up to 4000. A power-law fit of the class count gives `A^0.35` (HYP-9213: `A^0.37`).

**Certified versus actual.**

* THM-4556 (v) certifies switching for a share 0.1199 of odd exponents by finite 2-adic templates (total `<= 20`). The script reproduces the shares through `K = 18`.
* A parallel kernel ([numba script](../../04-computation/experiments/mersenne_certified_density_numba_20261006.py), [`.out`](mersenne_certified_density_numba_20261006.out)) extends this to `K = 27`: `0.1199, 0.1255, 0.1308, 0.1360, 0.1411, 0.1460, 0.1508, 0.1556` for `K = 20, …, 27` (FINITE-EXACT).
  * The increments decay only like `K^(−0.6)` (a two-point fit, NUMERICAL), so the certified share is still rising steadily.
* For the 543 switching odd `a <= 1200`, the template total at the actual least-lag merge has median 247; only 0.110 of odd `a` have it `<= 20`.
* So almost all switches happen late, as long orbits coalesce. But every one is still a template, i.e. a collision (THM-4555), only of large total.
* **Reformulation (structural, conditional on uniformity).** Let `S ⊂ Z_2` be the open set of odd 2-adic exponents whose post-run words (of `2·3^(a−1) − 1` and `2·3^(a−1−D) − 1`) collide at some finite total, for some odd `D`. An integer `a` can only use templates of total `O(a)`, because its orbit has `≈ 4.8a` steps. So the switching share near `A` approximates `μ(S ∩ {total <= cA})`, and HYP-9213 is essentially the statement `μ_2(S) = 1`.
  * The two-point fit would put the certified share near 0.95 at totals of a few hundred, the same order as the observed median of 247. Typed HEURISTIC.
* The smallest template is `a ≡ 95 (mod 128)`: post-run words `(2,6,1)` and `(4,1,1,3)` with common value `125/256`. Their root-collision normal forms are `(8,1) ~ (4,1,1,3)`, THM-4555's smallest sporadic collision. This was checked directly for `a = 95 + 128t`, `t <= 7` (script §H).

**HEURISTIC reading** (HYP-9214's model, applied to the Mersenne line).

* After its run, the orbit of `2^a − 1` behaves like a random walk in `log_2` with drift `log_2(3/4)` per odd step.
* At equal times the walks of exponents `a` and `a − D` start about `2D` bits apart. One-dimensional walks recur, so a fixed lag meets with probability tending to 1 as the orbit length `≈ 4.8a` grows.
* Chains through larger exponents make class minima rarer still. This is consistent with a class count growing slower than `sqrt(A)`.
* This is a heuristic, not a proof; the per-visit merge probability is HYP-9214's open point.

## 5. The mirror clocks are CRT-independent (PROVED at finite depth + FINITE-EXACT)

**Proposition 5.** Fix a template depth `K` and a backward depth `k`. Among odd exponents `a`:

* the event "`2^a − 1` has a certified trailing-ones switch of total `<= K`" is periodic modulo `2^(K−2)` (THM-4556 (iv));
* the event "`2^a − 1` has a smaller ancestor via an inverse word of length `<= k`" is periodic modulo `2·3^(k−1)` (THM-4554 (v)).

Their joint natural density is the product of their densities.

*Proof.* On odd `a`, the residues `a mod 2^(K−2)` and `a mod 3^(k−1)` are jointly equidistributed (CRT), and each event is a union of classes of one of them. ∎

**Census (script §F, depth 30, which reproduces THM-4554 (vi)'s 35 of 59).** Using the actual (unbounded-depth) switching status:

| odd `a` | `P(switch)` | `P(backward-minimal)` | `P(both)` | product |
|---|---|---|---|---|
| `[3, 121]` | 0.617 | 0.583 | 0.367 | 0.360 |
| `[3, 401]` | 0.825 | 0.590 | 0.495 | 0.487 |
| `[123, 401]` | 0.914 | 0.593 | 0.550 | 0.542 |

The class roots are backward-minimal at the base rate: 19 of 35.

**Reading.** The Mersenne exponent carries the forward and the backward information on independent axes: its 2-adic digits (forward, the clock of 3) and its 3-adic digits (backward, the clock of 2). They interact only through the size of `2^a − 1`, which is what the proof programme would need to exploit (THM-4554, remarks).

## 6. Collisions at −1 are hash collisions (FINITE-EXACT)

The reduced words (first letter `>= 2`) of total `A` number `M = 2^(A−2)`, and each has `f_u(−1) = N/2^(A−1)` with `N` odd (THM-4555 (ii)). Census (script §G):

| `A` | `M` | values `V` | `V/M` | colliding pairs | sporadic pairs | `H_2 − log_2 M` |
|---|---|---|---|---|---|---|
| 8 | 64 | 32 | 0.500 | 32 | 0 | −1.000 |
| 12 | 1024 | 484 | 0.473 | 624 | 28 | −1.150 |
| 15 | 8192 | 3765 | 0.460 | 5424 | 332 | −1.217 |
| 18 | 65536 | 29442 | 0.449 | 46180 | 3353 | −1.269 |

**Reading.**

* Up to total 8 the map is exactly 2-to-1: the root collision `(2, c, ·) ~ (c+2, ·)` is a bijection between the words that start with 2 and those that start with `>= 3`.
* Sporadic pairs, which join words in different root classes, first appear at `A = 9` (`(8,1) ~ (4,1,1,3)`). Their count grows by about 2.1 per unit of total, slightly faster than `M`.
* So the image falls below half, and the collision (Rényi-2) entropy falls short of `log_2 M` by a slowly growing deficit.
* This is the information-compression form of THM-4555's sporadic collisions. Uniform rewrites need such collisions, and they become relatively more common with total length, but only by a small factor per unit (compare D-ε of the mod18/19 note).

## 7. openai/math: the multiplication paper and the attached preprints (ANALOGY / DICTIONARY; claims reported as claims)

### 7.1 *Integer multiplication below n log n* (OpenAI, September 23, 2026)

**Claim.** A deterministic multitape Turing machine multiplies `n`-bit integers in `O(n (lg n)^(1−κ))`, with `κ = 2^(−182)`. This would refute an `Ω(n log n)` lower bound in that model. We read the LaTeX source; nothing was executed.

**Mechanism.**

* **Starting point.** Harvey–van der Hoeven's Gaussian resampling, then Nussbaumer's synthetic transforms. One coordinate is a polynomial modulo `y^r + 1`, so roots of unity are signed shifts.
* **The obstruction.** Scanning `Θ(n)` bits at `Θ(log n)` levels already costs `n log n`. Both data movement and butterflies must be cheapened.
* **Two tape procedures.**
  * Address-field interchange in `O(V u^τ)`. Its intermediate steps store XOR combinations, but its final effect is an exact permutation.
  * Simultaneous butterfly layers `H_0(u,v) = ((u+v)/2, (u−v)/2)` in `O(V d^(λ'))`.
* **The source of the saving.** A fixed finite linear network on `W` wires exchanges two banks and restores arbitrary scratch values (the transparent/catalytic construction of Buhrman–Cleve–Koucký–Loff–Speelman).
  * Each wire carries an `m × m` "frame" matrix. The total rank of frame changes is `s = Σ_e rank(M_head − M_tail)`.
  * The construction achieves `s < Wm` with `h = 100`, three-element subsets of `[h]`, `v = C(100,3)` and `m = h^3`.
  * Over `F_2` the identity on triples is written as the incidence Gram matrix (`|S ∩ T| mod 2`, rank `<= h`) plus "side wires" that cancel the intersection-1 pairs.
  * The bank swap is three shears, `(x,y) -> (x, y+x) -> (−y, y+x) -> (−y, x)`.
  * A rank-`a` edge costs `a` smaller interchanges, via a Bruhat-type factorization `A = E_1 Π E_2`. Each unit of `Π` is one shear-swap-shear.
  * Recursion: parameter `mf` costs `s` calls at `f` on volume `V/W`, and `s/W < m` gives a power saving.
* **Ordering by CRT.** Coefficients are placed on prime cyclic axes by a triangular CRT map, `b_i = a_i + μ_i Σ_(j<i) a_j P_j (mod s_i)`, executed as prefix-controlled cyclic shifts.

**Dictionary to the Collatz compiler.**

| multiplication paper | Collatz (this repo) | type |
|---|---|---|
| bank swap = three shears (XOR swap) | reset switch = `A_1` turned into `R` by the run, then `R` erased (Lemma 1) | ANALOGY (both are "exchange by composition of unipotent moves") |
| restored scratch: the network's net effect is independent of the scratch contents | `U ∘ R = U`: Collatz erases appended `01` blocks without trace, while runs are transcoded into powers of 3 | ANALOGY |
| rank saving `s < Wm`, iterated, gives `(lg n)^(1−κ)` | equal-time classes: `o(A)` certificates for `A` Mersenne numbers (HYP-9213, `~A^0.35`); finite templates certify 0.156 at total 27 | ANALOGY: the same shape of saving, not the same mechanism |
| triangular CRT `b_i = a_i + μ_i Σ_(j<i) a_j P_j` as prefix-controlled shifts | endpoint progression `z ≡ f_u(−1) (mod 3^j)`: the 2-adic prefix `u` fixes the 3-adic offset; `n + 1 = 2^A (z − f_u(−1))/3^j` (THM-4555 (i)) | EQUATION-LOOKALIKE with a real common shape (triangular change between two axes) |
| roots of unity as signed shifts in `Z[y]/(y^r + 1)` | trailing-ones lag `D` as the rotation `2^(−D)` in `Z/(2^a − 1)`: `m_D + 1 ≡ 2^(−D)(n + 1)` | DICTIONARY (this is why "shift distance" is the right name) |
| Gaussian resampling between prime-length and power-of-two transforms | the 2 ↔ 3 mirror: Moran duality (THM-4554 (iv)), CRT independence (Proposition 5) | ANALOGY |

**What might transfer.**

* The paper's real lesson is structural. A power saving can come from a fixed finite gadget whose "rank" is strictly below its "volume", applied recursively. The gadget can be hard to find but cheap to verify.
* The Collatz counterpart would be a finite family of collision templates whose iteration certifies a family with a strictly smaller number of independent certificates at each scale.
* The Mersenne plateaus show such a saving empirically. The certified part (0.156 by templates of total `<= 27`) does not yet show it.
* A provable version would need the late coalescence of section 4, i.e. HYP-9214's per-visit merge probability. That is where the open problem sits.

### 7.2 The three attached preprints

* **Unrestricted pro-modularity at the prime two** (October 4, 2026).
  * *Claim:* every continuous, odd, absolutely irreducible 2-adic representation of `G_Q`, unramified outside finitely many primes, occurs in the full completed Hecke algebra `T_2(N) = lim_k Z_2 ⊗ T_(<=k)(N)` at some odd tame level. There is no de Rham hypothesis, and the claim is occurrence, not classical modularity.
  * *Mechanism:* fixed-determinant pseudodeformation families, local block equivalences for `GL_2(Q_2)`, and patching.
  * *Bridge:* DICTIONARY only. The completed algebra interpolates Hecke eigenvalues across weights 2-adically, through Kummer-type congruences.
  * The Mersenne family is interpolated the same way. The post-run word of `2^a − 1` to total `K` depends only on `a mod 2^(K−2)` (THM-4556 (iv)), because `a -> 3^a` is 2-adically continuous with `ord_(2^K)(3) = 2^(K−2)`.
  * So the exponent `a` plays the part of a 2-adic weight. Switch templates are locally constant on that weight space, and THM-4556 (v)'s certified classes are open sets in it.
  * Nothing in the proof transfers.
* **Weil classes and Hodge classes on abelian powers** (September 30, 2026).
  * *Claim:* the rational Hodge conjecture for every self-power of an abelian sixfold with an imaginary-quadratic action whose Hermitian form is hyperbolic of signature `(3,3)`, and for abelian varieties of dimension `<= 5` with such an action. It relies on a companion CM theorem.
  * *Bridge:* CONTEXT. The S17 field `Q(sqrt(−7))` acts on `E = 49a1` and on the Klein quartic's Jacobian. For example, `E^6` with `Q(sqrt(−7))` acting through one embedding on three factors and through the conjugate embedding on the other three has signature `(3,3)`. Being a CM power, it lies in the companion CM case.
  * No Collatz content.
* **Weak mixing of triangular billiards with an irrational angle** (October 5, 2026).
  * *Claim:* the billiard flow in every triangle with an angle irrational relative to `π` is weakly mixing.
  * *Mechanism:* a measurable eigenfunction `Xf = iλf` on the flat double is expanded in angular Fourier modes. Energy identities telescope in `λ` to `||Xf||^2 + ||Yf||^2 <= λ^2 ||f||^2`, which forces `Yf = 0`, and rigidity then excludes it.
  * *Bridge:* ANALOGY. "No eigenfunction" is the billiard version of "no hidden clock".
  * In the Collatz notes the only clocks found are the 2-adic clock of 3 and the 3-adic clock of 2 (THM-4554–4556), and Proposition 5 shows that they are independent. The 2-adic map alone is Bernoulli, so its weak mixing is trivial; the open coupling is with size.
  * The telescoping-in-the-spectral-parameter step resembles the character identity of the procgen thread's three-mirrors note (sums of multiplicative moments = increments of `ρ_n(±1)`). EQUATION-LOOKALIKE.

### 7.3 Survey of the repository

A scout pass covered all 722 abstracts and family descriptions. About 25 LaTeX sources were read in detail; nothing was executed. **No preprint addresses Collatz, `|2^A − 3^l|`, linear forms in logarithms, or an irrationality measure of `log_2 3`.** The bridges below are ranked; all results are the preprints' claims.

1. **The circulant Hadamard conjecture** (family 179; Lean claimed).
   * *Claim:* real circulant Hadamard matrices exist only in orders 1 and 4. Hence Barker sequences of length `> 1` exist exactly for lengths 2, 3, 4, 5, 7, 11, 13.
   * *Mechanism:* Turyn's descent `n = 4u^2`. The alternating product `Δ(x) = Π_S x_S^((−1)^|S|)` over subsets of the primes of `u` is an algebraic integer whose conjugates all have absolute value 1, hence a root of unity (Kronecker).
   * *Bridge (verified here):*
     * Barker 7 (`+++−−+−`) has minus set `{3,4,6} = {1,2,4} + 2`: the trivial cycle's real code, the `(7,3,1)` Fano set.
     * Barker 11's plus set is `QR_11 + 8`, the Paley `(11,5,2)` set.
     * Barker 13's minus set `{5,6,9,11}` is a Singer `(13,4,1)` set.
   * So the `{7, 21}` object is also the Barker sequence of length 7. The Paley bridge note's Proposition 3 (a single doubling orbit is a difference set only at length 3) is the Collatz-side finiteness statement. TECHNIQUE candidate: the alternating-product rigidity for Gauss periods of codes (S17 used Gauss sums).
2. **The 9/4 bound for ω** (family 107).
   * *Claim:* `ω <= 9/4`, from the shared-leg entropy inequality `λ(T) >= e^(p_X H(q)) Π λ(T_i)^(q_i)`.
   * *Mechanism:* the proof runs on the variational identity `max_q e^(H(q)) Π A_i^(q_i) = Σ A_i` (checked numerically here).
   * *Bridge:* EQUATION-LOOKALIKE, and genuinely the same convex duality, as in the refuel-bill Kraft sums and THM-4554 (iv)'s Chernoff generating function `(3/2) ρ(θ)^d`.
3. **Finite tensor savings and exact Fourier circuits** (family 130), and integer multiplication (109).
   * *Mechanism:* one finite saving, amplified by tensor powers, with CRT factorization and uncomputation of borrowed workspace. The existence of a finite saving is proved by duality: otherwise a contradictory "price" functional would exist.
   * *Bridge:* TECHNIQUE (shape). Section 8, step 4 asks for exactly this shape. A "finite saving gadget versus price functional" duality resembles the repo's Bellman brackets and Kraft LP duality.
4. **Prefix instructions and incompressible flows; incompressible box transport** (family 376).
   * *Mechanism:* Turing-universal Moore-type generalized shifts `x -> B^k x + c` on cylinder boxes. A determinant-one history recorder makes them measure-preserving, and they are realized by divergence-free Navier–Stokes forcing ("incompressible" means divergence-free, not Kolmogorov).
   * *Bridge:* CONTEXT, and precise. Collatz lies outside this universal radix-power class, because its multiplier 3 is coprime to the radix 2. Its determinant-one completion is the `×2×3` solenoid natural extension (the coalescence synthesis, axis 1).
5. **Artin primitive roots for every admissible base** (family 029).
   * *Claim:* at least `c_a x/(log x)^2` primes `p` in `(x, 2x)` have `a` as a primitive root, for `a = 2, 3` among others; `{2, 3}` simultaneously only conditionally.
   * *Bridges:* the repo's Artin coordinate (HYP-9174, and the Artin audit of the procgen thread); the clock tower's primes `2·3^k + 1` (mod18/19 note, D-δ); and the 2026-10-06 S17 note's Proposition 4 denominators `m = 7p` (2 primitive mod `p`, `p ≢ 1 (mod 3)`) that carry `τ_7` for shortcut codes.
   * Infinitely many such `m` would follow from the claim in the progression `p ≡ 2 (mod 3)`. That form was not checked: CONDITIONAL.
6. **The irrationality exponent of π is 2** (family 017; Lean claimed).
   * *Mechanism:* a Roth-style interpolation determinant on `Y = e^X`, with a "collision saving".
   * It is the collection's only irrationality-exponent method, and it does not transfer to `log_2 3` as written (transcendental Taylor coefficients).
   * Even `μ(log_2 3) = 2` would give only `|2^A − 3^l| >= 3^l l^(−1−ε)` ineffectively, which is useless for finite cycle exclusion. CONTEXT for the expense-records note (`collatz_expense_diophantine`, 2026-10-05).
7. **The entropy-rate dimension formula** (family 148).
   * *Claim:* `dim μ = min{1, h_RW/χ}` with exact overlaps allowed.
   * *Bridge:* CONTEXT, and a useful distinction. Collatz block maps never coincide as maps, so the repo's collisions are coincidences of one orbit point (−1), not overlaps of maps; and `h* = H(log_3 2)` is already Eggleston/Moran.
8. **Context only:**
   * Steinitz–Bergström prefix discrepancy `Θ(sqrt d)` (097): LRC, and running-excess balance.
   * Bernoulli convolutions (153): the case `λ = 2/3` stays open (base-3/2 dictionary).
   * Catalan's constant (005): Frobenius-at-primes finite certificates.
   * Generalized star height `<= 3` (134): episodes form a prefix code with affine actions, a formal cousin of the rewrite compiler; but the Collatz first-descent language is not regular.
   * Thorp shuffle (238): butterfly switch layers.
   * Thompson's group F (248), and the dyadic Erdős similarity case (084).
   * Deligne–Drinfeld (008): pentagon and hexagon relations, cf. the carry cocycle's coherence.
   * Mean-payoff games (104): Bellman brackets.

Direct consequences for a Collatz approach: none found.

## 8. Toward proofs: what moved, and the next steps

**Moved.**

* The odd shift distance is a theorem for every source (Theorem 2), not a census fact. It is explained by one identity: the run turns "append 1" into "append 01" (Lemma 1).
  * A rewrite compiler can therefore restrict trailing-ones searches to good lags and obtain each paired lag for free.
  * The debt of a switching reset-2 source passes to the smaller reset-2 source `m_(D+1)`.
* The `{7, 21}` question has an exact answer: the two append maps generate `F_21 = Aut(P_7)`, and only at 7 (Proposition 3).
* The 37-of-60 census splits into two parts:
  * a certified 2-adic part, the template classes of THM-4556 (v), which give proved infinite families (`a ≡ 95 (mod 128)` is the first);
  * a late-coalescence part, which is the bulk (median template total 247) and is open: HYP-9213 and HYP-9214.
* Proposition 5 is a barrier statement in the repo's sense. At every finite depth the forward 2-adic data and the backward 3-adic data of the exponent are independent. A proof cannot come from combining the two clocks alone; it must use size (THM-4554, remarks).

**Not moved.** Collatz; HYP-9213 (`o(A)` classes); HYP-9214.

**Next steps (ranked).**

1. **A local merge lemma for the Mersenne line.** Two orbits that are within `c` bits of each other at equal times merge within `O(1)` steps with probability bounded below, in the 2-adic measure on the current digits. This is the per-visit probability that HYP-9214's model leaves open. With it, the coalescing-walk heuristic of section 4 would become a conditional theorem.
2. **Push the certified density past `K = 27`.** The numba kernel doubles its cost per level (7 s at 27). Reaching `K ≈ 40` needs a class recursion that refines only unresolved classes. This would show whether the certified share tends to the empirical 0.95–0.98 or saturates, i.e. whether `μ_2(S) = 1`.
3. **Wire Theorem 2 into the compiler** (`checked_switch_phase19`): good-lag search, paired lags, and debt hand-off to `m_(D+1)`. Then measure the residual seeds.
4. **A certificate-rank inequality in the multiplication paper's shape.** For the Mersenne family, `N(2A)/N(A) ≈ 1.25 < 2` empirically. A provable "strictly fewer independent certificates per doubling" for a template-defined subfamily would be the honest analogue of `s < Wm`.

## 9. Hostile checks and numerology

* **The parity law does not depend on 7.** Lemma 1 and Theorem 2 hold for every odd source. Proposition 3 only says that the two maps of Lemma 1 happen to generate the full Paley symmetry at `p = 7`.
  * The relation is "the same two maps", not "7 causes odd `D`".
  * `21 = |F_21| = R^3(0)` is the `k = 3` solution of `2^k + 1 = 3k` (EXACT, and coincidental: two solutions).
* **The shift dictionary cannot see the parity.** In `Z/(2^a − 1)` with `a` odd, the lags `D` and `D + a` give the same rotation `2^(−D)`. So the parity law is not a Mersenne-ring statement. It is the 2-adic congruence `3 ≡ −1 (mod 4)` in Lemma 1.
* **Exponents 7 and 21 are both class roots** (non-switching). So are all odd `a <= 27` except 13: a short-orbit effect. NUMEROLOGY.
* **"37 of 60" depends on the window.** At larger `a` the share is 0.95–0.98, with orbits stopped at 1. Every least-lag merge for `a <= 2000` is at a value `>= 911`, so no count rests on simultaneous arrival at 1.
* **The `a ≡ 95 (mod 128)` family:** `128 = 2^(K−2)` at `K = 9`, and `125 = 2^7 − 3`. The 7s are the clock's modulus, not Fano's. NUMEROLOGY guard.
* **The multiplication paper's triples are not designs:** all three-element subsets of `[100]`, with intersections 0–3. There is no Fano plane in it. Any `{7, 21}` reading of the gadget is EQUATION-LOOKALIKE at best.
* **Proposition 5 is an exact density statement at finite depth only.** The infinite-depth switching set is open in `Z_2`, and its density need not exist.

## 10. Audit

Pending.

## 11. Reproduction and references

**Reproduction.**

* `python3 04-computation/experiments/mersenne_switch_parity_f21_compression_20261006.py` runs in a few minutes and ends `ALL CHECKS PASSED`; its output is in [`.out`](mersenne_switch_parity_f21_compression_20261006.out).
* `python3 04-computation/experiments/mersenne_certified_density_numba_20261006.py 9 27` takes under 20 s (needs numba); its output is in [`.out`](mersenne_certified_density_numba_20261006.out).

**Repo inputs.**

* `checked_switch_phase19_20261004.md`, especially (3)–(5);
* [THM-4553](../../01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md), [THM-4554](../../01-canon/theorems/THM-4554-backward-sieve-keeps-a-positive-proportion-moran-duality.md), [THM-4555](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md), [THM-4556](../../01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md), [THM-4557](../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md);
* HYP-9213 and HYP-9214;
* `mod18_mod19_seven_sixtythree_fractal_20261006.md`;
* `collatz_paley_bridge_20261001.md`;
* `ramanujan_heegner7_trivial_cycle_20261006.md` (S17).

**Literature.**

* M. M. Ahmed, *The intricate labyrinth of Collatz sequences*, arXiv:1602.01617 (2016), Theorem 2.1: the root case of the reset switch.
* A. Schönhage and V. Strassen (1971); D. Harvey and J. van der Hoeven, *Integer multiplication in time O(n log n)*, Ann. of Math. (2021); H. Nussbaumer (1980): for the multiplication paper's lineage.
* H. Buhrman, R. Cleve, M. Koucký, B. Loff and F. Speelman, *Computing with a full memory: catalytic space* (STOC 2014): the restoration construction the multiplication paper cites.
* R. H. Barker (1953) and R. J. Turyn: Barker sequences and Turyn's descent.

**The openai/math preprints** (https://github.com/openai/math): unrefereed, model-written, partly formalized. They are cited by directory under `preprints/`:

* *Integer-multiplication-below-n-log-n-September-23-2026* (LaTeX source read);
* *The-circulant-Hadamard-conjecture-September-23-2026*;
* *Matrix-Multiplication-Nine-Fourths-October-2-2026*;
* *Finite-tensor-savings-and-exact-Fourier-circuits-September-25-2026*;
* *Prefix-Instructions-and-Incompressible-Flows-September-27-2026* and *Incompressible-Box-Transport-and-Finite-Computation-September-27-2026*;
* *Primitive-roots-for-every-admissible-integer-base-October-4-2026*;
* *The-irrationality-exponent-of-pi-is-2-September-24-2026*;
* *The-entropy-rate-dimension-formula-for-self-similar-measures-on-the-line-September-24-2026*;
* the three attached preprints (Weil classes on abelian powers; pro-modularity at 2; weak mixing of triangular billiards).
