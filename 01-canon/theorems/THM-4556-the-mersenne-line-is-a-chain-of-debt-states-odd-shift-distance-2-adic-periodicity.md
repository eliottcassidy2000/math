---
id: THM-4556
title: "The Mersenne line is a chain of debt states: U^a(2^a - 1) = oddpart((3^a - 1)/2) (binary repunit -> ternary repunit); the reset pairing glues {2^(2k-1) - 1, 2^(2k) - 1}; for every odd exponent with an equal-time Mersenne partner the least shift D is odd and the nearest partner is 2^(2k) - 1 = 3 (4^k - 1)/3 (the partner set is a union of reset pairs; the 63 = 3*21 pattern); switches are periodic in the exponent on the 2-adic clock of 3 (period 2^(K-2)); a finite class search certifies lower density >= 15719/131072 > 0.1199 of switching exponents among odd a"
status: "PROVED (i)-(iv); FINITE-EXACT (v), (vi); INDEPENDENTLY AUDITED (2026-10-06; corrections in MISTAKE-572; record in section 9 of the results note seven_twentyone_mersenne_openai_math_20261006.md)"
session: mac-mini-2026-10-06-mod1819
source: 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md
scripts:
  - 04-computation/experiments/sevens_20261006_mersenne_plateaus.py (+ .out, ALL CHECKS PASSED)
related:
  - 01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md (collisions, the switch m = (n+1)/2^D - 1)
  - 01-canon/theorems/THM-4554-backward-sieve-keeps-a-positive-proportion-moran-duality.md ((v): the 3-adic clock of Mersenne backward-minimality)
  - 05-knowledge/results/checked_switch_phase19_20261004.md (the reset switch (3) and the debt states (4)-(5))
  - OEIS A390816 (odd steps of 2^n - 1), A193688 (all steps of 2^n - 1; a(2n) = 1 + a(2n-1), Herrera 2023)
---

# THM-4556 — the Mersenne line is a chain of debt states

**Setting.**
* `U(x) = oddpart(3x+1)`, and `M_a = 2^a − 1`.
* `R_a = (3^a − 1)/2 = 111…1` in base 3 (a ternary repunit), so `R_0 = 0` and `R_a = 3R_(a−1) + 1`. This is the orbit of 0 under the undivided map `x -> 3x + 1`.
* `σ(n)` is the number of `U`-steps from `n` to 1. Orbits are stopped at 1 when counting meetings.
* Two odd numbers *merge at equal time* if `U^j(x) = U^j(y) ≠ 1` for some `j`.
* `a_k = (4^k − 1)/3` (`1, 5, 21, 85, …`) are the direct predecessors of 1.

**(i) Repunit transit (PROVED).** For `a >= 2`, `U^a(M_a) = oddpart(3^a − 1) = oddpart(R_a)`.

**(ii) The Mersenne line is a chain of debt states (PROVED).**
* For even `a = 2k >= 4`, `R_(a−1)` is odd and `U(R_(a−1)) = oddpart(R_a)`. Hence `U^a(M_a) = U^a(M_(a−1))`, so `M_(2k)` and `M_(2k−1)` merge at equal time and `σ(M_(2k)) = σ(M_(2k−1))`. This is the reset switch of `checked_switch_phase19` (3), going back to Ahmed 2016, Theorem 2.1. In standard steps it is A193688's `a(2n) = 1 + a(2n−1)` for `n >= 2` (A. Herrera, 2023).
* For odd `a >= 3`, `R_a = 3·2^v·oddpart(R_(a−1)) + 1` with `v = 1 + v_2(a−1)`. This is a debt state `(e, h, M) = (1, v, oddpart(R_(a−1)))` of `checked_switch_phase19` (4).
* Hence, for odd `a`, `σ(M_a) = σ(M_(a−1))` iff `σ(R_a) = σ(oddpart(R_(a−1))) − 1` (lag 1 in stopping time).
* If the debt state's orbit meets that of `oddpart(R_(a−1))` one step later at a value `≠ 1`, then `M_a` and `M_(a−1)` merge at equal time and the pair `{M_(a−1), M_a}` continues a plateau of `σ`.
* The converse (equal `σ` implies a merge before 1) is not proved. It holds for all odd `a <= 6000` (FINITE-EXACT, 2575 cases), where the last value before 1 is always 5, 85 or 341.

**(iii) Odd least shift (PROVED).** Let `a` be odd. If `M_a` merges at equal time with `M_b` for some odd `b < a`, it also merges with `M_(b+1)`.
* Hence the least shift `D = a − b` over all equal-time Mersenne partners is odd.
* The nearest partner `M_b` (least shift, i.e. largest `b`; `b` even) is a multiple of 3, namely `M_b = 3·a_(b/2)`. It is a leaf of the backward tree, and its cofactor `a_(b/2)` maps to 1 in one step.
* Every even partner `b >= 4` also brings `b − 1`, so the partner set is a union of reset pairs and its smallest member is odd. Other switches need not land on a multiple of 3: `M_31` merges with `M_30` and `M_29`, and later with `M_26` and `M_25`.

**(iv) 2-adic periodicity (PROVED).** Let `K >= 3` and `a >= 2`, and continue words past 1 with `U(1) = 1` and exponent 2.
* The 2-adic word of `M_a` after its run of `a − 1` ones, truncated at total `K`, depends only on `a mod 2^(K−2)`.
* Consequently, suppose the reduced words of `M_(a0)` and `M_(a0−D)` have colliding prefixes `u ~ u'` of total `K` with `|u'| = |u| + D`.
* Then for every `a ≡ a0 (mod 2^(K−2))` with `a > D` and `3^(a−1+|u|) > 2^K`, `M_a` and `M_(a−D)` merge at equal time by step `(a − 1) + |u|`. The merge happens exactly then if `u` is the shortest colliding prefix.

**(v) Certified density (FINITE-EXACT).** A complete search over the odd classes `c mod 2^(K−2)` looked for a partner class `c − D` (`D` odd, `<= 31`) whose word prefixes collide (equal total, length difference `D`, equal numerator).
* The certified share of odd exponents is `0.0156, 0.0312, 0.0469, 0.0605, 0.0713, 0.0801, 0.0879, 0.0950, 0.1016, 0.1080, 0.1141, 0.1199` for `K = 9, …, 20`.
* By (iv), every odd `a > 31` in a certified class has a uniform Mersenne switch: `a − 1 >= 31 >= D`, and the merge value exceeds `3^31/2^20 > 1`. So the set of such exponents has lower density `>= 15719/131072 = 0.11992…` among odd exponents.
* Allowing every `D <= 31`, even ones included, certifies exactly the same classes, and the least `D` is always odd.
* The first certified classes are `95 mod 2^7` (`K = 9`); `159` and `207 mod 2^8` (`K = 10`); and `31, 335, 367, 455 mod 2^9` (`K = 11`). All have `D = 1`.
* The independent audit reproduced the exact fractions `1/64, 4/128, …, 15719/131072`. All 366 odd `a <= 6000` in certified classes do merge with `M_(a−D)`.

**(vi) Census (FINITE-EXACT).**
* `σ(M_a)` takes 23, 37, 58 and 104 distinct values for `a <= 100, 400, 1200, 6000`.
* Every level set `{a >= 3 : σ(M_a) = s}` is a union of reset pairs (PROVED, by (ii)). The largest level set for `a <= 1200` has 118 members (`σ = 4316`).
* Odd `a` with an equal-time partner `M_b`, `b < a`: 37 of the 60 in `[3, 121]` and 164 of the 199 in `[3, 400]`.
* Of the 2500 odd `a` in `[10^3, 6·10^3]`, 2450 (0.980) share `σ` with a smaller exponent, and exactly these have an equal-time partner. 2204 (0.882) share it with `a − 1`.
* The audit reproduced all 6000 terms of the OEIS A193688 b-file from these `σ` values plus the halving counts.

## Proofs

(i) For the Terras map `T(n) = (3n+1)/2` (`n` odd), `T^j(2^a − 1) = 3^j 2^(a−j) − 1` for `j <= a`.
* So `T^a(M_a) = 3^a − 1`, reached after `a` odd steps and no extra halvings.
* The `a`-th odd step of `U` then strips the remaining factors of 2, giving `U^a(M_a) = oddpart(3^a − 1)`. ∎

(ii) `R_(a−1)` is a sum of `a − 1` odd terms, so it is odd for even `a`. Then `3R_(a−1) + 1 = R_a`, and `U(R_(a−1)) = oddpart(R_a)`.
* By (i), `U^(a−1)(M_(a−1)) = oddpart(R_(a−1)) = R_(a−1)`, so `U^a(M_(a−1)) = oddpart(R_a) = U^a(M_a)`. This value is not 1 for `a >= 4`.
* For odd `a`, `v_2(R_(a−1)) = v_2(3^(a−1) − 1) − 1 = 1 + v_2(a−1)` by lifting the exponent. Hence `R_a = 3R_(a−1) + 1 = 3·2^v·oddpart(R_(a−1)) + 1`.
* Write `t = σ(oddpart(R_(a−1)))`. Then `σ(M_(a−1)) = (a − 1) + t` and `σ(M_a) = a + σ(R_a)`. These agree iff `σ(R_a) = t − 1`.
* An equal-time merge forces equal stopping times. The converse needs equal last values before 1, which is checked, not proved. ∎

(iii)
* By (ii), `M_(b+1)` and `M_b` merge at time `b + 1 <= a − 1`, because `b + 1` is even and `>= 4`. (For `b = 1`, `M_1 = 1` has no merges.)
* The merge time `j` of `M_a` with `M_b` satisfies `j >= a`. For `j <= a − 1`, `U^j(M_a) + 1 = (3/2)^j 2^a`. Meanwhile `U(x) + 1 <= (3/2)(x + 1)` gives `U^j(M_b) + 1 <= (3/2)^j 2^b`, which is smaller.
* At time `j` both orbits of `M_b` and `M_(b+1)` sit at the same value, so `M_a` meets `M_(b+1)` at time `j` too.
* Hence any partner with odd `b` comes with a partner `b + 1` at shift `D − 1`, and the least `D` is odd. Then `b` is even, and `2^b − 1 = (2^(b/2) − 1)(2^(b/2) + 1) = 3·a_(b/2)`. ∎

(iv)
* `U^(a−1)(M_a) = 2·3^(a−1) − 1` (Terras form of (i)).
* The `U`-word of an odd `x`, truncated at total `K`, depends only on `x mod 2^(K+1)`.
* `2·3^(a−1) − 1 mod 2^(K+1)` depends only on `3^(a−1) mod 2^K`, and `3` has order `2^(K−2)` modulo `2^K` for `K >= 3`.
* The same holds for `M_(a−D)` with `a − D` in place of `a`.
* A collision is an identity of words, so THM-4555 (iv) applies to every such `a` with run `a − 1 >= D`, and the meeting time is `(a − 1) + |u|`.
* The common value is at least `3^(a−1+|u|)/2^K > 1`, so the meeting is before 1. An earlier merge would have the same total difference `D`, and so would itself be a shorter colliding prefix. ∎

(v), (vi) Script `sevens_20261006_mersenne_plateaus.py` (and the class search to `K = 20`, recorded in the results note). ∎

## Remarks

* **Seven and twenty-one.**
  * `7 = 111_2` and `21 = 10101_2` approximate two distinguished 2-adic points (`7 ≡ −1 mod 8`, `21 ≡ −1/3 mod 64`). One is `−1`, the repelling fixed point of `(3x+1)/2` and the limit of `M_a`. The other is `−1/3`, where `v_2(3x+1)` is infinite, and the limit of `a_k`.
  * `3·(−1/3) = −1` becomes `2^(2k) − 1 = 3a_k` at finite level, with `63 = 3·21`.
  * (iii) says the least-shift switch from an odd exponent lands on such a number. Other switches need not (`M_31` also merges with `M_29`).
  * ANALOGY: `F_21 = {x ↦ 2^i x + b mod 7}` is doubling plus translation modulo `M_3`. Nothing in the theorem uses it.
* **Mirror clocks.** Forward switches of `M_a` run on the 2-adic clock of 3 (`a mod 2^(K−2)`, (iv)). Backward-minimality of `M_a` to depth `k` runs on the 3-adic clock of 2 (`a mod 2·3^(k−1)`, THM-4554 (v)).
* **Concurrent extension.** The opus-2026-10-06-S18 note `05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md` (Theorem 2, parity law) proves that trailing-ones lags come in pairs `{D, D+1}` for every source, so the least lag is odd for every reset-2 source. (iii) is its Mersenne case `t = 1`. The same note certifies 0.1556 at template total 27, extending (v).
* **Scope.** These are statements about the Mersenne family and its rewrite rules. They transport certificates; they do not prove that any `M_a` reaches 1. Collatz is OPEN.

---

## Update (2026-10-07, mac-mini-2026-10-07-oaimath3): the switch density rises to 0.385

(v)'s certified density `15719/131072 = 0.1199` (exhaustive, `K = 20`; S19: 0.1740 at `K = 31`) is superseded by a computer-assisted bound.
* In the Terras clock, the lag-1 Mersenne pair merges with probability at least **0.3853**: float64 value iteration from below on the box `B(9,100)`, reproduced by an independent implementation on smaller boxes; THM-4569 (5).
* Each merge is decided by `a` modulo a power of 2, so the lower density of odd `a` with `σ(2^a − 1) = σ(2^(a−1) − 1)` is at least 0.3853.

## Update 2 (2026-10-07, same session): the density is 1 (THM-4581)

* The lag-1 Mersenne pair chain is absorbed almost surely (THM-4581 (3)). So the odd `a` with `σ(2^a − 1) = σ(2^(a−1) − 1)` have natural density 1, and the switching set has `μ_2(S) = 1`.
* With (ii) and S18's Proposition 6, this proves HYP-9213: there are `o(A)` distinct `σ`-levels.
* The density of odd `a` without a switch of template total `≤ K` is between `c K^(−1/2)` and `C K^(−1/2) (log K)^2` (THM-4581 (4)). The certificates (v) are finite-`K` instances of this.
