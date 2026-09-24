---
id: THM-4471
title: "The fixed-point 'proof of the Collatz conjecture' in arXiv:2502.20642 is invalid: its Theorems 2.1(5)-2.3(5) are false, its coefficients equally 'prove' 3n-1 reaches 1, and no contraction argument in |x-y| can prove Collatz"
status: >
  PROVED + INDEPENDENTLY AUDITED (a refutation; Collatz remains OPEN).
  T. Kawasaki, "A proof of the Collatz conjecture", arXiv:2502.20642v1
  [math.GM] (28 Feb 2025), defines weighted generalized pseudocontractions
  (WGP). It claims (Theorem 2.3(5)) that a WGP whose coefficients satisfy,
  at each pair, one of two alternatives has a fixed point attracting every
  orbit, and it applies this to the Collatz shortcut map on (N, |x-y|).
  (A) The theorem is false. The successor map x -> x+1 on (N, |x-y|)
  satisfies every hypothesis with the paper's own constants
  (lambda = 0, A = 1/2, B = 2, M = 2) and has no fixed point. The least
  counterexample has 3 points (a 3-cycle).
  (B) The proof of Theorem 2.1(5) needs alternative 1 at (p, Tp) or
  alternative 2 at the SWAPPED pair (Tp, p). Condition (5) only offers an
  alternative at each pair separately. For the paper's Collatz coefficients
  the uncovered case occurs at every odd p >= 3. The claimed decay
  d(T^n x, T^(n+1) x)^2 <= 2^-n d(x,Tx)^2 fails at x = 5, n = 1.
  (C) SHEET control: the paper's coefficient table, copied verbatim,
  satisfies every hypothesis for the 3n-1 map, whose orbit
  5 -> 7 -> 10 -> 5 never reaches 1.
  (D) More generally, no contraction-type argument in |x-y| can prove
  Collatz. Every odd step p >= 3 stretches the consecutive distance by at
  least 4/3, and up to 3/2, and Mersenne starts 2^m - 1 defeat every
  uniform orbit-radius bound. The surviving metric principles are
  reformulations: Caristi, and a uniformly discrete contraction metric,
  are each equivalent to Collatz; Bessaga's contraction metric is
  equivalent to its cycle half.
  The paper's Lemma 2.1 (the owner's triangle-inequality sandwich) is
  correct: it is the triangle inequality followed by AM-QM, and it
  discards the cross term +-2ab, which is where the sign law lives.
source: collatz-procgen-20260922 session, fixed-points lane (2026-09-24), prompted by the owner's pointer to arXiv 2502.20642 and its Lemma 2.1; gap and 3-point counterexample found by the orchestrator, successor counterexample and SHEET control by the lane; audited and promoted by the session orchestrator 2026-09-24
depends_on: []
related:
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (transport theorem; any proof must be side-aware)
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
  - 01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md
script: 04-computation/experiments/collatz_procgen_20260924_fixedpt_kawasaki.py
script_audit: 04-computation/experiments/collatz_procgen_20260924_fixedpt_orchestrator_check.py
output: 05-knowledge/results/collatz_procgen_20260924_fixedpt.out
output_audit: 05-knowledge/results/collatz_procgen_20260924_fixedpt_orchestrator_check.out
script_sha256: 2cc3fdd47e8dac126f4aad22a46d7a6076bb17e3d6f571122ca92fdf974f8437
script_audit_sha256: 6be831d6fab0bceddf7b0159ffe1cffd96b78cdcb4bb6faa343586ec865fbec0
output_sha256: b8095963e2b4b63e6c6139333dddf6e65db7d6278815c8b7f69c2c39e34f66d7
output_audit_sha256: 727315e5c97959186ed36200e1e6c17fbc9be6a2afaa4ef20b80090bb7ebac82
hash_basis: raw LF bytes
audit: >
  Two independent transcriptions of the paper's Section 3 coefficient
  table (the orchestrator's and the lane's) agree. The orchestrator's
  code (written before and independently of the lane) establishes:
  * Theorem 3.1's inequality and condition (5) hold at every pair < 400;
  * the proof-needed alternative fails exactly at the odd p >= 3 below 5000;
  * the x = 5 failure;
  * the 3-point counterexample (M = 4);
  * the successor counterexample with the paper's constants;
  * the verbatim-table 3n-1 control.
  The lane's pipeline was re-run, and its output is byte-identical (peak
  memory 481 MB).
---

# THM-4471 -- arXiv:2502.20642's fixed-point proof of Collatz is invalid

**PROVED + INDEPENDENTLY AUDITED.** Full audit: [collatz_procgen_20260924_fixed_points_kawasaki_audit](../../05-knowledge/results/collatz_procgen_20260924_fixed_points_kawasaki_audit.md).

## 1. The paper's framework

* **WGP.** A map `T` is a *WGP* with coefficient functions `alpha, ..., zeta : X x X -> R` if

  `alpha d(Tx,Ty)^2 + beta d(x,Ty)^2 + gamma d(Tx,y)^2 + delta d(x,y)^2 + eps d(x,Tx)^2 + zeta d(y,Ty)^2 <= 0`.
* **Condition (5).** Write `s1 = alpha + zeta + 2min(beta,0)`, `r1 = delta + eps + 2min(beta,0)`, `s2 = alpha + eps + 2min(gamma,0)` and `r2 = delta + zeta + 2min(gamma,0)`. Condition (5) of Theorems 2.1–2.3 asks, **for each pair** `(x,y)`, for
  * **alt1:** `s1 > 0` and `-r1 <= A s1`; or
  * **alt2:** `s2 > 0` and `-r2 <= A s2`.

  Theorem 2.3 adds `alpha + beta + zeta >= B` to alt1, `alpha + gamma + eps >= B` to alt2, and `|coefficients| <= M`.
* **Section 3.** The paper takes `T(1) := 1`, `x/2` (x even) and `(3x+1)/2` (x odd, `x >= 3`) on `(N, |x-y|)`, gives a piecewise table with entries in `{-2, ..., 2}`, and concludes Collatz with `lambda = 0, A = 1/2, B = 2, M = 2`.

## 2. The theorem

**(A) Theorem 2.3(5) is false.**
* Let `S(x) = x + 1` on `(N, |x-y|)`.
* Use the coefficients `(1, 1, -2, 0, 0, 0)` when `x >= y` and `(1, -2, 1, 0, 0, 0)` when `x < y`.
* With `t = |x - y|` the WGP left side is `-6t - 1 < 0` in both cases.
* alt1 holds when `x >= y`, and alt2 when `x < y`: `s = 1`, `r = 0`, and the `B`-sum is `2`.
* Every hypothesis holds with the paper's constants, yet `S` has no fixed point.
* The least counterexample has 3 points: on `{0,1,2}` the 3-cycle works with `M = 4`. Any 2-cycle violates (5).

**(B) The swapped quantifier.**
* Substituting `(x,y) = (p, Tp)` into the `beta`-eliminated inequality gives `s1(p,Tp) d(Tp,T^2 p)^2 + r1(p,Tp) d(p,Tp)^2 <= 0`.
* Substituting `(x,y) = (Tp, p)` into the `gamma`-eliminated one gives `s2(Tp,p) d(Tp,T^2 p)^2 + r2(Tp,p) d(p,Tp)^2 <= 0`.
* So the contraction step needs **alt1(p,Tp) or alt2(Tp,p)**. Condition (5) supplies "alt1 or alt2" at `(p,Tp)` and, separately, at `(Tp,p)`.
* For the paper's table the uncovered case, alt2 only at `(p,Tp)` and alt1 only at `(Tp,p)`, occurs at every odd `p >= 3`: 2499 of the steps with `p < 5000`.
* The decay claim is false at `x = 5`: `d(8,4)^2 = 16 > (1/2) * 9`.
* The paper's Collatz-specific algebra (its Theorem 3.1) is correct, up to one harmless slip. The error is entirely in its Section 2.

**(C) SHEET control.** The paper's table, copied verbatim, satisfies every hypothesis for `T_-(x) = x/2`, `(3x-1)/2`, which fixes 1 (checked for all pairs `< 400`). So the argument would equally "prove" that every positive integer reaches 1 under `3n-1`, but `5 -> 7 -> 10 -> 5`. By the transport theorem (mod-192 note, Theorem 6) any valid proof must use the sign. The paper's inequalities, which are even in `x -> -x` apart from the special point 1, cannot tell the sheets apart.

**(D) No contraction in `|x-y|`.**
* `sup_p d(Tp,T^2 p)/d(p,Tp) = 3/2`, attained at every `p = 3 (mod 4)`, and every odd `p >= 3` stretches by at least `4/3`.
* For `x = 2^m - 1`, `T^j x = 3^j 2^(m-j) - 1` for `j <= m`.
* Hence every hypothesis that forces one-step decay fails: Banach, Kannan, Chatterjea, Reich, and the corrected Theorem 2.1. So does every uniform orbit-radius bound.
* **The surviving principles are reformulations.**
  * Caristi (`d(x,Tx) <= phi(x) - phi(Tx)`, `phi >= 0`) is equivalent to Collatz, with `phi` the orbit's total variation.
  * A uniformly discrete complete contraction metric exists iff Collatz holds.
  * Some complete contraction metric exists iff there is no nontrivial cycle (Bessaga's converse, CITED), which is the cycle half only.

## 3. The owner's sandwich

The paper's Lemma 2.1 comes from the owner's sandwich

`|d(x,z) - d(z,y)| <= d(x,y) <= d(x,z) + d(z,y)`.

* **Lemma 2.1 is correct.** It is the triangle inequality followed by AM–QM, `(a+b)^2 <= 2(a^2+b^2)`, and it discards the cross term `+-2ab`.
* **The sign law is the two equality cases.** For a parity word `w`, `2^p |T^p x| = 3^a|x| + c_w` (upper equality) for `x > 0`, and `|3^a|x| - c_w|` (lower equality) for `x < 0`, with the middle point `z = 0`.
* **The paper's argument breaks where the cross term is `+2ab`:** the up-steps, where positive orbits realise the upper equality and consecutive distances grow.

## 4. What the correct fixed-point theorem gives

* The inverse branches `D(x) = 2x` and `E(x) = (2x-1)/3` are Banach contractions of `Z_2` (factor 1/2).
* Every parity word has exactly one fixed point, the cycle gate `c_w/(2^p - 3^a)`, so there are `2^p` rational periodic points per period.
* For `p <= 24` exactly 18 of them are integers: the five known cycles.
* Collatz's cycle half is the statement that this Banach fixed-point chain meets `Z_(>0)` only in `{1,2}`. Banach provides existence for free, and the whole problem is integrality at the gate.
