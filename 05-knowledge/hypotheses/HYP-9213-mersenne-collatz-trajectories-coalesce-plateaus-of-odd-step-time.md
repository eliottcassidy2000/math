---
id: HYP-9213
title: "Mersenne coalescence: the odd-step stopping time sigma(2^a - 1) takes o(A) distinct values for a <= A (empirically about A^0.37: 23, 37, 58, 104 values for A = 100, 400, 1200, 6000), so almost every reset pair {2^(2k-1) - 1, 2^(2k) - 1} shares sigma with a smaller Mersenne number (for a <= 6000 every such coincidence is an equal-time merge)"
status: >
  RESOLVED 2026-10-07: PROVED (THM-4581 (6d)). The lag-1 Mersenne switch is the pair chain from (1, 2) on (p, q_1), absorbed
  almost surely, so mu_2(S) = 1, and S18 Proposition 6 (with THM-4556 (ii), (iv)) gives o(A) distinct values. The finer count
  ~ A^0.37 (exponent 1 - alpha) remains NUMERICAL/HEURISTIC; THM-4581 (6f) gives alpha >= 1/2 if alpha exists. Earlier:
  FINITE-EXACT census to a = 6000 (audited, MISTAKE-572); certified lower densities 0.1199 (THM-4556 (v)) and 0.3853 (THM-4569).
source: mac-mini-2026-10-06-mod1819, 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md, section 1
related:
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
  - 01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md
  - OEIS A390816, A193688; Ohira-Watanabe arXiv:1104.2804 (growth ~ 13.45 a for all steps)
---

# HYP-9213 — the Collatz trajectories of 2^a − 1 coalesce

**Conjecture.** Let `σ(n)` be the number of odd steps from `n` to 1. Then `#{σ(2^a − 1) : a <= A} = o(A)`. Equivalently, the reset pairs that start a new level set have density 0.

**Evidence (FINITE-EXACT, `a <= 6000`).**
* The count of distinct values is 23, 31, 37, 50, 63, 80, 104 for `A = 100, 200, 400, 800, 1600, 3200, 6000`.
* The maximal runs of constant `σ` number 426. The longest has length 282 (`a = 3243..3524`). A least-squares fit of the counts gives `A^0.362`.
* Near `a = 6000` two values, 27728 and 28762, alternate in runs.

**Mechanism (PROVED pieces, THM-4556).**
* By (ii), consecutive pairs share `σ` exactly when the debt state `R_a = 3·2^v·oddpart(R_(a−1)) + 1` resolves with lag 1 in stopping time.
* By (iv)–(v), each collision-witnessed merge of total `T` holds on the whole class `a mod 2^(T−2)`. Non-resolution is never certified.
* Certified density is `>= 0.1199` among odd `a`. Empirically, on `[10^3, 6·10^3]`, 0.980 of odd `a` have a partner and 0.882 merge with `a − 1`.
* For `a <= 6000`, every `σ`-coincidence is an equal-time merge: the 102 pairs that start a new level are exactly the 102 non-merging pairs.

**Remark.** Conjecture HYP-9214 (debt resolution tends to 1 for generic large sources) would make this plausible but does not imply it, because the Mersenne family is thin.

---

## Update (2026-10-07, mac-mini-2026-10-07-oaimath3): PROVED (THM-4581)

* **The switch is a chain absorption.** S19's lag-1 switch `U^i(2X − 1) = U^(i+1)(2X/3 − 1)` (with `X = 3^(a−1)`) is absorption at `(0, 0)` of the Terras-clock pair chain from `(1, 2)`.
* **Theorem H (THM-4581 (3)).** That chain is absorbed almost surely. `X` is Haar on `1 + 8Z_2` as `a` runs over the odd 2-adic integers, and the four forced parities of `q_1` do not matter.
* **So the switching set has full measure:** `μ_2(S) = 1`. S18's Proposition 6 then gives the conjecture: the odd `a` with `σ(M_a) ≠ σ(M_(a−1))` have density 0, and the number of `σ`-levels is `o(A)`.
* **Still open.** The empirical count `A^0.37` corresponds to an any-lag exponent `α ≈ 0.63`. THM-4581 gives only `α ≥ 1/2`.
