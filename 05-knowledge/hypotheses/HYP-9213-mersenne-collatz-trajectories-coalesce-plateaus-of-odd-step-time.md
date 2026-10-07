---
id: HYP-9213
title: "Mersenne coalescence: the odd-step stopping time sigma(2^a - 1) takes o(A) distinct values for a <= A (empirically about A^0.37: 23, 37, 58, 104 values for A = 100, 400, 1200, 6000), so almost every reset pair {2^(2k-1) - 1, 2^(2k) - 1} shares sigma with a smaller Mersenne number (for a <= 6000 every such coincidence is an equal-time merge)"
status: >
  OPEN. FINITE-EXACT census to a = 6000 (independently audited 2026-10-06; corrections in MISTAKE-572); certified lower density
  15719/131072 of switching odd exponents (THM-4556 (v)); the pairing and odd least shift are PROVED (THM-4556 (ii)-(iii)). No proof of o(A).
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
