---
id: HYP-9213
title: "Mersenne coalescence: the odd-step stopping time sigma(2^a - 1) takes o(A) distinct values for a <= A (empirically about A^0.37: 23, 37, 58, 104 values for A = 100, 400, 1200, 6000), so almost every reset pair {2^(2k-1) - 1, 2^(2k) - 1} merges at equal time with a smaller Mersenne number"
status: >
  OPEN. FINITE-EXACT census to a = 6000; certified lower density 0.1199 of switching odd exponents (THM-4556 (v));
  the pairing and odd least shift are PROVED (THM-4556 (ii)-(iii)). No proof of o(A).
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
* The maximal runs of constant `σ` number 426, the longest of length 282.
* Near `a = 6000` two values, 27728 and 28762, alternate in runs.

**Mechanism (PROVED pieces, THM-4556).**
* By (ii), consecutive pairs continue a plateau exactly when the debt state `R_a = 3·2^v·oddpart(R_(a−1)) + 1` resolves with lag 1.
* By (iv)–(v), resolution is decided by finite 2-adic data on the exponent. The certified classes already have density `>= 0.1199` at `K = 20`, and the empirical rate is about 0.95.

**Remark.** Conjecture HYP-9214 (debt resolution tends to 1 for generic large sources) would make this plausible but does not imply it, because the Mersenne family is thin.
