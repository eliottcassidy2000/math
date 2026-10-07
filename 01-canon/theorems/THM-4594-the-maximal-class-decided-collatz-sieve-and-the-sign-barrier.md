---
id: THM-4594
title: "The maximal class-decided Collatz sieve and the sign barrier: certifying n mod 2^K by every smaller number reachable backward from any orbit point decided by the class (not only by descent, translation joins [Roosendaal/Barina] or depth-1 path merging and the odd-even-even rule [Angeltveit 2026]) leaves 6,915,181 classes mod 2^30 (0.644%; 0.957x Angeltveit's 2-adic rules) and 0.493% mod 2^34 (descent: 0.884%); hence a minimal Collatz counterexample n_0 (> 2^71) lies in those 6,915,181 classes, with n_0 mod 9 in {0,1,3,6,7} and n_0 not = 11247, 12191, 12799, 23743 mod 2^15 (the first exclusions beyond published rules); for unrefined 2-adic classes the class of the fixed point -1 is never certified and those of -5 and -17 stay uncertified to depth 120, and no weight log n + g(n mod M) and no pathwise pair-chain weight can decrease, so a proof using only 2-adic class-decided certificates must use the sign of n (3-adic refinement certifies all three classes)"
status: >
  PROVED: maximality of the certificate family (every certificate given by one affine m(n) on the class is a backward word;
  ratio 1 = translation join); the residue constraints on a minimal counterexample, given Barina 2025's verification to 2^71
  (every certificate at depth <= 30 applies for n > 2^6.76; exact per-certificate thresholds); -1 never certified BY UNREFINED
  2-ADIC CLASSES (3-adic refinement certifies it); the bounded-memory and pathwise-weight barriers. FINITE-EXACT:
  - exact class counts to K = 34; validation: descent = OEIS A076227 to K = 36; descent + joins gives Roosendaal's published 1720 at 2^16;
    branches agree with an independent big-integer brute force for K <= 11;
  - joins subsumed by branches for K <= 30;
  - -5 and -17 uncertified (unrefined classes) at every depth <= 120 (exhaustive; audit A2).
  NUMERICAL: the decay rate C K^(-3/2) 2^(-0.0500K), whose constant falls from about 5.6 (descent) to about 3.1 (maximal).
  CONJECTURE: the survivor lemma (HYP-9231). KNOWN prior art: Roosendaal's sieve (translation joins), Barina 2020/2025 (implementation,
  credited to Roosendaal; 2^34 descent+joins and depth-0 3^k sieves), Angeltveit arXiv:2602.10466 (depth-1 path merging AND the
  odd-even-even rule; deeper path merging described but not implemented). The full branch closure is new relative to these, but its
  gain over the published fixed-depth rules is about 4% (0.957x at K = 30), not the 18% first claimed (corrected after audit A2, MISTAKE-585).
  Independently audited (A2): counts re-derived exactly to K = 34; descent = A076227; thresholds peak at 108.01 = 2^6.755.
session: mac-mini-2026-10-07-golden
source: 05-knowledge/results/golden_collatz_resonance_20261007.md
scripts:
  - 04-computation/experiments/golden_20261007_readers/core/ (sieve.c, sieve_bounds.c, out_*.txt, verify_branch.py, classify_small.py, survcollide.c, path_cert.py, cert_depth.py + .out)
related:
  - THM-4581, THM-4590 (density-one coalescence; the sieve is its finite, all-relations, integer-certified version)
  - HYP-9231 (the survivor lemma), Barina 2020/2025, Roosendaal (ericr.nl/wondrous), Angeltveit 2026
---

# THM-4594 — the maximal class-decided sieve and the sign barrier

## Setting

* A **class** is a parity word of length `K`, i.e. a residue `n mod 2^K`, optionally refined by `n mod 3^J`.
* On a class, the orbit points are affine: `x_s = (3^(a_s) n + c_s)/2^s` for `s ≤ K`.
* A **certificate** is an affine `m(n) < n`, valid for all large `n` in the class, whose orbit meets `n`'s. By strong induction, a certified class is Collatz-safe once all smaller numbers are.
* **Descent:** `x_s` itself, when `3^(a_s) < 2^s`.
* **Join:** a translate `n − d` that merges with `n` (pair chain from `(0, −d)`).
* **Branch:** `m = (2^i x_s − c(u))/3^b` for a backward word `u` with `b` odd steps and ratio `3^(a_s − b)·2^(i−s) < 1`.
* Angeltveit's published 2-adic rules are a sub-family of branches: depth-1 path merging, and the odd-even-even rule (`m = T^k(n)` followed by `ℓ ≥ 1` odd steps and then 2 even steps joins `(m−1)/2`).

## Statements

1. **Maximality (PROVED).** Every certificate of the form "one affine `m(n)` on the class" is a branch. Ratio exactly 1 is a translation, i.e. a join. So branches ∪ joins is the maximal class-decided family.
2. **Exact counts (FINITE-EXACT).** Uncertified classes mod `2^K`:

   | `K` | descent | + joins (Roosendaal) | + depth-1 path merging | depth-1 + joins | Angeltveit's 2-adic rules (descent + path merging + odd-even-even; audit A2) | **maximal** |
   |---|---|---|---|---|---|---|
   | 16 | 2114 | 1720 | 1856 | 1562 | 1370 | **1363** |
   | 24 | 286581 | 234156 | 244392 | 206402 | 178304 | **172868** |
   | 30 | 12771274 | 10446423 | 9910223 | 8385079 | 7227826 | **6915181** |
   | 32 | 41347483 | 33880411 | 34512882 | — | — | **23797887** |
   | 34 | 151917636 | — | — | — | — | **84720794** (0.4931%) |

   * At `K = 30` the maximal sieve leaves 0.541× descent, 0.662× joins and **0.957× Angeltveit's 2-adic rules**. The two agree exactly for `K ≤ 14`.
   * Joins add nothing beyond branches for `K ≤ 30` (FINITE-EXACT; HYP-9231 is the general form).
   * With `n mod 9`: 0.409% remain at `K = 28`. The comparison figure 0.490% omitted the odd-even-even rule and is withdrawn (audit A2).
3. **A minimal counterexample (PROVED given Barina 2025, Collatz verified to `2^71`).**
   * Every certificate at depth `≤ 30` applies for `n > 2^6.76`, with exact per-certificate thresholds.
   * So a least non-convergent `n_0` satisfies `n_0 mod 2^30 ∈ U_30`, with `|U_30| = 6,915,181` (0.644%), and `n_0 mod 9 ∈ {0, 1, 3, 6, 7}`.
   * `n_0 ≢ 539, 615 (mod 1024)`, via `m = (27n+7)/32` and `(27n+3)/32`. These are already excluded by Angeltveit's odd-even-even rule. (Corrected after audit A2: they were first presented as new.)
   * The first exclusions beyond descent, joins and all of Angeltveit's 2-adic rules are at `2^15`: `n_0 ≢ 11247, 12191, 12799, 23743 (mod 2^15)`. For example, 12799 has `m = (6561n + 1953)/8192`, the depth-2 "m ≡ 8 mod 9" variant that Angeltveit describes but did not implement. Verified on 300 lifts (audit A2).
4. **The 2-adic sign barrier (PROVED for unrefined classes, `J = 0`; FINITE-EXACT depths).**
   * (a) For unrefined classes, the class of the 2-adic fixed point −1 (`n ≡ −1 mod 2^K`) is certified at **no** depth: on it, every backward word with `b ≤ a` has ratio `≥ 1`.
   * (b) For unrefined classes, those of −5 and −17 stay uncertified at every depth `≤ 120` (exhaustive; audit A2), while every other point of their cycles is certified by depth 9.
   * **With 3-adic refinement the barrier disappears** (audit A2). All three classes are certified at depth 0 by one backward turn around their cycle, and each certificate's fixed point is the cycle point. The ratio `2^L/3^a < 1` exactly because the cycle is negative:
     * −1: `(2n−1)/3` on `n ≡ 2 (mod 3)`;
     * −5: `(8n−5)/9` on `n ≡ 4 (mod 9)`;
     * −17: `(2048n − 2363)/2187` on `n ≡ 2170 (mod 2187)`.
   * (c) No weight `W(n) = log n + g(n mod M)` decreases along orbits over bounded windows. Write `M = 2^v·M_odd` and take `n ≡ −1 (mod 2^L·M_odd)` with `L > v`.
     * −1 is a fixed point of `T`. For `i ≤ L − v`, `T^i(n)` is odd and satisfies `T^i(n) ≡ −1 (mod 2^(L−i))` and `(mod M_odd)`, because `(3(−1)+1)/2 = −1` and 2 is invertible mod `M_odd`.
     * So `T^i(n) ≡ −1 (mod M)` and `g(T^i n)` is constant, while `log T^i(n)` increases by `log(3/2) + O(1/n)` per step.
     * `L` is arbitrary, so `W` increases over arbitrarily long windows. (Proof corrected after audit A2: the earlier period-`P` sketch controlled only multiples of `P`.)
   * (d) No function of the pair-chain state decreases along every transition. The positive pair `(2, 1)` produces the chain cycle `(0, 1) → (−1, 1/3) → (0, 1)` under the bits `(10)^∞`, and the Mersenne states `(h, 3^h − 1)` are fixed under all-odd bits.
   * *Remark (e).* A criterion stated only for rational or eventually periodic 2-adic points would also have to handle −1, −5 and −17.
   * So a proof that uses only 2-adic class-decided certificates must use the sign: the 2-adic tail `…000` of positive integers against `…111`. 3-adic data, or archimedean thresholds, are the other ways around the barrier.
5. **Certificates are deep (FINITE-EXACT up to the search budget).**
   * For 37 stopping-time record holders, the certificate depth `K(n)` is 0.77–1.00 × the stopping time `σ(n)`, and up to 12.4 × `log_2 n` (`n = 27`: `K = σ = 59`). Five were re-checked with 50× the budget, unchanged.
   * So certificates cannot be confined to the first `log_2 n` bits, which is where density methods stop.

## The remaining conjecture

**Collatz ⟺ every positive integer `n > 2^71` has a class-decided certificate, at some finite depth, whose threshold is below `n`.**

* Let `U_∞` be the closed, Haar-null set of 2-adic integers whose unrefined class is never certified. **Collatz ⟹ `U_∞ ∩ Z_{>0} = ∅`, but the converse fails for cycles.** A positive cycle's minimum `n_0 = c/(2^L − 3^a)` lies in a descent-certified class whose threshold equals `n_0` exactly (compare 1 in the class `1 mod 4`, threshold 1). The threshold form above is the correct equivalent. (Corrected after audit A2.)
* By 4, `U_∞` contains −1 and, to depth 120, −5 and −17 (unrefined).
* Measure methods (THM-4581/4590) cannot see the sign or the thresholds.

## Reading

* This is the integer-certified, all-relations version of THM-4590's coalescence.
* It improves verification sieves by a constant factor: about 1.85× over descent, but only about 4.5% over Angeltveit's published 2-adic rules at `K = 30`. It does not change the exponential rate `2^(−0.0500K)`, which is fixed by the never-descending walk excursions.

**Audit (2026-10-07, independent audit A2).**
* CONFIRMED: the counts to `K = 34` (independent search), descent = A076227, the thresholds (max 108.01 = `2^6.755` at Ellison's `(27, 17)`), and HYP-9231 to `s ≤ 28`.
* Corrected above (MISTAKE-585):
  * the 539/615 novelty (already Angeltveit's odd-even-even rule), and the gain over published rules (0.957×);
  * the sign barrier, restricted to unrefined 2-adic classes;
  * the false `U_∞` equivalence;
  * the proof of (c).
