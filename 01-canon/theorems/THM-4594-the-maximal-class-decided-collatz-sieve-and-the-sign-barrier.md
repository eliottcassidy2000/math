---
id: THM-4594
title: "The maximal class-decided Collatz sieve and the sign barrier: certifying n mod 2^K by every smaller number reachable backward from any orbit point decided by the class (not only by descent, translation joins [Roosendaal/Barina] or depth-1 predecessors [Angeltveit 2026]) leaves 6,915,181 classes mod 2^30 (0.644%) and 0.493% mod 2^34 (descent: 0.884%); hence a minimal Collatz counterexample n_0 (> 2^71) lies in those 6,915,181 classes, with n_0 mod 9 in {0,1,3,6,7} and e.g. n_0 not = 539, 615 mod 1024; but the class of the 2-adic fixed point -1 is never certified, the classes of -5 and -17 stay uncertified to depths 93 and 88, no weight log n + g(n mod M) and no pathwise pair-chain weight can decrease, so any proof must use the sign of n"
status: >
  PROVED: maximality of the certificate family (every certificate given by one affine m(n) on the class is a backward word;
  ratio 1 = translation join); the residue constraints on a minimal counterexample, given Barina 2025's verification to 2^71
  (every certificate at depth <= 30 applies for n > 2^6.76; exact per-certificate thresholds); -1 never certified; the
  bounded-memory and pathwise-weight barriers. FINITE-EXACT:
  - exact class counts to K = 34; validation: descent = OEIS A076227 to K = 36; descent + joins gives Roosendaal's published 1720 at 2^16;
    branches agree with an independent big-integer brute force for K <= 11;
  - joins subsumed by branches for K <= 30;
  - -5 and -17 uncertified to depths 93 and 88 (exhaustive at each depth).
  NUMERICAL: the decay rate C K^(-3/2) 2^(-0.0500K), whose constant falls from about 5.6 (descent) to about 3.1 (maximal).
  CONJECTURE: the survivor lemma (HYP-9231). KNOWN prior art: Roosendaal's sieve (translation joins), Barina 2020 (implementation,
  credited to Roosendaal), Angeltveit arXiv:2602.10466 (depth-1 predecessor joins; deeper versions not implemented).
  Found by the session's core reader; the 539 mod 1024 exclusion was re-checked independently here (500 lifts).
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

## Statements

1. **Maximality (PROVED).** Every certificate of the form "one affine `m(n)` on the class" is a branch. Ratio exactly 1 is a translation, i.e. a join. So branches ∪ joins is the maximal class-decided family.
2. **Exact counts (FINITE-EXACT).** Uncertified classes mod `2^K`:

   | `K` | descent | + joins (Roosendaal) | + depth-1 (Angeltveit-type) | depth-1 + joins | **maximal** |
   |---|---|---|---|---|---|
   | 16 | 2114 | 1720 | 1856 | 1562 | **1363** |
   | 24 | 286581 | 234156 | 244392 | 206402 | **172868** |
   | 30 | 12771274 | 10446423 | 9910223 | 8385079 | **6915181** |
   | 32 | 41347483 | 33880411 | 34512882 | — | **23797887** |
   | 34 | 151917636 | — | — | — | **84720794** (0.4931%) |

   * At `K = 30` the maximal sieve leaves 0.541× descent, 0.662× joins and 0.825× (depth-1 + joins).
   * Joins add nothing beyond branches for `K ≤ 30` (FINITE-EXACT; HYP-9231 is the general form).
   * With `n mod 9`: 0.409% remain at `K = 28`, against 0.490% for the best classical combination.
3. **A minimal counterexample (PROVED given Barina 2025, Collatz verified to `2^71`).**
   * Every certificate at depth `≤ 30` applies for `n > 2^6.76`, with exact per-certificate thresholds.
   * So a least non-convergent `n_0` satisfies `n_0 mod 2^30 ∈ U_30`, with `|U_30| = 6,915,181` (0.644%), and `n_0 mod 9 ∈ {0, 1, 3, 6, 7}`.
   * Exclusions beyond the classical sieves start at `2^10`: `n_0 ≢ 539, 615 (mod 1024)`, via the depth-5 branches `m = (27n+7)/32` and `(27n+3)/32`. Example: `455 → 683 → 1025 → 1538 → 769 → 1154 = T^10(539)`.
4. **The sign barrier (PROVED; FINITE-EXACT depths).**
   * (a) The class of the 2-adic fixed point −1 (`n ≡ −1 mod 2^K`) is certified at **no** depth: on it, every backward word with `b ≤ a` has ratio `≥ 1`.
   * (b) The classes of −5 and −17 stay uncertified to depths 93 and 88, while every other point of their cycles is certified by depth 9. The sieve isolates exactly the grand-orbit minima on `Z_{<0}`, as it must for a minimal counterexample.
   * (c) No weight `W(n) = log n + g(n mod M)` decreases along orbits over bounded windows. Take `n ≡ −1 mod 2^(L+P)` on a period-`P` orbit of `(3r+1)/2 mod M_odd`: the `g`-terms telescope while `log n` gains `P log(3/2) > 0`.
   * (d) No function of the pair-chain state decreases along every transition. The positive pair `(2, 1)` produces the chain cycle `(0, 1) → (−1, 1/3) → (0, 1)` under the bits `(10)^∞`, and the Mersenne states `(h, 3^h − 1)` are fixed under all-odd bits.
   * (e) Any criterion stated only for rational or eventually periodic 2-adic points would also cover −1, −5 and −17.
   * So every argument blind to the sign (the 2-adic tail `…000` of positive integers against `…111`) fails.
5. **Certificates are deep (FINITE-EXACT up to the search budget).**
   * For 37 stopping-time record holders, the certificate depth `K(n)` is 0.77–1.00 × the stopping time `σ(n)`, and up to 12.4 × `log_2 n` (`n = 27`: `K = σ = 59`). Five were re-checked with 50× the budget, unchanged.
   * So certificates cannot be confined to the first `log_2 n` bits, which is where density methods stop.

## The remaining conjecture

**Collatz ⟺ every positive integer `n > 2^71` has a class-decided certificate, at some finite depth, whose threshold is below `n`.** Equivalently, `U_∞ ∩ Z_{>0} = ∅`, where `U_∞` is the closed, Haar-null set of 2-adic integers whose class is never certified.

* By 4, `U_∞` contains −1 and (to depth ~90) −5 and −17.
* So the conjecture is precisely a statement that the positive integers avoid a null closed set that contains the negative cycles' minima.
* Measure methods (THM-4581/4590) cannot see the sign.

## Reading

* This is the integer-certified, all-relations version of THM-4590's coalescence.
* It improves verification sieves by a constant factor (about 1.8× over descent) without changing the exponential rate `2^(−0.0500K)`. That rate is fixed by the never-descending walk excursions.
