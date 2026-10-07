# Audit A — Collatz two-anchor theorems (THM-4600, THM-4601, HYP-9240; results note sections 0-7, 10-11)

Independent adversarial audit, 2026-10-07, for session mac-mini-2026-10-07-twoanchor.
* Scope: the committed text at `2ba93dab07` and the session's later working-tree additions (`heads_plus1.py`, `short_types.py`, `barrier_extend.py`).
* All verification used the auditor's own exact-arithmetic code in this directory; the session's scripts were read, not run.
* Saved by the main session from the auditor's final message, because subagents could not write report files.

## Verdicts

| # | Claim | Verdict |
|---|---|---|
| A1 | THM-4600 (1) run transparency | CONFIRMED |
| A2 | THM-4600 (2) periodic anchors | CONFIRMED, sharp, no off-by-one; CORRECTION: which cycle point |
| A3 | THM-4600 (3) conjugation formula | CONFIRMED |
| A4 | THM-4600 Reading on positivity; results §7 "No re-anchoring" | CORRECTION: false for actual chains |
| A5 | "a = 2 is new" | mild OVERCLAIM |
| B1 | Compiler rules are D-chain absorptions | CONFIRMED |
| B2 | (i) lag-one form | CONFIRMED |
| B3 | (ii) +1 barrier | CONFIRMED; three wording/typing fixes |
| B4 | (iii) universal state | CONFIRMED |
| B5 | (iv) compiler | identity CONFIRMED; CORRECTION: depths, class count; OVERCLAIM: "complete grammar", "Σu ≤ 22" |
| B6 | (v) "state depends only on (type, D)" | CORRECTION: false for D ≥ 5 |
| B7 | (vi) Mersenne shells | CONFIRMED |
| B8 | numbers | values CONFIRMED; CORRECTION: they are unconstrained-bit quantities |
| B9 | "(ii) is sharp", "any K-uniform rule", orphans | OVERCLAIM |
| C1–C4 | HYP-9240 | statement CONFIRMED; the reduction's rounding, the evidence ranges and "every K-uniform certificate" need fixes |

## Corrections

**A2. Which cycle point.** In THM-4600 (2), "c a point of that cycle" should be `c = c_w`, the point at which w is read. From another cycle point one reads the corresponding rotation of w (2348/2348 tests). The condition `v2(v−c) > Σw` is sharp: all 348 boundary cases fail.

**A4. Re-anchoring does occur.** Take an r = 1 source with `Y_end ≡ −1 mod 2^L`.
* The chain goes from (3, 1−27) to (−5, 3^(−5) − 1), anchored at −1, in 12 Terras steps; the limit pair is (−53, −1).
* The rest of the ones-run is then transparent.
* Verified on integer sources with K = 10, 50, 400 at post-run times 12..39.

**A5. Novelty.** Say "new in the repo (elementary; no literature search)" and cite the finite precursor in `collatz_collision_dp_20261007.md` (pure-two prefixes of lengths 1..16).

**B3. Wording of (ii).**
* "y_D* tends to 0" should read "reaches 0 (Terras time 1, resp. 3)".
* "before Terras time j + 1, i.e. before the end of the two-run" should read "at any Terras time ≤ j after x; in particular not within the two-run, which ends at 2J = j − r".
* `k_inf/D → 1` should be typed NUMERICAL: [0.75, 1.47] for D ≤ 3000; [0.971, 1.045] for 2000..3000; [0.966, 1.022] for 3001..8000.

**B5. Compiler.**
* *Depths.* The heads absorb at post-run Terras depth Σu + 1, i.e. 12 (r = 1) and 11 (r = 2), which is Terras time j + 11 and j + 9 after x. This holds on 320/320 sources, and these are the earliest possible absorptions.
* *Class count.* The r = 2 class count is 32 of 8192, not 16; c = 1 is valid.
* *Grammar off-by-one.* `heads_plus1.py` builds totals ≤ SMAX − 1, so its 6237 pairs have Σu ≤ 21.
* *Grammar not complete.* There are also:
  * child ladders i ≥ 2: `Σu = Σv + 2i`, `F_v(1) = 4^i F_u(1) + (4^i−1)/3`, final letters (c, c+2i); the first is (1,2,9)~(1,1,1,1,1,1), i = 3, depth 13;
  * source ladders: `Σv = Σu + 2i`, `F_u(1) = 4^i F_v(1) + (4^i−1)/3`, final letters (c+2i, c); the first is (12,1,1)~(3,1,1,3,2,6), depth 17.
* *Masses.* The exact absorbed mass by depth 22 is 0.016943 (r = 1) and 0.032581 (r = 2); the i = 1 grammar carries 82.0% and 90.4% of it. The r-conditioned masses at depth 20 are 0.012386 and 0.026512, so the comparison with the unconstrained BFS value 0.0241 is invalid.

**B6. Statement (v).** The end-of-run state is a function of (D, J), the same for r = 1, 2, and constant only for J ≥ J_0(D).
* J_0 = 3, 3, 6, 6, 8, 8, 10, 10, 15, 15 for D = 3..12, and up to 2965 for D ≤ 3000.
* For D = 5, 6 the states are (6,−296), (6,−404), (6,−485) for J = 3, 4, 5, then (6,−728).
* r enters only through the forced post-run prefix: (1,1) for r = 1, (1,0,0) for r = 2.

**B8. Realisable bits.** The fair-bit table and "shortest pattern 9" are unconstrained-bit quantities; the length-9 pattern starts with an even bit, which cannot occur after a two-run. Realisable values:

| s | 30 | 100 | 200 | 300 | 400 | 1000 |
|---|---|---|---|---|---|---|
| ρ_1(s) | 0.038 | 0.172 | 0.269 | 0.328 | 0.375 | 0.519 |
| ρ_2(s) | 0.053 | 0.184 | 0.278 | 0.338 | 0.384 | 0.529 |

**B9. Scope and sharpness.**
* "(ii) is sharp" should read "sharp up to an additive constant".
* "Any K-uniform rule has depth ≥ 2J" is false once 3-adic conditions are allowed: `m = (2n−1)/3 = 2^(K+1)(t/3) − 1` is a depth-0 certificate when 3 | t (394/1200 sources).
* Orphans are "2% of the odd exponents (50 of the 2500 odd a in [10^3, 6·10^3])".
* The maximum N_D is 2527 for D ≤ 8000.

**C2–C4. HYP-9240.**
* *Exact debt.* `k_entry = ⌈(s0 + σ_T(N_D))/2⌉ − odd(N_D)`. The sufficient condition is `excess(N_D) < s0/2`; the version with `⌈s0/2⌉` fails when s0 and σ_T are both odd.
* *Evidence ranges.* Debts lie in [0.75D, 1.47D] for D ≤ 3000 and in [0.90D, 1.10D] for 163 ≤ D ≤ 3000. The excess is ≤ 9 for N ≤ 880 (attained at N = 871). s0 ∈ [1.79D, 2.22D] for D ≥ 100.
* *Scope.* Restrict to K-uniform collision certificates on unrefined 2-adic classes. These are exactly the generalized deletions `3^a(n+1)/2^b − 1`, i.e. (b−a)-chains.

## Confirmed numbers

* **THM-4600 (1)–(3):** 138,502 exact checks, including negative k, plus an exhaustive check mod 2^14.
* **(i):** checked on 3000 sources.
* **(ii) limit chains:**
  * D ≤ 3000: none absorbs; only D = 1, 2 go to 0; in phase / out of phase 2288 / 710; minimum debt 3 at D = 3, 4; minimum out-of-phase debt 13; max N_D = 880 (at D = 253).
  * 3001 ≤ D ≤ 8000: none absorbs; 3610 / 1390; max N_D = 2527.
  * The finite-j reduction is exact on 8947 actual chains.
* **(iii):** 840 sources, K = 4..1500, J = 3..150; validity exhaustive mod 27·2^12.
* **(iv):** the identity holds on 13,272 (K, J, c) triples; 320 sources certified with exact words.
* **(vi):** the shell bijection checked for m = 4..10; both large examples verified by direct orbits (common odd endpoints of 2983 and 8310 bits); the shell fractions reproduced exactly.
* **Monte Carlo** from (3, −26) with fair bits: P(100) = 0.178, P(1000) = 0.523, P(4000) = 0.722. The BFS counts for lengths 9–20 match.
* **Haar-weighted branch coverage:** 0.646 ± 0.014.

## Prior art

Neither of these states the +1-anchored relation:
* "corresponding stems" for consecutive-integer pairs (arXiv:1511.09141);
* coalescence (arXiv:2005.09456).

## Scripts

All in this directory, each with a matching `.out`:
* `a1_transparency.py`
* `a2_barrier.py`
* `a2b_barrier_extend.py`
* `a2c_debt_ratios.py`
* `a3_universal_state.py`
* `a4_compiler.py`
* `a5_types_mersenne.py`
* `a6_numerics.py`
* `a6b_grammar.py` (outputs at depth 19 and `_23`)
* `a6c_grammar_depths.py`
* `a7_misc.py`
* `a8_branch_coverage.py`
