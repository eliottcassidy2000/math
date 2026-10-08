# Audit E of THM-4606 (coalescence needs contraction)

Independent adversarial audit, run 2026-10-08 as a session subagent. The subagent could not write report files in its environment, so the main session saved its final message here verbatim.

**The main theorem holds and I found no major errors.** For p ≥ 5, q_p(e) < 1 for every nonzero integer e and q_p(e) → 0 as |e| → ∞. For p = 1 and p = 3, translates merge almost surely. The proofs are complete. What needs fixing:
* one sentence that is false as written (statement 5);
* one numerical value that is about three times too high (q_13);
* several scope and wording problems.

Everything was re-derived by hand and re-implemented from scratch in Python 3.10 standard library. No existing file was edited and nothing was committed.

## Verdicts

| # | Statement | Verdict | Findings |
|---|---|---|---|
| 1 | Step table | CORRECT | E6 (nit) |
| 2 | Escape lemma | CORRECT | E7 (nit) |
| 3 | Integer translations | CORRECT WITH FIXES | E4 (minor); E9, E10 (nit) |
| 4 | Dichotomy in p | CORRECT WITH FIXES | E3 (minor); E12 (nit) |
| 5 | Integers (density transfer) | CORRECT WITH FIXES | E1 (minor, false as written); E8 (nit) |
| 6 | Sheet-blindness | CORRECT WITH FIXES | E5 (minor) |
| 7 | Numerics | CORRECT WITH FIXES | E2 (minor); E11 (nit) |
| — | Reading section | CORRECT WITH FIXES | E13, E14 (nit); E3 repeats here |

## Findings

**E1 (minor; statement 5 is false as written).** The sentence "For p = 5 this density is below 0.0087 for every K" holds only for e = ±1.
* Exact densities at K = 14: e = ±3 gives 0.0444, e = 6 gives 0.0230, e = 5 gives 0.00995.
* Witness: the whole class n ≡ 16 (mod 32) merges with n + 3 at t = 5, for example T⁵(16) = T⁵(19) = 3. The chain word from (0,3) is 00001, so q_5(3) ≥ 1/32.
* Fix: "For p = 5 and e = ±1 …". For e = ±1, q_5 = 0.00833 ± 0.00005, so the 0.0087 bound is safe there.

**E2 (minor; statement 7, status line, results note §0 and the results INDEX).** The value q_13 ≈ 5·10⁻⁶ comes from a single event and is about three times too high.
* My direct-orbit run found 8 merges in 4,000,000 pairs: 2.0·10⁻⁶, with 95% interval [0.86, 3.94]·10⁻⁶. The author's 5·10⁻⁶ lies outside that interval.
* Exact enumeration gives a rigorous floor: P(absorbed by t = 33) = 9099/2³³ = 1.06·10⁻⁶.
* Adding the Monte Carlo tail after t = 33 gives q_13 ≈ (1.6 ± 0.4)·10⁻⁶. The tail rests on only 2 events, so the ± is rough.
* Optionally sharpen q_11 to 2.43·10⁻⁵ ± 0.04·10⁻⁵ (rigorous floor 101183/2³² = 2.356·10⁻⁵).

**E3 (minor; statement 4, Reading, HYP-9244).** "p = 3 … with tail ≍ T^(−1/2)" overstates THM-4581. Its rate (statement 4) is only at sketch level, and the upper bound carries a factor (log T)². Fix: write c·T^(−1/2) ≤ P(no merge by T) ≤ C·T^(−1/2)(log T)², at sketch level.

**E4 (minor; Proof 3).** The sentence "the path either stays at k = −1 with |f| growing geometrically, or returns" is wrong.
* With E = pe, the run at k = −1 is the affine map E ↦ (pE + p − 1)/2. Its length is exactly v_2((p−2)·pe_1 + p − 1).
* The run is infinite only at the fixed point pe_1 = −(p−1)/(p−2), and there f stays constant rather than growing. Example: for p = 5, starting from e = −1/3 the path sits at (−1, −4/15) forever.
* For integer e this never happens when p ≥ 5, since (p−1)/(p−2) is not an integer. So the theorem's conclusion stands; replace the sentence with that argument.

**E5 (minor; statement 6).** "Statements 1–5 hold verbatim for px + r with e replaced by re" needs two qualifications.
* In the px + r map's own coordinates the step table picks up factors of r, so |η| ≤ |r|/2, which is not < 1 once |r| ≥ 3. The statements are verbatim only after conjugating back (e ↦ e/r).
* Statements 3 and 5 then cover only translations divisible by r. An integer translation not divisible by r becomes a non-integer e, where the growth lemma's integrality fails. For example, 5x + 3 with translation −1 corresponds to e = −1/3, exactly the stuck start in E4.
* My breadth-first search found escape paths for e = ±1/3, ±2/3 (p = 5, 7) and ±1/7, ±2/7 (p = 5). So the extension is very likely true, but it is not proved as written. Either restrict the statement to translations re with e a nonzero integer, or add the extension.

**E6 (nit; statement 1).** The sharp bound is |η| ≤ 1/2, attained. The stated "1/2 + 1/(2p) < 1" fails at p = 1, which the Setting allows, and at p = 1 the two multipliers coincide. State |η| ≤ 1/2 and take p ≥ 3 for the multiplier sentences.

**E7 (nit; statement 2).** The escape lemma is purely qualitative. With η = 1/2 the constant A_p(1/2) is about 10^981,058 for p = 5, 10^216,011 for p = 7, 10^135,166 for p = 9 and 10^93,571 for p = 13. Say so, so the lemma is not read as giving a numerical bound.

**E8 (nit; Proof 5).** Equal-time integer merges with unequal odd counts can only occur for at most one integer per (t, parity word), so they have density zero and the claimed equality holds. The proof should say this. I found none for n ≤ 2¹⁶ and t ≤ 40.

**E9 (nit; Proof 3, unequal times).** The clause "a merge at k_n ≠ 0 …" repeats the case "equal times with unequal odd counts", because with k_0 = 0, k_n is the odd-count difference. Merge the two.

**E10 (nit; Proof 3).** The last algebraic step for e < 0 is omitted. (p/2)(|e|/2 + 1/(2p) − 1/(p−2)) ≥ (p/4)|e| − 1 needs p ≥ 10/3: it holds for p ≥ 5 and fails at p = 3 (constant −5/4). Add one line.

**E11 (nit; status, statement 7, results note §3).**
* The first-merge times are labelled "Exact", but shortest_merge.py discards states with |f| > 10⁴. That cut is not a valid exclusion beyond depth about 14.
* My search prunes only states that provably cannot be absorbed in time (|f_t| > 2^(L−t) − 1, valid because |f'| + 1 ≥ (|f| + 1)/2). It confirms every value: 3, 11, 4, 11, 18, 23, 5, 13, 8, 18 and 6 for p = 3, 5, 7, 9, 11, 13, 15, 17, 21, 23, 31, with word counts 1, 1, 1, 2, 3, 3, 1, 1, 1, 1, 1.
* For p = 19, 25, 27 and 29 there is no merge up to depth 26, rigorously; for p = 19 none up to 30.
* The "200,000 Haar pairs" are simulated chain paths with random coins (equal in law), not actual orbit pairs.
* "Merges for p = 5 occurred as late as Terras time 188 … none later" is specific to that sample. I saw merges at t = 202–352 in 1.9 million samples, and none after 400.
* The "exact lower bound" column can be replaced by much stronger rigorous bounds:

  | p | rigorous lower bound |
  |---|---|
  | 5 | 18401/2²² = 0.004387 |
  | 7 | 449933/2²² = 0.10727 |
  | 9 | 36571/2²⁴ = 0.002180 |
  | 11 | 2.356·10⁻⁵ |
  | 13 | 1.059·10⁻⁶ |

**E12 (nit; statement 4).** "Among the maps px + 1 (p odd) … iff p < 4" should say p ≥ 1 odd; for negative p the natural criterion would be |p| < 4. Worth adding that p < 4 is exactly the Matthews–Watts contraction threshold m_0·m_1 < d^d = 4.

**E13 (nit; Reading).**
* "Two orbits coalesce exactly when that common dynamics contracts" is true only within the px + 1 family. The same results note has a rank-3 map on Z_5 that contracts yet does not coalesce (numerically); qualify the sentence.
* "|f|^(−θ) is the escaping quantity" is never proved; at departures it grows by 2^θ for both coins. Label it heuristic.

**E14 (nit; prior art).** The file makes no novelty claim, which is appropriate. My brief searches found no prior statement of this theorem. Suggested context to cite:
* Kontorovich–Lagarias (arXiv:0910.1944): the 2-adic 3x+1 and 5x+1 maps are measure-theoretically conjugate. So THM-4606 separates them by an additive property that the conjugacy does not preserve.
* Matthews–Watts (Acta Arithmetica 43, 1984): their contraction criterion has the same threshold.
* The two-point-motion picture for random iterated function systems: synchronization exactly when the Lyapunov exponent is negative (Homburg–Kalle, Advances in Mathematics 482, 2025).

## Evidence

**Task A: the table and Proof 1.**
* The hand derivation for general odd p reproduces the table. At every flip k' = k + 1 − 2β, whatever the sign of k.
* Table against direct 2-adic steps: 0 mismatches on 60,000 states. These included p = 1, |k| ≤ 6, and e with denominators beyond powers of p.
* Every displayed f' formula as an exact rational identity: 0 mismatches on 200,000 states with |k| ≤ 8. Off departures the two coins always give one 1/2 and one p/2; at departures both give 1/2.
* Chain against actual big-integer orbits: 0 failures across 3,000 pairs (308 absorptions), plus 4,000 pairs with k_0 < 0 done 2-adically.
* Flip coins read from actual orbits: P(up) between 0.492 and 0.508 in every cell, and lag-1 correlation at most 0.003. The skeleton is a simple random walk for every odd p.

**Task B: the escape lemma.** Every step checks:
* ξ_n is i.i.d. fair;
* ln λ ≥ ξ − ln p at departures;
* the number of departures is at most the number of visits to 0;
* C(2m,m)/4^m ≥ 1/(2√m) for m ≤ 10⁶;
* the exact visit law for n ≤ 400 respects exp(−(m−1)/(2√n)), with ratio at most 1;
* both δ_1 and δ_2 tend to 0.

The induction constants are consistent, with slack. The two per-step inequalities had 0 violations on 866,535 steps. The bound is uniform in k_0 and holds for every rational 2-adic e_0.

**Task C: the growth lemma.**
* 300,000 starts (odd p from 5 to 63, |e| ≤ 5000) passed everything with 0 failures:
  * the closed form e* = (p/2)^(j+1)·e_1 + A_j;
  * 0 < A_j < (p/2)^(j+1)/(p−2);
  * both sign-wise bounds;
  * |e*| ≥ (p/4)|e| − 1;
  * e* a nonzero integer;
  * no visit to (0,0);
  * the run-length formula.
* The thresholds x* = 4/(p−4) are right.
* My own breadth-first search reproduces the author's escape words for p = 5 and 7.

**Task D: statements 3–6.**
* The unequal-time argument is correct.
* Density transfer: in 9 (p, e, K) configurations with K up to 16, every residue was tested. Chain absorption agreed with the literal integer merge every time, and the verdict did not depend on the lift.
* The conjugation identity had 0 failures, and residue counts match across conjugate maps.
* p = 1 is correct, and p = 3 is correct as an almost-sure statement.

**Task E: numerics.**
* About 10.3 million actual orbit pairs (512–3000-bit y), plus my own chain runs, exact enumerations, and a path-by-path check that the author's update matches direct orbits (0 disagreements).
* There were no unequal-time meetings in any direct sample.

| p | author | direct orbits (mine) | hybrid: exact mass + Monte Carlo tail |
|---|---|---|---|
| 5 | 0.00834 ± 0.00020 | 0.008283 ± 0.000117 (600k); 0.00788 ± 0.00028 (100k, 3000-bit); 0.008412 ± 0.000144 (400k) | 0.008330 ± 0.000051 |
| 7 | 0.1084 ± 0.0007 | 0.10870 ± 0.00040 | 0.10857 ± 0.00005 |
| 9 | 0.00226 ± 0.00011 | 0.002282 ± 0.000062 | 0.002253 ± 0.000011 |
| 11 | 2.5·10⁻⁵ | 2.73·10⁻⁵ (109 of 4M) | 2.43·10⁻⁵ ± 0.04·10⁻⁵ |
| 13 | 5·10⁻⁶ | 2.0·10⁻⁶ (8 of 4M) | 1.6·10⁻⁶ ± 0.4·10⁻⁶ |

* Exact per-time counts match the direct histograms (|z| ≤ 1.2). The author's 400-step cutoff is harmless.
* The reset word for p = 2^a − 1 absorbs at time a + 1 for a = 2 to 7 (p = 3 to 127), so q_p ≥ 2^−(a+1) holds.

**Task F: the Reading and prior art.** Covered by E3, E13 and E14. The section's caveats (no individual orbit decided; nothing about divergence) are accurate.

## Not verified
* THM-4581 itself, which statement 4 relies on for p = 3.
* Any upper bound on q_p beyond statistical ones: theory gives only the qualitative escape lemma.
* q_13 beyond t = 33, which rests on 2 Monte Carlo events.
* The prior-art search was brief (about 8 queries) and cannot rule out an existing statement.
* The algebra is checked by hand plus exact identities at many points, not by computer algebra.

## Files
Everything is in `04-computation/experiments/reframes_20261007/audit_E/`:
* `a_table_check.py` / `.out` and `a4_reseed.py` / `.out` (Task A)
* `b_escape_check.py` / `.out` (Task B)
* `c_growth_bfs.py` / `.out` (Task C)
* `d_density_transfer.py` / `.out` (Task D)
* `e_direct_orbits.py` (run by `run_direct.sh`) → `e_direct_orbits.out`
* `e_shortest_exact.py` (run by `run_exact.sh` and `run_exact2.sh`) → `e_shortest_exact.out`, `e_shortest_exact2.out`
* `e_chain_tail.py` (run by `run_tail.sh` and `run_p5_confirm.sh`) → `e_chain_tail.out`, `e_p5_confirm.out`

## Sources
- Kontorovich–Lagarias, Stochastic models for the 3x+1 and 5x+1 problems: https://arxiv.org/abs/0910.1944
- Matthews map and the Matthews–Watts conjectures (MathWorld): https://mathworld.wolfram.com/MatthewsMap.html
- Homburg–Kalle, Iterated function systems of affine expanding and contracting maps on the unit interval: https://dare.uva.nl/id/77a6d642-a097-4f02-9538-08fc39e560f6
- Burson, On the distribution of the first point of coalescence for some Collatz trajectories: https://arxiv.org/pdf/2005.09456
- Lagarias, The 3x+1 problem: an annotated bibliography: https://arxiv.org/pdf/math/0309224
