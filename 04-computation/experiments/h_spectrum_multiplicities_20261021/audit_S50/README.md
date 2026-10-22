# Independent audit of the S50 packet (THM-4744, THM-4745, HYP-9355, HYP-9356)

These are the independent auditor's own scripts and outputs, kept verbatim as audit evidence. They import none of the session's code. `alib.py` holds the shared helpers, and `a0_selftest.py` checks them against networkx and brute force first.

Scripts marked ✓ print ALL CHECKS PASSED. The others are measurement or replication scripts, whose verdicts are read off their tallies.

| script | output | content |
|---|---|---|
| `a0_selftest.py` ✓ | `a0_selftest.out` | helper self-test |
| `a1_product_law.py` ✓ | `a1_product_law.out` | the chain product law for 1001 compositions; `NT_s` counts for `s = 3..17`; path-by-path restriction bijection for 46 chains (4430 paths); the exchange inequality |
| `a2_recognition.py` ✓ | `a2_recognition.out` | scores, `c3 = m − 2`, rigidity (also by networkx), spine recovery, non-isomorphism, converse, palindromes, for `m ≤ 12` |
| `a3_ear_census.py` ✓ | `a3_ear_census.out` | an independent ear census from `C_3`, giving strong-class counts through order 8, the inventories below 375 at orders 9–11, `E(m)` through 12, and `inventory_375.json` |
| `a4_gluings.py` | `a4_gluings.out` | the exact-gluing census and the failure modes, plus 8 extra pairs |
| `a5_inventory_checks.py` ✓ | `a5_inventory_checks.out` | all 184 `P_d` below 375; the `N(n,h)` table; the true cubic laws for `h = 45, 75`; the second-floor counterexamples; the order-7 holes |
| `a5b_m3_predictions.py` | `a5b_m3_predictions.out` | shows that the first version's (M3) predictions were the truncated series |
| `a6_misc_claims.py` ✓ | `a6_misc_claims.out` | the exact HYP-9355 tables (`n = 7` by all `2^21` labelled tournaments), Paley values, the `q = 27` normalization, other numerical claims |
| `a7*_residue_replication.py`, `a7c_pool_n11.py` | `.out` | independent residue sampling at `n = 11, 12`, and the pooled verdict on P2 |
| `a8_repunit_primes.py` | `a8_repunit_primes.out` | the repunit-versus-prime comparison and the residues mod 16 |
| `b1_oneway_gap.py` ✓ (with `b1_lib.py`) | `b1_oneway_gap.out` | second pass: the 6,688 configurations (chains with `N ≤ 8`) that the first written proof of THM-4745 (5) did not reach, certificates for the 2,238 with `|T| ≤ 13`, all completions for 574, the Step 5 identities for 237 chains, the Init/Term lemma on 400 random tournaments |
| `b2_patch_witnesses.py` ✓ | `b2_patch_witnesses.out` | second pass: the explicit (M), (E), (E′) witness for every one of the 6,688 configurations (THM-4745 Step 4) |
| `b3_ear_census_625.py` ✓ | `b3_ear_census_625.out` | second pass: an independent ear census below 625 (13,439 / 8,769 / 3,536 / 658 strong classes at orders 9–12, none at order 13, `E(13) = 1`) |
| `b4_steps12.py` ✓ | `b4_steps12.out` | second pass: Steps 1 and 2 of THM-4745 (5) (all 196 back-arc choices for `F`, `F′`, the 78 one-way combinations, the 35 interior configurations) |

The audit's verdicts and the session's corrections are recorded in MISTAKE-652 and in the Audit sections of THM-4744, THM-4745 and the note `05-knowledge/results/h_spectrum_multiplicities_20261021.md`.
