---
id: THM-4603
title: "Run-compressed pair chains: (1) drift lemma - if u - c = 3^k (v - c') with c, c' T-periodic points (e.g. cycle points of U-words w, w', or the halving point 0) and v2(v - c') >= j*Lambda (Lambda = lcm of the Terras periods), both orbits run their periodic patterns in lockstep for j common periods and the debt drifts by the difference of the odd densities per Terras step (THM-4600 is the zero-drift case); (2) ladder completeness - in an admissible pair-chain state (3^max(0,-k) e in Z), every equal-time merge in which both orbits have an odd value before the merge is a ladder collision: the last distinct odd values satisfy z' = 4^i z + (4^i - 1)/3, i >= 1 (child or source ladder), so the ladder compiler of THM-4601 (iv), with all i and both directions, is complete, and together with the concurrent signed-gap decoder it enumerates all class-decided common futures from a marked state; (3) conditional on eventual periodicity of the rational limit orbit, a long run of any periodic pattern in either orbit sends any state to absorption, re-anchoring, zero drift without anchoring, or a linear debt drift; from the universal state (3, 1-27) a child ones-run re-anchors at -1 with debt -5 and a source ones-run drifts at +1/2 (a mechanism for unbounded post-run debt)"
status: >
  PROVED: (1); (2) for admissible states in which both orbits have an odd value before the merge (the 4z+1 ladder and the odd predecessors
  (2^a m - 1)/3 are classical; (1) and (2) are elementary bookkeeping on classical facts); (3) the finite-L reduction (that of THM-4601 (ii))
  and the classification conditional on the eventual periodicity of the rational limit orbit (open in general: a 3x+d periodicity question,
  which contains the Collatz problem on negative integers when c' = -1). FINITE-EXACT: (1) on 720 session tests (k in [0,6]) and 2,748 exact
  tests by audit C (k in [-9,9], boundary precision, integral cycle points, the halving point 0); (2) on 1651 session merges and 2,948 audit
  merges from 10 named and 1,500 random admissible states; (3) all 532 table lines re-derived by audit C, 550 limit pairs all periodic,
  9 transitions checked on actual integers. Independently audited 2026-10-07 (audit C); corrections applied (scope of (2), conditional
  classification, realisable lists, SHIFT label, heavy-tail wording); MISTAKE-587. The same scope correction to (3) was made independently and
  concurrently by codex-complement (05-knowledge/results/collatz_run_classification_scope_20261007e.md): without an authenticated
  eventual-periodicity witness the outcome is UNRESOLVED; its exact control is the universal-state child run on (2,4), which lands on
  distinct cycles of periods 12 and 6, both of odd density 1/3.
session: mac-mini-2026-10-07-twoanchor (continuation)
source: 05-knowledge/results/runcompress_orphans_cayley_20261007.md
scripts:
  - 04-computation/experiments/runcompress_20261007/ladder_and_drift.py, drift_check_actual.py (ALL CHECKS PASSED), transitions.py (+ .out)
  - 04-computation/experiments/runcompress_20261007/audit_C/ (c1_drift_lemma.py, c2_ladder_completeness.py, c3_transitions.py, c3_heavy_tail.py, REPORT.md)
related:
  - THM-4600 (run transparency: the zero-drift case); THM-4601 (the residual branch; (iv) the ladder compiler, made complete here)
  - THM-4581 (the pair chain); THM-4594 (class-decided certificates)
  - 05-knowledge/results/collatz_twoanchor_head_decoder_20261007b.md and collatz_complement_ladders_20261007c.md (concurrent signed-gap decoder)
---

# THM-4603 — run-compressed pair chains

## Setting

* T is the Terras map on `Q ∩ Z_(2)`.
* A pair-chain state is the relation `u = A(v) = 3^k v + e` at aligned Terras times. The debt k changes by `par(u) − par(v)` per step.
* The state is **admissible** if `3^max(0,−k) e ∈ Z`.
* A T-periodic point c has a Terras period P_c and an odd density. Examples:
  * the cycle point `c_w = B_w/(2^(Σw) − 3^|w|)` of a U-word w, with density |w|/Σw;
  * the halving point 0, with density 0.

## Statements

**(1) Drift lemma.**
* Let c, c′ be T-periodic points and `Λ = lcm(P_c, P_(c′))`.
* Suppose `u − c = 3^k (v − c′)` (k ∈ Z) and `v2(v − c′) ≥ jΛ`. Then:
  * for jΛ Terras steps both orbits follow the periodic parity patterns of c and c′;
  * `T^(jΛ)(u) − c = 3^(k′)(T^(jΛ)(v) − c′)` with `k′ = k + jΛ(dens(c) − dens(c′))`.
* Strict inequality `>` is needed only for an exact last U-letter.
* Both deviations lose one bit of 2-adic precision per Terras step, so the runs end together.
* THM-4600 is the case c = c′.
* The case of the source on +1 (density 1/2) and the child near 0 (density 0) gives drift +1/2. This is the collapse of the reset children D ≤ 2 in THM-4601 (v), whose debt grows like J.

**(2) Ladder completeness.**
* Let the state at time 0 be admissible. Suppose the two orbits first coincide at Terras time s with zero debt, and each orbit has an odd value in [0, s). After a two-run both orbits start odd, so this holds there.
* Let Z be the first common odd value, and z_u, z_v the last odd values before s. Then `z_u ≠ z_v`, both map to Z under U, and for some i ≥ 1 either:
  * `z_v = 4^i z_u + (4^i−1)/3` (child ladder, letters c and c+2i); or
  * the same relation with u and v exchanged (source ladder).
* The heads `h_u`, `h_v` before z_u, z_v satisfy:
  * `|h_v| = |h_u| + k_0`;
  * `Σh_u = Σh_v ± 2i` (equal time; for an even start, Σh includes the initial halvings);
  * merge time `s = max(Σh_u, Σh_v) + 1`, independent of the final letter;
  * at an anchored state with anchor c*, `F_(h_v)(c*) = 4^i F_(h_u)(c*) + (4^i−1)/3`, or the mirror for source ladders.
* Both hypotheses are needed:
  * the reset state (1,1), where u = 3v + 1 with v odd, merges at s = 1 with no odd value before s;
  * inadmissible starts give half-ladders `z′ = 2^b z + (2^b−1)/3` with b odd (audit C: b = 1, 3).
* **Consequences.**
  * The ladder compiler of THM-4601 (iv) (anchor +1, debt 3), with every i and both directions, is complete.
  * The concurrent signed-gap decoder (collatz_twoanchor_head_decoder_20261007b.md, collatz_complement_ladders_20261007c.md) is complete for the identity `F_v(Y) = S^r(F_u(X))`, `S(z) = 4z + 1`, at a fixed marked state. It states that it is "not for arbitrary common-future diagrams". (2) supplies exactly that step: every such merge has this form, with r ≠ 0 at the last distinct odd values.

**(3) Long-run transitions** (conditional).
* Let one orbit run near the periodic orbit of c′ for L Terras steps. The other is then near the rational point `A(c′)` (child run) or `A^(−1)(c′)` (source run), and for s < L the chain follows the exact limit pair.
* If the limit pair is eventually periodic, and an authenticated witness is supplied (true in all 550 computed cases; open in general), then after a transient τ independent of L exactly one of the following holds. Without a witness the outcome is UNRESOLVED (collatz_run_classification_scope_20261007e.md).
  * **absorption** at a fixed time;
  * **re-anchoring**: same cycle, in phase, giving a fixed anchored state depending only on L mod the period;
  * **zero drift, not anchored**: same cycle out of phase, or a different cycle of equal odd density;
  * **drift**: different densities, debt `k_τ + ρ(L − τ) + O(1)` with ρ the density difference.
* After a two-run the first post-run letter is 1 or ≥ 3, so runs beginning with the letter 2, and runs near 0, cannot start there.
* From the universal state (3, 1−27), the realisable transitions among the tabulated cycles (one phase each; `transitions.out`) are:

  | outcome | child runs | source runs |
  |---|---|---|
  | re-anchoring | letters 3, 4, 6 (debts 4, 9, 30); the ones-run re-anchors at −1 with debt −5 after 12 steps (limit pair (−53, −1)) | letters 4, 6 (debts 3, −30) |
  | positive drift, +0.167 … +0.5 per step | 5, 7, 8, (1,6), (3,5) | 1 (+1/2; debt 6 at step 8; child limit 25/27 → trivial cycle), (1,2), (1,1,2), (1,1,1,2), (1,1,2,3), (1,2,1,3), (1,1,1,3), (1,1,2,2) |
  | negative drift, −0.176 … −0.433 per step (debt paid down) | (1,1,2), (1,1,3), (1,2,2), (1,1,1,2), (1,1,1,3), (1,1,2,2), (1,1,2,3) | 5, 7, 8, (1,5), (1,6), (3,4), (3,5) |
  | absorption at a fixed time | (1,1,4) at 23, (1,1,4,2) at 37 | (1,2,5) at 20, (1,5,2) at 31 |

## Proofs

**(1).** Iterate `T_β(u) − T_β(c) = λ_β(u − c)` with λ ∈ {3/2, 1/2} along both periodic patterns. The ratio of the accumulated factors over Λ steps is `3^(Λ·(dens(c) − dens(c′)))`.

**(2).**
* Strictly before s the orbits differ.
* Their last odd values are distinct and lead, via halving chains, to the common odd value Z.
* For integers the ladder is the classical description of the odd U-preimages of Z.
* In general, the equal-time and debt conditions give `3^(|h_u|+1) e_0 + 3·2^(g_u) B_u + 2^(t_u) = 3·2^(g_v) B_v + 2^(t_v)`, where t is the Terras time of the last odd value and g counts initial halvings. By admissibility the first term is ≡ 0 mod 3. So `2^(t_u) ≡ 2^(t_v)` mod 3, the letter difference is even, and the relation is a ladder.

**(3).** Exact rational iteration of the limit pair until (u, v) repeats; the drift over the period is then read off. ∎

## Reading

* The pair chain is a debt walk whose long-run increments are differences of odd densities. Re-anchoring keeps the debt bounded; drift moves it linearly in the run length.
* **Positive-drift runs explain why no bounded post-run depth suffices.** A source ones-run of length L raises the debt by about L/2: a new Mersenne-like excursion.
* **They do not cause the s^(−1/2) tail.** Such runs have Haar probability 2^(−L); the tail is the diffusive return time of the zero-drift debt walk (THM-4581). From the universal state, survival to s = 3000 is 0.316 overall and 0.312 among paths with no source ones-run ≥ 14 (audit C).
