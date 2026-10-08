---
id: THM-4605
title: "Mersenne-line exits by pair-chain absorption and the stopping-time invariant: an exit M_E ~> M_(E-D) is certified by the absorption of THM-4581's Mersenne-switch chain run with the source x = 2*3^(E-1) - 1 as reference, which needs only 3^(E-1) mod 2^W (divide-and-conquer parity vectors, the recursion of Elsenhans 2025); without absorption the two orbits share no value in the window (no-meeting lemma); certificates with merge value >= 3 preserve (sigma_T(x_K), odd count of M_K) and conversely, so for exponents whose orbits reach 1 the deletion clusters are its level sets, and no route of such certificates from any E > 61821 reaches any K <= 12800; the supplied t = 23 source and its 24 children form a closed fan whose bottom M_99708993677 escapes at Terras time 19,000,765 with 1916 deletions (D = 1..1910, 1929..1932, 1935, 1936; deepest M_99708991741), a certified class of 1949 exponents that one orbit computation would ground; the 2084-deletion level-2 cluster is not absorbed by 2^29 - 9; 1999 of the first 2000 members of the unpaid progression t = 13847 mod 27648 have certified exits"
status: "PROVED: statements 1, 2, 3 (soundness, no-meeting lemma, fast parity vectors) and 8 (stopping-time invariant and its corollaries). FINITE-EXACT: statements 4-7 and 9 (exact integer pair chains; coins from exact residues; the escape cross-checked from affine maps mod 2^(2^20)). The chain itself is KNOWN (THM-4581 (6d)); density-one exits for each D follow from THM-4581 (3) (note section 4.4); the exponent 1/2 of the no-exit tail for finite lag sets is THM-4593. Independently audited after checkpoint 023a49bcf, and again (focused) after the corrections (MISTAKE-588). Session opus-2026-10-07-S21."
session: opus-2026-10-07-S21 (owner directive: codex-grounding's status update; 'the next decisive target is a terminating, source-preserving rule for one of those actual children')
source: 05-knowledge/results/mersenne_line_barriers_20261007.md
scripts:
  - 04-computation/experiments/mersenne_line_barriers_20261007.py (+ .out; ALL CHECKS PASSED) (A), (A2), (B)-(F), (M)
  - 04-computation/experiments/mersenne_fan_escape_20261007.py (+ .out; ALL CHECKS PASSED) (P), (G), (H), (I)
  - 04-computation/experiments/mersenne_fan_block_20261007.py (+ .out; ALL CHECKS PASSED) (J), (K), (N)
  - 04-computation/experiments/mersenne_line_coalescence_20261007.py (+ .out; ALL CHECKS PASSED) (L), (Q), (Q2), (H2)
  - 04-computation/experiments/mersenne_fan_level2_20261007.py (+ .out; ALL CHECKS PASSED) level 2 to 2^29
related:
  - THM-4581 (mac-mini: the pair chain, its Haar absorption, (6d) the Mersenne-switch chain); THM-4593 (partner coalescence: the no-exit tail for finite lag sets); THM-4601 (i) (the reset rule) and (ii) (the +1 barrier); THM-4556 (ii), (iv)
  - 05-knowledge/results/collatz_residual_t23_20261007f.md and collatz_grounding_progress_20261007f.md (codex-grounding: the t = 23 exit and the 24 children)
  - 05-knowledge/results/collatz_parameter_complement_synthesis_20261007e.md (the family E(t), the unpaid progression)
  - HYP-9242 (orphan law; its partner criterion is the invariant of statement 8), HYP-9243 (diffusive coalescence)
---

# THM-4605 — Mersenne exits by pair-chain absorption; the stopping-time invariant; the t = 23 fan escapes

## Setting

* `T` is the Terras map and `M_E = 2^E − 1`. Since `T^j(M_E) = 3^j 2^(E−j) − 1` for `j ≤ E`, after `E − 1` steps the source is at `x = 2·3^(E−1) − 1`.
* The child `M_(E−D)` is at `y = 2·3^(E−D−1) − 1`, with `y + 1 = 3^(−D)(x + 1)`.
* **Source-reference pair chain.** This is THM-4581's table with x as reference. THM-4581 (6d) runs the same chain with the child as reference, starting at `(D, 3^D − 1)`.
  * The coin is `β_n = T^n(x) mod 2`.
  * The integer state is `(k, a)` with `y_n = 3^k x_n + e` and `a = 3^max(0,−k) e`, starting at `(−D, 1 − 3^D)`.
  * Updates, with `σ = a mod 2` and `m = −k`:

    | σ | β | update |
    |---|---|---|
    | 0 | 0 | `a/2` |
    | 0 | 1 | `(3a + 1 − 3^k)/2` if k ≥ 0; `(3a + 3^m − 1)/2` if k < 0 |
    | 1 | 0 | `(3a + 1)/2` if k ≥ 0, `(a + 3^(m−1))/2` if k < 0; then `k + 1` |
    | 1 | 1 | `(a − 3^(k−1))/2` if k ≥ 1, `(3a − 1)/2` if k ≤ 0; then `k − 1` |

* `s(K) = σ_T(x_K)` is the Terras time from `x_K` to 1 (∞ if never). When it is finite, `o(K)` is the number of odd steps of `M_K`'s orbit to 1.
* A negative "within W" means no absorption at any Terras time ≤ W − 9, or ≤ W − 8 where the final state is also checked.

## Statements

1. **Soundness (PROVED).**
   * If the chain is at `(0, 0)` at time n, then `T^n(x) = T^n(y)`. Hence `M_E` and `M_(E−D)` have a common future, and `M_E` reaches 1 iff `M_(E−D)` does.
   * The coins up to time n depend only on `x mod 2^(n+1)`, i.e. on `3^(E−1) mod 2^n`.
   * Trivial-cycle value coincidences are not claimed.
2. **No-meeting lemma (PROVED).** Suppose `E ≥ (1 + log₃4)·N + D + 2`.
   * If the chain is not absorbed at any time ≤ N, then `T^m(x) ≠ T^(m')(y)` for all `m, m' ≤ N`.
   * If it is first absorbed at `n ≤ N`, every meeting with `m, m' ≤ N` has `m = m' ≥ n`.
   * Proof: compare the coefficients of x in `(3^A x + C)/2^m = (3^B y + C')/2^(m')`; the constants are below `3^(N+D+1) 4^N < x`.
   * Every negative below satisfies the hypothesis.
3. **Fast parity vectors (PROVED).**
   * For x known mod `2^n`, `par(x, n) = (V, a, c)` with `T^n(x') = (3^a x' + c)/2^n` for all `x' ≡ x mod 2^n`.
   * Halves combine by `x_mid = (3^(a1) x + c1) >> n1`, `a = a1 + a2`, `c = 3^(a2) c1 + 2^(n1) c2`. The cost is `O(M(n) log n)`.
   * This is the recursion of Elsenhans (arXiv:2502.16743), with the affine form of Terras, Everett and Lagarias.
4. **The supplied residual (FINITE-EXACT).**
   * `E = 99708993705` (t = 23 of `E(t) = 924745897 + 2^32 t`) has exits exactly `D ∈ {3, 4, 7, …, 28}` (W = 7200, D ≤ 40), all absorbed at Terras time 6553. This is codex-grounding's set of 24 children.
   * The original parent `M_99708993713` (W = 8000, D ≤ 64) exits:
     * at D = 1–4 at time 19;
     * at D = 7, 8 at 35 (D = 8 is the source);
     * at D = 5, 6, 11, 12, 15–36 at 6553. These include the 24 children, and D = 36 is the fan bottom `M_99708993677`.
5. **The fan is closed (FINITE-EXACT).**
   * Within W = 16384 and D ≤ 128, the exits of each child are exactly the fan members below it.
   * The bottom `M_99708993677` has none.
6. **The escape (FINITE-EXACT).**
   * For the bottom, no chain with D ≤ 256 is absorbed before Terras time **19,000,765**.
   * Run with collapse, all chains D ≤ 4000 behave as follows:
     * the group of D = 1 reaches its final membership `{1, …, 1910} ∪ {1929, 1930, 1931, 1932, 1935, 1936}` (1916 members) at step 2,741,933;
     * that group is absorbed at Terras time 19,000,765.
   * Hence `M_99708993677 ⇝ M_(99708993677 − D)` for these 1916 values; the deepest is `M_99708991741`.
   * The parent, the seven exponents between it and the source, the source, the 24 children and these 1916 exponents (**1949 exponents**) have a certified common future.
   * Cross-check:
     * `T^n(x) ≡ T^n(y_D) mod 2^(2^20)` at n = 19,000,765 for D = 1 and 256;
     * the odd-step counts are 9500295, 9500296, 9500551;
     * the residues differ at n − 1.
7. **Level 2 (FINITE-EXACT negative).**
   * The other 2084 chains with D ≤ 4000 form two interleaved clusters (D ranges 1911–2902 and 2895–4000), which merge at step 4,474,989.
   * The merged cluster is not absorbed at any Terras time ≤ `2^29 − 9`.
8. **Stopping-time invariant (PROVED).** Let `1 ≤ D ≤ E − 2`.
   * (a) An absorption with merge value ≥ 3 gives `s(E) = s(E − D)` and, when finite, `o(E) = o(E − D)`.
   * (b) Conversely, equal finite invariants force absorption by time `s(E) − 3`, at merge value ≥ 8. HYP-9242's corrected status (mac-mini, audit C) states the same.
   * (c) An absorption at a value ≤ 2 occurs at value 2, with one orbit already at 1 and the other one step from 1. It happens iff `s(E) ≠ s(E − D)` and `2o − s` agree; no pair with K ≤ 12800 qualifies.
   * **Corollaries.**
     * Among exponents whose orbits reach 1, deletion clusters are exactly the level sets of `(σ_T(M_K) − K, o(K))`. This is HYP-9242's partner criterion, for every such K. Part (a) also holds when the orbits never reach 1, which the next corollary uses.
     * Since `s(E) > (E − 1)·log₂3` and `max_(K ≤ 12800) s(K) = 97982`, no route of certificates with merge values ≥ 3 from any `E > 61821` reaches any `K ≤ 12800`.
     * A route reaching a verified exponent from such an E must contain a complete orbit computation.
     * Along a deletion route the time offset is exactly the deleted D: one Terras step per bit. An orbit needs at least `1 + log₂3` steps per bit, typically about 8.64.
     * One orbit computation for any member of the class of statement 6 grounds all 1949 exponents.
   * Checked: all 919 absorptions with 4 ≤ E ≤ 160 (the biconditional over all 12,560 pairs), agreement of `s(K)` with mac-mini's exact table for K ≤ 2000, and the absence of cycle-absorption pairs for K ≤ 12800.
9. **The unpaid progression (FINITE-EXACT).**
   * For `t ≡ 13847 mod 27648` (source `M_E(t)` fixed), 1999 of the first 2000 members have certified exits: 1970 within W ≤ 16384 (D ≤ 64), and 29 within `2^22` (D ≤ 256).
   * The member `t = 3525143` has no exit with D ≤ 256 within `2^26` (no absorption at any time ≤ `2^26 − 8`).

## Not claimed

* Grounding of any of the 24 children. Their certified class is ungrounded, and by statement 8 deletion routes cannot ground it.
* Any density gain: each pointwise certificate covers one 2-adic cell of mass about `2^(−cost+34)`.
* That the level-2 cluster is never absorbed. Chains with D > 4000 were not run.
* Universal Collatz.
