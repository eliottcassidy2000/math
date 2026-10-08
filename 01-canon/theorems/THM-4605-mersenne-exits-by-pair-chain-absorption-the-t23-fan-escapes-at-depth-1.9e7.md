---
id: THM-4605
title: "Mersenne-line exits by pair-chain absorption: an exit M_E ~> M_(E-D) is certified by the absorption of the source-reference pair chain of x = 2*3^(E-1) - 1 (THM-4581's table), which needs only 3^(E-1) mod 2^W, and divide-and-conquer parity vectors reach depth 2^25 in seconds; the supplied t = 23 source M_99708993705 and its 24 children form a closed fan whose bottom M_99708993677 has no exit with D <= 256 before Terras time 19,000,765, where all 256 deletion chains (collapsed into one at step 2,741,933) are absorbed, so the fan has a certified common future with M_99708991767 (D = 1910); the next cluster is not absorbed within 2^27; 1999 of the first 2000 members of the unpaid progression t = 13847 mod 27648 have certified exits"
status: "PROVED: soundness (statement 1), correctness of the fast parity vectors (statement 2). FINITE-EXACT: statements 3-7 (exact integer pair chains; coins from exact residues; the escape independently cross-checked from affine maps mod 2^(2^20)). Session opus-2026-10-07-S21."
session: opus-2026-10-07-S21 (owner directive: codex-grounding's status update; 'the next decisive target is a terminating, source-preserving rule for one of those actual children')
source: 05-knowledge/results/mersenne_line_barriers_20261007.md
scripts:
  - 04-computation/experiments/mersenne_line_barriers_20261007.py (+ .out; ALL CHECKS PASSED)
  - 04-computation/experiments/mersenne_fan_escape_20261007.py (+ .out; ALL CHECKS PASSED)
  - 04-computation/experiments/mersenne_line_coalescence_20261007.py (+ .out; ALL CHECKS PASSED)
related:
  - THM-4581 (mac-mini: the pair chain and its Haar absorption); THM-4601 (i) (the reset rule: even E ~> E-1 at Terras time 3)
  - 05-knowledge/results/collatz_residual_t23_20261007f.md and collatz_grounding_progress_20261007f.md (codex-grounding: the t = 23 exit and the 24 children)
  - 05-knowledge/results/collatz_parameter_complement_synthesis_20261007e.md (the family E(t), the unpaid progression)
  - HYP-9242 (orphan law), HYP-9243 (diffusive coalescence)
---

# THM-4605 — Mersenne exits by pair-chain absorption; the t = 23 fan escapes

## Setting

* `T` is the Terras map, and `M_E = 2^E − 1`. Since `T^j(M_E) = 3^j 2^(E−j) − 1` for `j ≤ E`, after `E − 1` steps the source is `x = 2·3^(E−1) − 1`.
* The child `M_(E−D)` is at `y = 2·3^(E−D−1) − 1`, with `y + 1 = 3^(−D)(x + 1)`.
* **Source-reference pair chain** (THM-4581's table with x as reference):
  * coin `β_n = T^n(x) mod 2`;
  * integer state `(k, a)` with `y_n = 3^k x_n + e`, `a = 3^max(0,−k) e`, starting at `(−D, 1 − 3^D)`.
  * Updates, with `σ = a mod 2` and `m = −k`:

    | σ | β | update |
    |---|---|---|
    | 0 | 0 | `a/2` |
    | 0 | 1 | `(3a + 1 − 3^k)/2` if k ≥ 0; `(3a + 3^m − 1)/2` if k < 0 |
    | 1 | 0 | `(3a + 1)/2` if k ≥ 0, `(a + 3^(m−1))/2` if k < 0; then `k + 1` |
    | 1 | 1 | `(a − 3^(k−1))/2` if k ≥ 1, `(3a − 1)/2` if k ≤ 0; then `k − 1` |

## Statements

1. **Soundness (PROVED).**
   * If the chain is at `(0, 0)` at time n, then `T^n(x) = T^n(y)`. Hence `M_E` and `M_(E−D)` have a common future, and `M_E` reaches 1 iff `M_(E−D)` does.
   * The coins up to time n depend only on `x mod 2^(n+1)`, i.e. on `3^(E−1) mod 2^(n+1)`.
   * Trivial-cycle value coincidences are not claimed.
2. **Fast parity vectors (PROVED).**
   * For x known mod `2^n`, `par(x, n) = (V, a, c)` with `T^n(x') = (3^a x' + c)/2^n` for all `x' ≡ x mod 2^n`.
   * Halves combine by `x_mid = (3^(a1) x + c1) >> n1`, `a = a1 + a2`, `c = 3^(a2) c1 + 2^(n1) c2`. Cost `O(M(n) log n)`.
3. **The supplied residual (FINITE-EXACT).**
   * `E = 99708993705` (t = 23 of `E(t) = 924745897 + 2^32 t`) has exits exactly `D ∈ {3, 4, 7, …, 28}` within W = 7200, D ≤ 40, all absorbed at Terras time 6553. This is codex-grounding's set of 24 children.
   * The original parent `M_99708993713` exits to `M_99708993677` (D = 36) at 6553.
4. **The fan is closed (FINITE-EXACT).** Within W = 16384 and D ≤ 128, every exit of every child lands on another fan member, and the bottom `M_99708993677` has none.
5. **The escape (FINITE-EXACT).**
   * For the bottom, the chains for D = 1…256 collapse into one chain at step 2,741,933. It is the first to absorb, at Terras time **19,000,765**, so no exit with D ≤ 256 occurs earlier.
   * The chains for D = 1500 and 1910 collapse into it before that time; D = 1911 does not.
   * Hence `M_99708993677 ⇝ M_99708991767`. So the source, its 24 children and the parent have a certified common future with `M_99708991767`.
   * Cross-check: `T^n(x) ≡ T^n(y_D) mod 2^(2^20)` at n = 19,000,765 for D = 1 and 256, with odd-step counts differing by D, and the residues differ at n − 1.
6. **Level 2 (FINITE-EXACT negative).** The chains for D = 1911, 1912, 2500 and 4000 of the bottom are not absorbed within `2^27` Terras steps.
7. **The unpaid progression (FINITE-EXACT).**
   * For `t ≡ 13847 mod 27648` (source `M_E(t)` fixed), 1999 of the first 2000 members have certified exits: 1970 within W ≤ 16384 (D ≤ 64), and 29 within `2^22` (D ≤ 256).
   * The member `t = 3525143` has no exit with D ≤ 256 within `2^26`.

## Not claimed

* Grounding of any of the 24 children: their certified descendant `M_99708991767` is itself ungrounded.
* Any density gain. Each pointwise certificate covers one 2-adic cell of mass about `2^(−cost+34)`.
* Universal Collatz.
