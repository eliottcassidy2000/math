---
id: THM-4601
title: "Two-anchor reduction of the residual first-reset-2 branch: for residual sources n = 2^K t - 1 (x = 2*3^(K-1) t - 1 = 1 + 2^j u, j >= 3, J = floor((j-1)/2) letters 2), (i) the branch is the lag-one pair (x, x-1) with x = 1 mod 8; (ii) +1 barrier: for no deletion child h_D = (n+1)/2^D - 1 with D <= 3000 does the class-decided (equal-time, zero-debt) collision complete before the end of the two-run (the 2-adic limit chains at x* = 1 never absorb); (iii) for J >= 3 the K -> K-3 and K -> K-4 chains end the two-run in the universal state (3, 1-27), X - 1 = 27(Y - 1), independent of K, t, J (clearing identity F_222(x) - 1 = 27(F_411(y) - 1)); (iv) head collisions at +1 (|v| = |u|+3, sum u = sum v + 2, F_v(1) = 4F_u(1)+1) give certificates 1^(K-1) 2^J u c ~ 1^(K-4) (4,1,1) 2^(J-3) v (c+2) uniform in K >= 4, J >= 3, c >= 1, e.g. (1,10)~(1,1,2,2,3) and (10)~(3,1,1,3); (v) the variable-depth rule: six run types, each a fixed pair-chain problem in the post-run bits, Haar coverage -> 1; (vi) on the Mersenne line the certified density in every shell v2(K-1) = m >= 4 is one number rho_r(s), and 2^1889-1 ~> 2^1886-1, 2^5249-1 ~> 2^5246-1 share one head pair"
status: >
  PROVED: (i); (ii) for D <= 3000 (exact rational limit chains; the reduction from finite j to the limit is proved); (iii); (iv) (identity of
  affine maps, child validity forced by the odd endpoint); (v) (the reduction; Haar coverage -> 1 is THM-4581); (vi) the shell law.
  FINITE-EXACT: the two explicit head pairs and 240 random sources with K <= 1001, J <= 150; the two Mersenne examples; the universal state on
  640 sources. NUMERICAL: Haar coverage from the universal state (0.18 at depth 100, 0.53 at 1000, 0.72 at 4000) and the J-table (reset child
  0.45 -> 0.003 as J = 1 -> 32, children 3/4 flat at 0.35-0.42 within 400 steps); Mersenne shells 0.26-0.32 at depth 200. CONJECTURE: (ii) for
  all D is HYP-9240. NOT CLAIMED: coverage of every positive integer (orphans exist, THM-4556 (vi)).
session: mac-mini-2026-10-07-twoanchor
source: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md
scripts:
  - 04-computation/experiments/twoanchor_20261007/twoanchor_core.py (parts B-E; ALL CHECKS PASSED, 119914 assertions)
  - 04-computation/experiments/twoanchor_20261007/child_compare.py (+ .out), mersenne_shells.py
related:
  - THM-4600 (run transparency at the anchors)
  - THM-4555, THM-4556 (collisions at -1; the Mersenne debt chain, orphans in (vi))
  - THM-4581 (Haar absorption of every pair chain)
  - THM-4594 (the sign barrier; (ii) is its positive-cycle twin)
  - 05-knowledge/results/collatz_terminal_lifts_20261007.md (the residual branch), collatz_reset2_rules_20261007.md (the compiler at -1)
  - HYP-9240
---

# THM-4601 — two-anchor reduction of the residual first-reset-2 branch

## Setting

* `U(x) = oddpart(3x+1)`, `T` is the Terras map, and `F_w` is the affine map of a word w.
* A **residual source** is `n = 2^K t − 1` with t odd, K ≥ 2, and first reset letter 2, i.e. `v2(3^K t − 1) = 1`.
* Its word begins with `1^(K−1)`, reaching `x = 2·3^(K−1) t − 1`. Write `x = 1 + 2^j u` with u odd, so that first reset 2 ⟺ j ≥ 3.
* The source then has `J = ⌊(j−1)/2⌋` letters 2, followed by letter 1 (`r = j − 2J = 1`) or a letter ≥ 3 (r = 2).
* On the Mersenne line (t = 1), `j = 3 + v2(K−1)`.
* **Deletion children:** `h_D = (n+1)/2^D − 1` (`1 ≤ D ≤ K−1`), with run ends `y_D = (x+1)/3^D − 1`.
* The D-chain is the pair chain of THM-4581 for `(x, y_D)`: state `(D, 3^D − 1)`, Terras shift D. An absorption at equal time is a certificate `n ⇝ h_D < n`.
* Every K-uniform collision rule of the tiling compiler (`collatz_reset2_rules_20261007.md`, §2) is of this form.

## Statements

**(i) Lag-one form.** `x − 1 = 3z + 1`, where `z = 2·3^(K−2) t − 1` is the last ones-run point of h_1.
* So the reset rule is the immediate merge of the translation pair (x, x−1) when x ≡ 5 mod 8 (`U(4w+1) = U(w)`).
* The residual branch is that pair with x ≡ 1 mod 8.

**(ii) The +1 barrier.** For `s < j`, the actual parities of x's and y_D's Terras orbits are those of the 2-adic limits `x* = 1` and `y_D* = (2 − 3^D)/3^D`. So the D-chain equals the *limit chain* up to time j. For every `D ≤ 3000` the limit chain never absorbs:
* D = 1, 2: `y_D*` tends to 0, and the debt grows by one every two steps.
* D ≥ 3: `y_D*` makes exactly D denominator-clearing odd steps, lands on a positive integer `N_D ≤ 880`, and enters the trivial cycle with debt `k_∞(D) ≥ 3`:
  * 2288 values in phase, 710 out of phase;
  * `min k_∞ = 3` at D = 3, 4;
  * `k_∞(D)/D → 1`.
* Hence for D ≤ 3000 no D-chain is absorbed before Terras time j + 1 after x, i.e. before the end of the two-run. No class-decided collision of n with h_D, which is necessarily equal-time with zero debt, can complete earlier. A sporadic orbit coincidence for one particular n is not excluded and is irrelevant for rules.

**(iii) Universal residual state.** If J ≥ 3 (j ≥ 7), then at the end of the two-run (Terras time 2J after x) both the D = 3 and D = 4 chains are in state `(3, 1 − 27)`:

    X − 1 = 27 (Y − 1),   X = T^(2J)(x),  Y = T^(2J)(y_3) = T^(2J)(y_4),

whatever K, t and J are.
* The proof uses the clearing identity `F_(2,2,2)(x) − 1 = 27(F_(4,1,1)(y_3) − 1)`, both sides equal to `(x+63)/64` for the child.
* y_3 has letters (4,1,1) iff x ≡ 1 mod 128; y_4 reaches the same value with letters (2,2,1,1).
* THM-4600 (a = 2) then carries the relation through the remaining J − 3 letters 2.
* Likewise `(6, 1 − 3^6)` for D = 5, 6 once J ≥ 6.

**(iv) Two-anchor compiler.** Let u, v be words with

    |v| = |u| + 3,   Σu = Σv + 2,   F_v(1) = 4 F_u(1) + 1.

Then for all K ≥ 4, J ≥ 3, c ≥ 1, as affine maps of n,

    F_(1^(K−4) (4,1,1) 2^(J−3) v (c+2)) ((n+1)/8 − 1)  =  F_(1^(K−1) 2^J u c)(n).

Consequences:
* Every positive odd n whose actual U-word begins `1^(K−1) 2^J u c` merges with `h_3 = (n+1)/8 − 1 < n`. That condition is a congruence for n mod `2^(K−1+2J+Σu+c+1)`.
* The displayed child word is the child's actual word.
* Explicit head pairs, with post-run Terras depths 11 and 9:
  * r = 1: `u = (1,10)`, `v = (1,1,2,2,3)`, `F_u(1) = 7/1024`, `F_v(1) = 263/256`;
  * r = 2: `u = (10)`, `v = (3,1,1,3)`, `F_u(1) = 1/256`, `F_v(1) = 65/64`.

**(v) The variable-depth rule.** Classify residual sources by run type: j = 3, 4, 5, 6, or j ≥ 7 with r ∈ {1, 2}.
* For every type and every D, the D-chain state at the end of the two-run depends only on (type, D):
  * exception: D ≤ 2 with j ≥ 7, where it also depends on J (debt ≈ J);
  * for D ∈ {3, 4} and j ≥ 7 it is the universal state of (iii).
* Certification of n by h_D is therefore the absorption of a fixed pair chain, driven by the post-run parity bits. It is uniform in K and t, and in J for D ≥ 3.
* Under Haar measure every such chain absorbs almost surely (THM-4581).
* So, under Haar, the deletion certificates cover a set of residual sources of measure 1. Coverage within post-run depth s tends to 1 as s → ∞.

**(vi) Mersenne shells.** For odd K with `m = v2(K−1) ≥ 4`, put `J = 1 + ⌊m/2⌋` and `r = 1 + (m mod 2)`.
* The post-run class `w′ = 3^(J−3)(3^(K−1) − 1)/2^(m+2) mod 2^L` runs bijectively over the odd classes as `(K−1)/2^m` does.
* Hence within each shell the natural density of K with `2^K − 1 ⇝ 2^(K−3) − 1` within post-run depth s is `ρ_r(s)`. This is the Haar measure of absorption from (iii), and it is the same for all m of the same parity.
* Examples with the same head pair `(10, c) ~ (3,1,1,3, c+2)`:
  * `2^1889 − 1 ⇝ 2^1886 − 1` (m = 5, J = 3, c = 3);
  * `2^5249 − 1 ⇝ 2^5246 − 1` (m = 7, J = 4, c = 1).

## Proofs

**(i).** Direct: x − 1 = `2·3^(K−1) t − 2 = 3z + 1`.

**(ii).**
* Since `x ≡ 1` and `y_D ≡ y_D* (mod 2^j)`, the Terras iterates agree modulo `2^(j−s)`. So the parity bits agree for s < j.
* The chain state at time s depends only on the initial state and the base's first s bits.
* For each D ≤ 3000 the limit chain is computed exactly in rationals:
  * a clearing phase on `a/3^d`;
  * then an integer orbit to {0} or to {1, 2}.
* In {1, 2}, the in-phase debt is constant and the out-of-phase debt oscillates between k and k+1, so absorption is decided at entry.

**(iii).**
* `F_4(y_3) = (x−17)/144`, `F_(4,1)(y_3) = (x+31)/96`, `F_(4,1,1)(y_3) = (x+63)/64`.
* `F_(2,2,2)(x) = (27x+37)/64`, so `F_(2,2,2)(x) − 1 = 27((x+63)/64 − 1)`.
* The letters are exact iff x ≡ 1 (mod 32, 64, 128) successively.
* For y_4: `F_(2,2,1,1)(y_4) = (x+63)/64` likewise.
* Then apply THM-4600.

**(iv).**
* Ones-runs give `x + 1 = 27(y + 1)`.
* Then (iii), then THM-4600, then the head identity. For anchored X, Y it reads `F_v(Y) − 4F_u(X) = F_v(1) − 4F_u(1) = 1`, from the slope relations.
* Finally `3(4Z+1) + 1 = 4(3Z+1)` gives the final letters c, c+2.
* Validity: if `F_word(m) = E` is an odd integer with m an odd integer, then every intermediate value is an odd integer with exactly the stated valuations.
  * A value of negative 2-adic valuation keeps a negative valuation under `(3·+1)/2^b`.
  * An even integer would make the next value non-integral.

**(v).**
* For s ≤ 2J < j the bits are those of the limits, which are determined by (j, D).
* After the run, the post-run point Y is Haar-distributed as t (or K) varies over the class.

**(vi).**
* `3^(2^m q) − 1 = 2^(m+2)·(odd)`, and `q ↦ (3^(2^m q) − 1)/2^(m+2)` is a 2-adic isometry on odd q.

∎

## Numbers (NUMERICAL)

**From the universal state** (20,000 fair-bit runs), P(absorbed within s) is:

| s | 30 | 100 | 300 | 1000 | 4000 |
|---|---|---|---|---|---|
| P | 0.054 | 0.18 | 0.34 | 0.53 | 0.72 |

* Shortest absorbing patterns have length 9.
* Counts of absorbing patterns for lengths 9–20: 1, 1, 3, 7, 15, 30, 67, 147, 301, 658, 1357, 2783.

**Child comparison** (Haar, 400-bit t, budget 400 post-run steps, 300 sources per row):

| J | reset child (D = 1, 2) | D = 3, 4 | any D ≤ 8 |
|---|---|---|---|
| 1 | 0.45 / 0.55 | | |
| 3 | 0.35 / 0.38 | 0.34–0.42 for all J ≥ 3 | 0.52–0.60 for all J ≥ 2 |
| 16 | 0.08 / 0.10 | | |
| 32 | 0.003 / 0.007 | | |

(reset-child cells: r = 1 / r = 2)

**Mersenne shells** (K < 12000, depth ≤ 200), m = 4..8: 0.283, 0.283, 0.255, 0.319, 0.304.

## Limits

* **Positive integers.** If both orbits reach 1 with unequal odd-step counts, the pair never merges at equal time. So for actual integers a single D succeeds only in the proportions above.
* **Orphans.** THM-4556 (vi): 2% of Mersenne exponents in [10^3, 6·10^3] have no deletion partner, and hence no K-uniform certificate. A complete rule must include certificates of depth ∝ K (descent).
* **(ii) is sharp.** The rule's depth 2J + O(1) is optimal among deletion rules (for D ≤ 3000).
