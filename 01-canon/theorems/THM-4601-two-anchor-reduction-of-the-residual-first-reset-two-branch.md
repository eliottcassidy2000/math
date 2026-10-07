---
id: THM-4601
title: "Two-anchor reduction of the residual first-reset-2 branch: for residual sources n = 2^K t - 1 (x = 2*3^(K-1) t - 1 = 1 + 2^j u, j >= 3, J = floor((j-1)/2) letters 2), (i) the branch is the lag-one pair (x, x-1) with x = 1 mod 8; (ii) +1 barrier: for D <= 8000 no D-chain (deletion child h_D = (n+1)/2^D - 1, or any generalized deletion child 3^a(n+1)/2^b - 1 with b - a = D) is absorbed at any Terras time <= j after x, because the 2-adic limit chains at x* = 1 never absorb; (iii) for J >= 3 the K -> K-3 and K -> K-4 chains end the two-run in the universal state (3, 1-27), X - 1 = 27(Y - 1), independent of K, t, J (clearing identity F_222(x) - 1 = 27(F_411(y) - 1)); (iv) ladder compiler at +1: head pairs with |v| = |u|+3 and either Sum u = Sum v + 2i, F_v(1) = 4^i F_u(1) + (4^i-1)/3 (child ladder) or Sum v = Sum u + 2i, F_u(1) = 4^i F_v(1) + (4^i-1)/3 (source ladder) give certificates uniform in K >= 4 and J >= 3, e.g. (10)~(3,1,1,3) and (1,10)~(1,1,2,2,3) (i = 1); (v) for every D the end-of-run state is a function of (D, J), constant for J >= J_0(D), so certification is a fixed pair-chain problem in the post-run bits, uniform in K, t and J >= J_0(D), with Haar coverage -> 1; (vi) on the Mersenne line the certified density in every shell v2(K-1) = m >= 4 is one number rho_r(s), and 2^1889-1 ~> 2^1886-1, 2^5249-1 ~> 2^5246-1 share one head pair"
status: >
  PROVED: (i); (ii) for D <= 8000 (exact rational limit chains; the finite-j reduction is proved and exact: on 8947 actual chains the states
  agree with the limit through time j and differ at j+1); (iii); (iv) (identity of affine maps; the child's word is forced by the odd endpoint);
  (v) (the reduction; Haar coverage -> 1 is THM-4581); (vi) the shell law. FINITE-EXACT: the explicit head pairs on 240 + 320 random sources
  (K <= 1001, J <= 150); the i = 1 head grammar (6237 pairs, Sum u <= 21); the three Mersenne examples (common odd endpoints of 2983, 8310 and 9706 bits);
  the universal state on 840 sources (K = 4..1500, J = 3..150). NUMERICAL: post-run absorption from the universal state, rho_1(s) =
  0.038, 0.172, 0.269, 0.328, 0.375, 0.519 and rho_2(s) = 0.053, 0.184, 0.278, 0.338, 0.384, 0.529 at s = 30, 100, 200, 300, 400, 1000;
  the J-table (reset child 0.45 -> 0.003 as J = 1 -> 32, D = 3, 4 flat at 0.3-0.42); Haar-weighted coverage of the branch by D <= 8 within
  400 post-run steps 0.646 +- 0.014. CONJECTURE: (ii) for all D is HYP-9240. NOT CLAIMED: coverage of every positive integer (THM-4556 (vi)).
  Independently audited 2026-10-07 (audit A): core mathematics CONFIRMED; corrections applied (post-run depths 12/11, r = 2 class count 32,
  (v) depends on (D, J) until J_0(D), grammar not complete (ladders), +1 -> -1 re-anchoring does occur, scope of the barrier); MISTAKE-586.
session: mac-mini-2026-10-07-twoanchor
source: 05-knowledge/results/twoanchor_reset2_friezes_20261007.md
scripts:
  - 04-computation/experiments/twoanchor_20261007/twoanchor_core.py (parts B-E; ALL CHECKS PASSED), barrier_extend.py (+ .out, D <= 8000)
  - 04-computation/experiments/twoanchor_20261007/child_compare.py, mersenne_shells.py, short_types.py, heads_plus1.py (+ .out)
  - 04-computation/experiments/twoanchor_20261007/audit_A/ (independent audit scripts and outputs, REPORT.md)
related:
  - THM-4600 (run transparency at the anchors)
  - THM-4555, THM-4556 (collisions at -1; the Mersenne debt chain, orphans in (vi))
  - THM-4581 (Haar absorption of every pair chain)
  - THM-4594 (the 2-adic sign barrier; (ii) is its positive-cycle twin for deletion certificates)
  - 05-knowledge/results/collatz_terminal_lifts_20261007.md (the residual branch); collatz_reset2_rules_20261007.md (compiler at -1);
    collatz_collision_dp_20261007.md (finite precursor of (ii): pure-two prefixes of lengths 1..16)
  - HYP-9240
---

# THM-4601 — two-anchor reduction of the residual first-reset-2 branch

## Setting

* `U(x) = oddpart(3x+1)`, `T` is the Terras map, and `F_w` is the affine map of a word w.
* A **residual source** is `n = 2^K t − 1` with t odd, K ≥ 2, and first reset letter 2, i.e. `v2(3^K t − 1) = 1`.
* Its word begins `1^(K−1)`, reaching `x = 2·3^(K−1) t − 1`. Write `x = 1 + 2^j u` with u odd, so that first reset 2 ⟺ j ≥ 3.
* The source then has `J = ⌊(j−1)/2⌋` letters 2 (the two-run, ending at Terras time `2J = j − r`). It is followed by letter 1 (`r = j − 2J = 1`) or by a letter ≥ 3 (r = 2).
* On the Mersenne line (t = 1), `j = 3 + v2(K−1)`.
* **Deletion children:** `h_D = (n+1)/2^D − 1` (`1 ≤ D ≤ K−1`), with run ends `y_D = (x+1)/3^D − 1`.
* The D-chain is the pair chain of THM-4581 for `(x, y_D)`: state `(D, 3^D − 1)`, Terras shift D. Its absorption is a class-decided certificate `n ⇝ h_D < n`.
* More generally, `m = 3^a (n+1)/2^b − 1` with `3^a < 2^b` has run end y with `x + 1 = 3^(b−a)(y + 1)`. Its collisions are absorptions of the (b−a)-chain.
* These generalized deletion children are exactly the K-uniform collision certificates on unrefined 2-adic classes, i.e. with no condition on t mod 3. A 3-adic refinement gives others: if 3 | t then `m = (2n−1)/3 = 2^(K+1)(t/3) − 1 < n` has `U(m) = n`, a depth-0 certificate.

## Statements

**(i) Lag-one form.** `x − 1 = 3z + 1`, where `z = 2·3^(K−2) t − 1` is the last ones-run point of h_1.
* So the reset rule is the immediate merge of the translation pair (x, x−1) when x ≡ 5 mod 8 (`U(4w+1) = U(w)`).
* The residual branch is that pair with x ≡ 1 mod 8.

**(ii) The +1 barrier.**
* For `s ≤ j`, the D-chain's state equals that of the *limit chain* for `x* = 1` and `y_D* = (2 − 3^D)/3^D`, because the parity bits agree for s < j.
* For every `D ≤ 8000` the limit chain never absorbs:
  * D = 1, 2: `y_D*` reaches 0 (at Terras time 1, resp. 3), and the debt then grows by one every two steps.
  * D ≥ 3: `y_D*` makes exactly D denominator-clearing odd steps, lands on a positive integer `N_D` (≤ 880 for D ≤ 3000, ≤ 2527 for D ≤ 8000), and enters the trivial cycle with debt `k_∞(D) ≥ 3`.
    * In phase / out of phase: 2288 / 710 for D ≤ 3000, and 3610 / 1390 for 3001 ≤ D ≤ 8000.
    * `min k_∞ = 3` at D = 3, 4.
    * NUMERICAL: in-phase `k_∞/D ∈ [0.75, 1.47]` for D ≤ 3000, `[0.971, 1.045]` for 2000 ≤ D ≤ 3000, `[0.966, 1.022]` for 3001 ≤ D ≤ 8000.
* Hence for D ≤ 8000 no D-chain is absorbed at any Terras time ≤ j after x, in particular not within the two-run. This holds for deletion and generalized deletion children alike.
* The bound is exact: the actual states always differ from the limit at time j + 1.
* A sporadic orbit coincidence for one particular n is not a class-decided certificate and is not excluded.

**(iii) Universal residual state.** If J ≥ 3 (j ≥ 7), then at the end of the two-run (Terras time 2J after x) both the D = 3 and D = 4 chains are in state `(3, 1 − 27)`:

    X − 1 = 27 (Y − 1),   X = T^(2J)(x),  Y = T^(2J)(y_3) = T^(2J)(y_4),

whatever K, t and J are.
* The proof uses the clearing identity `F_(2,2,2)(x) − 1 = 27(F_(4,1,1)(y_3) − 1)`, both sides equal to `(x+63)/64` for the child.
* y_3 has letters (4,1,1) iff x ≡ 1 mod 128; y_4 reaches the same value with letters (2,2,1,1).
* THM-4600 (a = 2) then carries the relation through the remaining J − 3 letters 2.

**(iv) Ladder compiler at +1.** Let u, v be words with `|v| = |u| + 3`, satisfying one of:
* **(child ladder, index i ≥ 1)** `Σu = Σv + 2i` and `F_v(1) = 4^i F_u(1) + (4^i − 1)/3`;
* **(source ladder, index i ≥ 1)** `Σv = Σu + 2i` and `F_u(1) = 4^i F_v(1) + (4^i − 1)/3`.

Then for all K ≥ 4, J ≥ 3, c ≥ 1 (child ladder), as affine maps of n,

    F_(1^(K−4) (4,1,1) 2^(J−3) v (c+2i)) ((n+1)/8 − 1)  =  F_(1^(K−1) 2^J u c)(n),

and symmetrically with final letters `(c+2i, c)` for source ladders.

Consequences:
* Every positive odd n whose actual U-word begins `1^(K−1) 2^J u c` merges with `h_3 = (n+1)/8 − 1 < n`. That condition is a congruence for n mod `2^(K−1+2J+Σu+c+1)`. The child's displayed word is its actual word.
* Explicit i = 1 head pairs, absorbing at post-run Terras depth `Σu + 1`, i.e. at Terras time j + 11 and j + 9 after x:
  * r = 1: `u = (1,10)`, `v = (1,1,2,2,3)`, `F_u(1) = 7/1024 → F_v(1) = 263/256`, post-run depth 12;
  * r = 2: `u = (10)`, `v = (3,1,1,3)`, `F_u(1) = 1/256 → F_v(1) = 65/64`, post-run depth 11.

  These are the earliest possible absorptions for each run type.
* Head grammar (FINITE-EXACT, heads_plus1.py; audit A):
  * there are 6237 child-ladder i = 1 pairs with `Σu ≤ 21`, and no source head is a prefix of another;
  * a pair's source class has measure `2^(r−Σu)` within its run type;
  * the union is 0.0139 (r = 1) and 0.0295 (r = 2).
* The i = 1 grammar is not complete. Exhaustive enumeration of absorbing post-run patterns also finds:
  * child ladders i ≥ 2, the first being `(1,2,9)~(1,1,1,1,1,1)` with i = 3, depth 13;
  * source ladders, the first being `(12,1,1)~(3,1,1,3,2,6)`, depth 17.

  The exact absorbed mass by post-run depth 22 is 0.016943 (r = 1) and 0.032581 (r = 2); the i = 1 grammar carries 82.0% and 90.4% of it.

**(v) The variable-depth rule.**
* For every D, the D-chain state at the end of the two-run is the limit chain's state at Terras time 2J. It is:
  * a function of (D, J) alone;
  * the same for r = 1 and r = 2;
  * constant for `J ≥ J_0(D) = ⌈s_e(D)/2⌉`, where `s_e(D)` is the limit's entry time into {1, 2}. J_0 = 3, 3, 6, 6, 8, 8, 10, 10, 15, 15 for D = 3..12, and up to 2965 for D ≤ 3000.
* The states (FINITE-EXACT; short_types.py, audit A):

  | D | J = 1 | J = 2 | J = 3 | J = 4 | J = 5 | J ≥ J_0(D) |
  |---|---|---|---|---|---|---|
  | 1, 2 | (D, (3^D−1)/2), anchored at −1/2 | (2, 1), anchored at −1/8 | (3, 1) | (4, 1) | (5, 1) | never constant: (J, 1), debt J |
  | 3, 4 | (D, (3^D−1)/2) | (4, 10) | | | | **(3, −26)** from J = 3 |
  | 5, 6 | (D, (3^D−1)/2) | (6, 91) | (6, −296) | (6, −404) | (6, −485) | (6, −728) from J = 6 |

* r enters only through the forced post-run prefix of the child's bits: (1, 1) if r = 1, (1, 0, 0) if r = 2.
* Certification of n by h_D is therefore the absorption of a fixed pair chain driven by the post-run bits. It is uniform in K, t, and J ≥ J_0(D).
* Under Haar measure every such chain absorbs almost surely (THM-4581). So the deletion certificates cover a set of residual sources of Haar measure 1, and coverage within post-run depth s tends to 1 as s → ∞.

**(vi) Mersenne shells.** For odd K with `m = v2(K−1) ≥ 4`, put `J = 1 + ⌊m/2⌋ ≥ 3` and `r = 1 + (m mod 2)`.
* The post-run class `w′ = 3^(J−3)(3^(K−1) − 1)/2^(m+2) mod 2^L` runs bijectively over the odd classes as `(K−1)/2^m` does.
* Hence within each shell the natural density of K with `2^K − 1 ⇝ 2^(K−3) − 1` within post-run depth s is `ρ_r(s)`, the same for all m of the same parity.
* Examples with the same head pair `(10, c) ~ (3,1,1,3, c+2)`:
  * `2^1889 − 1 ⇝ 2^1886 − 1` (m = 5, J = 3, c = 3);
  * `2^5249 − 1 ⇝ 2^5246 − 1` (m = 7, J = 4, c = 1).
* The r = 1 pair `(1,10, c) ~ (1,1,2,2,3, c+2)` also occurs: `2^6129 − 1 ⇝ 2^6126 − 1` (m = 4, J = 3, c = 1; common odd endpoint of 9706 bits).

## Proofs

**(i).** Direct: x − 1 = `2·3^(K−1) t − 2 = 3z + 1`.

**(ii).**
* Since `x ≡ 1` and `y_D ≡ y_D* (mod 2^j)`, the Terras iterates agree modulo `2^(j−s)`, so the parity bits agree for s < j.
* The chain state at time s depends only on the initial state and the base's first s bits.
* For each D ≤ 8000 the limit chain is computed exactly in rationals:
  * a clearing phase on `a/3^d`;
  * then an integer orbit to {0} or to {1, 2}.
* In {1, 2} the in-phase debt is constant and the out-of-phase debt oscillates between k and k+1, so absorption is decided at entry.

**(iii).**
* `F_4(y_3) = (x−17)/144`, `F_(4,1)(y_3) = (x+31)/96`, `F_(4,1,1)(y_3) = (x+63)/64`.
* `F_(2,2,2)(x) = (27x+37)/64`.
* The letters are exact iff x ≡ 1 (mod 32, 64, 128) successively.
* For y_4: `F_(2,2,1,1)(y_4) = (x+63)/64`.
* Then apply THM-4600.

**(iv).**
* Ones-runs give `x + 1 = 27(y + 1)`.
* Then (iii) and THM-4600.
* For anchored X, Y the slope relations give `F_v(Y) = 4^i F_u(X) + (4^i−1)/3` (child ladder).
* `z ↦ 4z + 1` iterated i times preserves U and adds 2i to the letter. So the endpoints agree with final letters c and c+2i.
* Validity: if `F_word(m) = E` is an odd integer with m an odd integer, then every intermediate value is an odd integer with exactly the stated valuations.
  * A value of negative 2-adic valuation keeps a negative valuation under `(3·+1)/2^b`.
  * An even integer would make the next value non-integral.

**(v).**
* For s ≤ 2J < j the bits are those of the limits, which depend only on D.
* After the run, the post-run point Y is Haar-distributed (given r) as t, or K within a shell, varies over the class.

**(vi).**
* `3^(2^m q) − 1 = 2^(m+2)·(odd)`, and `q ↦ (3^(2^m q) − 1)/2^(m+2)` is a 2-adic isometry on odd q.

∎

## Numbers (NUMERICAL)

**From the universal state with the realisable post-run bits:**

| s | 30 | 100 | 200 | 300 | 400 | 1000 |
|---|---|---|---|---|---|---|
| ρ_1(s) | 0.038 | 0.172 | 0.269 | 0.328 | 0.375 | 0.519 |
| ρ_2(s) | 0.053 | 0.184 | 0.278 | 0.338 | 0.384 | 0.529 |

* With unconstrained fair bits from (3, −26): P(absorbed within s) = 0.18, 0.53, 0.72 at s = 100, 1000, 4000. The shortest such pattern has length 9, but it starts with an even bit and is not realisable after a two-run.

**Child comparison** (Haar, 400-bit t, budget 400 post-run steps, 300 sources per row):

| J | reset child (D = 1, 2) | D = 3, 4 | any D ≤ 8 |
|---|---|---|---|
| 1 | 0.45 / 0.55 | | |
| 3 | 0.35 / 0.38 | 0.3–0.42 for all J ≥ 3 | 0.52–0.60 for J ≥ 2 |
| 16 | 0.08 / 0.10 | | |
| 32 | 0.003 / 0.007 | | |

(reset-child cells: r = 1 / r = 2)

* Haar-weighted coverage of the whole branch by D ≤ 8 within 400 steps: 0.646 ± 0.014 (audit A).
* **Mersenne shells** (K < 12000, depth ≤ 200), m = 4..8: 0.283, 0.283, 0.255, 0.319, 0.304.

## Limits

* **Positive integers.** If both orbits reach 1 with unequal odd-step counts, the pair never merges at equal time. So for actual integers a single D succeeds only in the proportions above.
* **Orphans.** THM-4556 (vi): 2% of the odd exponents (50 of the 2500 odd a in [10^3, 6·10^3]) have no equal-time Mersenne partner. Such sources need other certificates: 3-adic ones (above), non-Mersenne children, or descent of depth ∝ K.
* **Sharpness.** (ii) is sharp up to an additive constant:
  * no generalized deletion certificate with D ≤ 8000 completes by Terras time j after x;
  * the explicit heads complete at j + 11 and j + 9, the earliest possible for D = 3, 4;
  * the depth beyond 2J is unbounded (heavy tail).
