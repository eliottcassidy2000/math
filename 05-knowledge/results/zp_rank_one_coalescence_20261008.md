# Z_7, Z_11, Z_13: rank-one coalescence for the base-p Collatz maps is proved, and it cannot by itself finish Collatz

Session opus-2026-10-08-S22.
* **Owner directive (verbatim):** "keep aiming at remaining collatz steps, think about Z_7, Z_11 and Z_13".
* **Reading.** The repo's live Collatz frontier is mac-mini's coalescence programme (HYP-9244; THM-4606–4609). It studies generalized (Matthews–Watts) Collatz maps on the `d`-adic integers `Z_d`.
  * `Z_2` is complete. `Z_3` and `Z_5` were tested numerically.
  * Its first open "least requirement" is rank one for `d ≥ 3` (`debt_rank_least_requirements_20261008.md` §4: "Not yet written").
  * `Z_7`, `Z_11`, `Z_13` are the next prime bases. On each, the canonical rank-one map is the **base-p Collatz map**

        C_p(x) = x/p if p | x,   C_p(x) = ⌈(p+1)x/p⌉ otherwise,

    which is `3x + 1` (in Terras form) for `p = 2`.
  * Other readings of the triple are in §6.
* **Scripts.** All are in `04-computation/experiments/zp_rank_one_20261008/`, each with a `.out`, each printing ALL CHECKS PASSED.

| script | contents | runtime |
|---|---|---|
| `zp_chain_checks.py` | (A) table vs rationals and direct orbits, (B) local lemmas, (C) descent, (D) level weight, (E) cycles, (F) dictionary | 6 s |
| `zp_basins_survey.py` | (S1) cycles and basin shares for `p ≤ 31`, (S2) basin boundaries to `10^7`, (S3) merges of consecutive integers | 40 s |
| `zp_repunit_lines.py` | repunit lines `p = 7, 11, 13`, `K ≤ 8000`: (R1) cycles, (R2) rivers, (R3) horizon, (R4) orphans; data `repunit_orbits_p{7,11,13}_K8000.txt` | 18 min |
| `zp_tails.py` | merge-time tails for `p = 2, 3, 7, 11, 13, 17, 23` | 3 min |

* **Canon:** THM-4610 (PROVED), HYP-9245 (OPEN).

**Status.**
* PROVED (THM-4610): rank-one Haar coalescence for `C_p`, every odd prime `p`, together with:
  * density-one integer forms;
  * repunit exits for almost every exponent;
  * a `T^(−1/2)` lower bound;
  * density-zero basin boundaries;
  * the no-go 5(f).
* FINITE-EXACT:
  * the local lemmas on boxes;
  * the descent path for `|e| ≤ 20000`;
  * the cycles of `C_p`;
  * rivers on repunit lines to `K = 8000`;
  * the dictionary facts of §6.
* NUMERICAL: basin shares and boundaries, merge tails, repunit orphan fractions.
* HEURISTIC: HYP-9245.
* DICTIONARY: §6.
* Collatz OPEN.

## 0. Results

1. **The rank-one row of HYP-9244 is proved for the base-p Collatz maps (THM-4610).**
   * For every odd prime `p`, two Haar `p`-adic orbits related by `u = (p+1)^k v + e` merge at equal time almost surely, for every admissible start. This includes `y` and `y + e` for every nonzero integer `e`.
   * The proof follows THM-4581 (`p = 2`), with three changes:
     * At every step, the normalized offset is multiplied by the branch factor (`1/p` or `(p+1)/p`) of a fresh uniform digit. The debt is an exactly fair but *lazy* random walk: it flips with probability `2/p`.
     * The step "no more flips ⟹ equal digits", true for `p = 2`, fails for `p ≥ 3`. It is replaced by Lévy's Borel–Cantelli lemma plus injectivity of the digit map.
     * Accessibility comes from an explicit descent: depart, divide, return. It shrinks `|e|` by a factor of about `(p+1)/p²`.
   * Consequences:
     * density-one integer merging;
     * almost-sure repunit exits `R_E ⇝ R_(E−D)`;
     * `P(no merge by T) ≥ c T^(−1/2)`;
     * basin boundaries of density 0.
2. **The no-go: coalescence cannot carry Collatz from "almost all" to "all" (THM-4610 5(f)).**
   * Every structural input of the proof also holds for `p = 3` and `p = 11`: rank one, contraction, the fair skeleton and accessibility. Yet `C_3` and `C_11` have a second cycle on the positive integers.
   * The second basin is not a sliver:
     * the cycle at 7 (length 9) for `p = 3` captures **96.7%** of `n ≤ 10^7`;
     * the cycle at 642 (length 57) for `p = 11` captures **1.8%**.
   * Up to `10^6` the same holds for `p = 17, 23, 29, 31`. Only the trivial cycle appears for `p = 5, 7, 13, 19`.
   * So **of the owner's three bases, `Z_7` and `Z_13` behave like `3x + 1` and `Z_11` does not.**
   * The remaining Collatz step is arithmetic. It must see the cycle and divergence structure of `3x + 1` itself, which no coalescence or rank argument can.
   * THM-4590 (mac-mini) already shows the `p = 2` face of this on the *negative* integers (three cycles). The base-p maps show it on the positive integers, within the same family as `3x + 1`.
3. **Basins (NUMERICAL; HYP-9245 (a)).**
   * THM-4610 proves that basin boundaries have density 0. The convergence runs on a `ln N` clock and is slow.
   * At `N = 10^7`, the boundary density is the fraction 0.577, 0.883, 0.948, 0.959 (`p = 3, 11, 17, 23`) of `2s(1−s)`, its value for independent landing.
   * These fractions track the Haar no-merge probability at the typical orbit length `T_N = ln N/|Λ_p|`: 0.505, 0.842, 0.902, 0.923.
4. **Repunit lines** (`R_K = p^K − 1`, the base-p analogue of the S21 Mersenne line; `K ≤ 8000`).
   * Every `R_K` ends in a known cycle. For `p = 7, 13` all end in the trivial cycle. For `p = 11`, 50 end in the 642-cycle, in just 5 partner classes.
   * Partner classes are rivers: level sets of the multiplication count (one cycle-entry exception, `p = 11`, `K = 3, 4`).
   * The residual horizon is `ln(p+1)/|Λ_p|` per digit, to 0.2%.
   * The orphan fraction settles near `0.56 K^(−1/2)` for `p = 7`. For `p = 11, 13` the exponent is not yet resolved at `K = 8000` (HYP-9245 (b)).
5. **Merge tails.** `√T q(T)` reaches 10.87 (`p = 2`) and plateaus at 9.75–9.77 (`p = 3`) by `T = 10^5`. For `p ≥ 7` it is still rising (25.8, 40.3, 49.9 for `p = 7, 11, 13`), as a lazier skeleton predicts.

## 1. The maps

* `C_p(x) = x/p` on `pZ_p`, and `(Px + p − δ(x))/p` otherwise, where `P = p + 1` and `δ(x) = x mod p`.
* **Place in the literature.** These are Carnielli's `T_d` (arXiv:0810.5169) and the multiplier-`(p+1)` members of Hasse's class and of Möller's family.
  * Möller conjectured eventual periodicity of every orbit for multipliers `m < d^(d/(d−1))`, a condition `m = p + 1` meets.
  * arXiv:2111.06170 extends Tao's almost-bounded-orbits theorem to a class containing them.
* **Structure.**
  * Translation-only: all multipliers are `≡ 1 mod p`.
  * Two-valued: multipliers 1 and `P`.
  * Contracting: `Λ_p < 0`.
  * Rank one: the debt lives in `⟨P⟩`.
  * Among two-valued translation-only maps with `m_0 = 1`, `P = p + 1` is the only contracting multiplier for `p ≥ 3`.
* **The pair chain** (THM-4610's table) tracks `u_n = P^(k_n) v_n + e_n` exactly, with integer state `A = P^max(0,−k) e`.
  * The digit of `u` is `(β + A) mod p`.
  * (A) checks the table against exact rationals on 300,000 random steps, and against direct big-integer orbits for 1200 starts (324 merges). The relation holds at every step, and no equality occurs before absorption.

## 2. THM-4610 in brief

* **Local lemmas (statement 1).**
  * Factor law `|f′| ≤ F|f| + (p−1)/p`, with `F` the branch factor of the reference digit.
  * Direction: toward 0 costs `P/p`; away pays `1/p`.
  * Flip law: `1/p`, `1/p`, `(p−2)/p` when `p ∤ A`; no flips when `p | A`.
  * Exact runs at `k = 0`.
  * Valuation runs at level `h` controlled by `μ_h = 1 + v_p(h)`, with exactly `p − 1` digits leaving.
  * (B) checks all of these exhaustively for `p = 2, 3, 5, 7, 11, 13`. The worst additive term is `(p−1)/p` exactly.
* **Drift (statement 2).**
  * `κ(θ) = (1/p)p^(−θ) + ((p−1)/p)(P/p)^θ` is strictly convex with `κ(0) = κ(1) = 1`. So any `θ ∈ (0, 1)` works.
  * The level weight `s < 1` solves `g(s) ≤ 1`.
  * (D): at `θ = 1/2`, `(κ, least s) = (0.9659, 0.8966)`, `(0.9623, 0.8539)`, `(0.9658, 0.8108)`, `(0.9703, 0.7878)`, `(0.9769, 0.7624)`, `(0.9793, 0.7544)` for `p = 2, 3, 5, 7, 11, 13`. `p = 2` reproduces THM-4581's `s = 0.8966`.
* **Theorem (statement 3)** follows the THM-4581 template:
  * recurrence, with the Borel–Cantelli fix;
  * tightness at returns;
  * accessibility;
  * Lévy's 0–1 law.
* **Descent (statement 4).**
  * (C): for `2 ≤ |e| ≤ 20000` the path strictly shrinks `|e|`, with worst ratio 0.445, 0.241, 0.164, 0.100, 0.084 for `p = 3, 5, 7, 11, 13`.
  * The two-step paths settle `e = ±1`.
  * Every `|e| ≤ 3000` descends to 0.
* **What is not claimed.**
  * The `T^(−1/2)` upper rate.
  * Any statement about particular integers.
  * The general two-valued translation-only case. The proof transfers given accessibility; THM-4610's last section sketches this.

## 3. The no-go and the basins

* **(S1) Cycles and basin shares** (starts `≤ 10^6`):

| p | cycles (minimum: length, basin share) |
|---|---|
| 3 | 1: 3, 0.0332; **7: 9, 0.9668** |
| 5 | 1: 5, 1.0 |
| **7** | **1: 7, 1.0** |
| **11** | **1: 11, 0.9820; 642: 57, 0.0180** |
| **13** | **1: 13, 1.0** |
| 17 | 1: 17, 0.7066; 79: 49, 0.2934 |
| 19 | 1: 19, 1.0 |
| 23 | 1: 23, 0.7664; 82: 72, 0.2336 |
| 29 | 1: 29, 0.9084; 111: 97, 0.0916 |
| 31 | 1: 31, 0.9722; 389: 108, 0.0278 |

* **Logic of the no-go.** THM-4610 holds for every odd prime. Its proof uses rank one, contraction, the fair skeleton and accessibility, and nothing about cycles. An argument from these inputs alone to "every positive integer reaches 1" would prove that `7` reaches 1 under `C_3`, which is false. This is PROVED given the explicit cycle (FINITE-EXACT).
* **What it means for `3x + 1`.**
  * Any proof of the Collatz conjecture must use input beyond HYP-9244's least requirements: input that distinguishes `p = 2` (no second cycle known) from `p = 3, 11, 17, …`.
  * The cycle half (Diophantine; finite checks per period, the S15 "LRC type") and the divergence half are exactly what coalescence does not see.
  * In the reading of the 2026-10-04 atlas's "two walls" (ANALOGY), coalescence sits on the DRIFT side and cannot cross the INTEGRAL wall.
* **(S2) Boundaries** (`n ≤ 10^7`):

| p | second-basin share | boundary density at `10^3 … 10^7` | `q_eff = boundary/(2s(1−s))` | Haar `q_p(T_N)` (`zp_tails`) |
|---|---|---|---|---|
| 3 | 0.967 | 0.134, 0.069, 0.047, 0.040, 0.0366 | 0.577 | 0.505 (`T_N = 92`) |
| 11 | 0.018 | 0.026, 0.035, 0.032, 0.032, 0.0315 | 0.883 | 0.842 (116) |
| 17 | 0.292 | 0.402, 0.425, 0.408, 0.401, 0.3925 | 0.948 | 0.902 (143) |
| 23 | 0.234 | 0.479, 0.360, 0.347, 0.348, 0.3439 | 0.959 | 0.923 (169) |

  * The boundary has density 0 (THM-4610 5(e)). At these sizes, though, it is still most of `2s(1−s)` for `p ≥ 11`: integers below `10^7` take only 92–169 steps to descend, and the lazy skeleton rarely merges a pair that fast.
  * The Haar tail at `T_N` predicts the size of `q_eff` well, slightly underestimating it. That is the expected direction, since smaller `n` have shorter orbits.
  * For `p = 3`, boundary × `√ln N` settles at 0.147–0.150 from `10^6` to `10^7`.
* **(S3)** For `p = 3`, the share of consecutive pairs `n ≤ 10^6` merging at equal time above the cycles is 0.459, against the Haar prediction 0.498 at the typical orbit length `T = 79` (40,000 Haar paths).
* **Relation to THM-4590.**
  * THM-4590 proves `o(x)` cuts for `3x + 1` and shows the negative-integer basins of −1, −5, −17. That is the `p = 2` instance of the same phenomenon, with three cycles.
  * Here the base-p family puts it on the positive integers. In `C_3` the "second" basin is the majority.

## 4. Repunit lines (the S21 analogue)

* `R_K = p^K − 1` (all digits `p − 1`) rises for `K − 1` steps to `x_K = p P^(K−1) − 1`, like `M_K` to `2·3^(K−1) − 1`.
* The deletion children `R_(K−D)` give admissible starts `(−D, 1 − P^D)`, absorbed almost surely (THM-4610 5(b)).
* **(R1) Terminal cycles.**
  * `p = 7` and `p = 13`: all 8000 end in the trivial cycle.
  * `p = 11`: 7950 end in the trivial cycle and 50 in the 642-cycle. Those 50 form only 5 partner classes (bottom, size): (16, 3), (109, 16), (191, 5), (1013, 5), (3605, 21).
  * Coalescence groups the exceptional repunits together.
* **(R2) Rivers.**
  * Partner classes (equal height `n − K`, multiplication count, entry and cycle) are level sets of the multiplication count within a cycle. This is the base-p form of mac-mini's identity that Mersenne rivers are level sets of the odd count.
  * The one exception in 24,000 is `p = 11`, `K = 3, 4`. These have the same height and count but enter the trivial cycle at 9 and 10: a coincidence on the cycle, not a merge above it.
* **(R3) Horizon.** The mean residual steps per digit for `K ≥ 2000` are 12.690, 17.889, 20.452, against `ln(p+1)/|Λ_p|` = 12.716, 17.891, 20.474.
* **(R4) Orphans** (orphan fraction × `√K` per dyadic window from [800, 1600) to [6400, 8000]):
  * `p = 7`: 0.589, 0.565, 0.547, 0.581;
  * `p = 11`: 0.799, 0.773, 0.673, 0.581;
  * `p = 13`: 0.883, 0.684, 0.652, 0.634.
  * Base 7 supports the exponent 1/2. Bases 11 and 13 are still drifting (local exponents 0.55–0.87).
  * HYP-9245 (b) states the universal law and how to refute it. For `p = 2` the same law is HYP-9242 (fitted 0.59, CI [0.42, 0.75]).

## 5. Merge tails (NUMERICAL; `zp_tails.py`)

The start `(0, 1)`, 4800 paths per base:

| p | `√T q(T)` at T = 10, 10², 10³, 10⁴, 10⁵ | `q(10^5)` |
|---|---|---|
| 2 | 2.05, 4.71, 8.31, 9.94, 10.87 | 0.0344 |
| 3 | 2.15, 4.95, 8.22, 9.77, 9.75 | 0.0308 |
| 7 | 2.85, 7.55, 16.28, 23.10, 25.76 | 0.0815 |
| 11 | 3.02, 8.51, 20.45, 34.35, 40.25 | 0.1273 |
| 13 | 3.05, 8.74, 21.88, 39.08, 49.87 | 0.1577 |
| 17 | 3.09, 9.14, 24.24, 49.00, 66.28 | 0.2096 |
| 23 | 3.12, 9.43, 26.25, 58.98, 89.14 | 0.2819 |

* `p = 2` approaches THM-4581's measured `√T q ≈ 11`, and `p = 3` has plateaued.
* For larger `p` the skeleton flips only at rate `2/p`, and the offset contracts at level 0 only through digit-0 steps. So the `T^(−1/2)` regime arrives later, with a larger constant.
* The lower bound `q ≥ c T^(−1/2)` is proved (THM-4610 5(c)). The matching upper bound is not claimed.

## 6. Other readings of Z_7, Z_11, Z_13 (DICTIONARY; checked in (F))

* **Period-4 cycles of `3x + 1`.**
  * `{7, 11, 13} = {|16 − 3^l| : l = 2, 3, 1}`.
  * The 16 fixed points of `T^4` (Terras) are `c_w/(16 − 3^l)`, one per class mod 16. Besides 0, −1 and the trivial cycle they form three non-integral 4-cycles: `{1, 8, 4, 2}/13`, `{5, 11, 20, 10}/7` and `−{19, 23, 29, 38}/11`.
  * The trivial cycle `{1, 2}` is the *integral* member of the `l = 2` (denominator-7) family. This is the 2-adic face of S17's "the trivial cycle's code is `{1,2,4}/7`".
  * In S6/S8 terms the 13- and 7-cycles are first-descent centres (contracting by 3/16 and 9/16 per 4 steps), and the −11-cycle is a rising centre (27/16).
* **Collatz multipliers acting on `Z_p`.**
  * `⟨2⟩ = QR_7`, `⟨3⟩ = QR_11` and `⟨3⟩ ∪ 2⟨3⟩` on `Z_13` are connection sets of circulant tournaments: the Paley tournaments `P_7`, `P_11`, and a cubic-residue circulant on 13 vertices.
  * Their affine automorphism groups have orders 21, 55, 39: `Z_7 ⋊ ⟨2⟩ = F_21` (S18's Prop. 3), `Z_11 ⋊ ⟨3⟩ = F_55` (cf. THM-4592's `BS(1,3) mod 11`), and `Z_13 ⋊ ⟨3⟩ = F_39`.
  * Since `3 ≡ 2^4 mod 13`, both Terras branches act on `Z_13` modulo `⟨3⟩` as the same rotation `2^(−1)`.
* **Further remarks.** `7 · 11 · 13 = 1001`. `7 = 111₂` maps to `13 = 111₃` along the Mersenne line through `11 = T(7)` (THM-4556 (i)).
* None of these readings carries a mechanism for the remaining step; they are recorded as dictionary.

## 7. Connections

* **HYP-9244** (mac-mini). THM-4610 proves its rank-one row for the canonical maps on every odd prime base. The no-go shows that HYP-9244's three requirements (accessibility, contraction, recurrence) characterize Haar coalescence, not integer convergence.
* **THM-4581** (`p = 2`). The proof is its architecture. THM-4610 isolates what was parity-specific: the digit-equality step, and the absence of lazy steps.
* **THM-4590.** Its `o(x)` cuts for `3x + 1`, and the negative-integer basins, are the `p = 2` instances of 5(e) and of §3.
* **THM-4606–4609** (mac-mini): the other rows (expanding maps, `Z_2` classification, rank ≥ 3 and rank 0).
* **HYP-9242 / HYP-9243 / THM-4605** (the Mersenne line). The repunit lines behave the same way: rivers, horizon constants, deletion exits for almost every exponent. Rivers carry the stopping-time invariant. The orphan law is HYP-9245 (b).
* **S17 / S18 / THM-4592.** See §6.

## 8. Status and next steps

| item | status |
|---|---|
| THM-4610 statements 1–5 (rank-one coalescence on `Z_p`; descent; density-one forms; repunit exits; lower rate; boundaries; no-go) | PROVED |
| local lemmas on boxes; descent to `|e| ≤ 20000`; cycles of `C_p`; rivers to `K = 8000`; §6 facts | FINITE-EXACT |
| basin shares and boundaries; merge tails; orphan fractions | NUMERICAL |
| boundary law and repunit orphan law (HYP-9245) | OPEN (HEURISTIC + NUMERICAL) |
| Collatz | OPEN |

**Next steps.**
1. **The general two-valued case.** Prove accessibility for two-valued translation-only maps, or classify its failures by HYP-9244's obstruction lemma. This would close the rank-one row entirely.
2. **The `T^(−1/2)` upper bound for `C_p`.** THM-4581 (4)'s sketch should carry over, with constants growing like `p`.
3. **The arithmetic step.** For `3x + 1`, the input coalescence cannot supply is the cycle and divergence structure. The base-p family (`p = 3, 11` false; `p = 5, 7, 13` apparently true) is a test bench. Any proposed proof mechanism must fail for `C_3` and `C_11`.
4. **HYP-9245.** Run the tests in its refutation section.
