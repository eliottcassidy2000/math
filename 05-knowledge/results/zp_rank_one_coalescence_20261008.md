# Z_7, Z_11, Z_13: rank-one coalescence for the base-p Collatz maps is proved, and it cannot by itself finish Collatz

Session opus-2026-10-08-S22.
* **Owner directive (verbatim):** "keep aiming at remaining collatz steps, think about Z_7, Z_11 and Z_13".
* **Reading.** The repo's live Collatz frontier is mac-mini's coalescence programme (HYP-9244; THM-4606–4609). It studies generalized (Matthews–Watts) Collatz maps on the `d`-adic integers `Z_d`.
  * `Z_2` is complete. `Z_3` and `Z_5` were tested numerically.
  * Its first open "least requirement" is rank one for `d ≥ 3` (`debt_rank_least_requirements_20261008.md` §6 item 1; at commit 97b640656 its §4 ended "Not yet written").
  * `Z_7`, `Z_11`, `Z_13` are the next prime bases. On each, the canonical rank-one map is the **base-p Collatz map**

        C_p(x) = x/p if p | x,   C_p(x) = ⌈(p+1)x/p⌉ otherwise,

    which is `3x + 1` (in Terras form) for `p = 2`.
  * Other readings of the triple are in §6.
* **Scripts.** All are in `04-computation/experiments/zp_rank_one_20261008/`, each with a `.out`, each printing ALL CHECKS PASSED.

| script | contents | runtime |
|---|---|---|
| `zp_chain_checks.py` | (A) table vs rationals and direct orbits, (B) local lemmas, (C) descent, (D) level weight, (E) cycles, (F) dictionary | 4 s |
| `zp_basins_survey.py` | (S1) cycles, basin shares and cycle shapes for `p ≤ 31`, compared with Carnielli's Table 1, (S2) basin boundaries to `10^7`, (S3) merges of consecutive integers, (S4) the negative side (the sign-conjugate maps `F_p`), (S5) the 1-cycle equation and the single candidate for `p = 2` | 47 s |
| `zp_boundary_ratio.py` | (added after the audit) boundary density against `2s(1−s) q_p(T_N)` at `N = 10^3 … 10^7`, with `q_p(T_N)` from `10^5` paths | 25 s |
| `zp_repunit_lines.py` | repunit lines `p = 7, 11, 13`, `K ≤ 8000`: (R1) cycles, (R2) rivers, (R3) horizon, (R4) orphans; data `repunit_orbits_p{7,11,13}_K8000.txt` | 18 min |
| `zp_orphan_fit.py` | (added after the audit) orphan exponents with likelihood intervals from the committed repunit data; validity of HYP-9245 (b)'s refutation rule | 8 s |
| `zp_tails.py` | merge-time tails for `p = 2, 3, 7, 11, 13, 17, 23`, with binomial errors | 2 min |

* **Canon:** THM-4610 (PROVED; independently audited, corrections applied, MISTAKE-649), HYP-9245 (OPEN; restated after the audit).

**Status.**
* PROVED (THM-4610): rank-one Haar coalescence for `C_p`, every odd prime `p` (the almost-sure part of HYP-9244's requirement 1 for these maps), together with:
  * density-one integer forms;
  * repunit exits for almost every exponent;
  * a `T^(−1/2)` lower bound;
  * density-zero basin boundaries.
* KNOWN:
  * the barrier of §3 (the `3x − 1` map, and the sign-conjugate `F_p` for every `p`; Akin 2004);
  * the cycles of `C_p` (Carnielli's Table 1, computed by K. R. Matthews).
* FINITE-EXACT:
  * the local lemmas on boxes;
  * the descent path for `|e| ≤ 20000`;
  * the cycles of `C_p`, re-verified with their shapes (positive starts `≤ 10^6`, negative starts `≥ −10^5`);
  * rivers on repunit lines to `K = 8000`;
  * the dictionary facts of §6.
* NUMERICAL: basin shares and boundaries, merge tails, repunit orphan fractions, all with error bars.
* HEURISTIC: HYP-9245.
* DICTIONARY: §6.
* Collatz OPEN.

## 0. Results

1. **The almost-sure part of HYP-9244's requirement 1 (rank one, `d ≥ 3`) is proved for the base-p Collatz maps (THM-4610).**
   * For every odd prime `p`, two Haar `p`-adic orbits related by `u = (p+1)^k v + e` merge at equal time almost surely, for every admissible start. This includes `y` and `y + e` for every nonzero integer `e`.
   * The proof follows THM-4581 (`p = 2`) along the route of HYP-9244's requirement-1 sketch, whose level-weight drift audit G checked (mac-mini), with three changes:
     * At every level `h ≥ 1`, the normalized offset is multiplied by the branch factor (`1/p` or `(p+1)/p`) of a fresh uniform reference digit, up to a bounded additive term. The debt's flip skeleton is an exactly fair simple random walk, but the debt itself is *lazy*: it moves with probability `2/p` when `p ∤ A`.
     * The step "no more flips ⟹ equal digits", true for `p = 2`, fails for `p ≥ 3`. It is replaced by Lévy's Borel–Cantelli lemma plus injectivity of the digit map.
     * Accessibility comes from an explicit descent: depart, divide, return. It shrinks `|e|` by a factor of about `(p+1)/p²`.
   * Consequences:
     * density-one integer merging;
     * almost-sure repunit exits `R_E ⇝ R_(E−D)`;
     * `P(no merge by T) ≥ c T^(−1/2)`;
     * basin boundaries of density 0.
2. **The barrier, seen inside the family (KNOWN; THM-4610's remark after statement 5).**
   * Coalescence cannot carry Collatz from "almost all" to "all", and this is known.
     * Every Haar statement is invariant under `x ↦ −x`, which conjugates `3x + 1` to `3x − 1`. The `3x − 1` map has three known cycles on the positive integers.
     * Akin (2004) makes the same point for the 2-adic conjugacy method: it is "too general", working for every odd multiplier, and `5x + 1` has extra cycles.
     * The same conjugation gives, for every `p`, a map `F_p` with all of THM-4610's properties and `p − 1` fixed points.
   * In the base-p family the unconjugated maps show it too (cycles KNOWN from Carnielli's Table 1, computed by K. R. Matthews; re-verified here). `C_p` has a second positive-integer cycle for `p = 3, 11, 17, 23, 29, 31`, and its basin is not a sliver:
     * the cycle at 7 (length 9) for `p = 3` captures **96.7%** of `n ≤ 10^7`;
     * the cycle at 642 (length 57) for `p = 11` captures **1.8%**.
   * Only the trivial cycle appears on the positive integers for `p = 5, 7, 13, 19` (starts `≤ 10^6`). On the negative integers, `C_p` has extra cycles for `p = 3, 5, 13, 19, 23` and none for `p = 7` (starts `≥ −10^5`, (S4)).
   * **How the owner's three bases compare:**
     * `Z_7` has no extra cycles of either sign. It is the only prime `3 ≤ p ≤ 31` without one in these ranges ((S5)); Carnielli notes `d = 7` as the first base with none.
     * `Z_13` has `3x + 1`'s pattern: only the trivial cycle among the positive integers, and an extra cycle among the negatives (least `|x|` 147, length 36).
     * `Z_11` has an extra positive cycle (642, length 57).
   * The remaining Collatz step is arithmetic. It must see the cycle and divergence structure of `3x + 1` itself, which no coalescence or rank argument can. THM-4590 (mac-mini) shows the `p = 2` face of this on the negative integers.
3. **Basins (NUMERICAL; HYP-9245 (a)).**
   * THM-4610 proves that basin boundaries have density 0. The convergence runs on a `ln N` clock and is slow.
   * At `N = 10^7`, the boundary density is the fraction 0.577, 0.883, 0.948, 0.959 (`p = 3, 11, 17, 23`) of `2s(1−s)`, its value for independent landing.
   * The Haar no-merge probability at the typical orbit length `T_N = ln N/|Λ_p|` is 0.484, 0.839, 0.896, 0.924 (± 0.002, from `10^5` paths). So the ratio is 1.19, 1.05, 1.06, 1.04.
   * Over `N = 10^3 … 10^7` the ratio stays between 1.04 and 1.65 and drifts down at every `p`; for `p = 3` it runs 1.65, 1.61, 1.40, 1.26, 1.19. Whether it tends to 1 is not resolved.
   * This is THM-4590's cut heuristic for `3x + 1`, transplanted to base `p`.
4. **Repunit lines** (`R_K = p^K − 1`, the base-p analogue of the S21 Mersenne line; `K ≤ 8000`).
   * Every `R_K` ends in a known cycle. For `p = 7, 13` all end in the trivial cycle. For `p = 11`, 50 end in the 642-cycle, in just 5 partner classes.
   * Partner classes are rivers: level sets of the multiplication count (one cycle-entry exception, `p = 11`, `K = 3, 4`).
   * The residual horizon is `ln(p+1)/|Λ_p|` per digit, to 0.2%.
   * Orphans (least members of partner classes) decay with fitted exponents 0.57 [0.21, 0.93], 0.69 [0.37, 1.01] and 0.71 [0.39, 1.03] for `p = 7, 11, 13` (`K ∈ [800, 8000]`, 95% likelihood intervals). All three intervals contain 1/2, and none resolves the exponent (HYP-9245 (b)).
5. **Merge tails.**
   * `√T q(T)` at `T = 10^4, 10^5` is 9.94 ± 0.43 and 10.87 ± 0.83 for `p = 2`, and 9.77 ± 0.43 and 9.75 ± 0.79 for `p = 3`. That is consistent with a `T^(−1/2)` plateau but does not establish one (4800 paths).
   * For `p ≥ 7` it is still rising: 25.8 ± 1.3, 40.3 ± 1.5 and 49.9 ± 1.7 for `p = 7, 11, 13` at `10^5`, as a lazier skeleton predicts.

## 1. The maps

* `C_p(x) = x/p` on `pZ_p`, and `(Px + p − δ(x))/p` otherwise, where `P = p + 1` and `δ(x) = x mod p`.
* **Place in the literature.**
  * These are Carnielli's `T_d` (Applied Mathematics E-Notes 15 (2015) 207–215). Its Table 1 lists their cycles for `d ≤ 150` and `|x| ≤ 50000`, computed by K. R. Matthews.
  * They are the multiplier-`(p+1)` members of Hasse's class and of Möller's family.
  * Bruschi (arXiv:0810.5169) studies the undivided form `p C_p`. Tremblay (arXiv:1907.13302) locates their cycles.
  * Möller conjectured eventual periodicity of every orbit for multipliers `m < d^(d/(d−1))`, a condition `m = p + 1` meets.
  * arXiv:2111.06170 (Gonçalves–Greenfeld–Madrid) extends Tao's almost-bounded-orbits theorem to a class containing them. Its Remark 1.6 quotes Matthews–Watts: such maps preserve Haar measure on `Z_p` and are strongly mixing.
* **Structure.**
  * Translation-only: all multipliers are `≡ 1 mod p`.
  * Two-valued: multipliers 1 and `P`.
  * Contracting: `Λ_p < 0`.
  * Rank one: the debt lives in `⟨P⟩`.
  * Among maps with `m_0 = 1` and one common multiplier `m ≡ 1 mod p` on all nonzero residues, `m = p + 1` is the only contracting one. Contraction means `m < p^(p/(p−1))`, and `p^(p/(p−1)) < 2p + 1`.
  * Two-valued maps that keep the multiplier 1 on some nonzero residues can contract with other multipliers, e.g. `Z_3` (1, 1, 7) and `Z_5` (1, 1, 1, 1, 26).
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
  * `κ(θ) = (1/p)p^(−θ) + ((p−1)/p)(P/p)^θ` is strictly convex with `κ(0) = κ(1) = 1`; `κ(1) = 1` because the multipliers sum to `p²`. So any `θ ∈ (0, 1)` works.
  * The level weight `s < 1` solves `g(s) ≤ 1`.
  * (D): at `θ = 1/2`, `(κ, least s) = (0.9659, 0.8966)`, `(0.9623, 0.8539)`, `(0.9658, 0.8108)`, `(0.9703, 0.7878)`, `(0.9769, 0.7624)`, `(0.9793, 0.7544)` for `p = 2, 3, 5, 7, 11, 13`. `p = 2` reproduces THM-4581's `s = 0.8966`.
* **Theorem (statement 3)** follows the THM-4581 template:
  * recurrence, with the Borel–Cantelli fix;
  * tightness at returns;
  * accessibility;
  * Lévy's 0–1 law.
* **Descent (statement 4).**
  * (C): for `2 ≤ |e| ≤ 20000` the path strictly shrinks `|e|`.
    * Over all these `e` the worst ratio is 0.5, at `|e| = 2`, for every `p`.
    * Over `|e| ≥ 1000` the worst ratio is 0.445, 0.241, 0.164, 0.100, 0.084 for `p = 3, 5, 7, 11, 13`.
  * The two-step paths settle `e = ±1`.
  * Every `|e| ≤ 3000` descends to 0.
* **What is not claimed.**
  * The `T^(−1/2)` upper rate.
  * Any statement about particular integers.
  * The general two-valued translation-only case. The proof transfers given accessibility; THM-4610's last section sketches this.

## 3. The barrier and the basins

* **(S1), (S4) Cycles, basin shares and shapes.**
  * Positive starts `≤ 10^6`; negative starts `≥ −10^5`.
  * The cycles agree with Carnielli's Table 1 (KNOWN; computed by K. R. Matthews for `|x| ≤ 50000`).
  * The shape `(k_1, k_2)` counts the division and multiplication steps in one period, and `λ = (1/p)^(k_1) ((p+1)/p)^(k_2)`.

| p | positive cycles (minimum: length, basin share; `(k_1, k_2)`, `λ`) | negative extra cycles (least `|x|`: length; `(k_1, k_2)`, `λ`) |
|---|---|---|
| 3 | 1: 3, 0.0332; **7: 9, 0.9668; (2, 7), 0.832** | −22: 5; (1, 4), 1.054 |
| 5 | 1: 5, 1.0 | −57: 10; (1, 9), 1.032 |
| **7** | **1: 7, 1.0** | **none** |
| **11** | **1: 11, 0.9820; 642: 57, 0.0180; (2, 55), 0.990** | **none** |
| **13** | **1: 13, 1.0** | **−147: 36; (1, 35), 1.029** |
| 17 | 1: 17, 0.7066; 79: 49, 0.2934; (1, 48), 0.914 | none |
| 19 | 1: 19, 1.0 | −954: 117; (2, 115), 1.010 |
| 23 | 1: 23, 0.7664; 82: 72, 0.2336; (1, 71), 0.893 | −800: 75; (1, 74), 1.014 |
| 29 | 1: 29, 0.9084; 111: 97, 0.0916; (1, 96), 0.893 | none |
| 31 | 1: 31, 0.9722; 389: 108, 0.0278; (1, 107), 0.964 | none |

  * The trivial cycle `1, 2, …, p` has shape `(1, p − 1)`. On the negative side, `−1, …, −(p − 1)` are fixed points.
  * **Signs.** Positive cycles have `λ < 1` and negative ones `λ > 1`. The same holds for `3x + 1`, where a positive cycle word needs `3^l < 2^A`.
  * **Run counts.** Every extra cycle has a single run of divisions, which makes it a *1-cycle* in the sense of Steiner and Simons–de Weger. The exception is the 642-cycle of `C_11`, which has two runs (a 2-cycle). (FINITE-EXACT.)
  * **One equation versus many ((S5), FINITE-EXACT; the reading is ANALOGY).**
    * *The equation.* Read a 1-cycle from the start `x` of its multiplication run. Then

          x (p^(k_1+k_2) − (p+1)^(k_2)) = c_w,    c_w = Σ_j (p+1)^(k_2−1−j) p^j (p − i_j),

      with `i_j ∈ {1, …, p − 1}` the digits along the run. All ten extra 1-cycles above satisfy it.
    * *`p = 2`.* The run's digits are all 1, so each shape `(k_1, k_2)` has a single candidate, `c = 3^(k_2) − 2^(k_2)`. Steiner (1977) proved that it is divisible only for the trivial cycle; checked here for `k_2 ≤ 2000`, where on the negative side only the cycles −1 and −5 occur. Steiner's proof, and its extensions to `m`-cycles with `m ≤ 91` (Simons–de Weger 2005; Hercher 2023), rest on linear forms in logarithms, because the near-balance `2^A ≈ 3^l` is rare.
    * *`p ≥ 3`.* Each shape has `(p − 1)^(k_2)` candidate numerators. The near-balance is also cheap: one multiplication moves `ln λ` by only `ln(1 + 1/p) ≈ 1/p`, so `1 − λ ≲ 1/p` is available for every `k_1`.
    * *The upshot.* Extra cycles appear for every prime `3 ≤ p ≤ 31` except 7, on one side or the other.
    * *Reading for `3x + 1`.* The cycle half is a statement about one equation per shape at rare near-balances. That is the regime where Baker-type bounds work, and the base-`p` maps leave it.

* **The barrier (KNOWN).**
  * Every statement of THM-4610 is a Haar statement, invariant under `x ↦ −x`. So it holds for the sign-conjugate `F_p(x) = −C_p(−x)` (`x/p` or `⌊(p+1)x/p⌋` on the positive integers), which fixes `1, …, p − 1`. For `p = 2` this is the `3x − 1` map, with three known positive cycles.
  * An argument from Haar properties alone to "every positive integer reaches 1" would therefore also prove it for `3x − 1`, which is false. This is the classical reason a proof must distinguish `3x + 1` from `3x − 1`. Akin (2004) makes the point for the 2-adic conjugacy, and Kontorovich–Lagarias (§1.1) list the three negative cycles of `3x + 1`.
  * The positive cycles of `C_3`, `C_11`, … show the same thing without conjugating.
* **What it means for `3x + 1`.**
  * Any proof of the Collatz conjecture must use input beyond HYP-9244's least requirements: input that tells `3x + 1` on the positive integers apart from `3x − 1`, and from `C_3`, `C_11`, ….
  * The cycle half (Diophantine; finite checks per period, the S15 "LRC type") and the divergence half are exactly what coalescence does not see. The shapes above suggest where the cycle half lives for `3x + 1`: in the rarity of the near-balances `2^A ≈ 3^l`, which the base-`p` maps lack (ANALOGY).
  * In the reading of the 2026-10-04 atlas's "two walls" (ANALOGY), coalescence sits on the DRIFT side and cannot cross the INTEGRAL wall.
* **(S2) Boundaries** (`n ≤ 10^7`). The ratio is `R(N) = boundary/(2s(1−s) q_p(T_N))`, from `zp_boundary_ratio.py`, with `q_p(T_N)` from `10^5` Haar paths (standard error ≤ 0.0016):

| p | second-basin share | boundary density at `10^3 … 10^7` | `q_eff = boundary/(2s(1−s))` at `10^7` | `q_p(T_N)` at `10^7` | `R(N)` at `10^3 … 10^7` |
|---|---|---|---|---|---|
| 3 | 0.967 | 0.134, 0.069, 0.047, 0.040, 0.0366 | 0.577 | 0.484 (`T_N = 92`) | 1.65, 1.61, 1.40, 1.26, 1.19 |
| 11 | 0.018 | 0.026, 0.035, 0.032, 0.032, 0.0315 | 0.883 | 0.839 (116) | 1.14, 1.12, 1.06, 1.06, 1.05 |
| 17 | 0.292 | 0.402, 0.425, 0.408, 0.401, 0.3925 | 0.948 | 0.896 (143) | 1.12, 1.13, 1.08, 1.07, 1.06 |
| 23 | 0.234 | 0.479, 0.360, 0.347, 0.348, 0.3439 | 0.959 | 0.924 (169) | 1.34, 1.12, 1.06, 1.04, 1.04 |

  * The boundary has density 0 (THM-4610 5(e)). At these sizes, though, it is still most of `2s(1−s)` for `p ≥ 11`: integers below `10^7` take only 92–169 steps to descend, and the lazy skeleton rarely merges a pair that fast.
  * **The ratio.**
    * `R(N)` is above 1 throughout and drifts down at every `p`. That is the expected direction, since smaller `n` have shorter orbits. Whether `R(N) → 1` is not resolved.
    * This is THM-4590's heuristic for `3x + 1`, transplanted to base `p`. There, cuts ≈ `(1 − Σ d_i²) q(4.8 log_2 x)` reproduces the decline, but the measured level is about 1.13 times the prediction.
  * **Densities.**
    * HYP-9245 (a) assumes that the basin densities have limits.
    * Gonçalves–Greenfeld–Madrid (Remark 1.7) observed numerically that the basins of the three cycles of `3x − 1` have densities of about 1/3 each. `C_3`'s shares of 0.033 and 0.967 show that such shares need not be equal.
  * For `p = 3`, boundary × `√ln N` is 0.150 at `10^6` and 0.147 at `10^7`. Two points establish no plateau.
* **(S3)** For `p = 3`, the share of consecutive pairs `n ≤ 10^6` merging at equal time above the cycles is 0.459. The Haar prediction at the typical orbit length `T = 79` is 0.498 ± 0.003 from 40,000 Haar paths, and 0.502 ± 0.002 from `zp_boundary_ratio.py`'s `10^5` paths.
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
* **(R4) Orphans.**
  * **Partner classes.** A partner class is a set of exponents whose orbits first reach a cycle at the same element and at the same height `n − K`. Equivalently, `R_K` and `R_K′` merge with shift `K − K′` no later than their entry into the cycle. For every `K ≤ 8000`, that pair determines the multiplication count and the cycle (`zp_orphan_fit.py` (O1)). Orphans are the least members of the classes.
  * **Window fractions.** Orphan fraction × `√K` per dyadic window, from [800, 1600) to [6400, 8000], ± one Poisson standard error:
    * `p = 7`: 0.589 ± 0.157, 0.565 ± 0.130, 0.547 ± 0.107, 0.581 ± 0.175;
    * `p = 11`: 0.799 ± 0.183, 0.773 ± 0.152, 0.673 ± 0.119, 0.581 ± 0.175;
    * `p = 13`: 0.883 ± 0.193, 0.684 ± 0.143, 0.652 ± 0.117, 0.634 ± 0.183.
  * **Fits.** Poisson maximum-likelihood fits of `c K^(−α)` on `K ∈ [800, 8000]`, treating orphan events as independent, give `α` = 0.57 [0.21, 0.93], 0.69 [0.37, 1.01] and 0.71 [0.39, 1.03] for `p = 7, 11, 13` (95% profile-likelihood intervals). Every interval contains 1/2, and none resolves the exponent. Local exponents between windows carry standard errors of 0.4–0.8.
  * **Refutation rule.** HYP-9245 (b) states the law and a valid refutation rule, checked under the null by simulation in (O4). The original rule rejected a true law about 40% of the time.
  * **The case `p = 2`.**
    * On the Mersenne line the `K^(−1/2)` law is HYP-9243's prediction. HYP-9242 (fitted 0.59, CI [0.42, 0.75]) leaves faster decay open.
    * Those two count orphans among odd `K` only, since Mersenne class bottoms with `K ≥ 4` are odd. Per exponent their fraction is half as large.

## 5. Merge tails (NUMERICAL; `zp_tails.py`)

The start `(0, 1)`, 4800 paths per base. One binomial standard error of `√T q` is 0.01–0.02, 0.03–0.07, 0.2, 0.4–0.7 and 0.8–2.1 at the five times:

| p | `√T q(T)` at T = 10, 10², 10³, 10⁴, 10⁵ | `q(10^5)` ± s.e. |
|---|---|---|
| 2 | 2.05, 4.71, 8.31, 9.94, 10.87 | 0.0344 ± 0.0026 |
| 3 | 2.15, 4.95, 8.22, 9.77, 9.75 | 0.0308 ± 0.0025 |
| 7 | 2.85, 7.55, 16.28, 23.10, 25.76 | 0.0815 ± 0.0040 |
| 11 | 3.02, 8.51, 20.45, 34.35, 40.25 | 0.1273 ± 0.0048 |
| 13 | 3.05, 8.74, 21.88, 39.08, 49.87 | 0.1577 ± 0.0053 |
| 17 | 3.09, 9.14, 24.24, 49.00, 66.28 | 0.2096 ± 0.0059 |
| 23 | 3.12, 9.43, 26.25, 58.98, 89.14 | 0.2819 ± 0.0065 |

* **`p = 2`** is consistent with THM-4581's measured `√T q ≈ 11`.
* **`p = 3`.** The last two values, 9.77 ± 0.43 and 9.75 ± 0.79, are consistent with a plateau but do not establish one. mac-mini's `C_3` row agrees: `Z_3` (1, 4, 4) in `rank1_two_valued.out`, 1500 chains, gives 7.85, 9.09, 8.11 at `T = 1024, 4096, 16384`.
* **`q_3` at `T_N`.** The `T_N` column of `zp_tails.out` (4800 paths) gives `q_3(92) = 0.505 ± 0.007`. With `10^5` paths it is 0.484 ± 0.002 (`zp_boundary_ratio.py`), the value used in §3.
* For larger `p` the skeleton flips only at rate `2/p`, and the offset contracts at level 0 only through digit-0 steps. So the `T^(−1/2)` regime arrives later, with a larger constant.
* The lower bound `q ≥ c T^(−1/2)` is proved (THM-4610 5(c)). The matching upper bound is not claimed.

## 6. Other readings of Z_7, Z_11, Z_13 (DICTIONARY; checked in (F))

* **Period-4 cycles of `3x + 1`.**
  * `{7, 11, 13} = {|16 − 3^l| : l = 2, 3, 1}`.
  * The 16 fixed points of `T^4` (Terras) are `c_w/(16 − 3^l)`, one per class mod 16. Besides 0, −1 and the trivial cycle they form three non-integral 4-cycles: `{1, 8, 4, 2}/13`, `{5, 11, 20, 10}/7` and `−{19, 23, 29, 38}/11`.
  * The trivial cycle `{1, 2}` is the *integral* member of the `l = 2` (denominator-7) family. This is the 2-adic face of S17's "the trivial cycle's code is `{1,2,4}/7`".
  * In S6/S8 terms the 13- and 7-cycles are first-descent centres (contracting by 3/16 and 9/16 per 4 steps), and the −11-cycle is a rising centre (27/16).
* **Collatz multipliers acting on `Z_p`.**
  * `⟨2⟩ = QR_7`, `⟨3⟩ = QR_11` and `⟨3⟩ ∪ 2⟨3⟩` on `Z_13` are connection sets of circulant tournaments.
    * The first two give the Paley tournaments `P_7` and `P_11`.
    * On 13 vertices the connection set is `C_0 ∪ C_1`, a union of two quartic cyclotomic classes. Here `⟨3⟩ = {1, 3, 9}` is the subgroup of index 4.
  * Their affine automorphism groups have orders 21, 55, 39: `Z_7 ⋊ ⟨2⟩ = F_21` (S18's Prop. 3), `Z_11 ⋊ ⟨3⟩ = F_55` (cf. THM-4592's `BS(1,3) mod 11`), and `Z_13 ⋊ ⟨3⟩ = F_39`.
  * Since `3 ≡ 2^4 mod 13`, both Terras branches act on `Z_13` modulo `⟨3⟩` as the same rotation `2^(−1)`.
* **Further remarks.** `7 · 11 · 13 = 1001`. `7 = 111₂` maps to `13 = 111₃` along the Mersenne line through `11 = T(7)` (THM-4556 (i)).
* None of these readings carries a mechanism for the remaining step; they are recorded as dictionary.

## 7. Connections

* **HYP-9244** (mac-mini).
  * THM-4610 proves the almost-sure part of its requirement 1 for the maps `C_p` on every odd prime base. The route is the requirement's own sketch, with audit G's drift check.
  * The barrier (§3; KNOWN) shows that HYP-9244's three requirements (accessibility, contraction, recurrence) characterize Haar coalescence, not integer convergence.
* **THM-4581** (`p = 2`). The proof is its architecture. THM-4610 isolates what was parity-specific: the digit-equality step, and the absence of lazy steps.
* **THM-4590.** Its `o(x)` cuts for `3x + 1`, and the negative-integer basins, are the `p = 2` instances of 5(e) and of §3. Its cut heuristic is the model for HYP-9245 (a).
* **THM-4606–4609, THM-4611** (mac-mini): the other rows (expanding maps, `Z_2` classification, rank ≥ 3 and rank 0, block certificates).
* **HYP-9242 / HYP-9243 / THM-4605** (the Mersenne line). The repunit lines behave the same way: rivers, horizon constants, deletion exits for almost every exponent. Rivers carry the stopping-time invariant. The orphan law is HYP-9245 (b).
* **S17 / S18 / THM-4592.** See §6.

## 8. Status and next steps

| item | status |
|---|---|
| THM-4610 statements 1–5 (rank-one coalescence on `Z_p`; descent; density-one forms; repunit exits; lower rate; boundaries) | PROVED (audited) |
| the barrier (`3x − 1`, `F_p`); the cycles of `C_p` (Carnielli/Matthews) | KNOWN |
| local lemmas on boxes; descent to `|e| ≤ 20000`; cycles of `C_p` re-verified with their shapes; rivers to `K = 8000`; §6 facts | FINITE-EXACT |
| basin shares and boundaries; boundary ratios; merge tails; orphan exponents | NUMERICAL (with error bars) |
| boundary law and repunit orphan law (HYP-9245) | OPEN (HEURISTIC + NUMERICAL) |
| Collatz | OPEN |

**Next steps.**
1. **The general two-valued case.** Prove accessibility for two-valued translation-only maps, or classify its failures by HYP-9244's obstruction lemma. This would close the rank-one row entirely.
2. **The `T^(−1/2)` upper bound for `C_p`.** THM-4581 (4)'s sketch should carry over, with constants growing like `p`.
3. **The arithmetic step.**
   * For `3x + 1`, the input coalescence cannot supply is the cycle and divergence structure.
   * The base-p family is a test bench, alongside `3x − 1`. On the positive integers it is false for `p = 3, 11` and apparently true for `p = 5, 7, 13`. Any proposed proof mechanism must fail for `C_3`, for `C_11` and for `3x − 1`.
   * For the cycle half, the shapes in §3 point somewhere specific. Every extra cycle found is a 1- or 2-cycle at a cheap near-balance, with many candidate numerators per shape. Those are exactly the cycle types that linear forms in logarithms exclude for `3x + 1`, where each shape has a single candidate (Steiner; Simons–de Weger; Hercher).
   * **A naive count makes this quantitative ((S5), HEURISTIC).** Suppose each candidate is divisible with probability `1/|D|`.
     * Each base `p = 3, …, 31` then expects 0.73–0.96 positive 1-cycles (trivial cycle included) and 0.10–1.45 negative ones, the same order as found.
     * For `p = 2`, beyond the three shapes with `|D| = 1` (the cycles 1, −1, −5), the expectation is 0.30 positive and 0.18 negative from shapes with `k_2 ≤ 5`, all checked exactly. It is 0.005 and 0.019 from `5 < k_2 ≤ 2000`.
     * So the heuristic places every `3x + 1` 1-cycle in shapes small enough to check, while for `p ≥ 3` it expects about one per base, spread over many shapes.
     * Turning that small expectation into "none" is what Steiner's theorem does for 1-cycles. The open part of the cycle half is the same statement for every `m`.
4. **HYP-9245.** Run the tests in its refutation section.
