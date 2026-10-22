---
id: HYP-9369
title: "The loop expansion beyond regular tournaments, tested at scale, and the joint remainder. On almost-regular even-order tournaments the weight-graded second-order predictor (Session30 decorated cores, w = excess + leaves) has RMSE(y - y2) of about 200/n^4 in y = n(1 - kappa/rho): 0.048, 0.020, 0.0096, 0.0053, 0.0030, 0.0010 at n = 8 (all 85 classes), 10 (all 13,333 classes), 12, 14, 16, 20 (samples). The remainder must be bounded jointly: after second order the cluster part log Phi(1) - log rho and the finite-size part log(kappa/Phi(1)) miss by -(234..260)/n^3 and +(235..261)/n^3, per-tournament correlation -0.998 to -0.999, while their sum misses by about 1/n^3 on average with per-tournament spread about 200/n^5 (all 13 doubly regular 23-tournaments: -218.5, +225.3, joint +6.85). A terminal leaf acts as a pendant loop with the vector 1 in place of the loop diagonal d (decoration rule; square and claw identities)."
status: >
  OPEN (opus-2026-10-22-S54; sections 3-5 independently audited 2026-10-22,
  04-computation/experiments/loop_nonregular_20261022/audit_S54/REPORT.md: numbers confirmed, claw identity re-verified
  with Session30's engine on nine tournaments including random and transitive ones, notation and attribution corrected;
  sections 1-2 were outside the audit scope and are reproducible from the committed files). FINITE-EXACT: H for every class (Held-Karp, ar_tool.c) and exact forest numbers for
  n <= 16; the class lists at n = 4, 6, 8, 10 are complete (labelled masses 24, 2640, 3,230,080, 48,251,508,480 = the
  numbers of labelled regular 5-, 7-, 9-, 11-tournaments, count_regular.py). Samples at n = 12..20 come from a
  score-preserving 3-cycle-reversal chain (uniform in the limit; mixing not certified). Predictions are double precision,
  no fitted parameters. The predictor and its 5 + 26 cores are Session30's
  (05-knowledge/results/session30_nonregular_loops_20261022, imported unchanged); its n = 8 census numbers are reproduced
  exactly on this lane's independently generated class list. The decoration rule is PROVED for the two-arm family
  (algebra of proved closed forms) and VERIFIED to machine precision for the three-arm (claw) family, on almost-regular and
  (audit) on uniform random and transitive tournaments; the general rule is conjectural.
session: opus-2026-10-22-S54
source: 05-knowledge/results/session30_nonregular_loops_20261022/SYNTHESIS.md
scripts:
  - 04-computation/experiments/loop_nonregular_20261022/ar_tool.c: classes (canonical form + |Aut|), samples, exact H, forest numbers, deficits
  - 04-computation/experiments/loop_nonregular_20261022/budget_cw.py, crosscheck_s30.py: Session30 predictor and the exact z = 1 budget
  - 04-computation/experiments/loop_nonregular_20261022/predict_nonregular.py, summarize_tables.py (+ tables.out): baselines and tables
  - 04-computation/experiments/loop_nonregular_20261022/leaf_engine.py, leaf_closed.py, closed_forms.py, decoration_rule.py: independent leaf engine, closed forms, decoration rule
  - 04-computation/experiments/loop_nonregular_20261022/count_regular.py, complete_regular.py, deficit_*.txt: completeness and the deficit census
related:
  - 05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md
  - 05-knowledge/hypotheses/HYP-9367-three-loop-law-the-loop-expansion-is-asymptotic.md
  - 05-knowledge/hypotheses/HYP-9370-zero-free-disc-of-the-forest-polynomial-for-almost-regular-tournaments.md
  - 01-canon/theorems/THM-4761-random-tournaments-the-forest-identity-is-the-walsh-chaos-decomposition-of-kappa.md
  - 05-knowledge/results/session30_nonregular_loops_20261022/session30-deficit-six/INTEGRATION.md
---

# HYP-9369 — the almost-regular loop predictor at scale, and the joint remainder

## Credit and scope

The task (extend the linear-forest loop expansion from regular to almost-regular tournaments) was worked concurrently by Session30. Their commits 2b9185a253, 27fbfce27a and c2f3434066 (2026-10-22, 09:21–09:52) are audited and closed. They contain:
* the decorated-core engine with the grading `w = e + ℓ` (excess plus terminal leaves);
* the 5 weight-one and 26 weight-two cores;
* the closed path, lollipop and Y-tree kernels;
* the first/second-order predictor;
* the complete order-8 census;
* the flat-family constants and the lollipop obstruction.

This lane derived the path-tree, lollipop and star closed forms independently before reading those commits. This session's file listing before the fetch (about 09:58) shows `leaf_engine_selftest.out` at 09:25 and `leaf_closed.out` at 09:52; a later rebase reset the file times. Session30's lollipop and Y-tree formulas are in their commit 27fbfce27a (09:33). Priority is theirs, and nothing of theirs is re-derived here.

**Inheritance pass.**
* Closest proved mechanism: THM-4757 (forest identity, Möbius clusters) plus Session30's coefficientwise remainder theorem (`[z^L](log Φ − X0 − A1 − A2) = O_L(n^(−3))`; their §8).
* Hostile families: transitive tournaments (`H = 1`) and uniform random tournaments, whose score scale is `√n` (THM-4761).
* Corrected near miss: the path-tree-only extension (mean offset −0.62 at n = 8), and this lane's own measured-order grouping (§3).
* Least-used sidecar: the exact forest numbers `c_L = (n − L)! N_L`. They split `log κ` exactly into a cluster part and a finite-size part (§2).

This HYP adds:
* tests at scale (§1);
* the exact budget at `z = 1` (§2), which is Session30's open obligation 1;
* two documented failures (§3);
* the decoration rule (§4);
* a merge with the deficit-six thread (§5).

## 1. Prediction tables

`y = n(1 − κ/ρ)` with `κ = 2^(n−1)H/n!`. "Before" means the regular two-loop formula applied blindly, and the same plus the path trees (mac-mini's `almost_regular_tree.py`). "After" means Session30's weight-one and weight-two predictors. Scatter is the standard deviation of `y − pred`; offset is its mean.

| n | tournaments | regular formula: corr, offset | + path trees: corr, offset | weight 1: corr, RMSE | weight 2: corr, Spearman, scatter, offset, RMSE | weight 2 improves |
|---|---|---|---|---|---|---|
| 8 | all 85 classes | 0.790, +0.564 | 0.860, −0.622 | 0.891, 0.096 | 0.935, 0.901, 0.0458, −0.0138, 0.0479 | 73/85 |
| 10 | all 13,333 classes | 0.844, +0.630 | 0.924, −0.504 | 0.950, 0.084 | 0.974, 0.968, 0.0178, −0.0092, 0.0201 | 13,203/13,333 |
| 12 | 1500 samples | 0.856, +0.682 | 0.952, −0.424 | 0.973, 0.078 | 0.990, 0.987, 0.0081, −0.0052, 0.0096 | 1500/1500 |
| 14 | 1500 samples | 0.839, +0.721 | 0.963, −0.367 | 0.983, 0.073 | 0.995, 0.994, 0.0043, −0.0032, 0.0053 | 1500/1500 |
| 16 | 600 samples | 0.819, +0.751 | 0.974, −0.323 | 0.990, 0.068 | 0.998, 0.997, 0.0024, −0.0017, 0.0030 | 600/600 |
| 20 | 300 samples | 0.843, +0.796 | 0.986, −0.261 | 0.9945, 0.058 | 0.9993, 0.9989, 0.0009, −0.0004, 0.0010 | 300/300 |

Full numbers, including MAE, are in `tables.out`. The n = 8 row reproduces Session30's census report exactly: MAE 0.0831/0.0384, RMSE 0.0964/0.0479, 73 improved. Their n = 20 held-out errors (0.0007–0.0016 on three samples) fit the 300-sample RMSE 0.0010.

Four observations:
* **The regular formula does not transfer.** Applied blindly, it keeps a mean offset of about +0.6 to +0.8 at every `n`.
* **Path trees overshoot.** Adding them reverses the sign; the offset shrinks only like `1/n`.
* **Weight 1 leaves a near-constant offset.** The weight-one predictor's offset stays near −0.06 to −0.08, i.e. a second-order effect in `κ/ρ`.
* **The weight-2 error falls like n⁻⁴.** `n⁴·RMSE(y − y2)` is 196, 201, 199, 204, 194, 154 for `n = 8, 10, 12, 14, 16, 20` (`scaling_table.out`). The per-tournament scatter decays at the same rate, and Spearman reaches 0.9989 at `n = 20`. In `κ/ρ` the residual is about `200/n⁵` per tournament.

The formal theory only claims `O(n^(−2))` in `y`. Where the extra orders come from:
* The per-tournament scatter decays like `n^(−4)` in `y` because the joint log remainder of §2 has spread `≈ 200/n⁵`.
* Its mean `≈ 1/n³` shifts all tournaments nearly equally. In `y` this is a mean offset of about `−0.9/n²` at `n = 8, 10`, and less beyond. The additive form `1 + Y1 + Y2` of the predictor differs from its log form at third order and partly compensates.
* Ranking is therefore governed by the `n^(−4)` scatter: Spearman 0.901 → 0.9989.

## 2. The exact budget at z = 1: only the joint remainder is small

Exact forest numbers split `log(κ/ρ)` into two parts:
* the cluster part `A = log Φ(1) − log ρ`;
* the finite-size part `B = log(κ/Φ(1)) = log E_φ[n^L/(n)_L]`.

The predictor approximates `A` by `W1 + W2`, the weight-one and weight-two core sums at `z = 1`. It approximates `B` by the log of Session30's inverse-factorial operators through second order: `B_f = p1/n + [D²W1 + (2c1 − 1)DW1]/(2n) + (p2 − p1²/2)/n²`.

| n | tournaments | `n³(A − W1 − W2)` | `n³(B − B_f)` | correlation of the two errors | `n³·mean(joint)` | `n⁵·sd(joint)` |
|---|---|---|---|---|---|---|
| 8 | all 85 classes | −234.2 | +235.0 | −0.9983 | +0.84 | 233.5 |
| 10 | all 13,333 classes | −240.8 | +241.9 | −0.9986 | +1.15 | 213.3 |
| 12 | 1500 samples | −247.0 | +248.1 | −0.9988 | +1.14 | 198.9 |
| 14 | 1500 samples | −254.2 | +255.3 | −0.9989 | +1.12 | 194.0 |
| 16 | 600 samples | −259.6 | +260.6 | −0.9990 | +1.02 | 186.5 |
| 23 | all 13 doubly regular classes | −218.5 | +225.3 | (classes agree to ±0.12) | +6.85 | — |

Here `joint = log(κ/ρ) − W1 − W2 − B_f` (`scaling_table.out`, `budget*.txt`).

The `n = 23` row is a regular family. It uses mac-mini's exact forest numbers for the 13 doubly regular 23-tournament classes (`budget_drt23.out`; HYP-9371 lane). The same cancellation appears there: separate parts `∓220/n³`, joint `+6.85/n³`. This joint remainder is the second-order predictor's actual residual: `y − y2 = −0.0117` for every class, with `P_23` the outlier at `+7.02`. It is about 3% of either part. For doubly regular `T` the cluster part `A − W1 − W2` is the four-loop vacuum term plus all higher loops. Since `n³X_3(P_23) = −363.5` (HYP-9371), the five-loop and higher terms contribute about `+145/n³` at `n = 23`. That is the asymptotic-series growth that HYP-9371 records (3.5, 29, 341).

Three facts stand out:
* **Each part is third order with a huge constant.** `n³(A − W1 − W2)` runs from −234 to −260 and `n³(B − B_f)` from +235 to +261 over `8 ≤ n ≤ 16`, drifting by about 3.2 per unit of `n`.
  * The finite-size side is identified (`finite_size_third.py`). Write `log(n^L/(n)_L) = a/n + b/n² + c/n³ + …` with `a = C(L,2)`, `b = Σ_(i<L) i²/2`, `c = a²/3`. The third-order term `k3(a)/6 + Cov(a,b) + E c` under the one-loop measure tracks `n³(B − B_f)` with per-tournament correlation 0.980–0.986 at every `n` (`finite_size_third.out`). Its mean falls as 457, 352, 300, 272, 254 for `n = 8, …, 16`, towards the flat-limit value below.
  * In the flat limit `L/2 ~ Poisson(1)` (`ρ(z) → e^(z²)`), the inverse-factorial constants are exactly `E a = 3`, `Var(a)/2 + E b = 64/3` and `k3(a)/6 + Cov(a,b) + E c = 761/6 + 161/2 + 14 = 664/3 ≈ 221`.
  * The cluster side is of the size of the four-loop vacuum term: mac-mini reports `q³X_3 → −341` on doubly regular tournaments (HYP-9371).
* **The cancellation happens tournament by tournament.** The correlation of the two errors is `−0.998` or closer to `−1` at every `n`.
* **The joint remainder is two orders of magnitude smaller.** Its mean satisfies `n³·mean ≈ 0.84–1.15`, about 250 times smaller than either part. Its per-tournament spread decays two powers faster: `n⁵·sd ≈ 233 → 187`.

The same pattern already holds one order down for doubly regular tournaments:
* The cluster constant 29 (`q²X_2`, HYP-9367) and the finite-size second-order constants combine to the `log κ` coefficient 2 of `κ = e(1 − 3/(2n) + 25/(8n²))`.
* For flat almost-regular families Session30 has 33 (= 29 + 4) on the cluster side, against the `log κ` coefficient `65/8 − 25/8 = 5`.

**Consequence for the open uniform remainder.** A proof that bounds the omitted cores and the falling-factorial tail separately would have to control errors about 250 times larger than the joint error at `n ≤ 16`. Its constants would have to absorb `±240/n³`. The natural object is `log κ` itself. The exact Gamma average `κ = E[Φ(n/T)]`, `T ~ Gamma(n+1)` (THM-4757 §5), is one form in which both expansions live together.

## 3. Documented failures and boundaries

* **Asymptotic orders are the wrong bookkeeping at accessible n.** On flat almost-regular tournaments the lollipop and Y-tree are demoted by one power of `n`: `Σδ = 0`; `class_scaling.out` gives `n²X(0,1) → 4.0` and `n²X(−1,3) → 0` (with `n³X(−1,3) ≈ −2.2`). Grouping classes by these measured orders gives a second-order predictor that is clearly worse at small `n`. At `n = 8` its RMSE is 0.210 against 0.048 for the weight grading `w = e + ℓ`. The leaves must be counted by weight, as Session30 does.
* **Exact finite-size weights do not resum the expansion.** Replacing the `1/n` expansion of `n^L/(n)_L` by its exact values, applied to the Taylor coefficients of `exp(X0 + W1 + W2)` truncated to `L ≤ n − 1`, is catastrophic (`budget_cw8.out`):
  * At `n = 8`, the mean error in `y` is +8.3 (`T0`), −8.9 (`T1`) and +5.1 (`T2`), against RMSE 0.048 for the formal predictor.
  * The exact forest measure is strongly suppressed at large `L` by vertex exclusion (`φ_6 = 0.00049` against 0.127 for the one-loop series).
  * The falling-factorial weights (13.0 at `L = 6`, 52.0 at `L = 7`) amplify exactly those coefficients.
  * The formal expansion works because its two halves are expanded consistently in `1/n` (§2).
* **Random tournaments are outside the almost-regular grading.** See THM-4761 §2. The path-tree class is of order one there and the second-order leaf classes worsen the prediction.

## 4. The decoration rule: a leaf is a loop with 1 in place of d

Let `d = diag(H − G)`, the diagonal of the loop propagator, and `u = d + 1`. Use Session30's `G = −Ŝ²(I + Ŝ²)^(−1)` and `Q̃ = Ŝ(I + Ŝ²)^(−1)` (Session30's `h = Q̃1`; their `Q` is `zŜ`). Summing a skeleton over its pendant decorations replaces the end vector by `u`:

* **Two arms (square identity, PROVED).**
```
dumbbell + lollipop + path tree = −dᵀGd − 2dᵀG1 − 1ᵀG1 = −uᵀGu.
```
The three closed forms are THM-4757 (4) and Session30's lollipop and path kernels.
* **Three arms (claw identity, VERIFIED).** The claws are a centre joined by three bridges to ends that are leaves or carry one loop: one core each in classes `(2,0)`, `(1,1)`, `(0,2)`, `(−1,3)`. Their sum is
```
Y(u) = −(2/3) Σ_w [ (Gu)_w³ + 3 (Gu)_w (Q̃u)_w² ],
```
the Y-tree formula with `1 → u`.
  * Agreement is to machine precision: absolute differences below `1e-17`, relative about `1e-15`, at `n = 8, 12, 16` (`decoration_rule.out`).
  * The audit re-verified it with Session30's engine on nine tournaments, including uniform random and transitive ones: max difference `5e-15`; two-arm `1e-15`.
  * The identities therefore appear to hold for all tournaments, not only almost-regular ones.
  * The lane's brute-force self-test reaches only `L ≤ 6`. The `(1,1)` and `(2,0)` claw cores first occur at `L = 8` and `10`, so this check rests on the engines, which were validated against brute force at lower orders.
* **Conjecture (general rule).** For every skeleton, the sum over leaf-or-pendant-loop decorations of its pendant ends equals the skeleton with end vector `u`.
  * The local reason is that the lollipop kernel factorises as tail kernel × loop insertion: `4x²y²/((1 − x⁴)(1 + y²))` (Session30) is the path kernel `y²/(1 + y²)` times `4x²/(1 − x⁴)`, whose sum is `2·diag(H − G)`.
  * ANALOGY (audit, no identical statement found): tadpole dressing of an external source, and the articulation-point topological reduction of Mayer expansions.
  * On flat almost-regular families `d = −(2/n)1 + O(n^(−2))`, so decoration renormalises every leaf by `1 − 2/n`. This is why `n²X(0,1) → 4`: `X(0,1) ≈ (4/n)·(−X(−1,2))` with `X(−1,2) ≈ −1/n`.

## 5. Merge with the deficit-six thread

The deficit of a self-converse tournament is `Δ = h − b`, where `h = H/|Aut|` and `b` is the number of blue Hamiltonian-word orbits. A word is blue when positional reversal is an anti-automorphism. Aut acts freely and preserves blueness, so `b = #blue/|Aut|`. `ar_tool deficit` computes `Δ` for every self-converse class of the complete regular and almost-regular lists (`deficit_ar8.txt`, `deficit_reg9.txt`, `deficit_ar10.txt`, `deficit_reg11.txt`; FINITE-EXACT):

| n | family | classes | self-converse | min Δ | max b | Δ = 6 |
|---|---|---|---|---|---|---|
| 5 | regular | 1 | 1 | 0 (`C_5`) | 3 | — |
| 6 | almost-regular | 5 | 3 | 6 (`C_3[K_2]`: h = 15, b = 9) | 9 | `C_3[K_2]` |
| 7 | regular | 3 | 3 | 6 (`P_7`: h = 9, b = 3) | 9 | `P_7` |
| 8 | almost-regular | 85 | 15 | 56 | 39 | none |
| 9 | regular | 15 | 15 | 28 | 39 | none |
| 10 | almost-regular | 13,333 | 229 | 1592 | 193 | none |
| 11 | regular | 1223 | 229 | 1692 | 193 | none |

Points of contact:
* **Known centers.** The two hits are exactly the known deficit-six centers of orders 6 and 7, `(6,15,9)` and `(7,9,3)`. The order-5 center `(5,9,3)` is not regular.
* **Session30's completion theorem, confirmed.** The self-converse counts 15 = 15 and 229 = 229 confirm that self-converse almost-regular classes at even `n` match self-converse regular classes at `n + 1`.
* **Consistency with the inherited classification.** Almost-regular tournaments are strong (Session30's lemma). The census therefore checks the inherited classification of deficit-six centers through order 11 inside this family. It does not touch the residual orders 15–28.

**A counting bound (PROVED, elementary).** A blue word is fixed by its first `⌊n/2⌋` letters, chosen from distinct orbits of its reversal involution ψ. Each anti-automorphism carries at most `2^⌊n/2⌋ ⌊n/2⌋!` blue words, and there are `|Aut|` anti-automorphisms. So
```
b ≤ 2^⌊n/2⌋ ⌊n/2⌋!.
```
The observed maxima (39, 193) are far below the bound (384, 3840). For regular and almost-regular `T`:
* `H ≥ n!/(2 + o(1))^n` (Cuckler's lower bound for Hamilton cycles of regular tournaments, applied to `T` or to its regular completion);
* `|Aut| ≤ 3^((n−1)/2)` (Dixon 1967, Canad. Math. Bull. 10; Moon, Topics on Tournaments, §28 Ex. 5).

Then `n!/(n/2)! ≥ (n/2)^(n/2)` gives `h/b ≥ (n/(48 + o(1)))^(n/2) → ∞`, and a bounded deficit occurs in this family for only finitely many orders.

This is weaker than inherited canon for `Δ = 6`: the proved center cutoff already bounds every strong deficit-six center by order 28 (`deficit_six_centers_20261021b`, sharpened to 28 in commit 4119d59435). The counting bound is an independent, elementary reason specific to the regular regime, and it quantifies how fast the deficit grows there. Its threshold is not explicit, because Cuckler's `o(1)` is not.

**Next question (not pursued).** For a fixed anti-automorphism `ψ`, blue words are half-paths on the `⌊n/2⌋` orbits of `ψ`, closed by the middle arc `v_m → ψ(v_m)`. Their count has its own forest identity and a twisted one-loop gain from the skew matrix `S P_ψ`. This would make `Δ = h − b` loop-predictable in the regular regime. The two threads are complementary:
* the loop expansion lives where `h` is huge and bounded deficits cannot occur;
* the deficit classification lives at the padded, near-transitive boundary, where `h` is small.

**Onset 19, independently replayed (conditional).** `deficit_six_onset.py` (sympy, written from INTEGRATION.md §3 without its checker) rebuilds the degree-three-node series under `C = 0`:
```
G = z⁵/2 (1/p² − 1/u) + R_6/u + z¹⁰(1 + z + z³)/u²,   p = 1 − z − z³,   u = p(z²).
```
* The reduced denominator is `p²u²` (degree 18). The numerator is `z⁵ − 2z⁷ − 2z⁸ + z⁹ + 3z¹⁰ − z¹⁴ + 2z¹⁵ − z¹⁶ + z¹⁷ − z¹⁸`, with gcd 1.
* The denominator recurrence fails at index 18 and holds from 19.
* `g_5..g_16 = 1, 2, 3, 4, 7, 11, 18, 29, 47, 71, 114, 175`.

This agrees with Session30 and with the user's "onset 19 rather than 41". It checks the algebra only, not the closure.

The user-reported closure of the six residual orders (`R_6 = z⁵ + z⁶ + z⁷`) and the composed Collatz bounds 520,738, 1,039,055, 3,848,613 are not in the repository as of this session (Session30 recorded the same). Within the regular and almost-regular family they agree with everything above. Restricted to that family, `R_6` is `z⁶ + z⁷` through order 11 unconditionally (the census), and at all orders conditionally on the reported closure. The counting bound is a genuinely uniform-in-`n` argument for this family. It is the kind of all-scale statement the deficit-six and Collatz lanes still need for their padded and exterior families. It is not a substitute for their certificates, and no Collatz implication is claimed.

## 6. Statement and refutation targets

**Conjecture.** For almost-regular tournaments of even order `n`:
* (i) `|y − y2| ≤ C/n²` uniformly. Empirically `RMSE(y − y2) ≈ 200/n⁴` on random almost-regular tournaments.
* (ii) The joint remainder `log(κ/ρ) − W1 − W2 − B_f` is `O(n^(−3))` with a small mean constant (0.84–1.15 at `n ≤ 16`) and per-tournament spread `O(n^(−5))`.
* (iii) Separately, `A − W1 − W2` and `B − B_f` are each of size `≈ 250/n³`, with opposite signs and per-tournament correlation tending to −1.

**Refutation targets.**
* An almost-regular family on which `n²|y − y2|` grows.
* The joint remainder growing faster than either separate remainder.
* A pendant-decorated skeleton violating the decoration rule.
