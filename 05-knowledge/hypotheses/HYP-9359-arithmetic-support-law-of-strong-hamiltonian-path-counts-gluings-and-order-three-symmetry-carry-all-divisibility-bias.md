---
id: HYP-9359
title: "The arithmetic-support law of strong Hamiltonian-path counts: the number S(d) of strong tournament classes with d Hamiltonian paths satisfies log S(d) = smooth(log d) + 0.22[3|d] + 0.06[9|d] + 0.05[5|d] + small 2-adic terms, with primality neutral; the divisibility bias is carried almost entirely by exact one-vertex gluings (H = H(A)H(B)) and by automorphisms of order 3 (|Aut| divides H), and the generic remainder keeps only a weak attraction to 1 (mod 8) that decays with d"
status: >
  NUMERICAL (complete strong inventory for every H <= 2100, 723,544 classes): the coefficients, their stability between
  [301, 1098] and [1101, 2099], and the mechanism decomposition; all reproduced exactly by the independent audit
  (MISTAKE-657). Eleven predictions were committed before their evaluation (R1-R7: 917fef912a; G1-G4: be7c76115f), and
  all eleven PASSED. Their blindness is partial (audit chronology, audit_S51/c7):
  - R1-R7 precede the completion of the order >= 10 test data by 9.6 s, on local timestamps only; the remote push came
    later;
  - the order-9 slice of the test range (16% of its classes) had been examined in exploratory form through the class
    census n = 9 (correlations with 3|h, 5|h and primality);
  - G1-G4 were written after the test-range S(d) and the R verdicts were known; they were blind only to the
    decomposition.
  OPEN: the limits of the coefficients as d -> infinity (expected to tend to 0 with the Session25/HYP-9355 equidistribution
  picture), a closed-form model of the coefficients from the gluing and symmetry counts, and the mechanism of the residual
  2-adic attraction. Session opus-2026-10-21-S51.
session: opus-2026-10-21-S51
source: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md
scripts:
  - 04-computation/experiments/repunits_hyp9357_20261021/repunit_rarity.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/class_residues.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/refuted_leads.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/outliers.py (+ .out)
  - 04-computation/experiments/repunits_hyp9357_20261021/PREREGISTRATION_out_of_sample.md
  - 04-computation/experiments/repunits_hyp9357_20261021/audit_S51/ (independent audit: c4, c5, c5b, c7, c8)
related:
  - THM-4748, THM-4749 (this session)
  - HYP-9355 (H mod M equidistribution; the 2-adic pull at small n)
  - THM-4744/4745 (strong inventories; exact gluings)
  - Session24 multiplicity packet (polynomial degree law)
---

# HYP-9359 — the arithmetic-support law

## The law (NUMERICAL)

* **Data.** `S(d)` is the number of strong classes of all orders with `H = d`. It is complete for `d ≤ 2100`, by ears from `C_3` with the audit's method.
* **The model** is

  `log S(d) = c_0 + c_1 log d + c_2 log² d + c_3 log³ d + β_3[3|d] + β_9[9|d] + β_5[5|d] + β_p[d prime] + β_b[d blocked] + β_7[d ≡ 7 mod 8]`,

  where "blocked" means composite with every factorization through 7 or 21.

| coefficient | training `[301, 1098]` (399 values, s.d. 0.100) | test `[1101, 2099]` (500 values, s.d. 0.047) |
|---|---|---|
| `β_3` | +0.239 ± 0.014 | +0.219 ± 0.006 |
| `β_9` | +0.088 ± 0.018 | +0.060 ± 0.008 |
| `β_5` | +0.065 ± 0.014 | +0.054 ± 0.006 |
| `β_p` (prime) | −0.002 ± 0.014 | +0.005 ± 0.006 |
| `β_b` (blocked) | +0.028 ± 0.023 | +0.008 ± 0.010 |
| `β_7` (`≡ 7 mod 8`) | −0.027 ± 0.012 | −0.015 ± 0.005 |

## Pre-registered tests (all PASSED; blindness partial)

* **R1–R7** (commit 917fef912a; committed 9.6 s, by local timestamps, before the order-≥10 classes of the test range existed; the order-9 slice had been looked at in exploratory form). The predictions, with the observed values:
  * `β_3 ∈ [0.15, 0.40]` and largest;
  * `|β_p| < 0.05`;
  * `β_7 ∈ [−0.05, 0.02]`;
  * `S(2047) = 1599 < S(2049) = 1968`;
  * the base-4 repunit 1365 lies above the smooth fit (+0.148);
  * the order-14 coefficient of `[3|d]` is +1.19 ≥ 0.4;
  * the sign rule for `N(n, 2^k ± 1)` at `n = 10, 11, 12`, `k = 4..10`.
* **G1–G4** (commit be7c76115f; committed after the test-range `S(d)` and the R verdicts were known, before the decomposition was run on `[1101, 2099]`). The predictions, with the observed values:
  * without the classes that have an exact cover, `β_5 = −0.005` and `β_9 = −0.009`, while `β_3 = +0.135`;
  * with the order-3-symmetric classes also removed, `β_3 = +0.010`;
  * in the generic part, `β(1 mod 8) = +0.017 ± 0.005 > 0` and `β_7 = −0.012`;
  * the share of classes with `3 | |Aut|` drops from 0.0514 to 0.0447.
* **Two discrepancies in the pre-registration text** (they do not change any verdict):
  * the G-section quotes a training baseline `β(1 mod 8) = +0.049 ± 0.016` from an exploratory fit with extra 2-adic terms; the committed fit gives +0.046 ± 0.012;
  * the G-section says "trivial Aut" where G2 and the script use `3 ∤ |Aut|`. G2 and G3 pass under both readings (audit `c5b`).

## The mechanism (NUMERICAL decomposition of a PROVED pair of facts)

* **(i) Exact one-vertex gluings.** If `T = A ∪_v B` with `A`, `B` strong and `H(T) = H(A)H(B)` (THM-4745 (1)), then `H(T)` inherits every factor of `H(A)` and `H(B)`.
  * The cheapest factors are `3 = H(C_3)`, `5 = H(R_4)` and `9 = H(IC(2,2))`.
  * Removing every class with an exact cover lowers `β_3` from 0.237 to 0.139 (training) and from 0.217 to 0.135 (test). These fits use the decomposition model (`3|d`, `9|d`, `5|d`, `1 mod 8`, `7 mod 8`), without the prime and blocked terms, hence 0.237/0.217 here against 0.239/0.219 in the table.
  * It removes `β_5` and `β_9` entirely.
* **(ii) Order-3 symmetry.** Aut acts freely on Hamiltonian paths, so `|Aut(T)|` divides `H(T)`.
  * Of the strong classes with `301 ≤ H ≤ 1098`, 5.1% have `3 | |Aut|`; for `1101 ≤ H ≤ 2099`, 4.5%.
  * Removing them removes the remaining `β_3 ≈ 0.13–0.15`.
* **What is left.** After (i) and (ii), `β_3 = −0.015 ± 0.012` (training) and `+0.010 ± 0.005` (test, about 2 s.e.). So the two mechanisms carry almost all of the bias, not provably all of it.
  * The main arithmetic signal left is a weak attraction to `1 (mod 8)`: +0.046 ± 0.012 on `[301, 1098]` and +0.017 ± 0.005 on `[1101, 2099]`.
  * Primality is neutral in the full model, `β_p` in the table. The decomposition fits have no prime term.
  * The class-level 3-adic excess at `n = 9` (66,320 / 62,129 / 63,087 classes at `H ≡ 0, 1, 2 mod 3`) is this law seen at the class level. The labelled distribution is within 2% of uniform (HYP-9355).

## Consequences for the repunits

* For odd `k`, `2^k + 1` is a multiple of 3, so it receives both supports, and `2^k − 1` receives neither. That gives the comparisons 31/33, 127/129, 511/513 and 2047/2049.
* 63 is a multiple of 9, but both supports need a gap: `63 = 3·21 = 9·7`, and the order-8 atoms are rigid (THM-344).
* 255, 1023, 1365 are supported, and common.
* After the full model refitted on `[25, 1098]`, 63 (−1.68) and 127 (−0.55) remain outliers (`outliers.py`). Their remaining rarity is a small-order phenomenon:
  * `S_7(63) = 0`, a hole of the order-7 strong spectrum;
  * at orders 6–8, `H ≡ 7 (mod 8)` is depleted to 0.47–0.78 of uniform in the second octave `2 ≤ H/f(m) < 4`;
  * with at most two reversed arcs, `7 (mod 8)` is reached only along the lines of THM-4748 (2). This does not explain the depletion: the `7 (mod 8)` share by feedback arc number is not consistently low across orders 6–9 (`refuted_leads.py` (c)).

## Open

1. **Limits.** Do `β_3, β_5, β_9 → 0` as `d → ∞`? The symmetric share decreases with order. The gluing share needs a model: is it about `S(d/3)·g/S(d)`, with `g` the mean number of exact `C_3` attachments per strong class?
2. **A closed form.** Predict `β_3` from (i) and (ii) alone, by counting exact `C_3` gluings and order-3 symmetric classes per value, without fitting.
3. **The residual 2-adic attraction.** Find a mechanism. One candidate: one-apex and near-transitive families produce `H ≡ 1 mod 2^j` with high multiplicity, since their single-arc coefficients `2^(s−1)` are high powers of two.
4. **The order-7 hole at 63.** No ear of any strong 6-tournament has `H = 63`; all labelled out-sets were checked. Is there a structural reason, or is it a small-number accident?
