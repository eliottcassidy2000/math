---
id: HYP-9311
title: "CM harmonic law for the four Ramanujan families: for A(k) in {C(2k,k)^3, C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k) C(3k,k) C(6k,3k)} with Sun's weights w_k = d/dk log A(k), at every CM parameter m and every prime p not dividing 6m, sum_{k<p} A(k) m^-k w_k == q_p(m) sum_{k<p} A(k) m^-k (mod p). For A1 this is Kibelbek-Long-Moss-Sheller-Yuan 2016 Corollary 3; for A2-A4 it is OPEN, and it is EQUIVALENT to the Lucas-block detector congruence (S43b m3). It is the constant-term analogue of Z.-W. Sun's Conjecture 29. Heuristically it is the first-order truncation of the fixed-point equation q_0^(p-1) = 1 of the excellent Frobenius lift in the canonical (mirror / Serre-Tate) coordinate, whose logarithm has power-series part sum A(k) w_k x^k / sum A(k) x^k"
status: >
  OPEN for A2-A4, strongly supported. For A1 it is PUBLISHED: Kibelbek-Long-Moss-Sheller-Yuan 2016 Corollary 3
  (S43b M1; THM-4686 gives an alternative proof).
  - Equivalence (S43b m3). The p-shift relation also holds at carry indices:
    A(p+j) == A(1) (A(j) + p [A(j) w_j]) (mod p^2) for all 0 <= j < p, with [A(j) w_j] the p-integral product (491
    checks per family, p <= 61). So (L_p^A) is EQUIVALENT to the detector congruence (B1') whenever p does not
    divide A(1) m (truth table: 2790 rows, 0 disagreements). The law and the detector of HYP-9310(c) are one statement.
  - Evidence, recalibrated. Script, output and HYP were committed together (7f331c46e5), so the docstring
    "pre-registration" is not verifiable from git. By the equivalence, the test at the 33 detector survivors was implied
    for 7 <= p <= 200 by the detector's own confirmation data; only the primes 211..300 were new.
  - All 33 CM parameters pass at every prime 7 <= p <= 300 with p not dividing 6m (about 59 primes each): 0 failures.
    The S43b referee extended this to p <= 400 (A1 600, A2 750, A3 749, A4 375 prime-instances, 0 failures).
  - Generic m (2, 3, 5, 7, 10, 100) fail at 50-57 of about 59 primes in every family.
  - Every term A(k) m^-k w_k with k < p was p-integral: the poles of w_k at carry indices cancel against the zeros
    of A(k).
  - For A = C(2k,k)^3 the law is HYP-9310 = KLMSY Corollary 3, since its terms with (p-1)/2 < k < p vanish mod p.
  - PROVED ingredient: the p-shift lemma below.
  - Part (O), PASSED. This is the only prediction with a verifiable order: the script was committed in eb53cfd571
    (10:45:22) before the output in 91c22f341d (10:47:18). It was informed by the 8-parameter probe committed in the
    same commit.
    - The Lucas-block congruence holds mod p^3 at every ordinary prime of every CM parameter (908 instances,
      0 failures; reproduced by the referee).
    - At supersingular primes it holds only mod p^2: 16 sporadic exceptions out of 1038, where the chance expectation
      (the sum of 1/p) is 17.6 (corrected from "about 21", S43b m4).
    - Prior context: KLMSY remark that their Theorem 2 congruences "appear to hold modulo p^3" at ordinary primes for
      rational CM lambda. Long-Ramakrishna (Adv. Math. 2016, arXiv:1403.5232, Theorem 3) prove the lambda = 1 case
      mod p^3. (O) is a different, parallel block statement.
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/cm_harmonic_law_four_families_20261018.py (+ .out, .json): the law, 7 <= p <= 300
  - 04-computation/experiments/pshift_lemma_four_families_20261018.py (+ .out, .json): the p-shift lemma, exact, p <= 61
  - 04-computation/experiments/cm_detector_lucas_block_scan_20261018.py (+ .out, .json): the CM parameters
  - 04-computation/experiments/cm_lucas_block_mod_p3_probe_20261018.py (+ .out, .json): mod p^3 probe (split vs inert)
  - 04-computation/experiments/cm_lucas_block_ordinary_mod_p3_20261018.py (+ .out, .json): part (O), pre-registered in eb53cfd571
related: [HYP-9310, THM-4685, THM-4686, Kibelbek-Long-Moss-Sheller-Yuan 2016 Cor 3, Long-Ramakrishna 2016 Thm 3, Sun Conjecture 29 (arXiv:2210.07238), Beukers-Vlasenko Dwork crystals III (arXiv:2105.14841), Kibelbek-Long-Moss-Sheller-Yuan (arXiv:1210.4489), Gukov-Vafa 2002]
---

# HYP-9311 — the CM harmonic law for the four Ramanujan families

## Statement

The families and their weights `w_k = d/dk log A(k)` (Gamma continuation) are:

| family | `A(k)` | `w_k` |
|---|---|---|
| A1 | `C(2k,k)^3` | `6(H_(2k) - H_k)` |
| A2 | `C(2k,k)^2 C(3k,k)` | `3H_(3k) + 2H_(2k) - 5H_k` |
| A3 | `C(2k,k)^2 C(4k,2k)` | `4(H_(4k) - H_k)` |
| A4 | `C(2k,k) C(3k,k) C(6k,3k)` | `3(2H_(6k) - H_(3k) - H_k)` |

These are exactly the harmonic weights in Sun's Conjecture 29(i)-(iv). Let `m` be a CM parameter of the family.
Then for every prime `p ∤ 6m`,
```
(L_p^A)   Σ_(k=0)^(p-1) A(k) m^(-k) w_k  ≡  q_p(m) Σ_(k=0)^(p-1) A(k) m^(-k)   (mod p),     q_p(m) = (m^(p-1) - 1)/p.
```
Each term on the left is `p`-integral.

**CM parameters**, as found by the Lucas-block detector over `0 < |m| <= 10^5`:

| family | parameters |
|---|---|
| A1 | `-512, -64, -8, 1, 16, 64, 256, 4096` |
| A2 | `-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375` |
| A3 | `-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736` |
| A4 | `-32768, -3375, 1728, 8000, 54000` (the rational CM `j`-invariants) |

## p-shift lemma (PROVED, elementary)

Let `j` and `a` be carry-free: every factorial argument `cj` and `ca` is `< p`. Then
```
A(ap + j) ≡ A(a) A(j) (1 + a p w_j)   (mod p^2).
```
Proof:
- Write `(c(ap+j))! = (cap + cj)!`.
- The multiples of `p` contribute `p^(ca) (ca)!`.
- Each complete block contributes `Π_(i=1)^(p-1) (bp + i) ≡ (p-1)! (mod p^2)`, because `H_(p-1) ≡ 0 (mod p)`.
- The final partial block contributes `(cj)! (1 + cap H_(cj))`.
- The ratio `A` is balanced (`Σ ±c = 0`), so the powers of `p` and of `(p-1)!` cancel. What remains is
  `A(a) A(j) Π (1 + c a p H_(cj))^(±1) ≡ A(a) A(j)(1 + ap w_j)`.

Checked exactly: 9843 cases, `7 <= p <= 61`, 0 failures.

Consequently, for carry-free blocks the law is the second-order Lucas factorization
```
Σ_j t_(p+j) ≡ t_1 Σ_j t_j   (mod p^2),     t_k = A(k) m^(-k),
```
as for A1 (HYP-9310; KLMSY Lemma 10).

At `a = 1` the carry-free hypothesis is unnecessary (S43b m3). At carry indices `A(p+j) ≡ A(1)(A(j) + p[A(j) w_j])`
still holds, since:
- the extra factorial factor is `≡ 1 + Σ ±c_i H_(e_i) (mod p)` when `v_p(A(j)) = 1`;
- both sides vanish mod `p^2` when `v_p(A(j)) >= 2`.

So `(L_p^A) ⟺ (B1')` in every family.

## Heuristic: the truncated canonical coordinate is Teichmüller at CM points

Frobenius' method for the hypergeometric equation, with a triple indicial root 0 and `w_0 = 0`, gives the logarithmic
solution
```
y_1 = y_0 log x + Σ_k A(k) w_k x^k,      y_0 = Σ_k A(k) x^k.
```
So `Φ(x) = Σ A(k) w_k x^k / Σ A(k) x^k` is the logarithm of `q/x`, where `q = x e^(Φ(x))` is the canonical (mirror,
Serre-Tate) coordinate.

**Reformulation.** At `x_0 = 1/m` the law says `x_0^(p-1) (1 + p Φ_p(x_0)) ≡ 1 (mod p^2)`, where `Φ_p` is the
truncation at `k < p`. This is the first-order, truncated form of `q_0^(p-1) = 1`.

**Fixed points.** Beukers-Vlasenko (*Dwork crystals III*, Definition 7.14) define the excellent Frobenius lift by
`q^σ = γ^(p-1) q^p` (`γ = 1` here). A fixed point `q_0^σ = q_0` therefore satisfies exactly `q_0^(p-1) = 1`, and their
Corollary 8.2 gives the mod-`p^2` Cartier congruence at fixed points.

The identification "the CM points are these fixed points" is this repository's HEURISTIC, motivated by Deuring
(canonical lifts) and Serre-Tate. Beukers-Vlasenko establish CM fixed points only for the `n = 2` hypercubic family,
where `t(q)` is modular. For `n = 3`, which is A1, they conjecture only `t^σ(1/8) = 1/8`, that is `m = 64` (attribution
corrected, S43b m5).

Two further caveats:
- In the hypercubic coordinate `m = 1/t^2`, so `t ∈ Q_p` only when `m` is a square in `Q_p`.
- The reformulation needs `F_p` to be a unit (ordinary `p`). At supersingular primes `Φ_p` is not `p`-integral,
  although the law still holds.

**Status.** This is a HEURISTIC derivation. Making it rigorous needs Dwork-type congruences for the logarithmic period
(the truncations of `Φ`) at `|x_0|_p = 1`. It explains:
- why the law holds exactly at CM points and fails generically;
- why the correction is a Fermat quotient (the `p`-adic log of `x_0`);
- why the same `6 <-> q_p(m)` pairing appears in Sun's Conjecture 29. Sun's `k`-weighted Ramanujan-type sums are the
  `θ`-derivative of the same structure, with `log|m|` on the archimedean side.

## Part (O): third order at ordinary primes (pre-registered, PASSED)

At every CM parameter, the Lucas-block congruence
```
Σ_(k=p)^(2p-1) t_k ≡ t_1 Σ_(k<p) t_k
```
holds **modulo `p^3` at every ordinary prime** (`F_p = Σ_(k<p) t_k` a `p`-adic unit). At supersingular primes
(`F_p ≡ 0 mod p`) it holds only modulo `p^2`.

**Evidence.**
- Probe: `cm_lucas_block_mod_p3_probe_20261018.out`, at 8 parameters. The mod-`p^3` survivors were exactly the
  primes split in the CM field, apart from one sporadic inert prime.
- Prediction (O)/(S) committed in eb53cfd571 (10:45:22), before the output was written (10:45:43):
  `cm_lucas_block_ordinary_mod_p3_20261018.py`.
- Result for all 33 parameters, `7 <= p <= 300`:
  - 908 ordinary prime-instances with 0 failures mod `p^3`;
  - 1038 supersingular instances, of which 16 hold mod `p^3` sporadically, at `p ∈ {7, 11, 13, 17, 23, 41, 109}`.
    That is consistent with a chance rate of about `1/p`: the expected count is `Σ 1/p = 17.6` (corrected, S43b
    m4).

**Reading.** Expanding `A(p+j)` to second order and `m^(-(p-1)) = (1 + p q_p(m))^(-1)` mod `p^3` turns (O) into a
mod-`p^2` refinement of `(L_p^A)`. The refinement involves second-order weights
`w_j^(2) = Σ ±c^2 H^(2)_(cj)` and the terms `q_p(m) w_j` and `q_p(m)^2`, in the style of Sun's mod-`p^3` conjectures.
In the heuristic above, this is the next order of "the truncated canonical coordinate is Teichmüller". At
supersingular primes there is no unit root, and the expansion stops at first order.

## Context

- **Relation to Sun.** Sun's Conjecture 29 concerns Ramanujan-type points and the weights `(ak + b) w_k + const`.
  `(L_p^A)` has no polynomial weight and holds at every CM point, including the divergent `m = 1` for A1 and
  non-Ramanujan points such as `m = 108` for A2.
- **ANALOGY (Gukov-Vafa 2002; precision per S43b n9).** A `T^2` sigma model is rational exactly when both the complex
  structure `τ` and the Kähler modulus `ρ` are CM points of the same imaginary quadratic field. The four families are
  periods of elliptic curves and their symmetric squares. So the Lucas-block congruence (and this law) tests whether
  the complex-structure modulus `τ` is CM, which is the necessary half of rationality.

## Falsifiers

- A prime `p ∤ 6m` with `(L_p^A)` false at a listed CM parameter.
- A non-CM parameter (supersingular density 0) satisfying `(L_p^A)` for all large `p`.
