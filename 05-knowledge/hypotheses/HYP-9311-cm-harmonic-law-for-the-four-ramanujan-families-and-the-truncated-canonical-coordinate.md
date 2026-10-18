---
id: HYP-9311
title: "CM harmonic law for the four Ramanujan families: for A(k) in {C(2k,k)^3, C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k) C(3k,k) C(6k,3k)} with Sun's weights w_k = d/dk log A(k), at every CM parameter m and every prime p not dividing 6m, sum_{k<p} A(k) m^-k w_k == q_p(m) sum_{k<p} A(k) m^-k (mod p). It is the constant-term analogue of Z.-W. Sun's Conjecture 29. Heuristically it is the first-order truncation of the fixed-point equation q_0^(p-1) = 1 of the excellent Frobenius lift in the canonical (mirror / Serre-Tate) coordinate, whose logarithm has power-series part sum A(k) w_k x^k / sum A(k) x^k"
status: >
  OPEN, strongly supported.
  - Pre-registered prediction PASSED (written in the script docstring before the run).
  - All 33 CM parameters found by the p-adic CM detector (HYP-9310(c), |m| <= 10^5, four families) pass at every
    prime 7 <= p <= 300 with p not dividing 6m (about 59 primes each): 0 failures.
  - Generic m (2, 3, 5, 7, 10, 100) fail at 50-57 of about 59 primes in every family.
  - Every term A(k) m^-k w_k with k < p was p-integral: the poles of w_k at carry indices cancel against the zeros
    of A(k).
  - For A = C(2k,k)^3 the law is HYP-9310, since its terms with (p-1)/2 < k < p vanish mod p.
  - PROVED ingredient: the p-shift lemma below.
  - Part (O), PASSED (pre-registered in commit eb53cfd571 before the run): the Lucas-block congruence holds mod p^3
    at every ordinary prime of every CM parameter (908 instances, 0 failures) and only mod p^2 at supersingular primes
    (16 sporadic exceptions out of 1038).
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/cm_harmonic_law_four_families_20261018.py (+ .out, .json): the law, 7 <= p <= 300
  - 04-computation/experiments/pshift_lemma_four_families_20261018.py (+ .out, .json): the p-shift lemma, exact, p <= 61
  - 04-computation/experiments/cm_detector_lucas_block_scan_20261018.py (+ .out, .json): the CM parameters
  - 04-computation/experiments/cm_lucas_block_mod_p3_probe_20261018.py (+ .out, .json): mod p^3 probe (split vs inert)
  - 04-computation/experiments/cm_lucas_block_ordinary_mod_p3_20261018.py (+ .out, .json): part (O), pre-registered in eb53cfd571
related: [HYP-9310, THM-4685, Sun Conjecture 29 (arXiv:2210.07238), Beukers-Vlasenko Dwork crystals III (arXiv:2105.14841), Kibelbek-Long-Moss-Sheller-Yuan (arXiv:1210.4489), Gukov-Vafa 2002]
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
as for A1 (HYP-9310). The carry indices contribute through the poles of `w_k`; for A2-A4 they are part of the
full-period statement.

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
`q^σ = γ^(p-1) q^p` (`γ = 1` here). A fixed point `q_0^σ = q_0` therefore satisfies exactly `q_0^(p-1) = 1`. The CM
points are these fixed points: by Deuring, a CM curve is the canonical lift of its ordinary reduction, and its
Serre-Tate coordinate is a root of unity. Their Corollary 8.2 gives the mod-`p^2` Cartier congruence at such points.

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
    That is consistent with a chance rate of about `1/p`, since about 21 would be expected.

**Reading.** Expanding `A(p+j)` to second order and `m^(-(p-1)) = (1 + p q_p(m))^(-1)` mod `p^3` turns (O) into a
mod-`p^2` refinement of `(L_p^A)`. The refinement involves second-order weights
`w_j^(2) = Σ ±c^2 H^(2)_(cj)` and the terms `q_p(m) w_j` and `q_p(m)^2`, in the style of Sun's mod-`p^3` conjectures.
In the heuristic above, this is the next order of "the truncated canonical coordinate is Teichmüller". At
supersingular primes there is no unit root, and the expansion stops at first order.

## Context

- **Relation to Sun.** Sun's Conjecture 29 concerns Ramanujan-type points and the weights `(ak + b) w_k + const`.
  `(L_p^A)` has no polynomial weight and holds at every CM point, including the divergent `m = 1` for A1 and
  non-Ramanujan points such as `m = 108` for A2.
- **ANALOGY (Gukov-Vafa 2002).** A sigma model on an elliptic curve or a K3 surface is rational exactly at CM moduli.
  The four families are periods of elliptic curves and their symmetric squares. So the Lucas-block congruence (and
  this law) is a `p`-adic test for the moduli of rational sigma models.

## Falsifiers

- A prime `p ∤ 6m` with `(L_p^A)` false at a listed CM parameter.
- A non-CM parameter (supersingular density 0) satisfying `(L_p^A)` for all large `p`.
