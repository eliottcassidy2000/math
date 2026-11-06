---
id: THM-4885
title: "Representation mass of the one-apex join families. (A) For integers a, c >= 1 and b >= 0, the readout a g(v) - b e(v) + c (1 + v) of the one-apex orbit takes at most X on (X/a) ln(1 + a/c) + O(sqrt X) values of v (extends Session 16/26's (1,1,1) law). (B) The one-apex x one-apex joins (1 + u)(1 + v) + g(u) g(v) - e(u) e(v), u <= v, number (pi^2/24) X ln X + kappa X + o(X) up to X, kappa = -0.2388732. (C) So a one-apex x bank family has mean multiplicity M = sum_y (2/y_h) ln(1 + y_h/y_R) per odd number (48.1976 for the order-<=7 bank, 211.236 for order <= 8; M_m >= 2 (ln 2)^2 (m - 1)), and the one-apex x one-apex family has (pi^2/12) ln N + O(1), gaining (pi^2/12) ln 2 = 0.5701 per doubling. (D) FINITE-EXACT, pre-registered (HYP-9478): the fixed order-<=7 bank alone attains every odd N <= 2^30 except 7 and 21. Chan-Pak OPEN"
status: >
  PROVED:
  (A) extends the audited one-apex word law W(X) = X ln 2 + O(sqrt X) (Session 16
  feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md (1), 2026-10-20; Session 26 (4)), its case
  (a, b, c) = (1, 1, 1), by the same binary-product coupling;
  (B) and (C) are new, with kappa and M_m >= 2 (ln 2)^2 (m - 1) from the independent audit (U1, U2).
  FINITE-EXACT: (D), pre-registered as HYP-9478 P1-P4 (commit 09b7a56a98), CONFIRMED and reproduced independently.
  P1 and P2 test the implementation of (A)-(C); P3 and P4 are the genuine out-of-sample tests.
  KNOWN: the one-apex mass 2 ln 2; H >= 1 + 2 oc(T) with equality iff intercyclic, and the intercyclic spectrum
  (Grinberg-Stanley; Session 16 intercyclic_rays.md).
  INHERITED: the one-arc join law (session29 ADDENDUM; THM-4866 (B)); the one-apex triple and closed form (session179
  spectrum/REPORT.md (1); session124's f; Session191 one_apex/PROOF.md; Session 16 proofs.md section 2; THM-4876 (6);
  THM-4875 (D) and its audit); the families F1, F2 and their coverage (THM-4875 (D); Session 193; Session 203's
  corrected recertification).
  NOT CLAIMED: Chan-Pak completeness, which is OPEN. A mass formula is a first moment; it implies no pointwise
  coverage. The circle-method programme of the source note, section 7, is HEURISTIC.
  AUDIT: PASS WITH CORRECTIONS (3 majors, 16 minors, 6 upgrades; all applied; MISTAKE-725);
  04-computation/experiments/one_apex_join_mass_20261106_audit_S82/AUDIT.md.
session: opus-2026-11-06-S82
source: 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md
related:
  - THM-4875
  - THM-4876
  - THM-4878
  - THM-4881
  - THM-4866
  - THM-1370
  - HYP-9478
  - HYP-9468
  - HYP-9469
  - HYP-9457
  - HYP-9474
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - session203/coverage-recertification/README.md
verification:
  - 04-computation/experiments/chanpak_join_multiplicity_20261106.c
  - 04-computation/experiments/chanpak_join_multiplicity_20261106.out
  - 04-computation/experiments/chanpak_f1f2_cover34_20261106.c
  - 04-computation/experiments/chanpak_f1f2_cover35_bucketed_20261106.c
  - 04-computation/experiments/chanpak_mass_exploratory_20261106/
  - 04-computation/experiments/one_apex_join_mass_20261106_audit_S82/AUDIT.md
---

# THM-4885 - the representation mass of the one-apex join families

**Notation.**
- g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), and e(v) = [v = 0 or v even]. The rooted triple of the
  one-apex tournament T(v) at its apex is (g(v), e(v), 1 + v).
- For v >= 1 with top bit n, let rho(v) = sum over the zero bits i < n of 2^-(i+1). Then
  1 <= g(v) - v rho(v) <= 1 + #10(v).

## (A) One readout (PROVED; extends Session 16/26)

Let a, c >= 1 and b >= 0 be integers and Phi(v) = a g(v) - b e(v) + c (1 + v). Then

    #{v >= 0 : Phi(v) <= X} = (X / a) ln(1 + a / c) + O(sqrt(X / c) + (a / c) log X + b / c + 1).

The case (a, b, c) = (1, 1, 1), where Phi(w) = g(2w), is Session 16's W(X) = X ln 2 + O(sqrt X).

## (B) One-apex x one-apex joins (PROVED, new)

    #{0 <= u <= v : (1 + u)(1 + v) + g(u) g(v) - e(u) e(v) <= X} = (pi^2 / 24) X ln X + kappa X + o(X),

with kappa = kappa_1 - c' = 0.24748886 - 0.48636207 = -0.2388732, where:
- kappa_1 = lim_U [sum over u < U of phi(t_u)/(1 + u) - (pi^2/12) ln U], with phi(t) = ln(1 + t)/t and
  t_u = g(u)/(1 + u);
- c' = integral over 0 < mu < lambda < 1 of P(lambda mu (1 + r s) <= 1), with r, s independent uniform on [0, 1].

The constant pi^2/12 is the integral of ln(1 + t)/t over [0, 1].

## (C) Mean multiplicities (PROVED, new)

For realisable triples every value is a Hamiltonian path count (THM-4875 (D)), hence odd. Therefore:
- **One-apex x bank.** For a fixed finite set Y of rooted triples y = (y_h, y_e, y_R), the representations
  y_h g(v) - y_e e(v) + y_R (1 + v) of the odd N <= X have mean M_Y + O_Y(X^(-1/2)), where
  M_Y = sum_y (2/y_h) ln(1 + y_h/y_R). The bank of distinct rooted triples of order <= m gives:

  | m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
  |---|---|---|---|---|---|---|---|---|
  | M | 1.386 | 2.197 | 3.106 | 5.134 | 8.380 | 16.954 | 48.198 | 211.236 |

  The m = 1 value 2 ln 2 is KNOWN (Session 16/26). In general M_m >= 2 (ln 2)^2 (m - 1), because the order-<=m bank
  contains the one-apex triples with v < 2^(m-1).
- **One-apex x one-apex.** The pairs u <= v have mean (pi^2/12) ln X + O(1) over the odd N <= X. Over the odd N in
  [2^k, 2^(k+1)) the mean is (pi^2/12)(k + 2) ln 2 + 2 kappa + o(1), gaining (pi^2/12) ln 2 = 0.570091 per doubling.

## (D) Finite checks (FINITE-EXACT; pre-registered HYP-9478 P1-P4, CONFIRMED)

These come from `chanpak_join_multiplicity_20261106` to 2^30; the audit's independent sieve reproduces every line.
- **One-apex x bank (order <= 7).** The block mean of r2 is 48.1976 at k = 24..29 (M = 48.19761).
- **One-apex x one-apex.** The block means of r1 are 16.6250 and 17.1951 at k = 28, 29. (B) predicts 16.624975 and
  17.195065.
- **F2 alone with the fixed order-<=7 bank** attains every odd N <= 2^30 except 7 and 21. The block minimum of r2 is 4
  at k = 29.

**Proofs.** Source note, sections 2-4.
