---
id: THM-4885
title: "Representation mass of the one-apex join families. (A) For integers a, c >= 1 and b >= 0, the readout a g(v) - b e(v) + c (1 + v) of the one-apex orbit takes at most X on (X/a) ln(1 + a/c) + O(sqrt X) values of v. (B) The one-apex x one-apex joins (1 + u)(1 + v) + g(u) g(v) - e(u) e(v), u <= v, number (pi^2/24) X ln X + kappa X + o(X) up to X. (C) So a one-apex x bank family has mean multiplicity M = sum_y (2/y_h) ln(1 + y_h/y_R) per odd number (48.1976 for the order-<=7 bank, 211.236 for order <= 8), and the one-apex x one-apex family has (pi^2/12) ln N + O(1), gaining (pi^2/12) ln 2 = 0.5701 per doubling. (D) FINITE-EXACT and pre-registered (HYP-9478): the mass laws hold to four decimals at 2^28 to 2^30, and F2 alone (fixed order-<=7 bank) attains every odd N <= 2^30 except 7 and 21. Chan-Pak OPEN"
status: >
  PROVED: (A), (B), (C) (elementary; equidistribution of the one-apex fraction rho(v) in dyadic blocks).
  FINITE-EXACT: (D), pre-registered as HYP-9478 P1-P4 (commit 09b7a56a98) and CONFIRMED.
  INHERITED: the one-arc join law (session29 ADDENDUM; THM-4866 (B)); the one-apex triple and closed form
  (session179 spectrum/REPORT.md (1); session124's f; Session191 one_apex/PROOF.md; THM-4876 (6); THM-4875 (D) and its
  audit); the families F1, F2 and their coverage (THM-4875 (D), Session 193, Session 203's corrected recertification).
  NOT CLAIMED: Chan-Pak completeness, which is OPEN. A mass formula is a first moment; it implies no pointwise
  coverage.
  AUDIT: pending.
session: opus-2026-11-06-S82
source: 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md
related:
  - THM-4875
  - THM-4876
  - THM-4866
  - THM-1370
  - HYP-9478
  - HYP-9468
  - session203/coverage-recertification/README.md
  - session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md
verification:
  - 04-computation/experiments/chanpak_join_multiplicity_20261106.c
  - 04-computation/experiments/chanpak_join_multiplicity_20261106.out
  - 04-computation/experiments/chanpak_f1f2_cover34_20261106.c
---

# THM-4885 - the representation mass of the one-apex join families

**Notation.**
- g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), and e(v) = [v = 0 or v even]. The rooted triple of the
  one-apex tournament T(v) at its apex is (g(v), e(v), 1 + v).
- For v >= 1 with top bit n, let rho(v) = sum over the zero bits i < n of 2^-(i+1).

## (A) One readout

Let a, c >= 1 and b >= 0 be integers and Phi(v) = a g(v) - b e(v) + c (1 + v). Then

    #{v >= 0 : Phi(v) <= X} = (X / a) ln(1 + a / c) + O_(a,b,c)(sqrt X).

The constant in the O is uniform in the form needed for (B): O(sqrt(X/c) + (a/c) log X).

## (B) One-apex x one-apex joins

    #{0 <= u <= v : (1 + u)(1 + v) + g(u) g(v) - e(u) e(v) <= X} = (pi^2 / 24) X ln X + kappa X + o(X)

for a constant kappa. The constant pi^2/12 is the integral of ln(1 + t)/t over [0, 1].

## (C) Mean multiplicities

All the values in (A) and (B) are Hamiltonian path counts (THM-4875 (D)), hence odd. Therefore:
- **One-apex x bank.** For a finite set Y of rooted triples y = (y_h, y_e, y_R), the representations
  y_h g(v) - y_e e(v) + y_R (1 + v) of the odd N <= X have mean M_Y + O(X^(-1/2)), with
  M_Y = sum_y (2/y_h) ln(1 + y_h/y_R). The bank of distinct rooted triples of order <= m gives:

  | m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
  |---|---|---|---|---|---|---|---|---|
  | M | 1.386 | 2.197 | 3.106 | 5.134 | 8.380 | 16.954 | 48.198 | 211.236 |

- **One-apex x one-apex.** The pairs u <= v have mean (pi^2/12) ln X + O(1) over the odd N <= X. The mean over the odd
  N in [2^k, 2^(k+1)) is (pi^2/12)(k + 2) ln 2 + 2 kappa + o(1); it gains (pi^2/12) ln 2 = 0.570091 per doubling.

## (D) Finite checks (FINITE-EXACT; pre-registered HYP-9478 P1-P4, CONFIRMED)

These come from `chanpak_join_multiplicity_20261106` to 2^30.
- **One-apex x bank (order <= 7).** The block mean of r2 is 48.1976 at k = 24..29 (M = 48.19761).
- **One-apex x one-apex.** The block increments of r1 are 0.5701 at k = 23..29.
- **F2 alone with the fixed order-<=7 bank** attains every odd N <= 2^30 except 7 and 21. The block minimum of r2 is
  4 at k = 29.

**Proofs.** Source note, sections 2 and 3.
