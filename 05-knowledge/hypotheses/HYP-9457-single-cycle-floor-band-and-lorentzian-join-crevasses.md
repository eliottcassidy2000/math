---
id: HYP-9457
title: "The single-cycle floor and the crevasses: (a) the cycle ceiling H(T) >= f(n) + 2(hc(T) - 1) for every strong tournament (equality exactly at the Busch floor); (b) every uniquely Hamiltonian tournament is a one-arc join (CLOSED: KNOWN, Shelor 1977); (c) Lorentzian joins H = x_h y_h + x_R y_R - x_e y_e cover every odd number except 7 and 21, with a 3-adic singular series; pre-registered tests at order 10 and for joins of orders <= 7"
status: OPEN for (a), (c) (opus-2026-11-04-S73); (b) CLOSED - KNOWN (Shelor 1977: a vertex of in-degree 1, a single-vertex one-arc cut). Pre-registered P2, P4, P5 CONFIRMED; P6 REFUTED (ratio 1.37); P1 and P3 were sanity checks (audit, MISTAKE-705)
session: opus-2026-11-04-S73
date: 2026-11-04
related: [THM-4866, THM-4744, THM-4745, THM-4655, THM-1370, THM-4863, THM-4865, THM-4862, THM-4846, HYP-9453, HYP-9454, session29 ADDENDUM, polarized_path_calculus_20261103]
note: 05-knowledge/results/single_cycle_floor_crevasses_lorentzian_joins_20261104.md
script: 04-computation/experiments/scf_crevasses_20261104.c
---

# HYP-9457 - the single-cycle floor, the bottom band and the Lorentzian joins

## Notation

- H(T) is the number of Hamiltonian paths of T, and hc(T) the number of directed Hamiltonian cycles. A tournament on
  n >= 3 vertices is strong iff hc >= 1 (Camion).
- f(n) is Busch's strong floor: 3, 5, 9, 15, 25, 45, 75, 125 for n = 3, ..., 10 (THM-1370, THM-4744).
- f_k(n) = min{H : hc = k}.
- The single-cycle band at order n is [f(n), min over k >= 2 of f_k(n)).
- A one-arc join J = J(A, a; B, b) is A => B with the arc a -> b reversed. Its value is
  H(J) = x_h y_h + x_R y_R - x_e y_e, where:
  - x = (h(A), R_a(A), e_a(A)) is the rooted triple, with R_a(A) = sum over X of p(X) s_a(A - X) (klein's root
    response B_a, THM-4865);
  - y = x(B^op, b).

## Exploratory evidence (before this commit)

Program `scf_crevasses_20261104.c`, run on modes census 9, chains 17, joins 5 5 and joins 6 6.

**Single-cycle floor.**
- At every order 3-9, every floor class has exactly one Hamiltonian cycle; the floor-class counts are
  1, 1, 2, 2, 1, 5, 3.
- Every interval chain of order <= 17 (1,601 chains) is uniquely Hamiltonian.

**Cycle ceiling.**
- H >= f + 2(hc - 1) holds for every strong class at orders 3-9.
- Equality holds only at the floor classes.
- The k-cycle floors lie far above the ceiling. For example f_2(n) = 13, 23, 41, 69, 115 at n = 5, ..., 9, and in
  each case f_2(n) < f(n + 1).

**Bottom band.**
- The uniquely Hamiltonian stratum has 1, 1, 3, 8, 21, 55, 144 classes at orders 3-9 (Garey's F(2n - 6)).
- All of these classes are one-arc joins.
- The band [f(n), f_2(n)) consists of uniquely Hamiltonian classes by definition.
- Inside the band the strong spectrum misses exactly:
  - 21 at order 6;
  - 47, 55, 67 at order 8;
  - 77, 79, 83 at order 9.

**Lorentzian joins.**
- The formula is verified by direct counting on 3,000 random joins.
- With both sides of order <= 6 (254 distinct triples), every odd N <= 2,393 except 7 and 21 is a join value.
- Over the windows [32, 2048), every N has at least 4 representations.
- The mean representation count by N mod 3 is 1.53 / 0.70 / 0.77, residue 0 first. Mod 5 deviations are <= 0.13,
  and mod 7 deviations are <= 0.05.

## Pre-registered predictions

Confirmatory runs: `scf_crevasses_20261104 census 10` and `scf_crevasses_20261104 joins 7 7 200000`.

- **(P1) Single-cycle floor at order 10.** f(10) = 125. It is attained by exactly one class (THM-4744's E(10) = 1),
  which has exactly one Hamiltonian cycle.
- **(P2) Cycle ceiling at order 10.** No strong class of order 10 violates H >= f + 2(hc - 1). Equality holds only at
  the floor class.
- **(P3) Bottom band at order 10.**
  - Exactly F(14) = 377 classes have exactly one Hamiltonian cycle. This is Garey's count, a sanity check rather than a
    prediction.
  - Every one of the 377 is a one-arc join (a cut crossed backward by exactly one arc).
- **(P4) Second floor at order 10.** f_2(10) = min{H : hc = 2} lies in [175, 215], below f(11) = 225.
- **(P5) Join coverage.** With both sides of order <= 7, every odd N in [23, 20000] is a Lorentzian join value. 7 and
  21 are never join values (forced by THM-1370).
- **(P6) Singular series.** On the window [25000, 50000), with both sides of order <= 7:
  - the mean of r(N) on N = 0 mod 3, divided by its mean on N != 0 mod 3, lies in [1.6, 3.0];
  - every residue-class mean mod 5 and mod 7 lies within 0.2 of the overall mean.

Every outcome is reported, with the counterexample in the case of P1-P4.

## Results (2026-11-04, after the pre-registration commit 4d5099bc95)

Runs: `census 10` (4.5 min) and `joins 7 7 200000`. Outputs:
- `04-computation/experiments/scf_crevasses_20261104_census.out`;
- `04-computation/experiments/scf_crevasses_20261104_joins.out`.

- **(P1) CONFIRMED.** f(10) = 125, attained by exactly one class, which has exactly one Hamiltonian cycle.
- **(P2) CONFIRMED.** No strong class of order 10 violates the ceiling. The minimum slack is 0, and equality holds only
  at the floor class.
- **(P3) CONFIRMED.** There are 377 uniquely Hamiltonian classes at order 10, and 377 of them are one-arc joins.
- **(P4) CONFIRMED.** f_2(10) = 205, inside [175, 215] and below f(11) = 225. The single-cycle band at order 10 is
  [125, 205).
- **(P5) CONFIRMED, with a larger reach than predicted.** Every odd N in [23, 84705] is a Lorentzian join value with
  both sides of order <= 7. The first gap above 21 is 84707, at the edge of the range. 7 and 21 are never represented.
- **(P6) REFUTED as stated.** On [25000, 50000):
  - the mean of r(N) on N = 0 mod 3, divided by the mean on N != 0 mod 3, is 1.220 / 0.8905 = 1.37, not in
    [1.6, 3.0];
  - the mod 5 deviations (<= 0.04) and mod 7 deviations (<= 0.01) are within the predicted 0.2.

  The 3-adic local factor is real but weaker than at orders <= 5 and <= 6, where the ratio was about 2.0-2.1.

## Open statements after the tests

- **(a) The cycle ceiling.** H(T) >= f(n) + 2(hc(T) - 1) for every strong tournament, with equality exactly at the
  floor. At odd order this is G >= f(n) - 2 for the deck reading of THM-4863. It holds through order 10.
- **(b) CLOSED, KNOWN.** Every uniquely Hamiltonian tournament has a one-arc cut, which can always be taken to be a single
  vertex of in-degree 1 (Shelor 1977's complete insertion construction; session103; operation_closed_synthesis section
  II; klein THM-4862 (b), (f)).
- **(c) Join completeness.** Every odd N other than 7 and 21 is a Lorentzian join value. This holds through 84,705 with
  pieces of order <= 7. As a Goldbach-type statement it would follow from a lower bound on representation counts of the
  kind the comet shows.
