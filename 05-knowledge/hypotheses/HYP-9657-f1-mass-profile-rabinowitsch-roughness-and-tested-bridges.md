---
id: HYP-9657
title: "Representation mass of one-arc joins of one-apex tournaments (the F1 family): (A) the class means of the representation count r(N) modulo an odd prime p favour N = 0 and N = -1 (mod p) and show no Kronecker-symbol dependence, i.e. the local densities are of additive-divisor (sum of two products) type in N and N + 1, not of class-number type; (B) consequently integers N with N(N + 1) free of small odd primes - in particular the Heegner numbers d >= 7, for which d and (d + 1)/4 are both prime - carry the least mass and are over-represented among the misses of F1-type families (r(7) = r(43) = r(67) = r(163) = 0); (C) bridges tested and REFUTED: H mod 4 or 8 is not a function of the F_2 Redei / Monsky data, and a positive Monsky nullity does not force 3 | H (first failure: order 5, H = 13)"
status: >
  (A), (B): NUMERICAL (in-sample: odd N in [2^18, 2^21), script f1_mass.py, output f1_mass_insample_21.out), with
  three PRE-REGISTERED out-of-sample tests Q1-Q3 on [2^21, 2^23) recorded below BEFORE the run (pre-registration
  commit on origin/main preceding the run; all three CONFIRMED, section "Outcomes"). The mechanism in (A) is a heuristic reading of the two
  sub-families s = 0 and s >= 1; no asymptotic formula is proved. (C): FINITE-EXACT refutations (cycle_gas.py, every
  labelled tournament of orders 3-6). The Chan-Pak conjecture (every odd number other than 7 and 21 is a
  Hamiltonian-path count) remains OPEN; nothing here is evidence for or against it beyond what the cited canon
  already records.
session: opus-2026-11-10-S90
source: 05-knowledge/results/s90_rigidity_redei_20261110/README.md
scripts:
  - 04-computation/experiments/s90_rigidity_redei_20261110/f1_mass.py (+ f1_mass_insample_21.out, f1_mass_prereg_23.out)
  - 04-computation/experiments/s90_rigidity_redei_20261110/cycle_gas.py (+ .out)
related:
  - 01-canon/theorems/THM-5124-*, THM-5125-* (this session)
  - 01-canon/theorems/THM-4875-* (D) and THM-4885-* (the F1 family and its mean mass (pi^2/12) ln N), THM-4888-* (mod-3 law)
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-4094-*, THM-4884-* (the 29 values missed by the O/Z/J readout)
  - 05-knowledge/results/bsd_selmer_growth_20261013a.md (Monsky matrix of the reciprocity tournament)
  - 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md (table of largest-anomaly problems)
---

# HYP-9657 - where the mass of a join family is small, and what that says about Heegner numbers

## Setting

For h >= 0 let psi(0) = 0, psi(2h) = psi(h) + h, psi(2h + 1) = psi(h) + floor(h/2) + 1; the one-apex tournament with
reduced word of value n = 2h has 2 psi(h) + 1 Hamiltonian paths (a restatement of Session 16's weighted-inversion
formula and THM-4876 (6); see the note). The F1 family (THM-4875 (D), THM-4885) is, in these coordinates,

    N = (2 psi(h1) + 1)(2 psi(h2) + 1) + 2^s (2 h1 + 1)(2 h2 + 1) - [s = 0],      h1, h2, s >= 0,

and r(N) is the number of triples (h1, h2, s) with value N. THM-4885 gives the mean of r; this hypothesis concerns
its arithmetic fluctuation.

## Statements

**(A) Additive-divisor local densities (NUMERICAL).** For every odd prime p, the mean of r over the residue class
N = c (mod p) is largest for c = 0 and c = -1, and the remaining p - 2 classes are equal up to sampling error. In
particular the two classes with Kronecker symbol (-N/p) = +1 and = -1 are not distinguished beyond that. In-sample
(odd N in [2^18, 2^21), overall mean 11.03): p = 3: 12.12 (c = 0), 11.49 (c = -1), 9.48 (c = 1); p = 5: 12.37, 11.87,
others 10.30; p = 7: 12.01, 12.05, others 10.48-10.91; p = 11: 11.99, 11.78, others 10.76-10.93; p = 17: 12.56, 12.00,
others 10.80-10.94; p = 13, 19, 23: the same order with a smaller excess (about 0.2).
Reading: the triples with s >= 1 solve N = x1 x2 + 2^s y1 y2 and those with s = 0 solve N + 1 = x1 x2 + y1 y2, with
(x_i, y_i) on the graph of the one-apex count; a sum of two products has about 1/p more solutions in the class 0
(mod p). The known mod-3 bias of the misses (THM-4875 (D), THM-4888) is the case p = 3 of this pattern
(class N = 1 mod 3 is the only class that is neither 0 nor -1).

**(B) Rabinowitsch roughness (NUMERICAL).** Call N rough if N(N + 1) has no odd prime factor below a bound. Rough
N carry the least mass: in-sample, N prime with (N + 1)/4 prime has mean r = 7.90 (ratio 0.72) and miss rate 3.5%
(ratio 2.6); N(N + 1) free of odd primes <= 29 has mean 7.98 and miss rate 3.6%; 15 | N(N + 1) has mean 12.83 and
miss rate 0.59%. A Heegner number d >= 7 is prime and (d + 1)/4 is prime (the value at 0 of Euler's polynomial
x^2 + x + (d + 1)/4; classical, CITED), so the Heegner numbers are rough in this sense, and
r(7) = r(43) = r(67) = r(163) = 0 while r(11) = 2, r(19) = 4. Among the 30 integers N = 4q - 1 < 3000 with N and q
prime, eight have r = 0: 7, 43, 67, 163, 211, 547, 1171, 1723.
This is a statement about F1 (a density-one-at-best family); it does not single out 7 and 21 in the full spectrum,
whose exclusion is the small-order argument of THM-1370, and 21 = 3 * 7 is not rough.

**(C) Bridges tested and refuted (FINITE-EXACT).** For a tournament with adjacency matrix A let
M = A + diag(row parities) (the Monsky / Redei-type matrix of the reciprocity tournament, convention of
bsd_selmer_growth_20261013a.md).
1. H mod 4 is not a function of (corank of the F_2-Laplacian, nullity(M^2 + M + I), nullity(I + A)): it is on
   orders 3 and 4, fails on 3 of 6 classes at order 5 and on 6 of 16 at order 6.
2. "nullity(M^2 + M + I) > 0 implies 3 | H" holds on orders 3 and 4 and fails at order 5 (out-neighbour masks
   [18, 8, 3, 5, 14], H = 13); with the pure adjacency matrix it holds through order 5 and fails at order 6
   (masks [32, 1, 3, 7, 15, 30], H = 17).
So beyond Redei's parity theorem no congruence for H is read off these F_2 invariants at small orders; the exact
relations that do hold are the specialisation identities of THM-5125.

## Pre-registered tests (recorded before the run)

Universe: odd N in [2^21, 2^23), r computed by `python3 f1_mass.py 23 21 prereg`.
- **Q1.** For each p in {3, 5, 7, 11, 17} the two residue classes with the largest mean of r are 0 and -1.
- **Q2.** On N prime with (N + 1)/4 prime: mean r at most 0.80 of the overall mean, and miss rate at least 1.8 times
  the overall miss rate.
- **Q3.** On 15 | N(N + 1): mean r at least 1.12 of the overall mean, and miss rate at most 0.6 of the overall.
The in-sample analogues were run first (same script, window [2^18, 2^21)): Q1 holds for all eight primes tested;
Q2 ratios 0.717 and 2.56; Q3 ratios 1.163 and 0.430.

## Outcomes

Run after the pre-registration commit (`f1_mass_prereg_23.out`; odd N in [2^21, 2^23), 3,145,728 integers, overall
mean r = 12.19, miss rate 0.82%).
- **Q1 CONFIRMED.** Classes 0 and -1 are the top two for p = 3, 5, 7, 11, 17, and also for 13, 19, 23
  (p = 3: 13.31, 12.93 against 10.33; p = 5: 13.63, 13.09 against 11.41; p = 17: 13.75, 13.22 against 12.02).
- **Q2 CONFIRMED.** N prime with (N + 1)/4 prime (9,636 integers): mean ratio 0.715, miss-rate ratio 3.31.
- **Q3 CONFIRMED.** 15 | N(N + 1): mean ratio 1.169, miss-rate ratio 0.361.
The mean ratios are stable between the two windows (0.717 / 0.715 and 1.163 / 1.169) while the miss-rate contrast
widens (2.56 to 3.31 and 0.43 to 0.36), as expected when the overall mass grows like ln N.

## Non-consequences

No completeness or density statement follows. The profile is specific to F1; the O/Z/J readout's 29 missed values
(THM-4884) contain 7, 43, 67, 163 as well, but that list was not re-analysed here and several of its members are
not rough. The comparison with Ramanujan's ternary form x^2 + y^2 + 10 z^2 (odd exceptions 3, 7, 21, 31, 33, 43, 67,
79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719; Ono-Soundararajan 1997 under GRH) is typed ANALOGY in the
note: there the mass is a class number plus a cusp-form coefficient, here no Kronecker-symbol dependence is seen.
