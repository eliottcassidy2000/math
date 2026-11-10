---
id: HYP-9657
title: "Representation mass of one-arc joins of one-apex tournaments (the F1 family) and tested bridges. (A) For the primes p = 3, 5, 7, 11, 13, 17, 19, 23 (and 41, 43) the class means of the representation count r(N) are largest at N = 0 and N = -1 (mod p), the first excess coming entirely from the padded sheet s >= 1 and the second from the sheet s = 0; the excess is absent at p = 29, 31, 37, 53, 61, 73, 89, 127 and is not of size 1/p. (B) Integers N with N(N + 1) free of small odd primes carry the least F1 mass; the Heegner numbers 7, 43, 67, 163 (d and (d + 1)/4 prime) are such integers and have r = 0. (C) Bridges tested and REFUTED: H mod 4 or 8 is not a function of the F_2 Redei / Monsky data; a positive Monsky nullity does not force 3 | H; the interaction of two active vertices is not bounded below by half the larger increment (fails at order 8). (D) Program statement: a circle-method route to every large N needs a family with more than X^2 members; k apexes over a chain pass from k = 6."
status: >
  (A), (B): NUMERICAL (in-sample: odd N in [2^18, 2^21), f1_mass_insample_21.out), with three PRE-REGISTERED
  out-of-sample tests Q1-Q3 on [2^21, 2^23) recorded BEFORE the run (commit cf65422aaa precedes the output) and all
  CONFIRMED, also by the independent audit with its own code. The statements were RESTRICTED after the audit: the
  first version said "for every odd prime", "the remaining classes are equal up to sampling error" and "no
  Kronecker-symbol dependence"; all three are false on the same windows (audit M3) and were withdrawn (MISTAKE-752).
  KNOWN before this session (S82 note one_apex_join_representation_mass_pi2_over_12_20261106.md section 5, S78 note):
  a bounded local factor at small odd primes with the class N = 0 over-represented, the mod-3 profile of r and of
  the F1 misses, reduced mass on primes, and the rotating mod-7 bias of the one-apex image (THM-4878 (5), HYP-9480
  (A)). New here: the class N = -1 at p >= 5, the attribution of the two excesses to the two sheets, the
  Rabinowitsch stratum. (B'): EXPLORATORY, in-sample. (C): FINITE-EXACT refutations (cycle_gas.py, every labelled
  tournament of orders 3-6; interaction_gap_refuted.py). (D): HEURISTIC threshold and an elementary bound, no
  theorem. AUDITED (05-knowledge/results/s90_rigidity_redei_20261110/audit/AUDIT.md): PASS WITH CORRECTIONS, all
  applied. The Chan-Pak conjecture remains OPEN; nothing here is evidence for or against it beyond what the cited
  canon already records.
session: opus-2026-11-10-S90
source: 05-knowledge/results/s90_rigidity_redei_20261110/README.md
scripts:
  - 04-computation/experiments/s90_rigidity_redei_20261110/f1_mass.py (+ f1_mass_insample_21.out, f1_mass_prereg_23.out), f1_formula_check.py (+ .out)
  - 04-computation/experiments/s90_rigidity_redei_20261110/cycle_gas.py, ozj_mass.py, interaction_gap_refuted.py, ear_weyl_identity.py, apex_sums.py (+ .out each)
  - 05-knowledge/results/s90_rigidity_redei_20261110/audit/audit_f1_classes.py, audit_f1_split.py, audit_gap_witness.py (+ outputs; the audit's own code)
related:
  - 01-canon/theorems/THM-5124-*, THM-5125-* (this session)
  - 05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md (S82: local factors, section 5; circle-method discussion), 01-canon/theorems/THM-4875-* (D), THM-4885-* (the F1 family and its mean mass), THM-4878-* (5), 05-knowledge/hypotheses/HYP-9480-*
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md, THM-4094-*, THM-4884-* (the 29 values missed by the O/Z/J readout)
  - 05-knowledge/results/bsd_selmer_growth_20261013a.md (Monsky matrix of the reciprocity tournament)
  - 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md (table of largest-anomaly problems)
---

# HYP-9657 - where the mass of a join family is small, and what that says about Heegner numbers

## Setting

For h >= 0 let psi(0) = 0, psi(2h) = psi(h) + h, psi(2h + 1) = psi(h) + floor(h/2) + 1; the one-apex tournament with
reduced word of value n = 2h has 2 psi(h) + 1 Hamiltonian paths (a restatement of Session 16's weighted-inversion
formula and THM-4876 (6); see the note). The F1 family (THM-4875 (D), THM-4885) is, in these coordinates,

    N = (2 psi(h1) + 1)(2 psi(h2) + 1) + 2^s (2 h1 + 1)(2 h2 + 1) - [s = 0],      h1, h2, s >= 0

(checked against 961 explicit joins, f1_formula_check.py), and r(N) is the number of ORDERED triples (h1, h2, s)
with value N. THM-4885 counts unordered pairs u <= v instead (mean 13.01 against 12.19 for r on [2^21, 2^23)); the
zero sets agree, the individual values differ (r(11) = 2 and r(19) = 4 are 1 and 3 as pairs), and with the pair
count the pattern of (A) already fails at p = 13 (audit minor 12). Everything below is for r as defined here.

## Statements

**(A) Two sheets, two favoured classes (NUMERICAL, restricted to the primes listed).**
- For p in {3, 5, 7, 11, 13, 17, 19, 23} the two residue classes with the largest mean of r are N = 0 and N = -1
  (mod p), on both windows; the audit found the same at p = 41, 43 and found it ABSENT at p = 29, 31, 37, 53, 61, 73,
  89, 127 (at p = 73 and 127 the two classes are below average). In-sample values (odd N in [2^18, 2^21), overall
  mean 11.03): p = 3: 12.12 (c = 0), 11.49 (c = -1), 9.48 (c = 1); p = 5: 12.37, 11.87, others 10.30; p = 7: 12.01,
  12.05, others 10.48-10.91; p = 11: 11.99, 11.78, others 10.76-10.93; p = 17: 12.56, 12.00, others 10.80-10.94;
  p = 13, 19, 23: the same order with an excess of about 0.2.
- *Sheet decomposition (audit, audit_f1_split_23_21.out).* The triples with s >= 1 solve N = x1 x2 + 2^s y1 y2 and
  those with s = 0 solve N + 1 = x1 x2 + y1 y2, with (x_i, y_i) on the graph of the one-apex count. The excess of the
  class N = 0 comes entirely from the sheet s >= 1 and the excess of the class N = -1 entirely from the sheet s = 0;
  the two excesses are nearly equal at every prime: relative to the sheet mean, 0.56 / 0.49 (p = 3), 0.35 / 0.33 (5),
  0.23 / 0.28 (7), 0.18 / 0.18 (11), 0.036 / 0.039 (13), 0.26 / 0.23 (17), 0.045 / 0.051 (19), 0.043 / 0.048 (23),
  about 0 (29, 31, 37), 0.030 / 0.035 (41), 0.098 / 0.099 (43).
- *What is NOT true* (withdrawn after the audit). The size of the excess is not 1/p and is not monotone in p (the
  generic count of solutions of a sum of two products would give 0.375, 0.208, 0.146, 0.092, 0.077, 0.059, ...).
  The other classes are not equal at p = 7 (10.48-10.91; the mod-7 bias of the one-apex image is KNOWN: THM-4878 (5),
  HYP-9480 (A), S82 note section 7). A split by the Kronecker symbol (-N/7) is visible (-0.25 in sample, -0.44 out of
  sample) but cannot be separated from the base-2 structure, because the squares mod 7 are the powers of 2; at
  p = 5 no split is seen and at the other primes it is not robust. So no class-number-type (L-function) structure
  is established or excluded by these data; what is established is the sheet mechanism at small primes.
- The mod-3 profile of r (class 1 lightest) was recorded in the S82 note section 5 (r by N mod 3 = 18.1, 12.5, 14.1
  in its normalisation) and the residue profile of the F1 misses in the S78 note (THM-4875 (D)).

**(B) Rabinowitsch roughness (NUMERICAL).** Call N rough if N(N + 1) has no odd prime factor below a bound. Rough
N carry the least mass: in-sample, N prime with (N + 1)/4 prime has mean r = 7.90 (ratio 0.72) and miss rate 3.5%
(ratio 2.6); N(N + 1) free of odd primes <= 29 has mean 7.98 and miss rate 3.6%; 15 | N(N + 1) has mean 12.83 and
miss rate 0.59%. (Reduced mass on primes was KNOWN: S82 note, 0.796 for F2.) A Heegner number d >= 7 is prime and
(d + 1)/4 is prime (the value at 0 of Euler's polynomial x^2 + x + (d + 1)/4; classical, CITED). Of the Heegner
numbers, 7, 43, 67 and 163 are rough in this sense (11 and 19 are not: 3 | 12, 5 | 20), and
r(7) = r(43) = r(67) = r(163) = 0 while r(11) = 2, r(19) = 4. Among the 30 integers N = 4q - 1 < 3000 with N and q
prime, eight have r = 0: 7, 43, 67, 163, 211, 547, 1171, 1723.
This is a statistical statement about F1 (a density-one-at-best family); it does not single out 7 and 21 in the
full spectrum, whose exclusion is the small-order argument of THM-1370, and 21 = 3 * 7 is not rough.

**(B') The O/Z/J system (EXPLORATORY, in-sample, not pre-registered; ozj_mass.py).** For the three-letter rooted
system S, B, F, with m(N) the number of states of count N reached after a B or F letter: on odd N in [3000, 12000]
the mean of m is 72.1; class means modulo 3 are 101.3 / 55.8 / 59.0 and modulo 5 are 94.8 / 66.9 / 64.4 / 66.6 /
67.6; primes have ratio 0.68 and Rabinowitsch-rough integers 0.65. So there the favoured class is N = 0 (products)
and the class N = -1 is nearly neutral: the role of N + 1 in (A) is specific to F1. The 29 odd N <= 12000 with
m(N) = 0 are the values of THM-4884 (FINITE-EXACT there through 50,001).

**(C) Bridges tested and refuted (FINITE-EXACT).** For a tournament with adjacency matrix A let
M = A + diag(row parities) (the Monsky / Redei-type matrix of the reciprocity tournament, convention of
bsd_selmer_growth_20261013a.md).
1. H mod 4 is not a function of (corank of the F_2-Laplacian, nullity(M^2 + M + I), nullity(I + A)): it is on
   orders 3 and 4, fails on 3 of 6 classes at order 5 and on 6 of 16 at order 6.
2. "nullity(M^2 + M + I) > 0 implies 3 | H" holds on orders 3 and 4 and fails at order 5 (out-neighbour masks
   [18, 8, 3, 5, 14], H = 13, the core of Session245's THM-5117); with the pure adjacency matrix it holds through
   order 5 and fails at order 6 (masks [32, 1, 3, 7, 15, 30], H = 17).
3. A second-order gap for THM-5124. For two vertices active over the rest, with increments d_a, d_b and interaction
   mu > 0, the bound 2 mu >= max(d_a, d_b) holds on every labelled tournament of orders 4-7 (pair {0, 1} over the
   full rest; exhaustive, confirmed by the audit on every pair; equality on 14 profiles, among them
   (c, d_a, d_b, mu) = (1, 4, 2, 2), (1, 4, 4, 2), (5, 4, 4, 2), (11, 8, 8, 4)) and FAILS at order 8: out-neighbour
   masks [246, 252, 248, 177, 160, 64, 152, 32], a = 1, b = 2, counts 15, 29, 29, 49, i.e. (15, 14, 14, 6)
   (interaction_gap_refuted.py, re-counted over permutations). The audit's witnesses at orders 8, 9, 10 have
   mu / max(d) = 3/7, 5/13, 9/25 = (2^j + 1)/(3 * 2^j + 1), j = 1, 2, 3, which suggests the limit 1/3 for the infimum
   (OPEN; no order-8 census was run).
So beyond Redei's parity theorem no congruence for H is read off these F_2 invariants at small orders; the exact
relations that do hold are the identities of THM-5125 (4).

**(D) Program statement (OPEN; HEURISTIC threshold, no theorem).** For a family F of tournaments with counts at
most X put r(N) = #{T in F : H(T) = N} and S(theta) = sum_{T in F} e(theta H(T)). A circle-method proof that every
odd N in the bulk of F is attained needs the contribution of the frequencies away from 0 and 1/2 to be o(|F| / X).
In every family measured the sum has square-root size at generic frequencies (apex_sums.py: at 16 random theta the
median of |S| is about 2 |F|^(1/2) for one apex, n <= 22, and for two apexes, n <= 11), i.e. r fluctuates at least
at Poisson size; with that, a pointwise bound of the required strength is possible only if |F|^(1/2) is at most
about |F| / X, i.e. |F| >= X^2. One-apex, F1 and the S, B, F system have |F| about X^1, X^1 (up to logarithms) and
X^1.58 and are below this threshold; k apexes over a transitive chain have 2^(kn) cut patterns and counts at most
k! (k+1)^n (PROVED, elementary), so |F| >= X^(k / log2(k+1)), which exceeds X^2 from k = 6. Known obstacle: at
rationals with Mersenne-type denominators the one-apex sum decays slowly (rate 0.93 per letter at 1/7, about 0.89 at
1/15, 1/31, 1/63; S82 note), so those points must be treated as major arcs with vanishing local factor.
No bound is proved; the exact differencing identity for one ear cube (ear_weyl_identity.py) is the available tool.
This item records where the exponent count first allows a statement about every large N rather than about almost
all N; it is not a claim that the route works.

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
The audit reproduced Q1-Q3 with its own code and then tested the primes the pre-registration did not name: the
pattern of Q1 is absent at p = 29, 31, 37, 53, 61, 73, 89, 127. The pre-registered statements were narrow enough to
pass; the general sentence written around them was not (see (A), "What is NOT true").

## Non-consequences

No completeness or density statement follows. The profile (A)-(B) is specific to F1; for the O/Z/J readout see
(B'), where only the class N = 0 is favoured. The comparison with Ramanujan's ternary form x^2 + y^2 + 10 z^2 (odd
exceptions 3, 7, 21, 31, 33, 43, 67, 79, 87, 133, 217, 219, 223, 253, 307, 391, 679, 2719; completeness of the list
under GRH: Ono-Soundararajan 1997) is typed ANALOGY in the note: there the mass is a class number plus a cusp-form
coefficient; here the data neither establish nor exclude a Kronecker-symbol dependence.
