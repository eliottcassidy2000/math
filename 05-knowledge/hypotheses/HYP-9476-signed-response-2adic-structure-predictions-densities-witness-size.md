---
id: HYP-9476
title: "Signed one-apex responses S(a,b) = {a n - b g(n)}: (A) pre-registered predictions P1-P6 for THM-4883 and their outcomes (P2-P5 confirmed; P1 structure confirmed, per-class attainment refuted; P6 partly refuted); (B) every one-side-bounded class of S(a,b) has natural density strictly between 0 and 1 on its unbounded side; (C) the shortest witness of a target r has O(log |r|) letters"
status: (A) P1 structure CONFIRMED and per-class attainment REFUTED, P2 CONFIRMED, P3 CONFIRMED, P4 CONFIRMED, P5 CONFIRMED, P6 PARTLY REFUTED (8/9 at j = 4 and j = 9); (B) OPEN (FINITE-EXACT window fractions 0.097-0.617 for 59 interior classes, 0.190-0.601 for 15 exterior slopes); (C) OPEN (FINITE-EXACT witness lengths); audited with THM-4883 (PASS WITH CORRECTIONS, MISTAKE-724)
session: kind-pasteur-2026-11-05-S12
date: 2026-11-05
related: [THM-4883, THM-4882, HYP-9475, session198/rational-response-coverage/PROOF.md]
script: 05-knowledge/results/signed_response_2adic_structure_20261105/verify.py; 05-knowledge/results/signed_response_2adic_structure_20261105/witness_length.c; 05-knowledge/results/signed_response_2adic_structure_20261105/PREREGISTRATION.md
---

# HYP-9476 - signed one-apex responses: predictions, densities, witness size

Label map: the frozen `PREREGISTRATION.md` calls the density statement "HYP-9476 (A)" and the witness-size
statement "HYP-9476 (B)"; in this file they are (B) and (C), and (A) holds the pre-registered outcomes.

## (A) Pre-registered predictions (commit 2c11d95d56, before the confirmatory run)

- P1 (structure census, 629 slopes, b <= 64): one full class iff b is not a power of two, a sign and a bound for
  every other class, zero violations below 2^22 - CONFIRMED. "Every bound is the extremum of its class" -
  REFUTED: 1,143 of 1,902 non-full classes miss their bound (first: 3/4, class 7 mod 8, bound -5, maximum -9); the
  bounds are extrema of groups of classes (all 759 groups attain), because a piece of scale 2^nu < 2^s spreads over
  several classes.
- P2 (fullness, b <= 48, |r| <= 2000) - CONFIRMED: 544,115 targets, zero failures.
- P3 (exact half-lines) - CONFIRMED on 186 slopes.
- P4 (C2 offsets) - CONFIRMED: 78,000 identities, six literal pairs with offsets -1 and +1.
- P5 (densities in (0.02, 0.98)) - CONFIRMED on the 50 classes decided in the run (0.097-0.617); the audit showed the
  9 classes reported as lower bounds are decided too (0.526-0.546), so 59 classes.
- P6 (witness size, 3/5 and 8/9) - 3/5 all clauses hold; 8/9 fails clause (b): maxlen over |r| <= 16 is 21 > 20, and
  (the audit's depth-35 search) -427, -425, -405, 277 need at least 36 > 35 letters: PARTLY REFUTED. The tie-break of
  clause (c) (first minimal witness in depth-first order) was fixed after the registration; margins 184/522 and
  160/433 against 1/2.

## (B) Density conjecture (OPEN)

For every reduced slope and every class C of P_s other than C*, S(a, b) n C has a natural density (relative to C on
its unbounded side) strictly between 0 and 1. Evidence: the 59 non-full classes of the 14 tested interior slopes have
window fractions 0.097-0.617 on 500 <= |r| <= 20,000 (all classes decided: for pieces of slope 1 or 1/2 every window
value is reached at n <= 4|r| + 3 before scaling); for 15 exterior slopes (2/1, 3/1, 3/2, 5/2, 4/3, 5/3, 5/4, 1/1,
0/1, 1/3, 1/4, 2/5, 1/2, 3/8, 9/8) every class has window fraction 0.190-0.601 (audit upgrade 3); all unchanged from
n < 2^22 to n < 2^26. A counting argument proves density < 1 only for pieces of slope far from 1 and 1/2; the
boundary-slope pieces (images of n - g(n) and n - 2g(n)) are the obstacle. (B) would rule out that a non-full class is
cofinite on its unbounded side. It would not by itself show that S(a, b) is not eventually periodic for even b (a
class can meet S in an intermediate density and still be periodic with a finer period). The statement that every
two-sided progression inside S lies in C* is already THM-4883 (T1), not part of (B).

## (C) Witness size (OPEN)

THM-4883 (T5) shows that for an interior slope the clock is the only source-independent move in the unit class, so a
short witness must come from a target-dependent prefix. Conjecture: for every admissible slope there are C, C' with a
witness of at most C log2(|r| + 2) + C' letters for every target r (equivalently n <= (|r| + 2)^C up to a constant).
Data (exhaustive searches of `witness_length.c`; the files are in the packet): 2/3 maxlen over |r| <= 2^j is 20, 24,
24, 27, 30 at j = 10..14 (about 2.1 j, every target found by length 34); 4/5 23, 29 at j = 8, 9, then one target
beyond length 33; 6/7 26 at j = 8, 9, then beyond 33; 3/5 (pre-registered) 9, 17, 22, 22, 32, 32 at j = 4..9; 8/9 21,
21, 21, 29, 32 at j = 4..8. The trailing clock counts have a geometric tail (2/3: 6962, 2855, 1103, 420, 159, 53, 13, 3
for t = 1..8, ratio about 0.4), as a random model with a fixed chance per prefix length predicts. THM-4882's audit
records that for denominator 11 small targets (6/11: 7 and -11) already need n > 2^24, so C grows with the clock
length. Session 198's pump constructions give O(|r| + 1) vertices for 2/3 and 4/5 (upper-replacement-fibers,
(16)-(17)).
