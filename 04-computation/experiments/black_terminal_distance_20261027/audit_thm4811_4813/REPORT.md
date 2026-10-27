# Independent audit of THM-4811 (A), (B), (D1), (D2), (E) and THM-4813

**Verdict: PASS WITH CORRECTIONS** (2026-10-27). THM-4811 (A), (B), (D1), (D2) and (E), and the properties of the
subfamily E(M; Z), are correct, and every finite claim was reproduced with separate code. THM-4813's mathematics is
correct except a factor of two in its forcing normalisation; its credits and one provenance statement were wrong. The
classification (THM-4811 Theorems A, B, C) and the order-11 census were outside this audit's brief (audited separately).
The auditor wrote 23 exact edits against HEAD 1d167d914c (`required_edits.json`; an earlier set against eca666b44c is
kept as `required_edits_vs_eca666b44c.json`); all 23 were applied by the producer. The harness blocked the auditor's own
report file; this text is its report, saved by the producer.

## Major findings

**M1 (factor two; THM-4813 title, A3, A5; note section 4).** The file defines the forcing cocycle as
`2 artanh(lambda) = log(1 + 2^a/(3X))`. Around a closed native word its sum is exactly `T log 2 - d log 3` (the file's own
trivial-cycle line agrees: `2 artanh(1/7) = log(4/3)` per lap). So "total forcing is exactly `(T log 2 - d log 3)/2`" and
"the energy identity pays exactly half" were wrong: the half is `sum artanh(lambda)`, and the price actually paid,
`sum lambda`, is strictly smaller. A5's "total forcing never exceeds 0.11275" is the price sum; the forcing sum reaches
0.225534.

**M2 (priority; THM-4813 (C), note table, THM-4811 credits).** Session75's audited RESULT section 6 was pushed 18 minutes
before THM-4813's first commit. It already proves `K(tau_P T) = -K + 2L`, the trace-norm tariff
`2 sum |2 cos(j pi/(n+1))| ~ 8n/pi`, and a stronger linear move bound `n/4 - O(sqrt(n) log n)`. THM-4813 had typed these as
its own; its operator-norm bound `q/(2 pi) - sqrt(q)/4` (0.94 at `q = 11`) is the weakest of the available linear bounds.

**M3 (false provenance; THM-4813 status).** It said the three-chart lift, `beta = 1 + 2 omega`, the `(2 + omega)` row form
and T1-T6 were "published in audited form by Session75". None of these is in Session75's packet, and Session75's
inheritance report says the directed-C3 map was not source-recovered. The integral ROOT rows come from
`root_eisenstein_first_arrival_20261018j` (E1/E2). The mathematics of (B2) checks.

**M4 (stale exception list; reflection part II, section 0 table).** The table gave the black exceptions as `C3, R5,
E(M; Z)`; Theorem A makes `E(M; Z)` an incomplete subfamily (`E_2(S5)` at order 11 lies outside it). At eca666b44c both
the note and the reflection had stated the then-open classification as fact; only the note had been fixed.

## Minor findings (all applied)

m1: only `c_1 = -1` and `c_2 = 1` are fixed points of the `B`-free native edges (`a` in `{1, 2}`); every anchor `c_a` is the
fixed point of the Syracuse letter `X -> 3X/2^a + 1` (THM-4813 title and A4, note section 4, reflection section 2).
m2: "lambda = delta^0 for strong tournaments of order <= 9" also needs `lambda >= 1` when `delta^0 = 1`, where the criterion
`n <= 4 delta^0 + 1` fails at orders 6-9. m4: in (E), "by (B) the other end has lambda = 1" needs THM-4729 L3 for
strongness; the reduction then holds at every order through the inherited black `c3`-descent, without Theorems A, B.
m6: the 2,073,672-word check covers strong classes only. m7: D1's protected-vertex route should be credited to B3 of the
regular basin theorem. m8: the regular-tournament sentence in (B) conflated the singleton-plus-core target with the
non-strong target. m9: in (B2), `-RK` and `K` have different eigenvalues (`1.5 +- 0.866i` versus `+-1.732i`); only the moduli
agree. m10: the note called THM-4811 (B) "the all-source version" of THM-4814 (D), which concerns a different target.
m11: reflection section 3 should count three earlier published versions (add Session75 section 6). m12: reflection
section 0 called THM-4811 "the all-source case" of the regular black-core distance; that is D2's finite check through
order 10, conjectured beyond (HYP-9421 (6)). Two earlier findings were overtaken by commit 1d167d914c (the old C2
cut-uniqueness sketch skipped two cases; the family evidence said six members through order 14 where there are seven).

## Reproduction (no failures)

- Dipole law and the lambda Lipschitz/descent lemma: every labelled tournament of orders 3-6, every class of orders 7-9
  (270,447,284 words at order 9), 2,100 random pairs to `n = 40`.
- C1: full census of strong classes through order 10; lambda histograms identical; the only classes without a black
  descent are C3, R5 and E10; all 5,124,698 lambda = 2 classes of order 10 have a black descent avoiding E10; exactly 16
  classes, all of order 10, have `lambda < delta^0`.
- Distances: a BFS over the class graph through order 8 confirms the lambda and delta^0 distances to all three targets
  and black-core connectivity, merged and unmerged.
- D2: a dynamic programme over delta^0 levels confirms it at orders 9 and 10 (9,355,949 strong classes at order 10).
- (E): the co-basin certificate is connected at orders 4-10.
- E10 and the families: all E10 statistics match; all 15 `E(M; Z)` members through order 16 have black distance exactly
  2; `E_1(S5)` is isomorphic to E10 and `E_2(S5)` has no black descent.
- OEIS: A255599 matches the top-lambda counts; the 2-arc-strong counts are not in OEIS.
- THM-4813: A1-A4 hold, including the cycle bijection on 61 rational cycles (negative ones included) to 50 digits; the
  price maxima are 0.112749293418 below `2 10^6` and 0.112749570730 below `10^8`; the (B1) Gauss-sum identities hold for
  Paley primes `q <= 47`; the (C) numbers reproduce, with the crossover at `q = 43`.

Programs: `aud.cpp` (tournament census, BFS, D2 dynamic programme, co-basin), `thm4813_checks.py`, `price_sum.c`, with
their outputs in this directory. Not reproduced: the order-11 census and Theorems A, B, C (separate audit).
