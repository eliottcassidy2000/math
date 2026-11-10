# Pre-registration: kind-pasteur-2026-11-09-S20

Pushed before the runs named under "Predictions". Exploratory results obtained before this file are listed under
"Already seen".

## Owner prompt (2026-11-09)

"spend another session on chan-pak and come up with creative infinite constructions that can go beyond our attempts at
finite constructions for proofs"

## Concurrent and inherited work, read first

- **Near-joins.** HYP-9453 and THM-4865 give h + 2^k B for any rooted base, with B = B_t(T) the root response. HYP-9453
  section 2 already proposes the response-interval induction.
- **One-arc join transfer.** The 3 x 3 law (h, e, R) -> M(G, u, z)(h, e, R) is SEAM.md, with closed entries in
  THM-4875 (B) and the descent lemma in THM-4875 (C).
- **Fixed-port tails.** THM-4888 and THM-5085; decoders in THM-5096; minimal memory in THM-5091.
- **Root changes.** HYP-9480 and HYP-9494 (two root changes reach every odd N > 3,780,967 to 10^9).
- **Graft content.** Session242 (THM-5107, pushed 20:54 during this session) uses Q = sum_X H(B[X]) H(B[V - X]), the
  source response of this packet, inherited from Session229's proper-graft row. This packet does not claim Q.

## Lane of this session

Power-of-two root responses ("dyadic ports") turn near-joins into pure cones h -> h + 2^K. The questions are whether
cones can be chained, which equation governs the chaining, and whether self-feeding families with bounded response
complexity can be infinite.

## Already seen (exploratory, before this file)

- **Gadget tables.**
  - Closed forms for every transitive gadget TT_m (entry position p, exit position q) were checked against brute force
    on 3,360 cases (24 random bases of order <= 6, m <= 7): 0 mismatches.
  - 3 x 3 gadget matrices were computed for all 26,705 labelled gadgets of order <= 5 and all 24,747 class gadgets of
    order <= 7: 0 linearity failures.
- **Dyadic-port census.** Among rooted tournaments of order <= 8, 69 of the 148 attained odd values <= 300 have a port
  with B = 2^j.
- **No persistence through the standard near-join.** 15 cases tested; none has a dyadic port afterwards.
- **Dyadic closure** (triples with B a power of two, closed under one-arc gadgets with the exit in the gadget and under
  sink padding, started at the single vertex):
  - gadgets of order <= 4: values {1, 9};
  - gadgets of order <= 5: values {1, 5, 9, 31, 37, 59, 63, 409, 1151} below 10^6.
- **Bounded odd-part closures** (the odd part of B at most r0; gadgets of order <= 5):
  - r0 = 15: 258 values, identical below 10^5 and 10^7;
  - r0 = 63: 2,103 values, identical below 10^6 and 10^7;
  - r0 = 255: 14,340 values below 10^6 and 14,392 below 10^7.
- **e = 0 sheet.** For h <= 65, at most 62% of the even R in [2h, 4h] are realised at order <= 8.

## Predictions

- **P1.** The dyadic closure under all 24,747 class gadgets of order <= 7, started at the single vertex, is finite: it
  has fewer than 500 distinct values below 10^10, and none in (10^9, 10^10].
- **P2.** With class gadgets of order <= 6 and r0 = 63, the bounded odd-part closure has exactly the same value set
  below 10^9 as below 10^8.
- **P3.** With gadgets of order <= 5 and r0 = 1023, the bounded odd-part closure has between 60,000 and 150,000 values
  below 10^9.
- **P4.** Take every dyadic base (h, e, 2^j) realised by a rooted tournament of order <= 8, and every transitive gadget
  TT_m with m <= 20, entry p >= 2 and exit q < p. Exhaustive search finds no base that admits two consecutive dyadic
  cone steps with both new ports inside the attached gadgets.

Outcomes will be appended below, unedited, after the runs.

## Outcomes

- **P1: CONFIRMED.**
  - Setup: the dyadic closure under all 22,243 distinct class-gadget matrices of order <= 7, started at the single
    vertex.
  - Result: 197 values, largest 1,995,807; there is no value above 10^7 (`out/p1_dyadic_closure_order7_1e10.out`).
- **P2: CONFIRMED.**
  - Setup: r0 = 63 with the 1,884 distinct class-gadget matrices of order <= 6.
  - Result: 12,352 values, all below 86,371,467, so the value sets below 10^8 and 10^9 coincide
    (`out/p2_r63_order6_1e9.out`).
- **P3: CONFIRMED.**
  - Setup: r0 = 1023 with the order-<= 5 gadgets.
  - Result: 93,159 values below 10^9, inside the predicted window [60,000, 150,000]. Only the largest value,
    340,931,297, lies above 10^8 (`out/p3_r1023_order5_1e9.out`).
- **P4: CONFIRMED.**
  - Setup: the 181 distinct dyadic base triples of order <= 8, with transitive gadgets of size m <= 20.
  - Result: 0 two-step chains. Only 3 base triples admit even one dyadic step: (1, 0, 2) -> 9, (1, 1, 1) -> 9 and
    (33, 4, 64) -> 289. The last uses an exit at the gadget source with e = 4 != 0 (`out/p4_two_step_chains.out`).

## Addendum (pushed before the runs it names)

**Exploratory, seen before this addendum.** Each closure was taken with the order-<= 5 gadgets below 10^9. For
r0 = 3, 7, 15, 31, 63, 127, 255, 511, 1023 the closures are finite, with these counts:

| quantity | r0 = 3 | 7 | 15 | 31 | 63 | 127 | 255 | 511 | 1023 |
|---|---|---|---|---|---|---|---|---|---|
| values | 35 | 83 | 258 | 724 | 2,103 | 5,580 | 14,392 | 36,765 | 93,159 |
| complete prefix | 13 | 35 | 35 | 67 | 127 | 431 | 2,303 | 5,153 | 20,921 |

The complete prefix is the least odd n, other than 7 and 21, that is missing from the closure.

- **P5.** With the order-<= 5 gadgets and r0 = 2047, the closure below 4 * 10^9 is finite. Its complete prefix lies in
  [40,000, 120,000].
- **P6.** In the same run, the number of values lies in [180,000, 300,000], and every value is below 2 * 10^9.

### Addendum outcomes

- **P5: CONFIRMED.**
  - Setup: r0 = 2047 with the order-<= 5 gadgets, below 4 * 10^9.
  - Result: the closure is finite, and its complete prefix is 55,343 (every odd n < 55,343 except 7 and 21), inside the
    predicted [40,000, 120,000]. The value sets below 10^9 and below 4 * 10^9 coincide
    (`out/p5_r2047_order5_4e9.out`).
- **P6: CONFIRMED.**
  - Setup: the same run.
  - Result: 233,115 values, inside [180,000, 300,000]. The largest is 340,931,297, below 2 * 10^9; it is also the
    largest value of the r0 = 1023 closure.

## Notes appended after the independent audit (outcome lines above unedited)

- **P2 wording (audit minor 8).** "all below 86,371,467" should read "at most 86,371,467": that number is the maximum
  value.
- **HYP-9656 (A) as worded (audit M5).** The closure C(r0) is a set of triples, and it is always infinite because
  (1, 0, 2^k) belongs to it for every k. The predictions P1-P6 concern the value set h(C(r0)). The programs bound
  both h and R, and that bound is exact because a gadget step gives h' >= R.
- **"Already seen" (audit minor 6).** "No persistence through the standard near-join" is true only of the 15 cases
  tested. P4's own 33 -> 289 is a near-join landing on a dyadic exit.
- **Reproduction.** The audit reproduced P1-P6 and the complete-prefix table with independent code. It also ran
  r0 = 4095, giving N(4095) = 145,433 and 580,505 values.

