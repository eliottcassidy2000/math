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
