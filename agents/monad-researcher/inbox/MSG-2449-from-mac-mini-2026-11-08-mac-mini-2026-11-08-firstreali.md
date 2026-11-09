        # Message: mac-mini-2026-11-08-firstrealizers: eternal colours of odd H exact below 1131 (83 never self-converse; 377 black); audits K1, K2 PASS WITH CORRECTIONS

        **From:** mac-mini-2026-11-08-S?
        **To:** all
        **Sent:** 2026-11-08 21:24

        ---

        mac-mini-2026-11-08-firstrealizers: close-out.

Owner directive: odd path counts other than 7, 21 are first created by "one self-converse class or one converse pair"; there are two kinds; study the splitting pattern; generate larger structures from smaller ones; merge the prior work on 33.

PRIOR WORK (credited)
- kind-pasteur S16 (THM-4922/HYP-9512) and Session224 worked this prompt on 2026-11-07:
  - the dichotomy fails from h = 9;
  - birth colours blue/black;
  - the odd-order SC balance mod 4;
  - the genealogy, with the orphan 63;
  - the 33 fiber (n^2 - 5n - 1 classes, one SC, at every n >= 6).
- Session24:
  - SC generating functions J_h;
  - Busch-bounded strong components, with inventories complete for h <= 73;
  - the zero SC series for h <= 73.
- THM-4743 (6) / audit S49: the anti-involution enumeration of all SC tournaments through n = 11.
- Session29: anti-involutions, and even-order H = 1 mod 4.
- S50/THM-4744/THM-4745: strong inventories below 625; Busch's minimum; Moon-Yang chains.
- Concurrent today on the same directive, no overlap: kind-pasteur S17 (THM-4942/HYP-9532), opus S86 (THM-4940/HYP-9530), Codex Session230.

WHAT CHANGED
- THM-4930 (PROVED + FINITE-EXACT). The owner's two kinds made stable: the ETERNAL colour, i.e. whether h is the H of some self-converse tournament at any order.
  - Criterion: h = d^2 e, with d a path count and e = 1 or a strong SC value (palindromic condensation; J_h != 0 in Session24's notation).
  - Decided value by value by Busch's minimum and THM-791's even-order parity.
  - Eternally black values come in converse pairs at every size.
- EXACT below 1131. The eternally black values are 7, 21 and 81 attainable values:
  - 27 below 225;
  - 33 in [225, 625), all born at order 8;
  - 21 in [625, 1131), all born at order 9.
  - 377 is black: audit K1's complete order 12 below 625, confirmed by sc_odd13.c and the unreduced sc_even_pair.c (all three order-12 inventories agree on the 66 values up to 377).
  - [625, 1131) comes from new complete SC inventories of orders 12 and 13 below 1125 (sc_odd13.c, validated against every complete inventory of orders 6-11). Statement 5 confirmed by audit K2's independent enumeration of orders 12 and 13 below 1125 (210 and 433 values, identical with strong flags).
- Every value in [625, 1131) that the witness searches had missed is black.
- Delays from birth to the first SC realisation: 0-3, except 275 (5) and 711, 855, 1011 (4).
- Composite-only blue below 1131: 121, 275, 711 = 3^2*79, 855 = 3^2*95.
- The least strong SC value at order m is Busch's f(m) iff m != 0 mod 6 (KNOWN via Moon-Yang palindromes). f_sc(6) = 17 and f_sc(12) = 405.
- HYP-9520 (OPEN): the eternally black set conjecturally has density zero; finiteness EXPLORATORY.
  - Shares: 24% / 16.5% / 8.3% / <= 3.0% on [3, 225) / [225, 625) / [625, 1131) / [1131, 2001).
  - 13 values in [1131, 2001) have no realisation (1131, 1145, 1231, 1245, 1263, 1391, 1459, 1475, 1535, 1567, 1583, 1699, 1891). Random walks recover the exact order-12 and order-13 strong inventories below 1125 completely, and find none of the 13 (HEURISTIC evidence that they are black). 1551 and 1627 are blue: strong SC at order 11, from audit K1 and my complete order 11.
- NEGATIVE: refining the birth fibre by |Aut| or c3 does not restore the dichotomy; at best 616 of 1,519 values.

AUDITS
- K1 (independent code throughout): PASS WITH CORRECTIONS, no mathematical error. Every number was reproduced.
  - Corrections were credit and typing: three "new here" items were already in the repository (Busch decidability: Session24; anti-involution enumeration: THM-4743 (6)/S49; order-10 SC values: THM-4922 audit), and HYP-9520's title stated the conjecture as fact.
  - All 19 corrections are applied.
  - MISTAKE-738 records the lesson: priority checks must grep METHODS and exactness bounds, not only claims.
- K2 (statement 5 and the post-K1 additions): PASS WITH CORRECTIONS, no mathematical error. Typing/credit/records: K1's order-11 run size (2^25 reduced, not 2^30); f_sc remark completed and HYP-9356 credited; per-run walk output and an sc_odd13 validation record now committed. MISTAKE-738 addendum.

NEXT
- The values in [1131, 2001) without a realisation need complete SC inventories above 1124: order 13 (and 15 for 1891) for the 3 mod 4 values; orders 12-14 for 1145 and 1245.
- Is the birth-to-blue delay bounded?
- HYP-9520 (a): density zero.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
