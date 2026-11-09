        # Message: opus-2026-11-09-S87: ternary Chan-Pak lanes + order-10 census (coverage to 29,035; strong spectrum misses only 10, 31, 64); packing-number metagraph lines (staircases, run-window theorem: tight-word laws proved at n != 0 mod 3); audit PASS WITH CORRECTIONS

        **From:** opus-2026-11-09-S?
        **To:** all
        **Sent:** 2026-11-09 09:17

        ---

        opus-2026-11-09-S87 close-out. Owner directive: prove ternary Chan-Pak and other proposed conjectures; metagraph black/blue lines between new node classifications.

WHAT CHANGED
- THM-4950 (ternary Chan-Pak).
  - PROVED: HYP-9530 (A) is equivalent to strong realisation of the irreducibles of M3 = {1 mod 3} minus {10, 31}: lanes p, pq, 10p, 31p, 31pq, plus 40, 250, 310, 9610, 29791.
  - FINITE-EXACT: the complete strong order-10 census has 10,070 H_3 values (max 31,549), with witnesses re-checked; the audit re-derived the set without canonical forms. Coverage goes from 5,398 to 29,035 from orders <= 10 alone (m* = 9679). Pre-registered Q1-Q3 CONFIRMED.
  - The strong base-3 spectrum misses only 10, 31, 64 below 28,981. This is exact below 337 (THM-4952 (E)'s method) and confirms HYP-9541 (4)'s strong-spectrum sentence.
- THM-4951 (packing-number nu lines; Session26's bound credited).
  - Tight lines: nu is read off the triple word; prefix and double staircases realise every split, the double ones in blue.
  - Run-window theorem: tight words = tight colourings, and blue = palindromic, at every n != 0 (mod 3). Block law at 3k + 2; one short block at 3k + 1, (k+1) 2^(k-2) words.
  - Census of all lines by nu-pair for n <= 9; tight lines 1,1,1,57,31,13,3673,1802,541,563609 (n = 3..12).
  - Blue H mod 4 lines: 1:2:1 to n = 9, not at 11; the (1,3) excess is -0.0009 / +0.0017 / -0.0001 at n = 13 / 15 / 17.
- HYP-9540: (A), (C) open only at n = 3k; (B) allowed pairs per colour; (D) sign law (-1)^(c+1).
- HYP-9530 (A) updated.

PRE-REGISTRATIONS
- Q1-Q3 CONFIRMED.
- R1 CONFIRMED (KNOWN mechanism); R2, R3 REFUTED.
- T1-T5 CONFIRMED.
- U1 REFUTED; U2-U4 CONFIRMED.

AUDIT: PASS WITH CORRECTIONS (2 majors in HYP-9540, 17 minors; all applied), MISTAKE-742. Follow-up on the post-audit (3e), (3f): PASS WITH CORRECTIONS. Both are correct (all 150,048 colourings with runs >= 2 at n <= 26, samples to n = 32). F1 is fixed: at n = 3k the run window fails exactly on an isolated 0 between two 1-runs of length >= 2, i.e. words containing 1110111; 0,1,3,8,20 at n = 6..18.

NEXT OBLIGATIONS
- HYP-9540 (A), (C) at n = 3k: the 1110111 words (isolated 0 between two 1-runs of length >= 2) need a construction beyond run windows; they are realised by other tight lines for n <= 12.
- HYP-9540 (D): explain the sign of |D_1| - |D_0|; test n = 19, 21.
- Ternary lanes P, Q beyond 29,038: an order-11 strong census or structured families.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
