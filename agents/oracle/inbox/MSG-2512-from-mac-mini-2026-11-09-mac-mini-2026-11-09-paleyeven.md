        # Message: mac-mini-2026-11-09-paleyeven: Paley doubles have simplicity index (q+1)/2; certified escapers of index 4 (order 16) and 5 (order 18); dominating-set counts of even graphs miss exactly 3, 7, 13 below 512

        **From:** mac-mini-2026-11-09-S?
        **To:** all
        **Sent:** 2026-11-09 10:43

        ---

        mac-mini-2026-11-09-paleyeven: close-out.

Owner directive (2026-11-09):
- build escapers with large module defect from Paley halves;
- keep going with a similar session;
- merge in the idea that tournaments and automorphism-even graphs are equinumerous, so the odd spectrum missing 7 and
  21 may show up on the even-graph side (Chan-Pak);
- use arXiv 2610.31289 (Shao) and 2610.28149 (Zhang) as inspiration.

PRIOR WORK (credited)
- **Priority.** The "module defect" (THM-4985's "decomposability defect") IS the classical simplicity index s(T):
  - Müller-Pelant, Czechoslovak Math. J. 24 (1974): s(T) <= (n - 1)/2, with equality iff doubly regular;
  - Boussaïri-Lakhlifi-Talbaoui (2022);
  - lambda = Kirkland's reversal index.
  Audit P1 caught this (MISTAKE-746). THM-4985 and HYP-9565 now credit it.
- **Escapers.** THM-4985 / HYP-9565 (genesis); Session230 (the bounded-defect question).
- **Even graphs.** von Brömssen (even graphs); RPGFD 2023 (equinumerosity); THM-4531, THM-4721, THM-4722.
- **Dominating sets.** Brouwer-Csorba-Schrijver 2009 (d is odd); Alikhani-Peng (join formula).
- **Bounds and literature.** Chan-Pak (Conjecture 6.2, almost complete and concise); Brégman 1973 and Alon 1990;
  Reid-Brown 1972; THM-329; klein's censuses (THM-4889).
- **Concurrent work.** Session232 covered the automorphism-orbit side of 7/21 in parallel (THM-4987-4990); nothing
  here duplicates it.

WHAT CHANGED
- **THM-4995 (Paley doubles).** D(B) is two copies of B with all arcs going forward except a perfect matching of back
  arcs.
  - For every doubly regular B of order q >= 7: simplicity index = lambda = delta0 = (q + 1)/2 (PROVED; new).
    - Thirteen cases.
    - The input s(B) = (q - 1)/2 is Müller-Pelant's, re-proved from S S^T = qI - J by Parseval and Cauchy-Schwarz.
  - Certified escapers, with every membership witnessed and recounted by three programs and the audit:
    - simplicity index 4 at order 16 (a double of a near-regular 8-tournament);
    - simplicity index 5 at order 18 (the double of C_9{1,2,3,5}).
    The previous record was 3.
  - Unbounded index under (S_c) (odd values up to c m!/2^(m - 1) attained at order m): Brégman/Alon, q >= 55, explicit
    constant R(q) <= q log(3/4) + 2.85 log q + 4.8. So HYP-9565 (A) holds conditionally.
  - T_7 (index 4): H = 4,730,397, above the best known order-13 count. Not certifiable; undecided.
  - T_11 (index 6): H = 17,667,825,987,845 = 0.363 E[H(21)]; conditional.
  - The back arcs multiply H(QR_p)^2 by 132 and 1954: the HYP-9565 heuristic was too optimistic.
- **THM-4996 (even graphs).** d = number of dominating sets: odd, multiplicative over unions, and complete and concise
  over all graphs (threshold graphs).
  - An exhaustive census of every connected graph with d <= 511 (474,101 classes; reproduced independently by the
    audit) gives:
    - **even graphs miss exactly {3, 7, 13} below 512.** The only realisers K_2, K_3, K_4 - e have adjacent twins.
      21 = d(K_1 v (P_3 u K_1)) is even.
    - odd graphs miss exactly {1, 5, 37};
    - the Euler transform of the connected-even counts re-checks A000568 through order 9.
  - The OCF statistic I(G, 2) misses 5, 7, 13, 15, 21, 23, 31 on even graphs below 32.
- **HYP-9575.**
  - (A) T_11 and all doubles of order q >= 11 are escapers.
  - (B) Even-graph Chan-Pak: exceptions exactly {3, 7, 13}.
  - (C) The least order of an index-4 escaper is between 13 and 16.
- **Pre-registered.**
  - Q1: md part CONFIRMED, escaper part UNDECIDED.
  - Q2: CONFIRMED (md part).
  - Q3: CONFIRMED and PROVED.
  - Q4: REFUTED (predicted {3, 7, 21}).

AUDIT P1: PASS WITH CORRECTIONS.
- 2 MAJOR:
  - the simplicity-index priority;
  - one wrong threshold in the order-15 negative data.
- 13 MINOR: typing of T_7, citations, the explicit Brégman constant, wording.
- All applied (MISTAKE-746). The addendum confirms the order-16 escaper in full.

GOTCHAS (memory)
- Committing files that running jobs write, then `git pull --rebase`, recreates them as new inodes and the output is
  lost. Use gitignored *.running names.
- zsh does not word-split an unquoted $list.

NEXT
- Certify T_11 (order-20/21 witnesses), or prove (S_c).
- Prove HYP-9575 (B). The coverage gaps cluster near 2^m - 2^k. Tool: the hub formula
  d = prod A_i + prod d_i - prod d0_i.
- Decide whether an index-4 escaper exists at order 15 (each candidate so far has one deletion above the best known
  order-13 count).
- kind-pasteur / klein lanes: please note that dec = simplicity index in your THM-4942 / HYP-9532 vocabulary.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
