        # Message: mac-mini-2026-11-02-triangles: x^2 - 29/16 vs the n = 4 metagraph (NUMEROLOGY: counts 8-4-3 agree, fibre shapes differ); residue pairs of primes until the first +-1 equal the Goldbach pairs when m has a prime neighbour; audited (MISTAKE-702)

        **From:** mac-mini-2026-11-02-S?
        **To:** all
        **Sent:** 2026-11-02 21:00

        ---

        mac-mini-2026-11-02-triangles close-out. Owner's pasted prompt of 2026-11-02: the hourglass atom; {10, 20, 39}; primes modulo even numbers until ±1; the jagged odd-pair triangle in FRACTRAN/TM/CA form; and x -> x^2 - 29/16 against the n = 4 metagraph's 8 tilings -> 4 classes -> 3 merged nodes. All work is on main and pushed.

WHAT CHANGED
- Note: 05-knowledge/results/quadratic_2916_odd_pair_triangles_n4_20261102.md.
- Scripts: 04-computation/experiments/circulant_fractran_goldbach_20261101/quadratic_triangle_n4.py
  (pre-registered) and residue_pairs_posthoc.py, with outputs; audit_triangles/.
- Ledger: MISTAKE-702.
- No new THM, no HYP.

STATUS
- CITED: THM-4139 (the complete rational preperiodic graph of x^2 - 29/16); Poonen; Walde-Russo.
- Concurrent work credited: session146 prime290, the PRIME-MRB packet, hourglass_triangular_dynamics.
- FINITE-EXACT: pre-registered checks Q1-Q2, M1-M3, P1, P3 all PASSED.
- PROVED (elementary): the residue-pair proposition.
- Independent audit: SOUND WITH CORRECTIONS, all applied.

DECISIVE EVIDENCE
- The eight bounded odd/4 points map 8 -> 4 -> 3 onto the 3-cycle.
  - The n = 4 metagraph's 8 -> 4 -> 3 has different fibre shapes ({5,1,1,1} vs {2,2,2,2}).
  - The second stage is forced for any map of 4 onto 3.
  - No bijection matches the class line multigraph (a star with a loop) to f's graph (a triangle with a pendant).
  - Verdict: NUMEROLOGY.
- For even m, the residue pairs {r, m - r} hit by the odd primes scanned before the first prime = ±1 (mod m)
  - contain the Goldbach partitions of m always;
  - equal them when m - 1 or m + 1 is prime (proved);
  - strictly contain them for all 5,820 other even m <= 20000. This direction is OPEN: it needs a Goldbach-type
    statement.
  - The prime 2 has residue pair (2, m - 2), joining m to the previous even number.
- {10, 20, 39} = the chart_clique_norm carrier dimensions at n = 3, 4, 5.
- PRIMEGAME's trial-division region and the odd-pair triangle are both first-hit selections from jagged triangles,
  i.e. nested counter loops in machine form.

NEXT OBLIGATIONS
- A proof of strict containment R(m) ⊋ G(m) for even m with no prime neighbour. It may follow from
  Goldbach-in-intervals results.
- The hourglass itself: see hourglass_triangular_dynamics_20261103.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
