        # Message: mac-mini-2026-11-01-lenses: THM-4857 (GF(2) automaton of a tournament classifies merged nodes exactly through n = 5, fails for all n >= 6; circulants = additive CA), THM-4858 (Goldbach as a 53-symbol radius-1 CA and a 568-fraction FRACTRAN program); audited (MISTAKE-701)

        **From:** mac-mini-2026-11-01-S?
        **To:** all
        **Sent:** 2026-11-01 10:11

        ---

        mac-mini-2026-11-01-lenses close-out. Owner directive of 2026-11-01: analogies and seamless translation between Turing machines, cellular automata, FRACTRAN and the tournament research, with the repository's open questions seen through these lenses. All work is on main and pushed.

CONCURRENCY
- Session135 (FRACTRAN/TM structure, checkpoint 2 at 09:22) is credited and not duplicated: TM realisations of
  tournament streams, square-sum and LRC certificate machines, a FRACTRAN<->TM compiler, a loss ledger.
- This session took the lanes it left: cellular automata, the GF(2) automaton of tournaments, and the Goldbach
  machines translated into each other.

WHAT CHANGED
- Canon: THM-4857 and THM-4858.
- Note: 05-knowledge/results/lenses_tm_ca_fractran_tournaments_20261101.md.
- Scripts: 04-computation/experiments/circulant_fractran_goldbach_20261101/: machine_translations (+ posthoc),
  tournament_automaton, automaton_similarity_types, automaton_similarity_n7; audit_lenses/.
- Ledger: MISTAKE-701.
- No HYP opened.

PROOF STATUS
- PROVED (elementary or standard constructions) and FINITE-EXACT.
- Pre-registered T2b and T3 FAILED as stated: the checker counted helper-loop returns as TM steps. Post-hoc runs
  confirm the program.
- Independent audit: SOUND WITH CORRECTIONS, all applied.

DECISIVE EVIDENCE
- THM-4858: strong Goldbach fails iff
  - Leng's 25-state TM halts, iff
  - its 53-symbol radius-1 cellular automaton writes HALT (lockstep 200,000 steps), iff
  - LENGGAME halts (568 fractions over 349 primes; two-stack encoding; stacks agree at the first 3,000 TM steps,
    with 24.4M FRACTRAN steps).

  GOLDBACHGAME (37 fractions) is a different algorithm. The sizes measure encodings, not the predicate.
- THM-4857: the GF(2) automaton x -> Ax of a tournament is converse-blind, so it lives on merged nodes.
  - Its similarity type classifies merged nodes exactly for n <= 5: 2, 3, 10 types; at n = 5 these are the owner's
    ten nodes, a small-n coincidence.
  - It fails for every n >= 6: 24/34, 51/272, and n >= 8 by counting (A070933).
  - At n = 5: four of the five mixed nodes give reversible automata, node #6 runs a single 31-cycle, and TT5's
    automaton is nilpotent.
  - Circulants give additive automata with rank q - deg gcd(a, x^q - 1).
  - The Paley/QR-code dichotomy (Hamming at 7, Golay at 23) is KNOWN: THM-481.
- FRACTRAN priority = a transitive tournament on fractions. A tournament-scheduled FRACTRAN is well defined iff its
  priority tournament has a source on every realisable co-applicable set. The clean per-state form holds for fully
  guarded programs, not PRIMEGAME.
- Open questions through the lenses:
  - Π1: Goldbach, Legendre, the H-spectrum (bounded realisers);
  - Σ1: a sixth Fermat prime;
  - Π2: Collatz, twin primes;
  - decidable: LRC(n) for fixed n;
  - proved: square-sum, pure-blue.

NEXT OBLIGATIONS
- Optional: compile GOLDBACHGAME into a TM and CA with session135's compiler.
- Optional: test whether finer automata (over GF(3), or nonlinear rules) see beyond n = 5.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
