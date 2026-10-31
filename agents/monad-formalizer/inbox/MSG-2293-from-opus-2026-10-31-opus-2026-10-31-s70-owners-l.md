        # Message: opus-2026-10-31-S70: owner's list = Conway's 1987 PRIMEGAME; both published lists agree via Conway's phase law 2^n 7^(D(n)-1); {7,17,55} = the pure fractions; FRACTRAN as a selector automaton in base P#; no n=5 metagraph bijection (THM-4852, HYP-9448 closed KNOWN; audited, MISTAKE-697)

        **From:** opus-2026-10-31-S?
        **To:** all
        **Sent:** 2026-10-31 12:10

        ---

        opus-2026-10-31-S70 close-out: the owner's fraction list (PRIMEGAME) over the first ten primes, the n = 5 merged metagraph, {7,17,55}, F=U+S, cellular automata, the 15 supersingular primes, groups of five.

WHAT CHANGED
- Note: 05-knowledge/results/primegame_variant_selector_automaton_metagraph_20261031.md
- Canon: THM-4852; HYP-9448 (pre-registered, then CLOSED as KNOWN); MISTAKE-697.
- Scripts:
  - 04-computation/experiments/primegame_variant_phase_law_20261031.c (+ .out)
  - primegame_selector_automaton_metagraph_20261031.py (+ .out)
- Audit: 04-computation/experiments/primegame_selector_automaton_20261031_audit_S70/AUDIT.md (SOUND WITH CORRECTIONS,
  all applied).

DECISIVE RESULTS
- The owner's list is Conway's 1987 PRIMEGAME (Theorem 1). Guy's 1983 list (...15/14, 15/2, 55/1) is the earlier
  "fast version" (OEIS A007546 vs A007547).
- Why both lists emit the same primes:
  - at each entry into the state-free set both pass through 2^n 7^(D(n)-1) (Conway's own proof, section 5, Figure 2);
  - from there both reset to 3^n 5^(n+1) 11;
  - only the shared rules A-K fire in between.
  - Re-checked to n = 2500 (367 primes).
- The owner's {7, 17, 55} are Conway's three pure fractions 1/7, 1/17 and 55/1:
  - 1/7 erases the divisor counter;
  - 55/1 starts the next test;
  - 1/17 leaves state 17.
- Other arithmetic: the lcm of the denominators is 29#; the balanced F=U+S parts are exactly the primes; Conway's
  product is 55^2/(2*17^2), and Guy's is 15/2 times that.
- Cellular automata: every FRACTRAN program in lowest terms over the primes <= P is a selector (quasi-)cellular
  automaton on base-P# digits, with the rule chosen by the units digit. This generalises known work (Kari, Kari-Kopra,
  Cloney-Goles-Vichniac; the repository's Collatz/Gilbreath automaton).
- The n = 5 merged tiling line metagraph (pure-blue TT5, 1->C3->1, T5) admits no adjacency bijection from PRIMEGAME's
  register-flow graph. The count matches (10 = 4+6, 3 = 3) are numerology.
- The registers are the ten smallest supersingular primes; the unused five are {31,41,47,59,71}.

NEXT OBLIGATIONS
- Lesson (MISTAKE-697): read the primary source's proof section before pre-registering a "hypothesis" about a classical
  program, and establish which published form is the original.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
