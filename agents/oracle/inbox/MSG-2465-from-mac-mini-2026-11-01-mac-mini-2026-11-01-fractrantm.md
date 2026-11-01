        # Message: mac-mini-2026-11-01-fractranTM: THM-4855 (circulants of order <= 5 = palindromic path words; explicit witnesses for every odd q >= 7), THM-4856 (GOLDBACHGAME/FERMATGAME/LEGENDREGAME on PRIMEGAME's divisor loop); audited (MISTAKE-700)

        **From:** mac-mini-2026-11-01-S?
        **To:** all
        **Sent:** 2026-11-01 09:04

        ---

        mac-mini-2026-11-01-fractranTM close-out. Owner directive of 2026-11-01: why the strong points of symmetry are exactly the circulants of order <= 5; FRACTRAN and Turing-machine connections with more repository questions; the 25-state Goldbach machine; the five prime groups vs the first 25 square-sum graphs. All work is on main and pushed.

CONCURRENCY
- Sessions 129-134 and opus-2026-10-31-S71 were read first. They encoded no conjecture as a halting program, did
  nothing with Turing machines, and did not explain the circulant cut.
- Session65's audit of Leng's 25-state Goldbach machine was the starting point for the machine. Session135
  (concurrent, PRIMEGAME basins and acceleration) is credited.

WHAT CHANGED
- Canon: THM-4855 and THM-4856.
- Note: 05-knowledge/results/fractran_turing_circulant_goldbach_20261101.md.
- Scripts: 04-computation/experiments/circulant_fractran_goldbach_20261101/ (7 scripts with outputs; audit/).
- No HYP opened.

PROOF STATUS
- THM-4855 and THM-4856 are PROVED (elementary) and FINITE-EXACT.
- Pre-registered G4 and L3 FAILED as stated (a measurement artifact and a mis-specified control); post-hoc anatomy is
  recorded.
- Independent audit: SOUND WITH CORRECTIONS, all applied (MISTAKE-700). Fixes: a proof citation (now Burnside +
  normaliser), family B credited as THM-4730's three-interval path, the closure claim scoped, the count corrected to
  131,068, the Turing-machine link retyped DICTIONARY, the cost corrected to about 9N^2.

DECISIVE EVIDENCE
- THM-4855:
  - A Hamiltonian path of Circ(Z_q,S) is a word over S with distinct partial sums (a path word).
  - It is palindromic iff segment sums of the word and its reversal lie in S together. For prime q this means the
    word is a palindrome.
  - q = 3, 5: every sequencing is a palindrome (T5: 1111, 2112, 2222 = its three tilings).
  - q >= 7: explicit witnesses A1/A2/B/D with named failing segments, checked on all 131,068 connection sets to
    q = 35. This gives an all-order arc-level proof along THM-4730 (B)'s original explicit route; the theorem itself
    was already proved by the audited prefix-rerouting route.
- THM-4856: GOLDBACHGAME (37 fractions, start 43), FERMATGAME (36), LEGENDREGAME (33) halt iff the respective
  conjecture fails. All are built verbatim on PRIMEGAME's first eleven fractions (the divisor loop).
  - GOLDBACHGAME: N <= 1000.
  - FERMATGAME: F_0..F_4; its inverted control reaches F_5.
  - LEGENDREGAME: 120 windows.
  - Literal = dispatch on every checked configuration.
- Leng's 25-state machine:
  - one SCC, 14 scanner states, 88% of steps in its unary trial-division states;
  - per-stage cost n^3.26, against GOLDBACHGAME's about 9N^2 per target;
  - it shares only the trial-division design with PRIMEGAME (DICTIONARY); the structural link is GOLDBACHGAME's
    literal embedding of PRIMEGAME's kernel. The merged metagraph relations are NUMEROLOGY or GENERIC.
- Prime groups:
  - reading 3 added (A -> Q_15..Q_26, no cycles; D and E -> cycles);
  - interleaving obstruction: B, C, D interleave in the prime order, so no size-monotone graph property separates
    them;
  - Monster levels are arithmetic, not size. NUMEROLOGY.

NEXT OBLIGATIONS
- Optional: minimise the three programs. Compile the 21-state DeepMind Goldbach machine and compare.
- The H-spectrum statement (every odd H except 7, 21) is Pi_2. A halting encoding would need a size bound for
  minimal realisations; that is open.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
