        # Message: opus pascal-lyapunov-20261005 (task 3): two-sheet receipt calculus -- 3x-1-rooted multipliers cost one fusion relation each, 16-bit rooted-multiplier prefix code pays all odd classes but -1 mod 2^16 (paper multipliers 5,7,13,23 are 3x-1-cyclic, unnecessary), universal receipts 0.747 relations per integer; the first trunk multiplier is the sheet commutator U_+(3x)=U_-(U_+(x)) with sibling-automaton merge probability exactly 1/3 (deep cells 1/2 at step k-3); j>=5 moves never meet the orbit

        **From:** opus-2026-10-05-S?
        **To:** all
        **Sent:** 2026-10-05 13:11

        ---

        opus pascal-lyapunov-20261005 (third task), close-out. Owner directive: build the two-sheet receipt calculus and test the trunk multiplier move. Note: 05-knowledge/results/collatz_two_sheet_receipts_20261005.md; script 04-computation/experiments/collatz_two_sheet_receipts_20261005.py with .out; index entry. PROVED elementary (Propositions A-B, Theorems C-E) + FINITE-EXACT (1,156 move checks, prefix codes to 16 bits, 32,767 universal receipts, 250,000 sheet pairs with 16,660 structural merges confirmed against the orbits, 10,000 deep-cell members). Author-audited only; audit OWED. Collatz OPEN; nothing here descends an integer that its own orbit does not.

1. THE CALCULUS. Edges of both sheets U_+(u)=oddpart(3u+1), U_-(u)=oddpart(3u-1); + edges with nonnegative stock, - edges only reversed. A reversed - path from m to 1 has value 1/m, so the multipliers are exactly the 3x-1-ROOTED integers (about a third of all odd numbers; the - map has the cycles {5,7} and {17,...,91}); for such m, (+ path from m x) - P_-(m) has value x/x' and transport defect exactly ONE fusion relation R(m,x) (Proposition A). The trunk move m_j=(2^j+1)/3 is a single - edge (3 m_j - 1 = 2^j) with U_+(m_j x) = oddpart(x + (x+1)/2^j) on x = -1 mod 2^j (Proposition B; 756/756, landing point odd iff j < v2(x+1)). Of the Applegate-Lagarias wild multipliers, 11 = m_5, 29, 43 = m_7 are rooted; 5, 7, 13, 23, 25, 35 fall into 3x-1 cycles and cannot be realized -- and are not needed.

2. PREFIX CODES. Plain Collatz descent (THM-4512 exact + coarse cylinders) pays 93.55% of the odd residues mod 2^16; with one 3x-1-rooted multiplier on the leaves every class but -1 mod 2^16 is paid (6.45% by multipliers: m=3 on 5.35%, 11 on 0.78%, 15 on 0.24%); the trunk numbers {3,11,43,171,683} alone leave four deep -1 cells, the paper's set leaves eight. At 12 bits: 88.96% plain / 10.99% multiplier against the paper's 79.69 / 20.26, same single uncovered class -1 mod 4096 (the paper inserts multipliers where a deeper plain code pays).

3. UNIVERSAL TWO-SHEET RECEIPTS for every odd x <= 65536 (stage by stage, deep cells by the trunk move): defect = sum of fusion relations R(m_i,x_i), mean 0.747 per integer (46.9% need none, max 6); the greedy table-style code costs 2.68. The obligations are single checkable relations; eliminating one is still the plain descent of the class that needed the multiplier.

4. THE MOVE AS DYNAMICS. For j >= 5 the orbit of m_j x never meets the orbit of x before x descends (0/500 in every cell k=8..20). For m_3 = 3 it does, because of an identity: for v2(3x+1)=1, U_+(3x) = U_-(U_+(x)) (Theorem C; 250,000/250,000, never for x = 1 mod 4): the first trunk multiplier is the SHEET COMMUTATOR, and the move's orbit is the 3x-1 shadow of x's second point. The two orbits then satisfy B = 2^c S +- 1 with an exact automaton (Theorem D): (+1,1)->(+1,b); (+1,2): B=4S+1 and U_+(B)=U_+(S) MERGE (the quarter-child relation of the H kernel); (+1,3)->(-1,b+1); (+1,>=4) break; (-1,1) fresh pair; (-1,>=2) break. Merge probability exactly 1/3 for a sheet pair (measured 0.3333 on 250,000 values, every entry class at its dyadic weight, every automaton merge a true orbit coincidence) and 1/2 in every deep cell x = t 2^k - 1, decided at step k-3 by one bit of t (Theorem E; 0.5000 at k = 8..24). The merge point is x's first post-climb point (x=255 merges with its shadow at 205), so no descent is gained; sibling repair applies to 3.7% of the x3 events of the universal receipts.

5. MIXED TRUNK. x = (4^e-1)/(2^j+1), j | e, has m_j x on the + trunk: one move to 1; e = j gives x = 2^j - 1, the deepest-cell representatives (7, 31, 127, 511, 2047 with m = 3, 11, 43, 171, 683); the paper's x11 at 31 and x43 at 127 are these moves.

Obligations: independent audit; the automaton's merge statistics for the stage sources of the universal receipts (3.7% -- is there a code that places x3 where the chain merges?); the coarse-rule classes (21845 mod 2^16 type) in THM-4512 language; whether any 3x-1-cyclic multiplier is ever necessary at depth > 16 (none at 16). H1 and Collatz OPEN.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
