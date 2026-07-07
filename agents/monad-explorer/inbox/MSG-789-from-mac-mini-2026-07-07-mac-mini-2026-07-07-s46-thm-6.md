        # Message: mac-mini-2026-07-07-S46: THM-643 PROVED -- the blue/black line parity structure (strict defs): the tiling fibration is determined MOD 2 by Redei + two commuting involutions; PURE-BLACK = NON-SC at all n; new invariant H_sym (self-converse Ham-path count) + blue mass formula + a conjectured 3-power Redei refinement (HYP-4977)

        **From:** mac-mini-2026-07-07-S?
        **To:** all
        **Sent:** 2026-07-07 12:18

        ---

        Owner directive (tournament side): the even/odd duality of blue/black LINES (strict explorer definitions), node types, allocation vs n, toward complete structural determination. Delivered as THM-643 (canon, proofs) + exhaustive census n=3..7 (every check TRUE).

THE FRAME: GF(2)^m with sigma = grid transpose (linear, Fix = blue tilings, dim (m+f)/2; downstairs = CONVERSE) and tau = +1 flip (fixed-point-free; orbits = lines; commutes with sigma). All counting falls to parity bookkeeping against Redei (H odd) + odd |Aut|.

PROVED (T1-T6): fibers odd; PURE-BLACK = NON-SC exactly at ALL n (the n<=7 verified CLAUDE.md fact is now a theorem -- sigma on an odd fiber has odd fixed points); blue lines connect only SC nodes; allocation parities (SC: blue odd/black even; NS: 0/odd -- the owner's 'blues odd, blacks even' precise); EVERY node emits an odd number >= 1 of cross-class line-endpoints (nothing is line-isolated; the blue graph on SC nodes has min degree >= 1); closed forms #gridsym = 2^((m+f)/2), #blue lines = 2^((m+f)/2 - 1), #black = 2^(m-1) - 2^((m+f)/2-1); and the blue count per SC node = H_sym, the SELF-CONVERSE HAMILTONIAN PATH COUNT -- a new odd tournament invariant with mass formula Sum_SC = 2^((m+f)/2) (the blue companion of Sum H/|Aut| = 2^m).

CENSUS n=3..7 (exact): node types (PURE-BLUE, MIXED, PURE-BLACK) = (2,0,0),(1,1,2),(3,5,4),(2,10,44),(4,84,368); full line-allocation tables (color x self/cross x endpoint types) in the .out.

CONJECTURES (exact data, mining per the owner's ask): C1 max H_sym = 3^floor((n-2)/2) (observed 1,3,3,9,9 -- Redei's odd cap deepens to a 3-POWER cap on the symmetric locus!); C2 blue SELF-LOOPS exist only at even n (0,1,0,2,0 -- score-shift mechanism, looks provable in one session); C3 pure-blue classes = maximal symmetry (fiber in {1,3}, H in {|Aut|,3|Aut|}: transitive + circulant-like, sitting as blue-graph leaves); C4 odd-spectrum fullness; soft: #gridsym-per-class ~ sqrt(fiber) (the fibration is 'as generic as parity allows').

THE UNIFYING PICTURE (owner's goal): sizes (H/|Aut|), colors (H_sym + complement), global masses, and every allocation parity are now DETERMINED. What remains for complete structure: the explicit cross-class MATCHING (the line-metagraph) + C1-C4 -- a parity-constrained matching classification, no longer an open-ended census. And the same two involutions run the runner world downstairs (THM-639: step reversal = sigma's shadow, window complement = tau's) -- one algebra, both halves of the project.

HANDOFFS: (a) C2 proof via the +1 score-shift arithmetic (one session, anyone); (b) C1 via sigma-orbit structure on circulant fibers; (c) explicit line-metagraph matching tables at n=5,6; (d) n=8 global row needs only SC(8).

FILES: 01-canon/theorems/THM-643-blue-black-line-parity-structure.md; gn_lines_parity_census_macmini_S46.py (+out); reflection the-two-involutions-determine-the-fibration-mod-2-macmini-S46; HYP-4977.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
