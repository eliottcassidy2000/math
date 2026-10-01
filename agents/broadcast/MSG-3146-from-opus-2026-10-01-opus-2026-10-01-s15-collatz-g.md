        # Message: opus-2026-10-01-S15: Collatz graph is rigid (THM-4523, audited); connectivity from rigidity blocked (transfer barrier); Gamma_1 dense in Z_2 x Z_3

        **From:** opus-2026-10-01-S?
        **To:** all
        **Sent:** 2026-10-01 10:52

        ---

        opus S15, seventeenth and eighteenth notes (2026-10-01). The owner asked for Collatz's "uniqueness amongst infinite directed graphs with outdegree 1" (with the Beckâ€“Everett harmful-relation summary and an Eckmannâ€“Hilton remark), then asked me to "try proving connectivity using the rigidity theorem".

1. THM-4523: the Collatz graph is rigid. PROVED and INDEPENDENTLY AUDITED (blind subagent, own code, SOUND; corrections recorded in MISTAKE-553).
   - The unfolded backward tree of n is the 3-adic tree e(n)/o(n).
   - A separation lemma (a crossed matching forces 14(aâˆ’b) â‰¡ 0) shows that depth-D truncations separate labels mod 3^(2+âŒŠ(Dâˆ’4)/5âŒ‹). The observed s(D) is 2..8 for D = 4..26.
   - Hence vertices prime to 3 have pairwise non-isomorphic backward trees and Aut(N,T) = 1. The same holds on Z and on R_6, where the vertex 0 is the unique fixed point.
   - Corollaries:
     - A graph isomorphic to (N,T) has exactly one isomorphism to it.
     - (N, 3n+c) determines c.
     - On Z the graphs of 3n+c and 3nâˆ’c coincide (negation), so the statement there is câ€² = Â±c.
   - Eckmannâ€“Hilton point: the two preimage maps 2x and (xâˆ’1)/3 coincide only at âˆ’1/5, the fixed point of 6x+1 (6x+2 fixes âˆ’2/5). That is the only place a fibre can be symmetric, and âˆ’1/5 is 2-adically odd.
   - Near balance: 2vâˆ’1 â‰¤ depth â‰¤ 5vâˆ’6 with v = v_3(5w+1); w = 106288 balances to depth 34. These are Beckâ€“Everett's harmful relations turned around.
   - Theorem R_a: for n/2, an+1 with a an odd prime, backward trees separate vertices iff 2 is a primitive root mod aÂ² (Artin + Wieferich). The audit adds that whole-graph rigidity holds for every non-Wieferich a, e.g. a = 7 (sketch).
   - Proposition M: every residue graph Gamma_M is connected, so there is no periodic invariant; this also holds for 3nâˆ’1, which has three components.
   - Proposition L, in the language {f} only (not a PA/ZFC barrier): the cycle half is first-order; the divergence half is not.
   - Theorem C: the trivial component is characterised among all functional graphs by its local rules and one branch three-cycle.

2. Eighteenth note: no connectivity proof from rigidity.
   - Theorem B (transfer barrier): rigidity holds for 3nâˆ’1 on N (three components) and for 5n+1 (three cycles), and positivity is invisible at every finite depth. So rigidity is sheet-blind and drift-blind.
   - Proposition 1 (new): Gamma_1 is dense in Z_2 Ã— Z_3. The trunk points (4^jâˆ’1)/3 cover Z/3^k bijectively in j, and the Terras reduction does the rest. So every acyclic finite view of every integer occurs in the trivial component: divergent components are locally invisible.
   - Proposition 2: n âˆˆ Gamma_1 iff Î´_k(n) = min{m âˆˆ Gamma_1 : m â‰¡ n mod 2Â·3^k} is bounded. On 3nâˆ’1, Î´_k(5) and Î´_k(17) grow like 3^k.
   - Proposition 3: the minimal counterexample is the minimum of its component (n_0 â‰¡ 3 mod 4, n_0 â‰¢ 2 mod 3, n_0 â‰¢ 4 mod 9). The backward sieve is weak: about 0.305 of classes survive at depth 22. The same trap is realised by 5 and 17 on 3nâˆ’1.
   - What remains is the synthesis's 2-versus-3 transversality, reached from the graph side.

Directions D57â€“D64 are in the notes. For the LRC sessions, nothing new beyond the sixteenth note's D53.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
