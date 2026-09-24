        # Message: collatz-procgen-20260924 wave 8: THM-4469 Mahler bridge, THM-4470 pairing ladder, residue types side-blind

        **From:** mac-mini-2026-09-24-S?
        **To:** all
        **Sent:** 2026-09-24 02:21

        ---

        collatz-procgen-20260924 (mac-mini), wave 8. The owner asked what Collatz's fundamental features represent, for "one proves the other" links to famous problems, and to study incongruency, {2,3,11} and odd-square brackets, the doubling/(2n-1)/3 tree mod 192, and the analogy "graceful tree : 3N+1 :: square-sum arrangements".

WHAT CHANGED (all integrated in the synthesis, section 2f; reflection "Wave 8"; MISTAKES 2026-09-24):
1. THM-4469 (PROVED + INDEPENDENTLY AUDITED), the Mahler bridge. No-divergence on an adjacent supercritical block pair is EQUIVALENT to a generalized Mahler Z-number statement. Smallest instance: no positive integer's 3x+1 parity vector is eventually made of 0111101110 and 1101100111 iff no xi>0 has every xi(2187/1024)^j in Z+[4726,4727]/1163. The interval needs 1.88/p, and FLP, Dubickas and Bugeaud stop at 1/p. So Collatz, T1 or PC imply new Mahler-type theorems, and a Mahler-side proof settles a HARD slice of no-divergence (HYP-9134).
2. THM-4470 (PROVED + AUDITED), the pairing ladder. T(2i-1)+T(2i)=(2i-1)+2i. The only arithmetic-mean-fair consecutive pairings are 3n+1 and 3n-1, and the two sheets are the two consecutive pairings. Halving and up edges are two perfect difference systems (graceful = ANALOGY: no subtree is literally graceful). The pairing family shares all of this, yet:
   - a density-zero flip set gives a divergent orbit;
   - single flips create cycles at exactly 24 fragile pairs, all <=2308, checked to 1e7 (HYP-9135);
   - no periodic pairing is provable by bounded lookahead;
   - the landing-down price is >=29% (HYP-9136).
3. The owner's mod-192 tree programme is exact (E D^2 = S E is the "4x recursion"; the automaton mod 3*2^k has growth 4/3). Residue types are PROVED side-blind: x->-x transports every residue, Haar, integrality, l-adic and |x| datum to 3n-1. So a proof must use the sign law, integrality at the gates 2^p-3^a, and pointwise avoidance of the null set Bad. The HYP-9121 16-point trap is the mirror image of the 3n-1 cycles.
4. Implication atlas (58 nodes, 83 edges): no famous conjecture is known to imply Collatz, T1, NC or PC. Also:
   - Lang-Waldschmidt gives cycle length >= N^(1/2-eps);
   - abc is dominated by Baker at the cycle gate;
   - CST together with T1 implies Collatz;
   - Kohl: Collatz <=> a transitive class-transposition group;
   - Collatz's 1932 permutation shares the cycle gates.
5. Brackets: {2,3,11} is the prime escape set on (25/13, 25/11], with a complete (k,p) list. Square-sum with brackets is REAL; square-sum with Collatz is NUMEROLOGY (Collatz's exceptions are Diophantine).

DECISIVE EVIDENCE: all three lane pipelines re-run byte-identical, apart from timing lines. Independent orchestrator code checked:
- the M' census (10,7):4, (16,11):5, (19,13):1, (20,14):8;
- the cylinder classes 990 and 187, and the finite-depth equivalence;
- the 24 and 12 fragile pairs to i=20000;
- the Perron root 4/3 and the five Z-cycles on [-1e5,1e5].
The orchestrator check is procgen_atlas_20260924_orchestrator_check.py.

NEXT OBLIGATIONS: HYP-9134 (beat 1/p by a constant factor for one ratio 3^a/2^L); HYP-9135 (fragile pairs via the Simons-de Weger/Hercher mutant cycle equation); HYP-9127 is still the smallest zero-entropy instance. HYP numbers 9137-9159 remain reserved for this session.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
