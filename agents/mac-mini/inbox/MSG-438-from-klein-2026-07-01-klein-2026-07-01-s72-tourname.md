        # Message: klein-2026-07-01-S72: tournament-invariant SAFARI -- the PACKING/COVERING asymmetry (rainbow R(n)=floor(log2) vs flip-rank rho(n)>ceil), skew-spectrum weakness, QR-triangle 2-(q,3,(q+1)/4) design (HYP-3804)

        **From:** klein-2026-07-01-S?
        **To:** all
        **Sent:** 2026-07-01 12:08

        ---

        TASK (owner): explore other creative tournament-related new properties, even if not directly LRC-related.

An invariant safari on the cube model (exhaustive n=3..6), building on last session's flip-rank. Three findings.

(A) THE PACKING/COVERING ASYMMETRY (the most novel). Define the RAINBOW NUMBER R(n) = the max dimension of a subcube (fix arcs, free k) whose 2^k completions are ALL in DISTINCT iso classes -- the PACKING dual of the flip-rank rho(n) (the COVERING minimum, S71). Bounds: R <= floor(log2|G_n|) <= ceil(log2|G_n|) <= rho.
 VERIFIED exhaustive:
   n:        3   4   5   6
   |G_n|:    2   4  12  56
   R(n):     1   2   3   5   = floor(log2|G_n|) EXACTLY
   rho(n):   1   2   4   7   ( > ceil(log2) starting at n=6 )
   gap rho-R: 0  0   1   2   (GROWS)
 So the iso classes PACK at the information floor but COVER above the ceiling: covering every class is strictly HARDER than packing distinct ones under the S_n action. The group's folding OBSTRUCTS covering (it forces collisions you must absorb -- you can't separate what the group glued) but AIDS packing (fewer distinct things, and freedom to route around collisions). This reframes the S71 n=6 transition: at n=6 the covering problem first outgrows its info floor, while packing never does. The gap rho-R is a cheap diagnostic for 'how much symmetry is in the way.'

(B) max |Aut(T)| over iso classes = 3, 3, 5, 9 for n=3..6 (orbit = n!/|Aut|).

(C) THE SKEW-SPECTRUM IS A VERY WEAK INVARIANT (a clean negative). The eigenvalues of S = T - T^T (skew-symmetric +-1, so +- i mu_k) give only 1, 2, 2, 6 DISTINCT spectra among |G_n| = 2, 4, 12, 56 classes -- almost all tournaments are cospectral. Recording as a dead-end invariant: the T-T^T spectrum sees little beyond coarse cycle data.

(D) THE QR-TRIANGLE DESIGN (tournament <-> design bridge). The directed 3-cycles of a tournament form a 3-uniform hypergraph. For the quadratic-residue (doubly-regular) tournament on q = prime = 3 mod 4 vertices (i->j iff j-i is a QR mod q), this hypergraph is a 2-(q, 3, (q+1)/4) DESIGN: every vertex-pair lies in exactly (q+1)/4 directed triangles. VERIFIED q = 3, 7, 11, 19, 23 (lambda = 1, 2, 3, 5, 6). The merely-rotational tournament {1,2,3} on 7 is NOT a design (pair-coverages 1,2,3). So the 2-design property CHARACTERIZES doubly-regularity. (This QR tournament is exactly the one on the LRC atoms, HYP-3802; and it maximizes the OCF.)

HONEST: R(n)=floor(log2) is verified for n<=6 and conjectured general; (D) is likely known (doubly-regular tournaments <-> skew-Hadamard matrices <-> 2-designs; cf. the repo's prior five_cycle_3design_s24e) but freshly verified with the exact lambda=(q+1)/4; (C) is a documented negative. Exploratory per the owner's invitation, not tied to an LRC proof. The most novel is (A) -- the rainbow/flip-rank duality quantifying how the S_n action makes covering harder than packing.

Files: 04-computation/tournament_invariant_safari_klein.py (+out); 05-knowledge/hypotheses/HYP-3804-tournament-invariant-safari.md; 07-reflections/packing-is-easy-covering-is-hard.md.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
