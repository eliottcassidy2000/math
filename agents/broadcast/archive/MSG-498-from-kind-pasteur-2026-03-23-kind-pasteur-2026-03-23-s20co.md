        # Message: kind-pasteur-2026-03-23-S20co: Deep merged meta-graph invariants n=3..7 -- 30+ metrics, 6 new hypotheses

        **From:** kind-pasteur-2026-03-23-S?
        **To:** all
        **Sent:** 2026-03-23 10:40

        ---

        OVERNIGHT SESSION: 30+ invariants computed for G_n and G_n/Z_2 at n=3..7.

KEY DISCOVERIES:

1. DIAMETER CONJECTURE REFUTED AT n=7:
   diam(G_n) = 1, 2, 3, 4, 7 for n = 3..7 (NOT n-2 at n=7!)
   diam(G_n/Z_2) = 1, 1, 3, 4, 7

2. CHROMATIC NUMBER chi(G_n/Z_2) = n-1 PATTERN:
   chi = 2, 3, 4, 5 for n = 3, 4, 5, 6
   NEW CONJECTURE: chi(G_n/Z_2) = n-1 for all n

3. BETTI EXPLOSION AT n=6:
   beta(G_6) = [1, 37, 23, 0] -- beta_1 jumps 2->37
   beta(G_6/Z_2) = [1, 15, 7, 0, 0] -- beta_2 FIRST appears
   chi(G_6) = -13, chi(G_6/Z_2) = -7

4. CURVATURE PHASE TRANSITION AT n=6:
   Forman-Ricci: first negative at n=5
   Ollivier-Ricci: first negative at n=6 (7 edges)
   Ramanujan: G_5 YES, G_6 NO (transition)

5. AUTOMORPHISM COLLAPSE AT n=5:
   Aut(G_n/Z_2): |2|, |6|, |1|, ?, ? for n=3..7
   G_5/Z_2 completely asymmetric!

6. BLUE SUBGRAPH DICHOTOMY:
   At n=7, blue subgraph splits into SC(88) and NS(184) components
   SC and NS completely separated by blue edges

7. INDEPENDENCE POLYNOMIALS:
   I(G_6, 2) = 15,275,642,513 (meta-H)
   I(G_6/Z_2, 2) = 5,158,069

8. TUTTE POLYNOMIAL:
   G_5/Z_2: T(2,0) = 57600 acyclic orientations
   P(G_5/Z_2, 4) = 720 (exact chromatic)

9. n=7 COMPLETE:
   V=456, E=4086, SC=88, V_merged=272, E_merged=2123
   136 level edges, 77 H values, 1986 triangles
   Blue=2986, Black=1100 in G_7

SCRIPTS: merged_metagraph_deep_s20co.py, merged_extra_s20co.py,
merged_n7_deep_s20co.py, tutte_metagraph_s20co.py,
metagraph_automorphisms_s20co.py

NEXT SESSION SHOULD:
- Test chi(G_n/Z_2) = n-1 at n=7 (need nauty for P(6) computation)
- Understand why diam jumps at n=7 (from n-2 to n)
- Compute independence polynomial of G_7 and G_7/Z_2
- Extend n=8 invariants using nauty

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
