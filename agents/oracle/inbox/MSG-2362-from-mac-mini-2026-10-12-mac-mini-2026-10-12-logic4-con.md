        # Message: mac-mini-2026-10-12-logic4 continuation: HYP-9255 RESOLVED (THM-4626: line fibres make dual fields locally nilpotent; equivalent to JC_2; audited D), kappa(13) record (alpha(G_13) >= 36, obstruction global, external 470-hour proof CITED), alpha(G_17) >= 57

        **From:** mac-mini-2026-10-12-S?
        **To:** all
        **Sent:** 2026-10-12 19:31

        ---

        Continuation of mac-mini-2026-10-12-logic4 (owner: keep going on HYP-9255 and the kappa(13) SAT run).

HYP-9255 RESOLVED (THM-4626, PROVED, audited D PASS WITH CORRECTIONS applied): for a Keller map in any dimension and any dual field V_j = d/dF_j, generic rank zero of the cokernel (every geometric component of the generic fibre of F_(i != j) is a line) <=> V_j locally nilpotent <=> K[x] = R_j[F_j] <=> V_j surjective <=> H^1_dR(K[x]/K[F_(i != j)]) = 0. Proof: THM-4620(4) puts K[x] inside L_j[F_j] with V_j = d/dF_j, so every polynomial is killed by a power of V_j. Hence 'a surjective dual field in dimension three forces invertibility' (HYP-9255) is EQUIVALENT to JC_2 (slice surface = A^2 by cancellation + Kambayashi; padding a planar counterexample gives the converse), and the torsion sub-question is settled (generic rank zero => cokernel zero). The dual-field hypothesis is essential: x d_x - y d_y + d_z is surjective, divergence-free, not nilpotent. Plane precursors: Stein 1989, Berson-van den Essen-Maubach 2001, Gurjar-Masuda-Miyanishi 2013. THM-4620's torsion caveat removed.

kappa(13) (note section 10): in-repo FINITE-EXACT facts - alpha(G_13) >= 36 (CP-SAT, 0.3 s; explicit 34-set), so chi_f <= 169/36 < 5 and no counting/ratio/Delsarte (= 52.4)/theta bound can decide 5-colourability; G_13 triangle-free; every union of 0 with at most ten of the twelve norm circles (<= 141 vertices) is 5-colourable, so the obstruction is global; cube timings (14 neighbours fixed: ~10 s per cube, ~4e5 cubes; 7 fixed: > 600 s) put a naive cube-and-conquer at ~1e3 CPU-hours. decalion89's g13_chi.md proves chi(G_13) = 6 with the largest-class split s = 34..36, lex-leader under all 4,732 automorphisms, cube-and-conquer, kissat 4.0.4, ~470 kissat-hours plus ~500 drat-trim-hours on dozens of cloud machines, DRAT-checked, proofs deleted, unrefereed: CITED. In-repo status stays kappa(13) in {5, 6}. Full-instance runs (kissat 4.0.4, CaDiCaL 1.9.5, CP-SAT symmetry level 2) and the maximal-independent-set enumeration for the largest-class decomposition were still undecided at close-out; see their .out files under 04-computation/experiments/logic4_20261012/ (RUN_STATUS lines; kissat and CaDiCaL 1.9.5 were left running single-threaded, stop with pkill -f kappa_full_solver.py). alpha(G_17) >= 57 FINITE-EXACT in-repo (CP-SAT 0.1 s; a 58-set search was UNKNOWN after 1200 s), so alpha(G_17) = 57 would give kappa(17) = 6.

Next obligations: a replication of chi(G_13) = 6 needs ~60 hours on eight cores (the largest-class decomposition scripts kappa13_mis_enum.py / kappa_largest_class.py are in place); kappa(17) (alpha >= 51, decalion89: 5-6); HYP-9256 (q > 256); HYP-9257 (cusp-curve constituent).


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
