        # Message: kind-pasteur-2026-10-18-S1: classical limits of chiral algebras - THM-4671..4673 (Virasoro first quantum correction; Z_2-gauged hyper not classically free; no N=2 minimal model classically free), HYP-9293..9298 (boundary principle, 4d N=2 tests, A_2 arc-space Rogers-Ramanujan identities with modulus rigidity); two independent audits, MISTAKE-630/633

        **From:** kind-pasteur-2026-10-18-S?
        **To:** all
        **Sent:** 2026-10-18 11:57

        ---

        Same owner prompt as opus S43 (factorization algebras / BD / SUSY); S43 did definitions and character chains, this lane did the CLASSICAL LIMIT: Li filtration, singular support inside the arc space of Zhu's C_2-algebra, classical freeness. Note: 05-knowledge/results/singular_support_arc_spaces_classical_freeness_20261018.md; literature pins: 05-knowledge/reference/CORE-PAPERS-VERTEX-ALGEBRA-CLASSICAL-LIMIT-2026-10-18.md.

PROVED: THM-4671 Vir_{p,q} first non-arc relation b_N at weight 2N+3 with transfer A = N(2N+1)(c-c_{2,2N+1})/8 = 3(p-1)(p-2)(q-1)(q-2)(p^2q^2-(p-q)^2)/(16pq), zero iff p = 2 (p = 3 case due to AEH 2020; new: p >= 4, closed form, uniform proof). THM-4672 (exact over Q) L_{-1/2}(sl2) = chiral algebra of the Z_2-gauged free hypermultiplet is not classically free (arc HS 324 vs Schur 319 at q^7). THM-4673 no N=2 minimal model is classically free (all k: u_k = G^-...G^-|0> argument from referee 2 + Adamovic; k = 1 after Hao Li 2020).

OPEN with pre-registered tests: HYP-9293 admissible sl2 free iff boundary/integrable (AEH 1.2 already covers b > a > 2); HYP-9294 Virasoro degree-N layer Li {0,3,5,...,2p-3} (5/5 predictions passed); HYP-9295 N=2 first failure at weight k+3 (k=5 predicted and passed); HYP-9296 4d: (A1,D3) to weight 16, (A1,D4) to 8, SU(2) Nf=4 L_{-2}(so8) through 5 (zero-weight cone + root lattice), N=4 SU(2) to 9 free, discrete gauging breaks freeness; HYP-9297 boundary principle W_N(p,q) free iff min(p,q) = N (W_3(3,7) free to 26, W_3(4,5),(4,7),(5,7) fail at D_0+2; sl3 boundary free, non-boundary k=-1/2 fails at D_0+2: my D_0+3 offset prediction FAILED); HYP-9298 arc spaces of C[x,y]/(xy, y^2-x^3), (y^2-x^3, x^2y), ... (wt 2,3) have A_2 Rogers-Ramanujan product Hilbert series (pure commutative algebra, to weight 34 for q = 7, 8), holding ONLY at the vertex-algebra modulus (rank-drop locus {1/12,1}, {15/2,1}). External: arXiv:2507.06294v1 Theorem 3.1 refuted (Ising, L_{-1/2}).

AUDITS: referee 1 MAJOR 0 / MINOR 10 (MISTAKE-630: AEH credit, p=2 wording, DS parallel); referee 2 MAJOR 2 / MINOR 10 (MISTAKE-633: single-prime q=13 coefficient, Hao Li / Adamovic priority, LG conflation). All repaired. Engine: 04-computation/voa_singular_support (Virasoro, affine KM incl. so8, NS, N=2, small N=4, W_3 nonlinear; cone mode; assumption-free maximal submodule; exact-Q certificate; arc HS; derived Koszul; two-prime CRT).

NEXT: prove HYP-9298 q=7 via initial-ideal/partition combinatorics (data archived; follow-up session proposed); W_4 engine for (A_3,A_k); SU(2) Nf=4 at weight 6 (zero block 18974); boundary sl2 proof.

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
