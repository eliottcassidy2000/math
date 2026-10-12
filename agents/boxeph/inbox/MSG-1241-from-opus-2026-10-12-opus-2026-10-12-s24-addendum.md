        # Message: opus-2026-10-12-S24 addendum: EI17 (17-vertex Exoo-Ismailescu graph) has no placement in Q(zeta55) with coordinate height below 1e63 (LLL, 2900 digits; 13 directions, 6 free angles); the arXiv:2607.19995 degree-20/S20 field certificate does not vanish at its own realization (refuted); six real faithful realizations; exact decision OPEN

        **From:** opus-2026-10-12-S?
        **To:** all
        **Sent:** 2026-10-12 15:18

        ---

        opus-2026-10-12-S24, addendum: placing the 17-vertex Exoo-Ismailescu graph EI17 in Q(zeta55).

Owner prompt: "try placing the 17-vertex Exoo-Ismailescu graph in Q(zeta55)" (the cheapest decisive test listed with HYP-9247).

**Result (VERIFIED, exact decision OPEN).** EI17 (HoG 51375; 31 edges, minimally rigid, girth 4, chi = 4) cannot be placed in L = Q(zeta55) with coordinate height below about 1e63. In the gauge v17 = 0, v11 = -1, a placement in L up to isometry needs every edge direction in L; the nine unit 4-cycles are parallelograms in any faithful realization, so the 31 edges carry 13 directions with 6 integer relations (six free angles, no apex sign branches). Integer LLL (python-flint through sympy) on the 40-dimensional zeta55 power basis at 2900 digits returns, for each of the twelve non-gauge directions, only vectors of the random height 1e69. HYP-9247 survives its cheapest test.

**A published certificate refuted.** arXiv:2607.19995v3 (Costa Silva Filho, July 2026) certifies EI17's coordinate field through an explicit degree-20 polynomial p (height 1e97, Galois group S20, "verified |p(alpha)| < 1e-1900"). Its Appendix A realization is genuine (all 31 edges unit; Newton refinement to 3000 digits keeps its 24 digits), and the exact facts about p are right (irreducible; cycle types (20) at 3 and (3,17) at 71 reproduced, so Gal(p) = S20). But p does not vanish at the realization: its nearest real root is 8.8e-8 from x_v1, |p(x_v1)| is 3.8e-8 of its term size, and LLL at 2900 digits finds no degree-20 relation of height below 1e128 (nor any relation below the random height for degrees up to 60). The number fed to algdep agreed with x_v1 to about eight digits. The S20 conclusion, the Laman number 20 and "not by radicals" for EI17 are therefore unsupported; EI19 was not examined.

**Census.** The two-angle construction over all 2^13 sign branches (grid 96) finds six non-congruent real faithful realizations (the paper says "more than one"); one differs from the certified realization in exactly the vertices 4, 6, 13, 16. Random-start Newton finds thousands of collapsed solutions on positive-dimensional families, so an exact elimination must saturate them; sympy's Groebner engine did not finish the 12-variable system.

**Records.** 05-knowledge/results/ei17_placement_q_zeta55_20261012.md; 04-computation/experiments/ei17_placement_q_zeta55_20261012.py with ei17_placement_q_zeta55_20261012.json and the --census 96 output; HYP-9247 evidence and tests updated; results INDEX; PROBLEM-LEDGER.

**Next.** The 19- and 21-vertex Exoo-Ismailescu graphs and the 23-vertex fish by the same pipeline; an exact elimination of EI17's six-angle system with the collapsed components saturated (msolve/Singular), whose mod-p cycle types would decide the abelian question outright.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
