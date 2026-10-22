        # Message: mac-mini-2026-10-22-sierp: loop expansion of Hamiltonian paths - one loop = Wormald gain, parameter-free 2-/3-loop laws, doubly regular kappa = e(1 - 3/(2n) + 25/(8n^2))

        **From:** mac-mini-2026-10-22-S?
        **To:** all
        **Sent:** 2026-10-22 05:12

        ---

        mac-mini-2026-10-22-sierp: close-out (tournament side of the owner's 3^k+1 / Sierpinski / Collatz prompt; Session29 has the Collatz cylinders).

WHAT CHANGED (all pushed to main; note 05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md, sections 0-12; scripts in 04-computation/experiments/sierpinski_ocf_20261022/)

1. THM-4757 (PROVED, audited twice): a loop expansion of Hamiltonian paths.
   - kappa = H 2^(n-1)/n! = sum_L N_L(S)/(n)_L, where N_L is the sum over S-weighted directed linear forests. The N_L are THM-059's W-coefficients.
   - Moebius clusters are graded by excess. Regular T has no tree level.
   - One loop = log of Wormald's blow-up gain det(nI-S)/|det(inI-S)|, exactly.
   - Two loops (figure-eight, theta, dumbbell) in closed form via G = -S^2(I+S^2)^-1, H = S^2(I-S^2)^-1, O = S(I-S^2)^-1 (S = S/n).
2. HYP-9365: the parameter-free two-loop law kappa = rho(1 - y/n).
   - corr 0.998 on all 208 circulants at n = 25; 0.9975 on random regular n = 21.
   - Finds the H-maximal circulant at every n <= 25.
   - The doubly regular limit is e(1 - 3/(2n)); for R_n, y -> tan1 tanh1 + tan^2(1)/2 - 1.
3. HYP-9367: the three-loop law.
   - A general loop engine (15 three-loop cores) is validated against brute force. The residual after the second order is O(n^-3) on Paley, rotational, towers, all circulants and random regular tournaments.
   - Universal doubly regular curve: kappa = e(1 - 3/(2n) + 25/(8n^2) + O(n^-3)), with q^2 X_2 -> 29.
   - Four loops seen directly: the 5-census of doubly regular 23-tournaments has one free parameter, and H is linear in the transitive-5-set count across 12 non-Paley classes, to 3.7e-7.
4. Doubly regular universality. All propagators lie in span{I, J, S} and three-loop cores have <= 4 vertices.
   - H(D_19)/H(P_19) = 1 + 3.2e-5; three D_23 classes are within 1.9e-5 of P_23.
   - The audit adds a 5-census-equal pair that agrees to 2.8e-8.
5. HYP-9366: R_n and the transitive tournament are cospectral twins, so they have the same rho. F(T,-1) is a tangent-number OCF. The n mod 4 wobble of H(R_n) decays like (2/pi)^n.
6. Exhaustive data:
   - H(R_33), H(R_35);
   - the max H over all 1,495,297 regular 13-tournaments is 3719893 (a converse pair, as in HYP-9357; R_13 ranks 781st). The two-loop predictor ranks it 7th; rho alone ranks it 1,476,889th.
7. Earlier the same day:
   - THM-4756 (proves THM-086 and THM-094: F(T,x) = H A_n(x) mod 3);
   - THM-4754 (3-adic Sierpinski tournament, Collatz arc-flip law);
   - THM-4755 (blow-up gains, tower product with (28/26)^3);
   - HYP-9364.

PROOF STATUS: THM-4754..4757 PROVED and audited (audit 1: THM-4754..4756, HYP-9364; audit 2: THM-4757, HYP-9365; no mathematical error; credits and wording corrected). HYP-9364..9367 are open asymptotic statements with FINITE-EXACT evidence. c_2 = 25/8 is DERIVED (exact limits; the three-loop limit 29 by a 7-digit fit).

NEXT OBLIGATIONS:
- a rigorous O(1/n) bound on the excess >= 1 clusters (a proof of the principle);
- a closed form for the three-loop kernels and an exact rational X_2 for doubly regular tournaments;
- a non-regular / even-n theory (path trees give -delta^T(I+S^2)^-1 delta; the second order is needed);
- whether 3719893 is the global maximum at 13 (OEIS A038375 has data only to n = 11).


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
