        # Message: mac-mini-2026-07-20-S147: COMPLEX RADIAL LAYER OF GMC(2) CLOSED (Cauchy transform + d-zbar, no monodromy; charge-0 layer done for general non-Hermitian P) and the CROSS-SHELL COUPLING IS A CHARGE DESCENT down the charge Newton polygon -- a sharply-scoped PATH to full GMC(2). THM-1695

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 16:06

        ---

        OWNER: 'work the complex radial and the cross shell descent.'

PART A -- COMPLEX RADIAL, CLOSED (cleaner than the real jump, and no monodromy).
THM-1675 closed L(p^m)=0 forall m => p=0 for REAL p via the real-axis jump. For COMPLEX p there is a cleaner argument:
  Psi==0 (THM-1665) gives h(t) = int_0^inf e^{-v}/(1-t p(v)) dv == 1. With z = 1/t, h = z C_mu(z), where C_mu(z) = int_0^inf e^{-v}/(z - p(v)) dv is the CAUCHY TRANSFORM of mu = p_*(e^{-v}dv). So C_mu(z) == 1/z = C_{delta_0}(z).
  C_mu and 1/z are L^1_loc on C and agree OFF the arc {1/p(v) : v>=0} -- a MEASURE-ZERO curve -- hence agree as DISTRIBUTIONS on C. Applying d_zbar (with d_zbar 1/(z-w) = pi delta(z-w)): pi mu = pi delta_0, so mu = delta_0. But for a nonconstant polynomial mu({0}) = meas{v>=0 : p(v)=0} = 0 != 1 -- contradiction. So p is constant, and L(p) = p = 0 gives p == 0. QED.
No monodromy needed (unlike DvdK's Theorem 2, which needs the critical-value/Nilsson-class argument for the general CT case). Verified: 0 nullcone members among 50610 complex polynomials (deg<=3, Gaussian-integer coeffs), and z*C_mu(z) != 1 measured directly on nonconstant samples.
=> HYP-8350 IS NOW FULLY CLOSED (real by THM-1675, complex here), so the charge-0/radial layer of GMC(2) is done for GENERAL, NON-HERMITIAN P.

PART B -- THE CROSS-SHELL COUPLING IS A CHARGE DESCENT.
E[P^m] = sum over charge-balanced m-tuples of L( s^{(sum|k_i|)/2} prod lambda_{k_i} ). After L, a tuple's factorial argument is sum_i ( |k_i|/2 + deg lambda_{k_i} ) = sum_i phi(k_i), phi(k) := |k|/2 + deg lambda_k.
THE DOMINANT SHELL as m -> inf maximises sum phi(k_i) over balanced tuples -- the TOP EDGE OF THE CHARGE NEWTON POLYGON. For a SYMMETRIC top (both +-K present, deg lambda_{+-K} = d) the max is m/2 copies of each, giving the dominant term
    C(m, m/2) (lead lambda_K * lead lambda_{-K})^{m/2} (m(K/2 + d))!,   NONZERO.
So E[P^m] = 0 for all large even m FORCES lead lambda_K * lead lambda_{-K} = 0 -- a DESCENT STEP that drops one top charge's leading coefficient. Verified: whenever both top charges are present (ab != 0), some E[P^m] != 0 for m <= 8.
ITERATING shrinks the charge range until either:
  - the top charge is ONE-SIDED (only +K or only -K) -- where no balanced tuple uses it at full weight, so E[P^m] = 0 TRIVIALLY, and such charge-one-sided P are MATHIEU-ZHAO-HARMLESS (for fixed Q, QP^m cannot reach charge 0 for m >> 0); OR
  - the support collapses to CHARGE 0 ONLY -- where E[P^m] = L(lambda_0^m) = 0 forces lambda_0 = 0 by the now-closed Part A.
CONSEQUENCE (conjectural, pending the residual): the GMC(2) NULLCONE IS EXACTLY THE ONE-SIDED-CHARGE POLYNOMIALS, which are MZ-harmless -- i.e. GMC(2) IS TRUE. The descent is the charge analogue of the TNC coefficient ladder (THM-1610), and the one-sided terminus is DvdK's degenerate case (THM-1630).

HANDOFF.
HYP-8350 is FULLY CLOSED -- the radial layer of GMC(2) is settled for general P.
HYP-8470 is now a SHARPLY-SCOPED PATH to full GMC(2), with THREE residuals:
  (i) rigorous TOP-SHELL DOMINANCE -- that the top-edge factorial strictly beats all lower shells, INCLUDING the case where several tuples tie at the top edge and their leading coefficients might cancel. This is a per-shell Watson estimate (THM-1665's method) and is the hardest of the three.
  (ii) the ASYMMETRIC-TOP LP (|K| != |K'|), where the balancing tuple is not simply +-K and the top edge of the charge Newton polygon has to be read off a linear program.
  (iii) the uniform-in-Q m >> 0 bound making 'one-sided => MZ-harmless' precise.
THESE THREE PLUS THE CLOSED RADIAL LAYER WOULD PROVE GMC(2). That is the natural full-attack target, now sharply scoped.
SCOPE: Part A is a full proof (closes HYP-8350). Part B is a mechanism and a path, NOT a proof of GMC(2), which remains open but is reduced to the three residuals plus the MZ-harmlessness bookkeeping.

Artifacts: THM-1695; 04-computation/gmc2_complex_radial_and_charge_descent_macmini_S147.py (+out); HYP-8350 fully closed, HYP-8470 -> charge-descent path.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
