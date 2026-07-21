        # Message: mac-mini-2026-07-20-S152: SINGLE-STRADDLE GMC(2) PROVED IN CLOSED FORM -- E[P^{j*m0}] = C(j*m0,j) alpha^j L(Q^j) reduces the whole moment tower to the already-closed radial Laplace layer (THM-1675/1695); moment-count bound SUBSUMED; sole residual for full GMC(2) is now the localisation lemma. THM-1760

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 19:47

        ---

        OWNER: work the closed-form uniform proof.

I found a genuine reduction, not more verification -- and it PROVES the single-straddle case.

THE TOWER IDENTITY (proved for all j).
A single straddle is P = alpha Z^p + W q(V), with V = ZW = |Z|^2 and q(V) = sum_i beta_i V^{a_i} (charge +p carried by one term, charge -1 carried by r terms at radial degrees a_i, since W V^{a_i} = Z^{a_i} W^{a_i+1} has charge -1). Balance in P^m forces j copies of alpha Z^p and jp copies of W q(V), so m = j(p+1) = j*m0 and moments live ONLY at multiples of m0 = p+1. The balanced part of P^{j*m0} is the single multinomial term C(j*m0, j)(alpha Z^p)^j (W q)^{jp}, and since Z^{jp} W^{jp} = V^{jp} and E[V^k] = k! = L(V^k):
    E[P^{j*m0}] = C(j*m0, j) * alpha^j * L( Q(V)^j ),   Q(V) := V^p q(V)^p.
Verified exactly j=1..3 for p=1,2 and r=1,2,3; the balanced-term derivation is general.

THE REDUCTION -- A PROOF.
On the nullcone with alpha != 0, the scalar C(j*m0,j) alpha^j is nonzero, so L(Q^j) = 0 for ALL j >= 1. By THM-1675/1695 -- the radial Laplace nullcone, L(Q^j)=0 for all j implies Q == 0, CLOSED for complex Q via the Cauchy transform -- we get Q == 0. Since Q = V^p q^p and V^p is not identically 0, q == 0, so every beta_i = 0.
=> SINGLE-STRADDLE GMC(2) IS PROVED IN CLOSED FORM. No Groebner, no per-pattern bound: the whole moment tower factors through one radial polynomial Q, and the already-closed radial layer finishes it.

THE MOMENT-COUNT BOUND IS SUBSUMED.
THM-1740's per-straddle level r*m0 is NOT a separate thing to bound. The identity shows E[P^{j*m0}] IS L(Q^j) up to a nonzero scalar, so the number of levels needed to force Q == 0 is exactly the radial certifying level for Q (a degree-p(1+max a_i) polynomial), supplied by THM-1675/1695. So HYP-8540's bound, for a single straddle, is a CONSEQUENCE of the closed radial layer, not an independent conjecture -- the 2, 3, 7-flavoured moment counts were the radial layer's certifying levels all along.

MULTI-STRADDLE LOCALISES.
Witness aZ^2 + bW + cW^3 (straddles m0 = 3 and 5): E[P^3] = 6ab^2, E[P^4] = 288a^2bc, E[P^5] = 7200a^3c^2, E[P^6] = 360a^2b^4. With a != 0: E[P^3]=0 => b=0, then E[P^5]=0 => c=0 -- bottom-up (klein THM-1700), straddles firing at their own levels without interfering, exactly as THM-1745's max law predicts.

HANDOFF -- GMC(2) IS DOWN TO ONE GAP.
HYP-8590, THE LOCALISATION LEMMA, is now the SOLE residual for full GMC(2) (given span-2 THM-1600, complex radial THM-1695, and single-straddle THM-1760 are all closed). To prove: at the dominant straddle's return level (max r*m0), its radial factor L(Q_dom^j) is ISOLATED -- E[P^{j*m0_dom}] = (nonzero scalar) * L(Q_dom^j) + (terms from strictly lower straddles that vanish independently) -- so L(Q_dom^j) = 0 for all j and THM-1760's single-straddle reduction applies. Then peel the dominant straddle (coefficients forced to 0) and induct downward: a straddle descent. The LRC resonance (THM-1745) points at the tool: the dominant straddle is the LAST to return, so at its level it is ALONE -- a first-return / covering argument.
SCOPE: single-straddle GMC(2) is PROVED (the tower identity is general, the reduction rigorous, the radial layer itself proved). Multi-straddle is verified on one two-straddle witness; the localisation lemma is stated, not proved. The radial-reduction trick is specific to the Gaussian moment functional L and does not obviously transport to TNC's constant-term functional (HYP-8505 stays separate).

Artifacts: THM-1760; 04-computation/gmc2_uniform_proof_reduction_macmini_S152.py (+out); HYP-8540 subsumed, HYP-8590 (localisation lemma) new.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
