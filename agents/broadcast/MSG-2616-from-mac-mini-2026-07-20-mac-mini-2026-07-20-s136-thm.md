        # Message: mac-mini-2026-07-20-S136: THM-1540 -- the 2D NULLCONE CONJECTURE stated, TWO-CHARGE CASE PROVED, GMC(2) reduced to one open case (two-sided with >=3 charges), and the search finally made CONTROLLABLE. Plus: n=2 IS Duistermaat-van der Kallen, Gaussian-averaged

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 11:46

        ---

        OWNER: 'aim to finish GMC(2) by finishing the stronger 2 dimensional nullcone conjecture.' That reframing was the right move for a reason beyond strength -- see the methodological point below.

THE CONJECTURE. With N(E) := {P : E[P^m] = 0 for all m >= 1} and charge c = deg_Z - deg_W:
    NULLCONE CONJECTURE (n=2):  N(E) = {P : all charges >= 1} u {P : all charges <= -1}
                                     = exactly the STRICTLY ONE-SIDED charge support.

WHY IT IS THE BETTER TARGET -- the methodological point, and it matters.
A direct GMC(2) search has NO POSSIBLE POSITIVE CONTROL: a positive control would BE the counterexample being sought. That is why I refused to lean on the S135 sweep. The NULLCONE conjecture makes a POSITIVE PREDICTION -- one-sided P ARE in the nullcone -- so the machinery can be validated BEFORE its negatives are trusted. All five controls pass. That converts an uncontrolled sweep into meaningful evidence, and it is the whole gain from the reframing.

(A) THE EXACT POLAR REDUCTION.
With H_r(u) := P(ru, r/u) -- a Laurent polynomial in u whose EXPONENTS ARE EXACTLY THE CHARGES and whose r-degree is the total degree:
        E[P^m] = L( v -> CT_u( H_sqrt(v)(u)^m ) ),     L(v^k) = k!.
Verified exactly. CT_u is the charge-0 projection; L is the Gaussian average over the radius (r^2 ~ Exp(1)). So GMC(2) IS DUISTERMAAT-VAN DER KALLEN'S CONSTANT-TERM PROBLEM, GAUSSIAN-AVERAGED. That names the right context: DvdK is the model theorem of the whole Mathieu-Zhao area, and this says the Gaussian nullcone is a weighted average of Laurent nullcones.

(B) THE TWO-CHARGE THEOREM -- PROVED.
Let P = P_C + P_{-B} have exactly two charges C > 0 > -B. Write P_C = Z^C q(V), P_{-B} = W^B s(V) with V = ZW; g = gcd(B,C), B' = B/g, C' = C/g. Charge-0 requires k_C*C = k_{-B}*B with k_C + k_{-B} = m, forcing k_C = tB', k_{-B} = tC', m = t(B'+C') -- UNIQUE where it exists, IMPOSSIBLE otherwise. Hence
    E[P^m] = 0 unless (B'+C') | m,   and otherwise   = C(m; tB', tC') * L(F^t),  F = v^{CB'} q^{B'} s^{C'} != 0,
which is nonzero for large t by THM-1520's saddle lemma. SO A TWO-CHARGE P IS NEVER IN THE NULLCONE. Verified on 6 explicit P -- nonzero exactly at the predicted multiples of B'+C' (e.g. W^2+Z^3 gives 0,0,0,0,7200,0,0,0 with B'+C' = 5).

(C) THE GENERAL REDUCTION -- what the >= 3-charge case needs.
E[P^m] = sum_d [charge-0, degree-d coefficient of P^m] * (d/2)!, and the factorials weight the TOP degree overwhelmingly. Plotting the support in (charge, degree) coordinates, the leading term comes from the TOP EDGE of the Newton polygon where it crosses charge 0 -- and a single edge's monomials form a ONE-VARIABLE LAURENT polynomial. Two inputs: (i) the 1-variable Laurent nullcone lemma (CT(h^m)=0 for all m => h one-sided) -- 10320 two-sided h tested, no exceptions; (ii) DOMINATION of the leading term -- NOT established.

(D) EVIDENCE, NOW CONTROLLED. Two-sided P tested at support sizes 2/3/4/5: 400 / 10000 / 118800 / 886800, with ZERO in the nullcone, and the positive control passing.

(E) GMC(2) IS A TWO-LINE COROLLARY. P in N(E) => strictly one-sided => every charge of P^m is >= m => QP^m has charge 0 only if m <= deg_W(Q) => E[QP^m] = 0 for every m > deg_W(Q). QED

STATE OF PLAY -- I want to be exact about this. GMC(2) IS NOT FINISHED. Finished: the one-sided case (THM-1520), the two-charge case (B), the reduction of the rest (C), and the corollary (E). OPEN: two-sided P with THREE OR MORE distinct charges. That is the whole remainder.

HANDOFF -- and there is a piece of structure here worth acting on:
(i) HYP-8375 is the last case, with its route already isolated (Newton-polygon top edge -> one-variable Laurent). THE PLEASANT STRUCTURE: its missing input (ii) DOMINATION and HYP-8350's saddle-lemma gap are THE SAME ANALYTIC OBSTACLE APPEARING IN TWO PLACES -- in both, consecutive degree levels contribute in a ratio that is O(1) rather than decaying, so a genuine Laplace/saddle treatment is required instead of a crude bound. SOLVE THAT ESTIMATE ONCE AND IT LIKELY CLOSES HYP-8350, THM-1540(B) unconditionally, (C), the nullcone conjecture, AND GMC(2) TOGETHER. That makes HYP-8350 the single highest-leverage item in this entire thread, and I would put a session on it before anything else here.
(ii) HYP-8380: is the polar reduction (A) a known bridge? DEZ (arXiv:1506.05192) do not state it and I found no source, but DvdK is well trodden and someone may have written it down. Worth 20 minutes of literature search before anyone builds on it as novel. If it IS novel it is independently interesting -- it suggests transporting DvdK's n-variable theorem to GMC(2n) WHOLESALE, i.e. a route to the whole conjecture family rather than one case.
(iii) The 1-variable Laurent nullcone lemma is presumably classical (the Newton-polytope criterion; DvdK 1998 is its n-variable Mathieu-subspace form). FIND THE CITATION rather than reproving it -- 10320 tests say it is not in doubt, only unsourced. No priority is claimed for it here.

Artifacts: THM-1540; 04-computation/gmc2_nullcone_macmini_S136.py (+out); HYP-8375, HYP-8380.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
