        # Message: mac-mini-2026-07-20-S133: THM-1480 -- the GMC(4) counterexample VERIFIED (256 readings + real coords), EXPLAINED (it is the alternating binomial identity), and IMPROVED TO FOUR TERMS: P=(1+W)(Wbar-|Z|^2), Q=W. Plus the reframe: it LOCALIZES to n=4, it is NOT a corollary of ~JC(3)

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 10:41

        ---

        OWNER supplied an outside-source counterexample disproving the Gaussian Moment Conjecture at n>=4 and asked (a) did we already have this or anything stronger, (b) improve it.

(0) DID WE HAVE IT? NO -- AND THE REPO SAYS SO IN TERMS.
THM-1435 states: 'THE WITNESS IS NOT PRODUCED.' The repo held the de Bondt transport machinery (control-validated), a proof that the BCW step cannot be skipped, and the bracket 5 <= vcwd <= ~20. No witness, for GMC or for VC.
ONE GENUINE CONVERGENCE WORTH RECORDING: THM-1435 SS D found INDEPENDENTLY, from the Alpoge map, that HOMOGENEITY is the load-bearing obstruction. The GMC counterexample is ALSO essentially inhomogeneous -- its degree-1 and degree-3 homogeneous parts, taken alone, are NOT counterexamples (they are eventually zero, hence they SATISFY GMC). Two independent routes landing on inhomogeneity.

(A) VERIFIED EXACTLY, TWICE, TAKING NOTHING ON TRUST -- INCLUDING THE NOTATION.
All 4^4 = 256 assignments of {Z1,Z2,W1,W2} to the letters {Z,Zbar,W,Wbar} of two independent standard complex Gaussians were tested. EXACTLY 8 SURVIVE (a single orbit under Z<->W and conjugation), all giving E[P^m]=0 and E[QP^m]=m! for m=1..9, and NO reading gives a different growth law. An independent real-coordinate expansion in u1..u4 ~ N(0,1/2) reproduces 1,2,6,24,120 with zero imaginary part. It genuinely lives in FOUR REAL GAUSSIANS.

(B) THE MECHANISM IS ONE BINOMIAL IDENTITY.
Using C(m,k)*k!*(m-k)! = m!, both sides collapse:
    E[P^m]   = m! * sum_j (-1)^j C(m,j) = 0
    E[QP^m]  = m! * C(m,0)              = m!
MULTIPLYING BY Q MERELY SHIFTS THE SUMMATION INDEX BY ONE and cashes the boundary term. (1-Z) manufactures the alternating signs; the Gaussian pairing manufactures the factorials that cancel the binomial denominators exactly; Q collects the residue. Convergent re-derivation (a second route in-session reached the same), and it mirrors DEZ's own Prop 1.2.

(C) STRICT IMPROVEMENT -- FOUR TERMS, NOT SIX.
    P = (1+W)(Wbar - Z Zbar) = (1+W)(Wbar - |Z|^2),   Q = W
Identical m! growth, still cubic. THE LONE Zbar IN THE PUBLISHED SECOND FACTOR IS REDUNDANT. Exhaustive over the degree-<=3 complex-monomial box: ZERO counterexamples at support <= 3, many at support 4, so 4 IS MINIMAL IN THAT BOX -- and the box contains the published example, so this search has a WORKING POSITIVE CONTROL. Degree-<=2 with support <= 5 finds nothing, so cubic appears necessary (bounded, not proof).

(D) AN INFINITE FAMILY. Q = Z2^r gives E[QP^m] = (-1)^{r+1} m! sum_{j<r} (-1)^j C(m,j) -- m! times a degree-(r-1) polynomial in m. Verified r=1..4: m!, m!(m-1), m!(1-m+C(m,2)), degree-3. Every r >= 1 breaks GMC; the published case is r=1.

(E) RIGIDITY. Deforming (1-Z) -> (1-lambda Z) gives E[P^m] = m!(1-lambda)^m, so LAMBDA = 1 IS FORCED EXACTLY. The construction sits at the unique zero of the alternating sum -- not part of a continuous family.

(F) THE OBVIOUS ROUTE TO n=3 IS STRUCTURALLY BLOCKED. |Z|^2 and a real X^2 both have mean 1, but E[|Z|^{2k}] = k! is EXACTLY what cancels C(m,k), while E[X^{2k}] = (2k-1)!! is not. The analogous sum evaluates to 0,1,0,9,0,225 -- it does not vanish.

(G) THE REFRAME THAT MATTERS -- AND I HAD IT WRONG AT SESSION START.
GMC is Derksen-van den Essen-Zhao, arXiv:1506.05192, Israel J. Math. 219 (2017) 917-928. Their Thm 1.6 is ONE-WAY and stated GLOBALLY: GMC(n) for all n ==> JC(n) for all n. So ~JC(3) yields only 'GMC(N) is false for SOME finite N', with NO BOUND ON N. The explicit example LOCALIZES THE FAILURE TO n=4, which the published implication chain cannot do. IT IS A NEW LOGICAL FACT, NOT A WITNESS FOR A COROLLARY.
Sharpness is genuinely OPEN: GMC(1) is TRUE (DEZ Prop 4.2); general GMC(2) is OPEN (proved only for homogeneous P, Cor 4.4); GMC(3) is unknown.

PROVENANCE. The example is written in DEZ'S OWN NOTATION (their Prop 3.2 defines exactly W_j = (X_j - iY_j)/sqrt2, Z_j = (X_j + iY_j)/sqrt2) and is structurally the Gaussian transplant of DEZ'S OWN Prop 1.2 -- their half-disk counterexample to the Moment Vanishing Conjecture (P = (x+iy)^2, Q = x+iy, integral of P^m = 0 for all m but integral of QP^m != 0 for all m): same shape, same 'fails for every m'. A literature sweep found NO prior instance of the Gaussian version. The only 2026 item is Zihan Zhang's expository page (20 July 2026), which derives only the non-localized 'GMC(N) false for at least one finite N'.

FAILED CONTROL, REPORTED AS SUCH. The n=3 search this session is INCONCLUSIVE AND CARRIES NO EVIDENTIAL WEIGHT: its n=4 POSITIVE CONTROL also found nothing, because the real-coordinate expansion of (1+W)(Wbar-|Z|^2) has more terms and non-unit coefficients and falls outside the box. A failed control invalidates the negative result. Only (F)'s structural obstruction stands. I am flagging this loudly because a negative-looking table is exactly the kind of thing that gets cited later as evidence.

TWO SELF-CORRECTIONS:
1. ATTRIBUTION, SECOND ITERATION. My S129 amendment to THM-1300 said 'co-credit is owed to Akhil Mathew'. That OVER-CORRECTS. The primary wording thanks Akhil for ASKING THE QUESTION and credits CLAUDE FABLE with producing the example; the identification of 'akhil' as Akhil Mathew comes from a secondary post, not a primary source. Accurate: the map is Alpoge's announcement, produced with Claude Fable, with Akhil credited for posing the question. I have now mis-stated this TWICE IN OPPOSITE DIRECTIONS; THM-1300 now carries a standing note that further edits must quote the primary post verbatim rather than paraphrase.
2. MY 'HEAT SHADOW' FRAMING WAS TOO STRONG. The DEZ functional IS the heat semigroup at the origin ((e^{Delta/2} w^a z^b)(0) = a! delta_ab, verified). But GMC's hypothesis is the WHOLE SERIES (e^{Delta/2}P^m)(0) = 0 while VC's is the SINGLE TERM Delta^m P^m = 0. They coincide only for deg P = 2; for the JC-relevant quartic, VC's hypothesis implies GMC's and NOT conversely, so GMC's hypothesis class is strictly LARGER. Neither DEZ nor Zhao states this relation -- it is unsourced folklore and is flagged as such in canon.

HANDOFF -- three:
(i) HYP-8330 (SHARPNESS). GMC(2)-general and GMC(3) are the open cases. The sharp question: IS THE ALTERNATING-BINOMIAL COLLAPSE THE ONLY MECHANISM producing E[P^m] = 0 for all m? Every known instance -- DEZ's Prop 1.2 half-disk example and this one -- is that same identity in disguise. If it is the only mechanism, its dependence on k! may force an EVEN number of real variables, which would make n=4 sharp and settle GMC(3) affirmatively. ANYONE RETRYING THE n=3 SEARCH MUST CARRY A WORKING POSITIVE CONTROL -- expand (1+W)(Wbar-|Z|^2) in real coordinates and confirm the box contains it before trusting any negative.
(ii) HYP-8335 (THE r-FAMILY). Higher r vanishes for an initial stretch of m (r=2 at m=1; r=3 at m=1,2; r=4 at m=1,2,3), so it probes the 'm >> 0' quantifier far better than r=1, which fails immediately. Is there a construction failing FIRST at arbitrarily large m? That would show the quantifier does real work rather than being cosmetic.
(iii) The m! here has NO established connection to THM-1460's (n-1)! arborescence count or any other factorial in this repo. Do not cross-link them.

Artifacts: THM-1480; THM-1300 attribution note (second iteration); 04-computation/gmc_counterexample_verify_macmini_S133.py, gmc_mechanism_and_family_macmini_S133.py, gmc_three_gaussian_search_macmini_S133.py (+outs); HYP-8330, HYP-8335.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
