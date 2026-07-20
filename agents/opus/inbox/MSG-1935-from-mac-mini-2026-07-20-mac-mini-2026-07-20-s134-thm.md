        # Message: mac-mini-2026-07-20-S134: THM-1500 THE GMC MASTER THEOREM -- both counterexamples are ONE construction P=(1+Z)(W-g(Z)U), the m! is FORCED by Lagrange inversion (f=1/(1+s) is unique), g=((1+s)^{2/d}-1)/s recovers 2+s at d=1 and 1 at d=2, and n=3 is MINIMAL with n=2 unreachable. GMC(N) false for all N>=3

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 10:54

        ---

        OWNER supplied a claimed GMC(3) counterexample (5-term quartic in 3 real Gaussians) and asked for the better next question on 'binomial collapse', aiming at GMC(2)/(3)/(4).

(A) THE GMC(3) COUNTEREXAMPLE IS CORRECT.
Z = (X+iY)/sqrt2, W = Zbar, U = T^2/2;  P = (1+Z)(W - (2+Z)U),  Q = Z.
Verified EXACTLY two independent ways: the complex formalism (E[Z^a W^b] = a! delta_ab, E[U^k] = (d/2)_k) for m = 1..10, and a full REAL-COORDINATE expansion in a,b ~ N(0,1/2), T ~ N(0,1) -- 10 monomials -- for m = 1..6, giving E[P^m] = 0 and E[QP^m] = 1,2,6,24,120,720 with ZERO imaginary part. It genuinely lives in three independent real standard Gaussians.

(B) THE CONTROLLING SERIES, AND WHY IT IS UNIQUELY FORCED.
Take the construction in general: P = (1+Z)(W - g(Z)U), Q = Z, with U >= 0 INDEPENDENT of the complex pair (Z,W). Using E[W^r F(Z)] = r! [s^r] F(s) and writing Phi for the MOMENT GENERATING FUNCTION of U, ONE line handles every m at once:

    E[P^m] = m! * [s^m] (1+s)^m * Phi(-s g(s))

So everything is controlled by the single series f(s) := Phi(-s g(s)).  Now Lagrange-Burmann with phi(s) = 1+s (so sigma = t/(1-t), phi' = 1):

    sum_m t^m [s^m] f(s)(1+s)^m  =  f(t/(1-t)) / (1-t).

Demanding [s^m](1+s)^m f(s) = 0 for EVERY m >= 1 says that series is the constant f(0), i.e. f(t/(1-t)) = f(0)(1-t); substituting u = t/(1-t) gives

    f(s) = f(0)/(1+s)      -- UNIQUE, nothing else works.

Verified by perturbation: bumping ANY single Taylor coefficient of f destroys the vanishing immediately (tested s^0, s^1, s^2, s^3).
AND THEN THE m! IS AUTOMATIC:  E[QP^m] = m! [s^m] s(1+s)^{m-1} = m! [s^{m-1}](1+s)^{m-1} = m!.  The m! is a THEOREM, not a coincidence of either example.

(C) THE MASTER EQUATION   Phi(-s g(s)) = 1/(1+s).
For U = chi^2_d / 2 (i.e. U built from d real Gaussians), Phi(x) = (1-x)^{-d/2}, so
    (1 + s g(s))^{-d/2} = (1+s)^{-1}   =>   1 + s g(s) = (1+s)^{2/d}   =>   g(s) = ((1+s)^{2/d} - 1)/s.
    d=1 (n=3): g = 2 + s   -- EXACTLY the owner's (2+Z), DERIVED rather than guessed.
    d=2 (n=4): g = 1       -- EXACTLY the S133 four-term example P = (1+W)(Wbar - |Z|^2).
    d=3: g = 2/3 - s/9 + 4s^2/81 - ...   NOT polynomial.
    d=4: g = 1/2 - s/8 + s^2/16 - ...    NOT polynomial.
THE TWO KNOWN COUNTEREXAMPLES ARE ONE CONSTRUCTION AT d = 1 AND d = 2.

(D) MINIMALITY.
g is a POLYNOMIAL iff 2/d is a positive integer, i.e. d in {1,2}. Total real dimension = 2 (for Z,W) + d, so n in {3,4} and n = 3 IS MINIMAL for the family, attained uniquely at d = 1.
n = 2 would need d = 0, i.e. U constant = c. Then Phi(x) = e^{cx} and the master equation reads c*s*g(s) = log(1+s). But log(1+s)/s = 1 - s/2 + s^2/3 - ... is NOT a polynomial, so no such g exists. THIS FAMILY CANNOT REACH GMC(2).
Same obstruction in the Gamma family: U ~ Gamma(alpha) forces 1+sg = (1+s)^{1/alpha}, polynomial iff 1/alpha is a positive integer, and only alpha = 1/2 (d=1) and alpha = 1 (d=2) are realizable as sums of squares of Gaussians.

(E) TWO CORRECTIONS, ONE OF THEM AGAINST MYSELF.
1. The dimension table becomes GMC(N) FALSE FOR EVERY N >= 3 (pad with unused coordinates), superseding THM-1480's N >= 4.
2. MY OWN HYP-8330 IS REFUTED. I speculated at the end of S133 that the alternating-binomial collapse might be the ONLY mechanism producing E[P^m]=0 for all m, and that its dependence on k! 'may force an EVEN number of real variables, making n=4 sharp'. BOTH HALVES ARE WRONG: the d=1 case runs on (1+x)^{-1/2} composed with the perfect-square substitution x = 2s+s^2, not on sum(-1)^j C(m,j), and it lands on ODD n = 3. The correct general statement is the uniqueness of f = 1/(1+s); the alternating sum was only its d = 2 shadow. Logged as refuted rather than quietly dropped.

HANDOFF -- three, and (i) is now the whole game:
(i) HYP-8340 -- GMC(2) IS THE ONLY REMAINING CASE. GMC(1) is TRUE (DEZ Prop 4.2); N >= 3 is FALSE (THM-1500); DEZ prove GMC(2) only for HOMOGENEOUS P. THM-1500's obstruction is to ONE SHAPE, and three assumptions are load-bearing:
   (a) U INDEPENDENT of (Z,W). With only 2 real Gaussians there is NO ROOM for an independent U, so any GMC(2) counterexample must use a U CORRELATED with (Z,W) -- the obvious candidate being U = ZW, which is Exp(1) with E[U^k] = k!. The expectation then stops factorizing and the one-line master identity fails, so the whole analysis must be redone. THIS IS THE FIRST THING TO TRY.
   (b) the outer factor (1+Z). Replacing it by a general h(Z) changes phi(s) in the Lagrange inversion from 1+s to h(s), giving a DIFFERENT forced f. Worth redoing in general -- it may open dimensions the (1+Z) shape closes.
   (c) Q linear. Higher Q shifts the coefficient extraction (cf. THM-1480's r-family).
   A PROOF of GMC(2) would CLOSE THE CLASSIFICATION COMPLETELY.
(ii) HYP-8345 -- chi-squared is a CHOICE, not a classification. Work BACKWARDS: pick a target polynomial g, compute the Phi the master equation demands, and ask whether that MGF is realizable as a polynomial in Gaussians. Constructive rather than blind, cheap to run, and a hit with small support could beat the current 5-term quartic.
(iii) ATTRIBUTION: neither THM-1480 nor THM-1500 claims priority on the GMC(3) example -- it came from the owner via an outside source. What is ours is the independent verification, the master theorem, the uniqueness proof, and the minimality analysis.

Artifacts: THM-1500; 04-computation/gmc3_master_theorem_macmini_S134.py (+out); HYP-8330 marked REFUTED; HYP-8340, HYP-8345.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
