        # Message: mac-mini-2026-07-20-S135: THM-1520 -- GMC(2)'s ONE-SIDED-CHARGE BRANCH IS CLOSED (charge telescoping collapses 2 vars to 1; the saddle limit exp(a_{D-1}/(D a_D)) is never zero), two-sidedness is NECESSARY for failure; plus the TELESCOPING PRINCIPLE with 4 other repo instances

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 11:12

        ---

        OWNER: 'work on the GMC(2) proof, think forbidding one variable telescopes, see how that pattern can abstractly apply to other problems the repo has touched.' The hint was exactly right and it closed a branch.

THE GRADING. Two real Gaussians = one standard complex Z with W = Zbar, E[Z^a W^b] = a! delta_ab. Grade monomials by CHARGE c = deg_Z - deg_W. Then E ANNIHILATES EVERY NONZERO CHARGE. Cleanly: in polar form Z = sqrt(V) e^{i theta} with V = ZW ~ Exp(1), V and theta are INDEPENDENT and charge is exactly the theta-Fourier index; E is the charge-0 projection followed by Z^a W^a -> a!.

(A) THE TELESCOPING LEMMA -- the owner's hint made precise.
Charges ADD under multiplication. So if P's charge support is ONE-SIDED (say all charges >= 0), the only combination summing to 0 takes the charge-0 part from EVERY factor:
    E[P^m] = L(p_0^m),   L(v^k) := k!,   p_0 = charge-0 part of P.
TWO VARIABLES COLLAPSE TO ONE. Verified exactly on random one-sided P; and it correctly FAILS for two-sided P (P = Z + W - ZW gives E[P^m] = -1, 4, -18, 108 against L(p_0^m) = -1, 2, -6, 24), which is the point, since two-sided is exactly where the n >= 3 counterexamples live.

(B) THE SADDLE LEMMA.
L(f) = int_0^infty f(v) e^{-v} dv. For p = a_D v^D + a_{D-1} v^{D-1} + ... of degree D >= 1,
    p(v)^m = a_D^m v^{Dm} (1 + (a_{D-1}/a_D)/v + ...)^m ~ a_D^m v^{Dm} e^{m a_{D-1}/(a_D v)},
and the saddle of v^{Dm} e^{-v} sits at v = Dm, where that exponent is a_{D-1}/(D a_D). Hence
    L(p^m) / (a_D^m (Dm)!)  ->  exp( a_{D-1} / (D a_D) )   -- A LIMIT THAT IS NEVER ZERO.
Confirmed numerically to 6 digits: 1/e = 0.36787944 for v-1 (measured 0.367879 at m=14), e^3 = 20.0855369 for v+3 (measured 20.085523), e^{-3/2} = 0.2231302 for v^2-3v+2, and exactly 1 for v and for 2v^3+v. So L(p^m) != 0 for all large m whenever deg p >= 1, and therefore L(p^m) = 0 for ALL m forces p = 0. Independently, symbolic elimination over m = 1..D+1 finds NO nonzero p of degree 1, 2 or 3.

(C) COROLLARY -- GMC(2) HOLDS ON THE ONE-SIDED BRANCH.
Let P have one-sided charge support with E[P^m] = 0 for all m. By (A), L(p_0^m) = 0 for all m; by (B), p_0 = 0. So every charge of P is >= 1, hence every charge of P^m is >= m, hence QP^m has charge 0 only if charge(Q) <= -m, i.e. only if m <= deg_W(Q). Therefore E[QP^m] = 0 for every m > deg_W(Q). QED

(D) TWO-SIDEDNESS IS NECESSARY FOR FAILURE.
The contrapositive of (C): ANY GMC(2) counterexample must have charges of BOTH SIGNS. A real structural constraint on the last open case, and consistent with the n >= 3 counterexamples, whose P = (1+Z)(W - g(Z)U) runs from charge -1 upward.

(E) THE ABSTRACT PATTERN -- new reflection .
THE TELESCOPING PRINCIPLE: if an invariant is computed by a GRADED functional (one annihilating every nonzero grade) and the object's support in that grading is ONE-SIDED, the invariant COLLAPSES to a function of the grade-0 part alone, in strictly fewer variables -- and failure REQUIRES two-sided support. Four other repo objects have exactly this shape:
  * THM-1460, ordinal sums. T1 (+) T2 forbids arcs one way, making L_in BLOCK TRIANGULAR, so spec = spec(L1) u (|T1| + spec(L2)) and the arborescence count telescopes into the factors. Grading = which block.
  * THM-1440, skew-Seidel parity. p(-x) = (-1)^n p(x) is a Z/2 grading; at odd n the grade-0 coefficient must vanish, forcing the zero eigenvalue. Same fact as 'sin vanishes at 0 and cos does not'.
  * THM-506 / HYP-2514, the master cycle-packing polynomial. FACES ARE GRADED PIECES: the spectrum is the signed face, H the unsigned odd-only face, and which one collapses to a determinant is exactly why the spectrum is computable and H is not.
  * THM-1405 / THM-1420, cut + cycle. Star flips are the grade-0 gauge directions; holonomies are what survives. The bicycle space vanishing at odd n is one-sidedness over F_2.
Moral for future sessions: LOOK FOR THE GRADING FIRST. Every time this repo has recorded something as vanishing 'for parity reasons' or 'because the charges do not balance', this is what was happening.

HANDOFF -- two, and the first is a single clean gap:
(i) HYP-8350 -- THE ONE GAP. THM-1520 (C) is CONDITIONAL on the saddle lemma, which is a Laplace-method derivation WITHOUT explicit error bounds. To make it rigorous: substitute v = Dm(1+u), expand log(p(v)^m e^{-v}) about u = 0, and bound the tails |u| > m^{-1/2+eps} against the Gaussian core. The delicate point is the O(1) phase correction m*arg(p(v)/v^D) -> a_{D-1}/(D a_D), delicate precisely because p has COMPLEX coefficients. Standard but not free. Payoff: THM-1520 (C) becomes unconditional and GMC(2) reduces EXACTLY to its two-sided branch with nothing else outstanding.
(ii) HYP-8355 -- the two-sided branch is now ALL that remains of GMC(2), hence of the entire GMC classification (GMC(1) true; N >= 3 false by THM-1500; one-sided N = 2 closed here). A WARNING I want on the record, learned the hard way in S133: a computational sweep of this branch has NO POSITIVE CONTROL AVAILABLE, because a positive control at n = 2 would BE the counterexample being sought. My ~950k-polynomial sweep found nothing and must NOT be cited as evidence that GMC(2) is true -- read it as 'the obvious small cases are clear' and nothing more. Two structural leads instead: (a) at n = 2 there is no room for a U independent of (Z,W), so the natural candidate is U = ZW ~ Exp(1) with E[U^k] = k!, which is CORRELATED with (Z,W) -- the expectation stops factorizing and THM-1500's one-line master identity fails, so find what replaces Phi(-s g(s)) = 1/(1+s); (b) replacing the outer factor (1+Z) by a general h(Z) changes phi(s) in the Lagrange inversion from 1+s to h(s), giving a DIFFERENT forced f. Unexplored, and it may open exactly the dimension the (1+Z) shape closes.

Artifacts: THM-1520; 07-reflections/the-telescoping-principle-macmini-S135.md; 04-computation/gmc2_charge_telescoping_macmini_S135.py, gmc2_limit_and_twosided_macmini_S135.py (+outs); HYP-8350, HYP-8355.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
