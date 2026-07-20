        # Message: mac-mini-2026-07-20-S139: THM-1580 -- the ARBORESCENCE COUNT STRICTLY DOMINATES H as a fingerprint (298 distinct vs 77 at n=7) AND is polynomial-time while H is not; my own S138 prime-2 framing REFUTED (it is crossing multiplicity); HYP-8315 extended to n=8 with its wording corrected

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 12:49

        ---

        OWNER: 'work the arborescence leads.'  All four, with one refutation of my own.

(A) HYP-8320 ANSWERED -- Sa AND H ARE INCOMPARABLE, AND Sa IS FAR FINER.
Over all iso classes:
    distinct H     =  3,  7, 19,  77   at n = 4,5,6,7
    distinct Sa    =  4, 11, 43, 298
    distinct PAIR  =  4, 12, 51, 417   (out of 456 classes at n = 7)
For n >= 5 NEITHER refines the other. This is the shape of THM-506's (char, perm) result BUT WITH THE COMPUTATIONAL ADVANTAGE RUNNING THE RIGHT WAY: sum_r a_r is a determinant (Tutte/Kirchhoff, polynomial time) while H is not. THE CHEAPER INVARIANT IS ALSO THE STRONGER ONE, and the pair is stronger still. That is the practically useful sentence in the whole thread.

(B) HYP-8390b REFUTED -- AND IT WAS MY OWN FRAMING FROM LAST SESSION.
I guessed the prime-2 split explains why log(Sa) gains a size-dependent shift while log H does not. WRONG, AND THE LAW ITSELF SAYS SO: Sa(T1 (+) T2) = Sa(T1) * det(|T1| I + L_in(T2)) holds for BOTH parities of |T1| (re-verified on 15 factor pairs) and never mentions 2.
THE REAL CAUSE IS CROSSING MULTIPLICITY. In T1 (+) T2 nothing in T2 beats anything in T1, so every Hamiltonian path runs through ALL of T1 and then ALL of T2 -- it CROSSES THE CUT EXACTLY ONCE, with NO choice about where, hence H is cleanly multiplicative with no interaction term. An arborescence rooted in T1 instead lets EACH of the |T2| vertices INDEPENDENTLY choose a parent, from the |T1| vertices across the cut or from inside T2 -- |T2| independent crossings, each with |T1| options, and that is exactly where p = |T1| enters the shift. EVENNESS OF Sa IS A CONSEQUENCE (the option count p can be even), NOT A CAUSE.

(C) HYP-8315 EXTENDED TO n=8, AND ITS WORDING CORRECTED.
    n     transitive=(n-1)!   min found   min=transitive?   max found   maximiser scores
    4                     6           6              yes          10   [1,1,2,2]
    5                    24          24              yes          55   [2,2,2,2,2]
    6                   120         120              yes         333   [2,2,2,3,3,3]
    7                   720         720              yes        2744   [3,3,3,3,3,3,3]
    8                  5040        5040              yes       23363   [3,3,3,3,4,4,4,4]
MINIMUM = the transitive tournament at exactly (n-1)!, now on SIX data points (exhaustive n <= 6, hill-climbed n = 7,8).
MAXIMUM: THE OLD PHRASING 'regular/Paley maximises' DOES NOT PARSE AT EVEN n, where no regular tournament exists. Corrected statement: the maximiser is REGULAR at odd n and NEAR-REGULAR at even n -- as balanced as the parity allows. (n=7's 2744 is the Paley value, matching THM-1460's closed form.) THE n=7,8 MAXIMA ARE HILL-CLIMBING LOWER BOUNDS, NOT CERTIFIED.

(D) HYP-8390c ANSWERED, and the answer is that the question was not deep.
v_2 is additive along the ordinal-sum factorisation, and for the transitive tower the shift is exactly p = n-1, so v_2(Sa(TT_n)) = v_2((n-1)!) = (n-1) - s_2(n-1) -- Legendre. Verified n = 3..11. No 2-adic structure beyond Legendre acting on the size factors.

HANDOFF -- HYP-8410 is new and cheap, and I would take it first:
The pair (H, Sa) separates 417 of 456 classes at n = 7, leaving 39 UNSEPARATED.
 (a) WHAT ARE THE 39? Extract the colliding groups and look for structure. Are they switching-equivalent? THM-1440's skew-Seidel spectrum IS switching-invariant while H and Sa are NOT, so it should cut differently -- and A PATTERN AMONG THE 39 WOULD BE WORTH MORE THAN THE FINGERPRINT ITSELF.
 (b) Does one more cheap invariant finish n = 7 -- the skew-Seidel spectrum, or the permanent (THM-506)? If a small polynomial-time package is complete at n = 7, that is a genuinely useful tournament fingerprint and directly comparable to (char, perm).
Also still open:
 * HYP-8315's MINIMUM. Sa = prod of the nonzero eigenvalues of L_in, whose trace is always C(n,2); for the transitive tournament L_in is triangular with in-degrees 0,1,...,n-1. Majorisation on the in-degree sequence is the obvious attempt. NOTE AM-GM IS NOT AVAILABLE -- L_in is non-symmetric with complex eigenvalues (Paley's are (q +- i sqrt q)/2), which is exactly why the naive bound fails.
 * HYP-8390a. Is the prime disjointness (2 purely additive, odd primes purely multiplicative) forced or an accident of H? Sa is NOT a counterexample -- it is simply not on the multiplicative side, being neither multiplicative nor all-odd.

CAVEATS I want carried: the n=7,8 maxima are not certified; (A)'s incomparability is about VALUES, not usefulness (H carries the Redei/OCF structure Sa does not, and nothing here suggests replacing it -- the claim is narrowly that as a SEPARATOR Sa is finer and cheaper); and (B)'s crossing account is an EXPLANATION of a proved law, not a new theorem -- what is established there is the negative.

Artifacts: THM-1580; 04-computation/arborescence_leads_macmini_S139.py (+out); HYP-8390 marked partially refuted / partially answered; HYP-8315 updated; HYP-8410 new.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
