        # Message: mac-mini-2026-07-20-S137: TWO CHARTS ON ONE TRICHOTOMY -- {-1,0,1} and {1,1/2,0} are the same object under x to 2x-1, the tie is the shared fixed point, and four of this session's results collapse into one statement (skew-Seidel = 2p-1; fiber fraction = a tie probability; THM-1500's ladder IS the trichotomy; the nullcone is its pure corners)

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 12:22

        ---

        OWNER: 'look back at the summand and multiplicand graph, and ... think about how {-1,0,1} and {1,1/2,0} are functionally equivalent but have shown up in this repo many times each. this is what i mean by even/odd vs positive/negative.'

They are functionally equivalent because x -> 2x-1 carries one onto the other. What makes it worth a file is that this repo has been computing in BOTH CHARTS for months without ever naming the change of variables -- and once named, four separate results from this session turn out to be the same statement.

THE DICTIONARY.
    {0, 1/2, 1}   ADDITIVE / probability chart.  Averaging natural.  Involution p -> 1-p.   = the owner's POSITIVE/NEGATIVE
    {-1, 0, +1}   MULTIPLICATIVE / sign chart.   (-1)^k natural.     Involution S -> -S.    = the owner's EVEN/ODD
x -> 2x-1 and x -> (1+x)/2 are inverse. This is the Fourier-Walsh <-> probability dictionary (E[X] = 2p-1), and it is the repo's own 'two arithmetics' one level DOWN -- on the value set rather than on the algebra.
THE WHOLE CONTENT IS THE FIXED POINT. Both charts have exactly one and it is the same point: p = 1/2 <-> S = 0. That is the TIE -- the diagonal, the draw, the balance.

FOUR INSTANCES, ALL VERIFIED:

(1) A TOURNAMENT IS A PROBABILITY MATRIX, AND ITS SKEW-SEIDEL MATRIX IS 2p-1.
Read p_ij in {0,1/2,1} -- 1 if i beats j, 0 if it loses, 1/2 on the diagonal. Then S = 2p-1 is EXACTLY the skew-Seidel matrix of THM-1440: lands in {-1,0,1}, skew, zero diagonal. Verified on 200 random tournaments. And COMPLEMENTATION p -> 1-p IS S -> -S, with the tie as the unique fixed point.
PALEY MAKES IT LITERAL: the skew-Seidel matrix of the Paley tournament IS THE LEGENDRE SYMBOL chi(j-i), whose value set is exactly {-1,0,+1}; its probability avatar (1+chi)/2 is the quadratic-residue indicator with 1/2 sitting at 0. Checked q = 7, 11, 19.

(2) THE FIBER FRACTION IS A TIE PROBABILITY.
CLAUDE.md records f(n) = (1/2)_{n-2}/(n-2)! with GF (1-x)^{-1/2}. That equals C(2k,k)/4^k -- THE PROBABILITY OF AN EXACT TIE IN 2k FAIR COIN FLIPS -- and also E[U^k]/k! for U = T^2/2. All three agree exactly for k = 0..8. So the repo's 1/2 is a FAIR COIN, its 'two-sheeted branched cover' is that coin's square root, and the tie it computes is the 0 of the sign chart.

(3) THM-1500's ADMISSIBLE LADDER IS THE TRICHOTOMY.
Phi(x) = (1-x)^{-d/2}, admissible iff 2/d is a positive integer:
    d=0  exponent 0    -> -1  n=2  g = log(1+s)/s   OBSTRUCTED (GMC(2) open)
    d=1  exponent 1/2  ->  0  n=3  g = 2+s          MINIMAL counterexample
    d=2  exponent 1    -> +1  n=4  g = 1            the four-term example
The exponent ladder {0, 1/2, 1} maps under x -> 2x-1 onto {-1, 0, +1}, and THE MINIMAL COUNTEREXAMPLE SITS AT THE FIXED POINT. The same 1/2 that is the fiber fraction's fair coin is the 1/2 that makes (1+s)^{2/d} a PERFECT SQUARE. Square root and square are one one-half.

(4) THE TELESCOPING PRINCIPLE IS A STATEMENT ABOUT THE SIGN CHART.
THM-1520/1540's nullcone = strictly one-sided charge support = THE TWO PURE CORNERS of the charge-sign trichotomy {-1,0,+1}; anything touching 0, or touching both signs, is outside. So 'forbidding one variable telescopes' (S135) and 'the nullcone is one-sided' (S136) are both the statement that THE NULLCONE IS THE BOUNDARY OF THE TRICHOTOMY AND THE HARD CASE IS ITS INTERIOR -- which is instance (1)'s fixed point, one level up.

CENSUS -- both charts are everywhere, never named as the same thing:
    sign chart:        Legendre 84 files, sgn 20, skew-Seidel 7, odd-function 6
    probability chart: merged metagraph 69, fiber fraction 29, (1/2)_k 15, (1-x)^{-1/2} 9

THE HONEST LIMIT, and I want this attached to any use of it. THIS IS A CHANGE OF COORDINATES, NOT A THEOREM, and it proves nothing on its own. Two specific cautions recorded in the reflection: (a) the instance-(4) trichotomy is the SIGN of a Z-grading, not a three-valued quantity -- the 'balance' statistic used to illustrate it is crude and mishandles charge-0 monomials, so quote THM-1540 rather than the reflection; (b) the recurrence of 1/2 is PARTLY structural -- instances (2) and (3) share a genuine generating function (1-x)^{-1/2} -- and partly just that 1/2 is the commonest number in mathematics. I checked the first three share an actual identity and I make NO CLAIM that every 1/2 in the repo is this one.

HANDOFF -- one real question, HYP-8385.
The minimal GMC counterexample sits at the FIXED POINT of the parameter involution. IS THAT FORCED OR COINCIDENTAL? If the minimal member of a construction always lands on the self-dual point of its parameter involution, that is a principle worth having -- it would PREDICT where the minimal case sits in other families instead of requiring a search. Cheap test: take the r-family of THM-1480, the chi^2_d ladder itself, and the Gamma-alpha family, and check whether the extremal/minimal member is the fixed point of the natural involution on the parameter. ONE counterexample kills it; THREE hits would make it worth stating. Any claim needs an actual identity, not a numerical coincidence.

The usable form of all this is a HABIT rather than a result: when a repo object takes three values, or an involution has a fixed point, ask which chart you are in and write the other one down. Sign coordinates make parity, characters and determinants easy; probability coordinates make averaging, measure and generating functions easy. The translation is free.

Artifacts: 07-reflections/two-charts-on-one-trichotomy-macmini-S137.md; 04-computation/two_coordinate_systems_macmini_S137.py (+out); HYP-8385.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
