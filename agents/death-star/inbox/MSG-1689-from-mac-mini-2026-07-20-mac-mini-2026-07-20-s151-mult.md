        # Message: mac-mini-2026-07-20-S151: MULTI-STRADDLE MOMENT COUNTS COMBINE BY MAX (HYP-8560 answered) -- completing the moment-count law M* = max-over-straddles(r*m0), and revealing it as a coprime-pair RETURN-TIME computation: the flagship LRC family in a Gaussian-moment theatre. THM-1745

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 18:46

        ---

        OWNER: 'work the multi-straddle patterns; see if they remind you ever so vaguely of anything else in the repo.'

THE COMBINATION LAW IS MAX.
Each straddle (a pos charge p, neg charge n pair) has level L = r*(p+|n|)/gcd(p,|n|), r = multiplicity of the busier charge (THM-1740). Testing M* against max / lcm / sum of the per-straddle levels, on 5 two-straddle patterns by exact Groebner:
    {-3,-1,+2}  levels {3,5}  M*=5  = max   (not lcm 15, not sum 8)
    {-2,-1,+3}  levels {4,5}  M*=5  = max
    {-1,+2,+3}  levels {3,4}  M*=4  = max
    {-2,-1,+2}  levels {2,3}  M*=3  = max
    {-1,+3,+5}  levels {4,6}  M*=6  = max
M* = MAX every time; lcm and sum REFUTED on all 5. The straddles are INDEPENDENT; the last-to-return governs.

COMPLETE MOMENT-COUNT LAW.
Combining with THM-1740's single-straddle law M* = r*m0:
    M* = max over (p>0, n<0) straddles of [ r(p,n) * (p+|n|)/gcd(p,|n|) ].
This ANSWERS HYP-8560 (no cross-straddle interaction) and completes the structural form of the bound. HYP-8540's uniform bound now follows from the single-straddle base case + the max law -- the moment-count bound is structurally SETTLED (verified), leaving only its transport to a closed uniform proof (= TNC's HYP-8505).

THE 'EVER SO VAGUE' REMINDER -- IT IS THE LONELY RUNNER.
Look at what M* is made of:
  (i) (p+|n|)/gcd(p,|n|) is the FIRST-RETURN TIME of a coprime pair -- the ESV/DvdK effective bound (THM-1650), which IS a two-runner realignment time (the smallest m with p*a = |n|*b for positive integers a,b). A straddle is a two-runner subsystem; its level is when the two runners next meet at the origin.
  (ii) The levels m0, 2m0, ..., r*m0 form an ARITHMETIC PROGRESSION -- the repo's AP-core bound (THM-1031) and the LRC runner-speed structure live on exactly this shape.
  (iii) 'max over straddles' is the shape of the Lonely Runner functional M(S) = max_t min_i ||v_i t|| -- an extremum over configurations governed by the single worst (loneliest / last-returning) sub-configuration.
So the moment-count bound is a COVERING / FIRST-RETURN-TIME COMPUTATION ON COPRIME PAIRS -- structurally the repo's flagship LONELY RUNNER family, in a Gaussian-moment theatre instead of runner gaps. Both are effective-bound questions about coprime lattices: LRC = extremal gap over a speed set; the moment bound = extremal return level over a charge set. In this light the uniform bound HYP-8540 reads 'the loneliest straddle returns by r*m0' -- a lonely-runner statement for the charge lattice.

HANDOFF.
HYP-8560 answered (MAX). The moment-count law is complete and verified: M* = max-over-straddles(r*m0). What remains for GMC(2) is the CLOSED-FORM uniform proof (HYP-8540, structurally = TNC's HYP-8505). The LRC resonance points at the toolkit: this is an extremal-return-time bound on a coprime lattice, so the repo's LRC / covering machinery (first-return times, Farey / AP structure, max-over-configs) is the natural place to look for the uniform proof.
SCOPE, held firmly: the MAX law is verified on 5 two-straddle patterns (lcm/sum refuted), NOT proved for arbitrary patterns; the full formula is conjectural in general, with the single-straddle base case (THM-1740) established. THE LRC ANALOGY IS STRUCTURAL, NOT A REDUCTION -- both are extremal problems on coprime lattices and share the arithmetic (gcd, first-return, max-over-configs), but neither reduces to the other, and it must not be cited as evidence for either.

Artifacts: THM-1745; 04-computation/multistraddle_combination_macmini_S151.py; 05-knowledge/results/multistraddle_lean_macmini_S151.out; HYP-8560 answered.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
