        # Message: collatz-procgen-20260922 wave 23: THM-4524 selfie tournaments (OPEN-Q-060 answered; all-odd arc parity first at N=6) + THM-4525 edim_m(Q6)=15 (arXiv:2608.09983 Open Problem 1 settled)

        **From:** mac-mini-2026-10-01-S?
        **To:** all
        **Sent:** 2026-10-01 13:24

        ---

        collatz-procgen-20260922, wave 23 (2026-10-01): the owner's selfie-tournament prompt, and the open problems of arXiv:2608.09983.

1. THM-4524 — selfie tournaments (PROVED + INDEPENDENTLY AUDITED).
   - Loops versus the path. The N optional loops are the switching gauge bits: (tiling, loop set) -> tournament is exactly 2-to-1, with L and its complement giving the same tournament. The N-1 base-path arcs are the loops' discrete derivative: arc i+1 -> i is reversed iff exactly one end is looped. So C(N-1,2) + N = C(N,2) + 1, and the extra bit is the global complement.
   - Loop-Walsh degree. In the loop spins, H has Walsh degree <= 2 floor(N/3), so sum_L (-1)^|L| H(switch_L T) = 0 for every T. For N <= 5, H restricted to a switching class is an Ising model.
   - Fixed-point OCF and selfie Redei. The selfie Redei count is odd and is H + 2|L| (mod 4).
   - Owner's 5/6 conjecture: REFUTED in both readings.
     * Every arc lies on some HP for a fraction of tournaments rising to 1 (0.998 at N = 10).
     * The true threshold is a parity one and is reversed. "Every arc on an odd number of HPs" is impossible for N = 3..5, and N = 6 has exactly one such class (QR_7 minus a vertex). N = 9 has none, N = 10 has exactly two, and N = 0, 3 mod 4 is excluded by parity.
     * Cayley tournaments of odd abelian groups have every arc count even.
   - Shaving.
     * c(e) even iff e can be shaved keeping H odd.
     * All-even tournaments need exactly 2 deletions to make H even.
     * Every non-all-odd class with N <= 7 shaves down to one HP with H odd at every step.
   - OPEN-Q-060: ANSWERED in Burnside form. A049313 = the number of unlabeled Euler graphs whose automorphisms reverse an even number of edges, i.e. Mallows-Sloane twisted by the tournament torsor's cocycle (Brauer's permutation lemma). Verified through n = 10. The bijective form is still open.
   - New hypotheses:
     * HYP-9167: Paley minus a vertex is Redei-rigid, and all-odd tournaments exist iff N = 2 mod 4. The antipodal orbit is PROVED; FINITE-EXACT for q = 7..27.
     * HYP-9168: the HP-blocking number equals the Hall bound (FINITE-EXACT, N <= 9).
   - Audit.
     * The orchestrator's own C engine and Python: 58 checks, plus E1 for every n <= 10. Covered: all labeled tournaments N <= 7, all classes N <= 9, QR_q minus a vertex for q <= 23, all circulants up to 13, both all-odd classes at N = 10.
     * The lane's runner was re-run: 133060 checks, identical up to timing.

2. THM-4525 — edim_m(Q_6) = 15 (FINITE-EXACT + INDEPENDENTLY AUDITED). This settles Open Problem 1 of Allikvere's paper, whose bounds were 6..15.
   - Q_6 is the n = 5 tiling cube. No edge-multiset resolving set of size <= 14 exists.
   - Three exhaustive searches up to Aut(Q_6) (order 46080), with different normal forms, all find nothing at k <= 14: the lane's max-imbalance and min-imbalance forms, and the orchestrator's plain |A| >= |B| layer split. The orchestrator's search enumerated 14,168,149,784 leaves for k <= 14 with none resolving; at k = 15 its 1678 resolving leaves reduce to exactly the same 229 orbits (882 s).
   - The resolving 15-sets form exactly 229 orbits, all with trivial stabilizer (10,552,320 sets).
   - Multisets do not refine (two histograms can merge when a landmark is added), so no class-splitting pruning was used.
   - New lemmas:
     * resolving sets are asymmetric;
     * antipodal reversal;
     * counting bound >= 7;
     * an entropy bound edim_m(Q_d) >= exp((0.6215 - o(1)) d^(1/3)): superpolynomial growth, the lower half of Open Problem 2.
   - Explicit sets: Q_7 <= 19 (the paper had 63), Q_8 <= 26, Q_9 <= 38, Q_10 <= 48, Q_11 <= 65, Q_12 <= 76. So density 1/2 is far from optimal (Open Problem 4).
   - HYP-9169: ln edim_m(Q_d) = Theta(d^(1/3)). Sparse union bounds give ln M_d / d^(1/3) ~ 2.77 for d = 11..32 (double precision, not interval-certified).

Next obligations:
- HYP-9167: the non-antipodal arc orbits of QR_q - 0; N = 13 and N = 14 are the first open cases.
- HYP-9168.
- The bijective form of OPEN-Q-060.
- edim_m(Q_7), which lies in [8, 19].
- An interval-certified sparse union bound for HYP-9169.
- Still flagged, not started: the Allikvere LRC(14)/(15) claim (arXiv:2609.02604); downloading its certificates needs the owner's permission.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
