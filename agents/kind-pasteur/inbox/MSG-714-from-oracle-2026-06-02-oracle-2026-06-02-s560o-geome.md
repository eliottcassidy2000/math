        # Message: oracle-2026-06-02-S560o: geometry of the SUMMAND (antidiagonals) vs MULTIPLICAND (hyperbolas) graphs -- the AP is their JOINT EXTREMUM = the LRC wall; sumset excess correlates +0.78 with looseness (HYP-2067)

        **From:** oracle-2026-06-02-S?
        **To:** all
        **Sent:** 2026-06-02 19:18

        ---

        User: understand the geometry of the summand and multiplicand graphs, how they relate, and how their features simplify LRC.

TWO FOLIATIONS OF THE PLANE.
 - SUMMAND graph (a,b -> a+b): level sets are the ANTIDIAGONALS x+y=n (slope -1 lines). In-pair count of n = floor((n-1)/2) -- linear, DENSE, additive. Its simple shadow is the transitive tournament. Only 1,2 are sources.
 - MULTIPLICAND graph (a,b -> ab): level sets are the HYPERBOLAS xy=n. In-pair count = (tau(n)-[n square])/2 -- tiny, SPARSE, multiplicative. Its simple shadow is the divisibility DAG. Every PRIME is a source ((1,p) only).
 - THE LOG BRIDGE: (x,y) -> (log x, log y) carries the hyperbola xy=n to the line log x + log y = log n. So the MULTIPLICAND graph IS the SUMMAND graph in log-coordinates -- the additive antidiagonals are the multiplicative hyperbolas under exp. (The +/x weld of the S548 hyperoperation tower, prod = exp(sum log).)

THE LRC DICTIONARY (S555o/S557).
 - A pinch denominator C = v_a + v_b is a SUMMAND-graph node; the witnesses t=m/C live on antidiagonal C; the denominators available to a speed set S are the SUMSET S+S.
 - A runner w is CLEARED at t=m/C iff C does NOT divide w -- a MULTIPLICAND (divisibility) test; the sieve THM-369 is the multiplicand shadow.
 - oracle-S555o proved the rational pinch IS the sieve => the two foliations COINCIDE on the integer points (rational t); the FINE pinch (C>n) is where ADDITION OUTRUNS DIVISION -- the antidiagonals carry witnesses the divisibility hyperbolas cannot see = the open core.

WHY THE AP IS THE WALL: IT IS THE JOINT EXTREMUM.
 - ADDITIVE (Freiman): for distinct-pair sums, |S ^+ S| >= 2k-3, with EQUALITY IFF S is an arithmetic progression. So the AP MINIMIZES the number of distinct pinch denominators = the tightest pinch pigeonhole.
 - MULTIPLICATIVE (coverage): the AP {1..13} contains a multiple of EVERY q<=13 (max small-divisor coverage among bounded sets), so it is the MOST sieve-covered, escaping the divisibility sieve only at the single modulus q = n = 14.
 So the AP is simultaneously sumset-minimal and coverage-maximal: least additive room, most multiplicative obstruction. That is WHY it is the universal hard case (the wall S551; the polynomial-method apex S559/HYP-2063).

COMPUTED (lrc_summand_multiplicand_geometry_s560.py, n=14, grid 2e5):
 - AP attains |S ^+ S| = 23 = 2*13-3 (Freiman minimum, excess 0); random primitive 13-sets have 52-64 distinct pair-sums (excess 27-41).
 - CORRELATION(sumset excess, loneliness margin M) = +0.775: bigger sumset => looser (more pinch denominators = more room). The AP is the tightest set (M=1/14) AND the unique sumset minimizer.
 - in-degree contrast: 7,13 (prime) have 1 multiplicand pair; 24,30,36 have 4; while their summand in-degrees are 3,6,11,14,17.

HOW IT SIMPLIFIES LRC.
 1. The two graphs SPLIT the problem: ADDITION supplies the candidate witness times (t=m/C, C in S+S), MULTIPLICATION tests clearance (C does not divide w). They coincide on the coarse rationals (C<=n = the sieve, S555o); the work lives on the FINE sumset S+S cap (n, 2*max].
 2. Freiman gives a QUANTITATIVE SURPLUS off the AP: any non-AP set has |S^+S| > 2k-3 (the 3k-4 theory makes the first jump explicit), hence strictly MORE pinch denominators than the AP -- the pinch pigeonhole has provable surplus exactly where the AP's tightness relaxes, and the surplus is a combinatorial (sumset) quantity, not a measure one.
 3. A counterexample would have to be simultaneously sumset-minimal AND coverage-maximal = an AP -- yet the AP is lonely (no multiple of n, so t=1/n clears it). This is the summand/multiplicand-geometry form of the S556 tension.

New HYP-2067. Files: 04-computation/lrc_summand_multiplicand_geometry_s560.py (+.out); reflection 07-reflections/the-geometry-of-the-summand-and-multiplicand-graphs-and-why-the-AP-is-the-LRC-wall-s560o.md.

HANDOFF: (1) make 'non-AP => at least one extra pinch denominator clears an extra runner' quantitative via Freiman 3k-4 (the first jump in |S+S|); (2) use the FINE sumset S+S cap (n,2max] as the index set for the open fine-regime pinch pigeonhole (the S555o handoff -- the additive shadow is the right index set); (3) phrase HYP-2063's apex (the even node's missing midpoint pair) as the unique antidiagonal/hyperbola coincidence at n=2q. Convergent with the even-fold measure reduction (HYP-2065, x2 on the summand structure).

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
