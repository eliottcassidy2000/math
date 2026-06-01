        # Message: kind-pasteur-2026-03-06-S18e: SC maximizer n=7, sigma* theorem, grid overlap, Paley deletion

        **From:** kind-pasteur-2026-03-06-S?
        **To:** all
        **Sent:** 2026-03-06 03:40

        ---

        SESSION SUMMARY:

1. SC MAXIMIZER CONFIRMED AT n=7 (exhaustive 2^21 = 2M tournaments):
All 15 self-complementary score classes have SC achieving max H. Only self-comp scores (s_i+s_{n-1-i}=n-1) can have SC tournaments. Key finding: score (4,4,3,3,3,3,1) has NO SC tournaments (not self-comp), so apparent 'counterexample' was vacuous.

2. ANTI-AUTOMORPHISM INVOLUTION THEOREM:
- sigma^2 is always an automorphism of T (proved)
- At least one involutory anti-aut always exists (verified exhaustive n=4,5,6, random n=7)
- At even n: involutory anti-aut is ALWAYS fixed-point-free (proved: fixed pt implies score=(n-1)/2, non-integer)
- At n=6: 10/12 SC classes have all anti-auts as involutions; 2 have order-6 anti-auts BUT still have involutory ones

3. sigma* THEOREM (corrected):
sigma induces sigma* on DIRECTED odd cycles (must reverse direction when mapping). sigma* is:
(a) Always an involution (proved)
(b) Always an automorphism of conflict graph Omega(T) (proved)
(c) Fixed-point-free at even n; has fixed points at odd n

4. DISJOINT PAIR MECHANISM:
SC tournament at n=6 score (3,3,3,2,2,2): 8 three-cycles selecting one from each sigma-orbit, forming 4 complementary vertex-disjoint pairs. alpha_2=4, I.P.=[1,14,4], H=45. NSC: only 1 disjoint pair, alpha_2=1, H=43. Two routes to max H: more cycles (alpha_1) OR more pairs (alpha_2).

5. GRID OVERLAP STRUCTURE (user insight):
n-grid = n overlapping (n-1)-grids. Sigma-orbit deletion at even n gives special inherited symmetry. Transversals = 3-cycles. Even/odd n+2 step preserves parity structure.

6. PALEY DELETION GIVES MAXIMIZER:
T_7 (H=189) -> delete any v -> T_6 (H=45) = global max at n=6! Conjecture: T_p deletion always gives maximizer at p-1.

NEXT STEPS:
- Prove SC maximizer algebraically using sigma* + orbit structure
- Test at n=8 (even n, computationally expensive)
- Verify Paley deletion conjecture at p=11 (need H computation for T_11-v)
- Connect grid overlap to inductive proof

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
