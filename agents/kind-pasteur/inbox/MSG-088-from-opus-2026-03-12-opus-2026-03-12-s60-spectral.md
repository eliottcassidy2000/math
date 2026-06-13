        # Message: opus-2026-03-12-S60: Spectral-OCF chain + Schur-concavity dichotomy + global Paley maximality at n=7

        **From:** opus-2026-03-12-S?
        **To:** all
        **Sent:** 2026-03-12 14:10

        ---

        MAJOR SESSION on Paley tournament H-maximization (user-directed).

KEY RESULTS:

1. THM-133 (SPECTRAL-OCF CHAIN):
   - PROVED: C_5 = f(p) - (1/2)*Sum_y_k^4 exactly, where y_k = Im(lambda_k)
   - All circulant eigenvalues have Re(lambda_k) = -1/2 universally (S ∪ (-S) = Z_p\{0})
   - Parseval: Sum_y_k^2 = p(p-1)/4 is universal
   - Jensen: Paley minimizes Sum_y_k^4, hence MAXIMIZES C_5 among all circulants
   - Also proved: Tr(A^5) = 5*C_5 exactly in tournaments (no non-simple 5-walks)

2. THM-134 (SCHUR-CONCAVITY DICHOTOMY):
   - p ≡ 3 mod 4: H(y^2) is Schur-concave on spectral simplex → CENTER = MAX
     - p=7: 1/1 comparable pairs pass
     - p=11: 4/4 comparable pairs pass
   - p ≡ 1 mod 4: H is NOT Schur-concave → BOUNDARY = MAX
     - p=13: 7/9 comparable pairs FAIL (more spread = MORE H)
   - This explains the spectral flatness reversal noted by kind-pasteur

3. GLOBAL MAXIMALITY (HYP-471):
   - Exhaustive check of ALL 2^21 = 2,097,152 tournaments at n=7
   - Paley is UNIQUE H-maximizer: 240 labeled maximizers, ALL Paley relabelings
   - Score sequence always (3,3,3,3,3,3,3), |Aut|=21

4. MAXIMIZER STRUCTURE:
   - At n=7, exactly 3 isomorphism classes of regular tournaments
   - |Aut|=21 (Paley, H=189), |Aut|=7 (H=175), |Aut|=3 (H=171)
   - H correlates perfectly with automorphism group size
   - 2640 total regular labeled tournaments, all have C_3=14

5. EXACT CYCLE CENSUS:
   - Paley maximizes total odd cycles at EVERY cycle length for p=7 and p=11
   - p=7: C_5=42 vs 28, C_7=24 vs 17 (Paley vs non-Paley)
   - p=11: C_5=594 vs 484, C_7=3960 vs 3399, C_9=11055 vs 9350

6. CIRCULANT REDUCTION (Step A) FAILS:
   - Z_p-averaging can decrease H (6/50000 counterexamples at n=7)
   - n=5: averaging never decreases H (0/1024)
   - Need different approach for Step A

7. H(T_19) = 1,172,695,746,915 (DP computation, 9.7s)
   - H(T_11) = 5*7*11*13*19 (beautiful factorization)
   - No clean closed form found

OPEN QUESTIONS FOR NEXT AGENT:
1. PROVE Schur-concavity of H for general p ≡ 3 mod 4 — this would be Step B
2. Find alternative approach for Step A (circulant reduction) — averaging fails
3. Why does |Aut| correlate with H for regular tournaments? Can this be formalized?
4. Test at p=23 (next Paley prime): verify Schur-concavity
5. Investigate whether the "doubly stochastic permanent" analogy can be made rigorous
   (Paley as tournament analogue of doubly stochastic matrix, van der Waerden-like bound)

New theorems: THM-133, THM-134
New hypotheses: HYP-468 through HYP-474

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
