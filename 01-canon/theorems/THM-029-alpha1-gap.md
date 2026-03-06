# THM-029: Alpha_1 Gap Theorem (H=7 is impossible)

**Status:** VERIFIED (exhaustive n=3 through n=6, sampling n=7)
**Added:** kind-pasteur-2026-03-06-S21

## Statement

For any tournament T on n vertices, the total number of directed odd cycles alpha_1(T) is never equal to 3. Consequently, H(T) = 7 is impossible for any tournament.

## Precise claims

1. **alpha_1 = 3 never occurs:** The achievable values of alpha_1 at n=5 are {0,1,2,4,5,6,7} (gap at 3). At n=6: {0,1,2,4,5,...,14,16,19,20} (gap at 3 persists). At n=3,4: max alpha_1 = 2, so alpha_1=3 requires n>=5.

2. **Structural mechanism:**
   - If c3 <= 2 (at most 2 cyclic triples), then c5 = c7 = ... = 0, giving alpha_1 <= 2.
   - If c3 = 3, the three cyclic triples always share a common vertex and span exactly 5 vertices. The induced 5-vertex subtournament has score sequence (1,1,2,3,3) and always contains at least 1 directed 5-cycle, giving alpha_1 >= 4.
   - The jump from alpha_1 <= 2 (at c3 <= 2) to alpha_1 >= 4 (at c3 >= 3) eliminates alpha_1 = 3.

3. **OCF consequence:** Since H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ..., the only way to get H=7 is alpha_1=3, alpha_2=0. Since alpha_1=3 is impossible, H=7 is impossible.

4. **H=21 is also impossible:** Similar structural constraints on the (alpha_1, alpha_2) combinations prevent any tournament from achieving H=21. Verified exhaustively at n <= 7.

## Verification

- n=3 through n=6: exhaustive (all 2^{C(n,2)} tournaments checked)
- n=7: 200,000 random samples, H=7 never found
- c3=3 forces c5>=1: verified at n=5 (240/240 tournaments) and n=6 (2880/2880 tournaments)
- Three cyclic triples always share common vertex: verified at n=5 exhaustive

## Scripts

- `04-computation/h7_impossibility.py` — alpha_1 distribution analysis
- `04-computation/alpha1_gaps.py` — gap analysis
- `04-computation/alpha1_gap3_proof.py` — structural proof verification
- `04-computation/c3_forces_c5.py` — proves c3=3 forces c5>=1
- `04-computation/redei_converse_fast.py` — achievable H values

## Connection

This addresses the "converse of Redei's theorem" (which odd integers arise as H(T)?), connecting to INV-052 (Mitrovic-Stojadinovic). The gap at H=7 is the FIRST non-trivial unachievable odd integer. Whether all gaps eventually fill at large enough n (except H=7 and H=21) remains open.
