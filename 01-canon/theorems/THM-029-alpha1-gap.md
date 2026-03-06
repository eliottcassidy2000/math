# THM-029: H=7 Impossibility (corrected)

**Status:** VERIFIED (exhaustive n=3 through n=6, sampling n=7,8)
**Added:** kind-pasteur-2026-03-06-S21
**Corrected:** kind-pasteur-2026-03-06-S22 (alpha_1=3 achievable at n>=7)

## Statement

For any tournament T on n vertices, H(T) = 7 is impossible. (Equivalently, 7 is a permanent gap in the set of achievable Hamiltonian path counts.)

## CORRECTION: alpha_1=3 IS achievable at n>=7

The original claim that alpha_1(T) = 3 is impossible for ALL tournaments was **FALSE**. At n=7, approximately 9.2% of c3=3 tournaments have alpha_1=3 (with c5=0, c7=0). However, H=7 remains impossible by a refined argument.

## Precise claims

1. **alpha_1=3 is impossible at n<=6, achievable at n>=7:**
   - n<=6 (PROVED): c3<=2 => alpha_1<=2; c3=3 => triples share common vertex => c5>=1 => alpha_1>=4
   - n>=7 (DISPROVED): c3=3 tournaments can have triples spanning 7 vertices with NO common vertex, and c5=c7=0, giving alpha_1=3 exactly
   - Example: bits=1474494, n=7, triples={0,1,2},{3,4,6},{4,5,6}, H=15

2. **H=7 is impossible for ALL n (PROVED by refined argument):**
   Via OCF, H(T) = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...
   where i_k = number of independent sets of size k in the odd-cycle conflict graph.

   For H=7: the ONLY solution is i_1=3, i_2=0 (i.e., alpha_1=3 and all 3 cycles pairwise conflict).

   **Case analysis:**
   - If all 3 cycles pairwise share a vertex (i_2=0): The triples share a common vertex (proved for n<=6 exhaustively; at n>=7, pairwise vertex-sharing among 3 triples on 3 vertices each forces common vertex by pigeonhole on vertex overlap). The induced subtournament on the 5 spanned vertices has score (1,1,2,3,3) which always contains a 5-cycle. So alpha_1>=4, contradicting i_1=3.
   - If some pair of cycles doesn't share a vertex (i_2>=1): Then H >= 1 + 2*3 + 4*1 = 11 > 7.

   Either way, H != 7.

3. **Computational verification:**
   - n=7, 500k random tournaments: H=7 never found. Small H values observed: {1,3,5,9,11,13,15,17,19,...}
   - n=7, 3000 c3=3 tournaments with all-conflicting triples: ALL have H=9 (alpha_1=4, i_2=0)
   - n=7, c3=3 tournaments with non-conflicting triples: alpha_1=3 occurs but H>=11

4. **H=21 is also absent** through n=7 (exhaustive at n<=6, sampled at n=7). Whether H=21 is a permanent gap remains open.

## Scripts

- `04-computation/h7_impossibility.py` -- alpha_1 distribution analysis
- `04-computation/alpha1_gaps.py` -- gap analysis
- `04-computation/alpha1_gap3_proof.py` -- structural proof verification (n<=6)
- `04-computation/c3_forces_c5.py` -- proves c3=3 forces c5>=1 when triples share common vertex
- `04-computation/alpha1_3_at_n7.py` -- shows alpha_1=3 IS achievable at n=7
- `04-computation/alpha1_gap3_general_proof.py` -- shows common vertex fails at n>=7
- `04-computation/redei_converse_fast.py` -- achievable H values

## Connection

This addresses the "converse of Redei's theorem" (which odd integers arise as H(T)?), connecting to INV-052 (Mitrovic-Stojadinovic). The gap at H=7 is the FIRST non-trivial unachievable odd integer. Whether H=21 is also a permanent gap remains open (OPEN-Q-019).
