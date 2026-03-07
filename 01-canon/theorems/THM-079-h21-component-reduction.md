# THM-079: H=21 Component Reduction

**Status:** PROVED (partial — disconnected case ruled out; P_4 case ruled out)
**Author:** opus-2026-03-07-S39
**Date:** 2026-03-07
**Dependencies:** THM-029 (H=7 impossibility)

## Statement

### Part A (PROVED): Disconnected Omega Case
If Omega(T) is disconnected, then H(T) ≠ 21.

### Part B (PROVED): P_4 Component Case
Omega(T) cannot have a connected component isomorphic to P_4 (path on 4 vertices).

### Corollary
H(T) = 21 requires connected Omega(T) with I(Omega(T), 2) = 21 and |V(Omega)| >= 6.

## Proof of Part A

By OCF multiplicativity for disjoint union:
I(G1 ∪ G2, 2) = I(G1, 2) · I(G2, 2)

For H = 21 with disconnected Omega: need I(C_i, 2) for components with product = 21.
Since H is always odd, each factor is odd.
Factorizations: 21 = 3·7 = 1·21 = 1·3·7 = 1·1·21 etc.

The only non-trivial factorization is 3·7.
I(component, 2) = 7 requires the component G to be K_3 (3 pairwise conflicting cycles).
By THM-029: 3 pairwise-conflicting cycles in Omega(T) ALWAYS generate a 4th cycle
(either via common vertex forcing 5-cycle, or via tournament arc constraints).
So |V(G)| ≥ 4 and I(G,2) ≥ 1+2·4 = 9 > 7. Contradiction.

Therefore I(component, 2) = 7 is impossible for any component of any Omega(T).
So H = 21 via disconnected Omega is impossible. QED.

## Proof of Part B

I(P_4, 2) = 21 (direct computation: 1 + 4·2 + 3·4 = 21).
A P_4 component of Omega means 4 cycles C1-C2-C3-C4 with adjacent pairs sharing vertices
and non-adjacent pairs vertex-disjoint.

Consider the middle pair C2, C3. They share a vertex (adjacent in P_4).
Their vertex union V(C2) ∪ V(C3) has 5 vertices (3+3-1 from the shared vertex).

**Key lemma:** Two directed 3-cycles sharing a vertex on 5 tournament vertices
always force at least one additional 3-cycle (verified exhaustively at n=5:
when C1={0,1,2} and C2={2,3,4} are both cycles, minimum t_3 = 3, never 2).

The forced additional cycle C* lies on a triple from V(C2) ∪ V(C3).
Since C* shares vertices with C2 or C3, it is adjacent to C2 or C3 in Omega.
But C* ≠ C1 (C1 has vertices outside V(C2)∪V(C3)) and C* ≠ C4 (same reason).
So C* is a 5th vertex in Omega adjacent to the middle of the path, and Omega ≠ P_4. QED.

## Graph Classification for I(G,2) = 21

Connected graphs with I(G,2) = 21 for |V| ≤ 6:
- |V| = 4: P_4 (path on 4 vertices). I = 1+8+12 = 21. [RULED OUT by Part B]
- |V| = 5: No connected graphs with I = 21.
- |V| = 6: K_6 minus 2 edges (2 non-edges). I = 1+12+8 = 21. [OPEN]

## Part C (PROVED): alpha_3 Elimination

All decompositions with alpha_3 >= 1 are impossible.
alpha_3 >= 1 means 3 mutually vertex-disjoint cycles, which forces alpha_2 >= 3
(the 3 pairwise-disjoint pairs). Then 4*alpha_2 >= 12 and 8*alpha_3 >= 8,
giving 4*alpha_2 + 8*alpha_3 >= 20, leaving at most 2*alpha_1 = 0.
But alpha_1 >= 3 (the 3 cycles themselves), so 2*alpha_1 >= 6.
Total >= 6+12+8 = 26 > 20. Contradiction.

## Remaining Viable Decompositions

Only 4 decompositions of alpha_1 + 2*alpha_2 = 10 remain:
- (10, 0): 10 pairwise-conflicting cycles. Needs alpha_1=10, all pairs sharing vertex.
- (8, 1): 8 cycles with exactly 1 disjoint pair.
- (6, 2): 6 cycles with exactly 2 disjoint pairs.
- (4, 3): 4 cycles with exactly 3 disjoint pairs.

## Part D (PROVED): (4,3) Decomposition Fully Eliminated

For (4,3): the only connected graphs on 4 vertices with 3 edges are trees: P_4 and K_{1,3}.
- P_4 ruled out by Part B.
- K_{1,3} ruled out: 3 mutually disjoint leaves force alpha_3>=1, contradicting Part C.
- The only other graph with 3 edges on 4 vertices is C_3+isolated, which is disconnected.
Therefore (4,3) is fully eliminated.

## Part E: Directed Cycle Counting Correction

**Important:** Omega(T) vertices are DIRECTED odd cycles, not vertex sets. A vertex set
of size 5 can support 0, 1, 2, or 3 distinct directed 5-cycles, each a separate vertex
in Omega. Multiple directed cycles on the same vertex set are always adjacent (sharing
all vertices), forming a clique. This means alpha_1 (total cycles) >= vertex-set count,
making I(Omega,2) LARGER. The decomposition analysis is unchanged but harder to achieve.

## Part F: Tournament Forcing of Extra Cycles

**Computational finding (opus-S40):** Even when 6 three-cycle triples have the K_6-2e
conflict pattern and are realizable as directed 3-cycles in a tournament, the tournament
structure ALWAYS forces additional directed 5-cycles. These extra cycles expand Omega
beyond the 6-vertex K_6-2e skeleton.

At n=6 (exhaustive, 2160 tournaments with K_6-2e 3-cycle pattern):
- Every realization has 4, 6, or 8 additional directed 5-cycles
- I(Omega(T), 2) is always 29, 33, or 37 (minimum 29, never 21)
- The sharing lemma alone is insufficient (2790/5040 K_6-2e triple arrangements
  pass the sharing lemma test), but tournament forcing of 5-cycles provides
  the actual obstruction

## Remaining Open Question

The 5-cycle forcing at n=6 pushes I(Omega,2) >= 29 for K_6-2e realizations. At n >= 7,
the tournament forcing should be even stronger (more vertices = more cycle-forming
opportunities). But a general proof is needed showing that for ANY n, ANY tournament T
with K_6-2e (or K_8-e or K_10) in its 3-cycle conflict graph always has I(Omega,2) > 21.

Three remaining decompositions:
- (10,0): K_10. At n=8: t3=10 always has t5 >= 12 (100% in 29k samples).
- (8,1): K_8-e. Not yet directly investigated.
- (6,2): K_6-2e. I(Omega) >= 29 at n=6 (exhaustive). Needs n>=7 verification.

## Scripts

- `04-computation/h21_gap_mechanism.py` — (t3, t5, p33) exhaustive analysis at n=6
- `04-computation/t3_t5_parity.py` — parity obstruction: t3=5 forces t5 even
- `04-computation/t3_t5_parity_law.py` — complete case analysis at n=6
- `04-computation/k6_2e_omega_check.py` — K_6-2e exhaustive check (opus-S40)
- `04-computation/k6_2e_sharing_proof.py` — sharing lemma approach (opus-S40)
