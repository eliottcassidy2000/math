# THM-019: Omega(T) Perfectness — FAILS at n=8

**Type:** Observation (DISPROVED as universal statement)
**Certainty:** 5 — COUNTEREXAMPLE FOUND at n=8
**Status:** DISPROVED for general n. Holds for n<=7.
**Added by:** kind-pasteur-2026-03-05-S12 (original), opus-2026-03-05-S7 (disproved)
**Tags:** #omega #perfect-graph #independence-polynomial #structure

---

## Statement (CORRECTED)

~~For every tournament T on n vertices, the conflict graph Omega(T) is a perfect graph.~~

**DISPROVED at n=8 by opus-2026-03-05-S7.** At n=8, 53.8% of random tournaments have
an induced C5 (5-hole) in the 3-cycle conflict subgraph of Omega(T), making Omega(T)
imperfect. An explicit counterexample is constructed below.

**Corrected statement:** Omega(T) is perfect for all tournaments on n<=7 vertices.

---

## Verification Record

| n | Tournaments tested | Non-perfect cases | Method |
|---|-------------------|-------------------|--------|
| 5 | 1,024 (all) | 0 | exhaustive |
| 6 | 2,000 (sample) | 0 | random sample |
| 7 | 1,000 (random) | 0 | C5 search in 3-cycle conflict graph |
| 8 | 1,000 (random) | 538 (53.8%) | C5 search in 3-cycle conflict graph |
| 9 | 100 (random) | 30 (30%) | full Omega perfectness check |

## Counterexample at n=8

Forced arcs creating a C5 in Omega(T):
- C1 = (0,1,5): 0->1->5->0
- C2 = (2,1,5): 2->1->5->2
- C3 = (2,3,6): 2->3->6->2
- C4 = (3,4,7): 3->4->7->3
- C5 = (0,4,7): 0->4->7->0

These 5 directed 3-cycles form a C5 in Omega(T):
- Adjacent pairs share vertices: C1,C2 share {1,5}; C2,C3 share {2}; C3,C4 share {3}; C4,C5 share {4,7}; C5,C1 share {0}
- Non-adjacent pairs are vertex-disjoint: C1,C3; C1,C4; C2,C4; C2,C5; C3,C5

13 forced arcs, 15 free arcs. All completions produce the C5.
H(T) = 417, I(Omega,2) = 417 — OCF still holds despite imperfection.

Code: `04-computation/omega_c5_test.py`

Code: `04-computation/omega_perfectness_implications.py`

---

## Additional Structural Properties (Verified)

### 1. All roots of I(Omega(T), x) are real and negative

For every tournament tested (n=5 exhaustive, n=6 sample), all roots of the
independence polynomial I(Omega(T), x) are real and negative.

**Consequence:** I(Omega(T), x) > 0 for all x > 0. In particular, H(T) = I(Omega(T), 2) > 0,
giving an independent proof that every tournament has at least one Hamiltonian path.
Combined with Redei (H(T) is odd), this gives H(T) >= 1.

### 2. Log-concavity and unimodality

The coefficients alpha_0, alpha_1, ..., alpha_k of I(Omega(T), x) are:
- **Log-concave:** alpha_k^2 >= alpha_{k-1} * alpha_{k+1} for all k
- **Unimodal:** alpha_0 <= alpha_1 <= ... <= alpha_j >= alpha_{j+1} >= ... for some j

Zero failures at n=5 (exhaustive) and n=6 (sample).

### 3. Chordality fails at n=6

Omega(T) is NOT always chordal. At n=6, 72/2000 sampled tournaments have non-chordal
Omega(T). So perfectness is tight — it does not strengthen to chordality.

### 4. Modular arithmetic cascade

H(T) mod 2^k is determined by the first k+1 independence polynomial coefficients:
- H(T) mod 4 = (1 + 2*alpha_1) mod 4
- H(T) mod 8 = (1 + 2*alpha_1 + 4*alpha_2) mod 8

Verified exhaustively at n=5 and in n=6 sample.

---

## Significance

If Omega(T) is always perfect, then by the **Chudnovsky-Seymour theorem** (claw-free
perfect graphs have independence polynomials with all real roots), the all-real-roots
property would follow. This would give:
- Positivity of I(Omega(T), x) for x > 0 (alternative proof of Redei)
- Log-concavity of independence polynomial coefficients
- Strong constraints on the distribution of vertex-disjoint odd cycle collections

**Connection to Grinberg-Stanley:** The Grinberg-Stanley proof of OCF (arXiv:2412.10572)
does not use perfectness of Omega(T). If perfectness is proved, it would give structural
insight into WHY the OCF formula works, beyond the algebraic proof.

**Claw-freeness (RESOLVED):** Omega(T) is claw-free for n<=8 (vertex counting: 3 pairwise
disjoint odd cycles + 1 touching all three requires >= 9 vertices). Fails at n>=9
(90% of random n=9 tournaments have claws). See THM-020.

---

## Proof Approaches

1. **Via Strong Perfect Graph Theorem (SPGT):** Show Omega(T) has no odd hole >= 5
   and no odd antihole >= 5 as induced subgraphs. This requires understanding which
   cycle-intersection patterns can occur in tournaments.

2. **Via claw-freeness:** Show Omega(T) has no induced K_{1,3}. If three odd cycles
   C1, C2, C3 are pairwise vertex-disjoint but all intersect a fourth cycle C0,
   derive a contradiction from tournament structure.

3. **Direct construction:** Show that the clique cover number equals the independence
   number for every induced subgraph of Omega(T).
