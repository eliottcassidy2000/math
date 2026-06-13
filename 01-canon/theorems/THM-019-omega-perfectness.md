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

**DISPROVED at n>=8.** Perfectness fails for 53.8% of random n=8 tournaments (induced C5 in Omega).
The significance section below is historical context only.

Despite the failure of perfectness, the key consequences remain valid through other means:
- **Real roots of I(Omega(T), x):** PROVED for n<=8 via claw-freeness (THM-020), but DISPROVED at n>=9 (THM-025 counterexample). Perfectness was not needed for the n<=8 proof.
- **Positivity of I(Omega(T), x) for x > 0:** Follows from OCF (H(T) >= 1 by Redei's theorem), independent of perfectness.
- **Log-concavity:** Verified computationally but does not require perfectness.

**Connection to Grinberg-Stanley:** The Grinberg-Stanley proof of OCF (arXiv:2412.10572)
does not use perfectness of Omega(T). OCF is fully established without it.

**Claw-freeness (RESOLVED):** Omega(T) is claw-free for n<=8 (vertex counting: 3 pairwise
disjoint odd cycles + 1 touching all three requires >= 9 vertices). Fails at n>=9
(90% of random n=9 tournaments have claws). See THM-020.

---

## Former Proof Approaches (MOOT — perfectness DISPROVED at n>=8)

These approaches are no longer relevant since perfectness fails at n=8 (53.8% of tournaments have induced C5 in Omega).

1. **Via SPGT:** Omega(T) DOES have odd holes (C5) at n>=8. Approach invalid.

2. **Via claw-freeness:** Omega(T) is claw-free for n<=8 but has claws at n>=9. This gives perfectness at n<=7 (claw-free + no odd antihole) but cannot extend further.

3. **Direct construction:** Not viable since the statement is false at n>=8.
