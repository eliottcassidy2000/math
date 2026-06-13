# THM-050: Corrected Consecutive-Position Formula for M[a,b]

**Status:** PROVED (combinatorial proof + exhaustive verification n=3,4,5, sampled n=6)
**Added:** opus-2026-03-06-S26
**Proved by:** opus-2026-03-06-S26 (correcting kind-pasteur-2026-03-06-S25c)

## Statement

**Theorem.** For any tournament T on n vertices and any a ≠ b:

    M[a,b] = sum_{j=0}^{n-2} (-1)^j * N(a,b,j)

where N(a,b,j) = #{Hamiltonian paths P with {a,b} at consecutive positions {j, j+1}}.

Equivalently, N(a,b,j) = C(a,b,j) + C(b,a,j) where C(a,b,j) = #{paths with a at position j and b at position j+1}.

## Proof

**Case 1: A[a][b] = 1.** For any partition S ∪ R = V\{a,b}, every pair of sub-paths (path ending at a on S∪{a}, path starting at b on R∪{b}) concatenates through edge a→b to form a valid Hamiltonian path. So C(a,b,j) = Σ_{|S|=j} E_a(S∪{a}) · B_b(R∪{b}) and the alternating sum equals M[a,b]. Also C(b,a,j) = 0 since A[b,a] = 0.

**Case 2: A[a][b] = 0.** Then A[b][a] = 1. By THM-030 symmetry (Corollary 3), M[a,b] = M[b,a]. Applying Case 1 to M[b,a] with edge b→a gives M[b,a] = Σ_j (-1)^j C(b,a,j). Since C(a,b,j) = 0, we get N(a,b,j) = C(b,a,j).

## Correction Note

The original claim by kind-pasteur used only C(a,b,j) (directed, non-symmetrized). This fails at n=3 (6/8 mismatches), n=4 (64/64), and n=5 (960/1024). The symmetrized version N = C(a,b,j) + C(b,a,j) is necessary and sufficient.

## Corollaries

**Corollary 1.** M[a,b] = 0 for all a ≠ b (i.e., M is scalar) iff for every vertex pair {a,b}, the alternating consecutive-position count vanishes.

**Corollary 2.** The total edge usage Σ_{a,b,j} C(a,b,j) = (n-1)·H.

## Verification

- n=3: 8/8 match
- n=4: 64/64 match
- n=5: 1024/1024 match
- n=6: 30/30 sampled match

## Source Files

- `04-computation/consec_formula_correct.py` — proof and verification
- `04-computation/consec_formula_validity.py` — failure analysis of original claim
