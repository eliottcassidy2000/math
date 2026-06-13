# THM-213: Pfaffian of Paley Skew-Adjacency Minor

**Status:** VERIFIED (computational, 4 cases)
**Found by:** opus-2026-03-14-S89c
**Verified in:** `04-computation/pi_ratio_89c.py`

## Statement

For the Paley tournament P_p (p ≡ 3 mod 4, p prime), let S be the skew-adjacency matrix (S_{ij} = 1 if i→j, -1 if j→i, 0 if i=j). Then:

**det(S_{00}) = p^{(p-3)/2}**

where S_{00} is the (p-1)×(p-1) principal minor obtained by deleting row 0 and column 0.

Equivalently, since p-1 is even:

**Pf(S_{00}) = ±p^{(p-3)/4}**

## Verified Cases

| p | det(S_{00}) | p^{(p-3)/2} |
|---|------------|-------------|
| 3 | 1 | 3⁰ = 1 ✓ |
| 7 | 49 | 7² = 49 ✓ |
| 11 | 14641 | 11⁴ = 14641 ✓ |
| 19 | 16983563041 | 19⁸ ✓ |

## Proof Sketch

The eigenvalues of S are: 0 (once), +i√p ((p-1)/2 times), -i√p ((p-1)/2 times).

For a circulant skew-symmetric matrix, all cofactors are equal. The cofactor of S_{00} is:

det(S_{00}) = product formula from eigenvalues...

Since S has eigenvalues 0 and ±i√p with known multiplicities, the cofactor equals:

det(S_{00}) = (i√p)^{(p-1)/2} × (-i√p)^{(p-1)/2} / p = (√p)^{p-1} / p = p^{(p-1)/2} / p = p^{(p-3)/2}

Wait: this uses the matrix-tree-like theorem for skew matrices. The cofactor of the zero eigenvalue direction gives the "complexity" analog. Since the all-ones vector is the null eigenvector, and the circulant symmetry means all diagonal cofactors are equal, we get:

det(S_{00}) = (1/p) × |∏ nonzero eigenvalues| = (1/p) × (√p)^{p-1} = p^{(p-3)/2} ✓

## Connection

- The eigenvalues ±i√p come from Gauss sums involving e^{2πi/p}
- The Pfaffian p^{(p-3)/4} is related to the half-sum of eigenvalue phases
- For large p: Pf ≈ p^{p/4}, growing as the fourth root of p^p

## Related

- THM-212: p | H(P_p) (Burnside on Hamiltonian paths)
- The det formula det(A(P_p)) = (p-1)/2 × ((p+1)/4)^{(p-1)/2} (different matrix)
