# THM-156: Energy-Disjointness Identity

**Status:** PROVED (kind-pasteur-2026-03-13-S60)
**Verified:** p = 7, 11, 13, 17, 19 (all orientations, zero error)

## Statement

For any circulant regular tournament on Z_p with connection set S:

    disj3(S) = (p/4) * E(S) + C(p)

where:
- disj3(S) = number of pairs of vertex-disjoint 3-cycle vertex sets
- E(S) = additive energy of S = |{(a,b,c,d) in S^4 : a+b = c+d mod p}|
- C(p) = A(p) - B(p)/2 with:
  - A(p) = (K(p) + c3(c3-1))/2
  - B(p) = m^5/5 + m^4/2 + m/20 + mp/4
  - K(p) = -3p(p^2-1)(p^2-9)/320 (from THM-155)
  - c3 = p(p^2-1)/24, m = (p-1)/2

## Proof

Chain of three identities:

**Step 1.** From THM-155 (disjoint 3-cycle identity):
disj3 = -T5/2 + A(p), where T5 = tr(A^5)/5 = # directed 5-cycles.

**Step 2.** Eigenvalue expansion of T5:
Each non-trivial eigenvalue is lambda_t = -1/2 + iD_t.
Re(z^5) = -1/32 + 5D^2/4 - 5D^4/2.
Let S4 = sum_{t=1}^{m} D_t^4. Then:
5T5 = m^5 - m/16 + 5mp/8 - 5S4.

**Step 3.** Parseval identity for additive energy:
E(S) = (1/p) sum_t |lambda_t|^4, and |lambda_t|^2 = 1/4 + D_t^2.
So: E*p = m^4 + m/8 + mp/4 + 2S4.

**Step 4.** Eliminate S4 between Steps 2 and 3:
From Step 3: S4 = (Ep - m^4 - m/8 - mp/4)/2.
Substitute into Step 2: T5 = B(p) - pE/2.
Therefore: disj3 = -T5/2 + A = pE/4 + (A - B/2).

## Verified values of C(p)

| p  | C(p)      |
|----|-----------|
| 7  | -19.25    |
| 11 | 316.25    |
| 13 | 1345.50   |
| 17 | 9962.00   |
| 19 | 21617.25  |
| 23 | 78366.75  |

## Significance

This identity bridges tournament combinatorics and additive combinatorics:
- **Paley minimizes E** (QR sets have uniform Fourier spectrum |S_hat(t)|^2 = (p+1)/4 for p = 3 mod 4, achieving Cauchy-Schwarz equality)
- **Paley minimizes disj3** (since slope p/4 > 0)
- **Paley maximizes H** (from BIBD analysis: fewer disjoint pairs => more alpha_1 => higher H via OCF)

The energy minimization of QR sets is a consequence of the Gauss sum identity |g(chi)|^2 = p and the multiplicative structure of the Legendre symbol.

## Related

- THM-155 (disjoint 3-cycle identity): prerequisite
- HYP-671: hypothesis entry
- HYP-672: E determines H uniquely
- Scripts: thm156_energy_disjointness.py, energy_disjointness_bridge.py, paley_energy_minimizer.py
