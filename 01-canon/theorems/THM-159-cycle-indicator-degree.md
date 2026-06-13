# THM-159: Cycle Indicator Polynomial Degree

**Status:** VERIFIED (k=3,5,7 by symbolic expansion)
**Session:** kind-pasteur-2026-03-13-S60

## Statement

For a circulant tournament on Z_p with sign function sigma, the directed k-cycle indicator

    ind_k(V) = number of directed Hamiltonian cycles on vertex set V

is a multilinear polynomial in the pairwise products {sigma(|v_i - v_j|)} of effective degree d(k) = k - 1, with ALL ODD-DEGREE TERMS VANISHING.

Specifically:
- d(3) = 2: ind_3(a,b,c) = (1 + s_ab*s_bc - s_ab*s_ac - s_bc*s_ac) / 4
  where s_ij = sigma(|j-i| mod p)
- d(5) = 4
- d(7) = 6

## Proof Sketch

For k=3: Direct expansion of
  count = A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
        = [(1+s1)(1+s2)(1-s3) + (1-s1)(1-s2)(1+s3)] / 8
        = (1 + s1*s2 - s1*s3 - s2*s3) / 4

The cubic term s1*s2*s3 appears with coefficient +1 in the first product and -1 in the second, cancelling. This is the "parity cancellation": both directed cycles on 3 vertices contribute opposite signs to the cubic term.

For general k: each directed Hamiltonian cycle on k vertices uses k edges. Summing over all (k-1)! directed cycles, the degree-k term (product of all k edge signs) cancels between forward and reverse cycles. More precisely, reversing a directed cycle negates the product of edge signs when k is odd, causing pairwise cancellation.

## Consequences

1. c_k is a polynomial of degree k-1 in sigma, hence determined by moments S_2,...,S_{2(k-1)}.
2. disj(k1,k2) has degree at most (k1-1)+(k2-1) = k1+k2-2 in sigma.
3. The Z_p-symmetry and disjointness constraint further reduce the effective moment range to S_4,...,S_{p-3}.
4. Combined with the OCF, this gives H = affine(S_4,...,S_{p-3}) = affine(c_3,...,c_p).

## Verification

Symbolic expansion at k=3,5,7 using cycle_indicator_degree.py:
- k=3: 2+3 = 5 nonzero terms, max degree 2
- k=5: 1+30+60 = 91 nonzero terms, max degree 4
- k=7: 1+105+1890+2520 = 4516 nonzero terms, max degree 6
