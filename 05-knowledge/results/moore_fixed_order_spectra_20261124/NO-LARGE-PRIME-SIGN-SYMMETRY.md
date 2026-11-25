# No large prime order in an admissible Gewirtz sign symmetry

Consider the [full Gewirtz transport hypothesis](GEWIRTZ-SEVEN-EXCLUSION.md). Its branch-type matrix \(C=D-I\) is symmetric with diagonal \(-1\), other entries in \(\{\pm1\}\), and \(C^2=0\) over \(\mathbf F_3\). The number \(d_{ij}\) of sign-swapping routes for each branch pair lies in \(10,13,\ldots,46\) by the sharp-set restrictions. This note rules out every prime-order symmetry of these **sign data up to gauge** at primes \(p\ge17\). In particular the sign pattern cannot have a transitive symmetry group. This does not restrict arbitrary automorphisms of a Moore graph unrelated to the hypothesized Gewirtz labeling, or rule out asymmetric sign patterns.

## 1. Gauge integration for odd order, and one- or two-cycle obstructions

Suppose a permutation \(g\) of odd prime order \(p\) takes \(C_{ij}\) to \(s_i s_j C_{ij}\). Iterating \(p\) times shows that the accumulated sign is one common value \(\epsilon\) at every coordinate, since all matrix entries are nonzero. Replacing \(s_i\) by \(\epsilon s_i\) leaves their pair products fixed and makes the products around all orbits 1 (including fixed sites). Choose signs \(w_i\) along each orbit by integration. After conjugating \(C\) by the diagonal signs \(w_i\), the matrix is genuinely invariant under \(g\). This can complement \(d_{ij}\) in 56, so it preserves the allowed range.

For \(p\ne3\), over a splitting field any nonconstant root of unity for a moving \(p\)-cycle has an \(m\)-dimensional frequency space, where \(m\) is the number of moving cycles. Fixed sites have no contribution at nonconstant frequencies. If \(m=1\) the scalar square-zero part must vanish. This makes all columns in that cycle identical, hence all rows identical, forcing \(d=1\) for an internal pair, impossible. It covers all primes above 28 (indeed an order-\(p\) permutation of 57 sites then moves just one cycle).

For \(m=2\) the two internal circulant polynomials \(f_1,f_2\) each have constant coefficient \(-1\). Trace zero at every nontrivial root gives \(\Phi_p\mid(f_1+f_2)\), so each lag's combined coefficient is \(-2\bmod3\). Since the sums \(-2,0,2\) are distinct modulo 3, both internal blocks consist entirely of \(-1\). The nonconstant frequency matrix has zero diagonal. If \(3\) has even order modulo \(p\), inversion of a primitive \(p\)-th root lies in its Frobenius orbit. Thus the two off-diagonal entries, which are conjugate by inversion, must have product zero only when both are zero. Again the columns within each cycle repeat. This applies for \(p=17,19\), since \(3^8=-1\bmod17\) and \(3^9=-1\bmod19\).

In the two large remaining cases with fixed sites there is also an elementary row bound. Gauge fixed-site rows so their sign on the first moving cycle is \(+1\); each row is constant on each moving cycle and there are at most \(2^{m-1}\) remaining patterns. For two fixed rows of the same pattern the number \(w\) of differing entries is at most the number \(F\) of fixed sites. The integer formula from \(D=C+I\) is \(d=w-2\) when \(C_{ij}=+1\), and \(d=w+1\) when \(C_{ij}=-1\).

- For \(p=17,m=3\), \(F=6\) and at most four patterns, so two compatible-pattern rows would have \(d\le7\).
- For \(p=23,m=2\), \(F=11\) and at most two patterns. In a class of at least six fixed sites, every pair must have sign \(-1\), since \(+1\) forces \(d\le9\). Now any two agree on all six or more of those fixed coordinates as well as both moving cycles, so \(w\le5\), again impossible.

## 2. Complete finite check for three moving 19-cycles

The remaining prime and orbit pattern is \(p=19,m=3\). The trace argument across the three internal circulants forces their signs at every lag to agree (three \(\pm1\) entries sum to 0 modulo 3 only if all equal). Let \(g(x)\) be that common length-19 sign block. Its constant coefficient is \(-1\), and its other coefficients are symmetric at lags \(\pm i\). If all its coefficients are \(-1\), each nonconstant-frequency block has zero diagonal; rank at most one and the Hermitian relation force it to vanish, giving repeated interior rows. So \(g\) is nonconstant amongst 511 possibilities.

The nonconstant frequency field is \(\mathbf F_{3^{18}}\) because \(3\) has order 18 modulo 19. Write reversal as \(f\mapsto\bar f\). The three-by-three block matrix is Hermitian square zero, hence rank one at nonconstant frequency. Name the upper off-diagonal polynomials \(f,h,e\), in positions \(12,13,23\). Its minors require
\[
 f\bar f=h\bar h=g\bar g,\qquad e=\bar f h/g
 \quad\bmod\Phi_{19}.
\]
For the first two identities, equality modulo \(\Phi_{19}\) is equivalently equality of their cyclic autocorrelations modulo 3: the difference would be a constant 19-term string, and its lag-zero difference is already zero. Use these 9 nontrivial symmetric lag digits as a bucket key.

Here is the exact finite enumeration:

1. Bucket all \(2^{19}-2\) nonconstant sign strings, taking each \(f,h\) only up to independent rotations and overall sign. These correspond to changing the origins and the gauges on the second and third cycles. Interchanging the two cycles allows \(h\ge f\).
2. For each common symmetric \(g\), compute the field residue \(\bar f h/g\). A nonzero residue modulo \(\Phi_{19}\) determines at most one actual sign string: two would differ by the same constant in every coefficient, which cannot take a nonconstant \(\pm1\) string to another \(\pm1\) string unless it is zero. Discard when no such \(e\) exists.
3. Check the constant-frequency part separately. With \(G=g(1)\) modulo 3, the corresponding rank-one identity is \(G e(1)=f(1)h(1)\); if \(G=0\) also require \(e(1)=0\), as the zero-diagonal symmetric square-zero matrix here has rank at most one. The autocorrelation tests already ensure \(f(1)=h(1)=0\) in that case.
4. Compute each relative row dot product over the integers, and hence \(d_{ij}=(56-(C^2_{ij}+3C_{ij}))/2\) for each orbit type and lag. The replay checks divisibility, and finds that for every feasible affine triple the minimum complementary swap count \(\min_{i\ne j}\min(d_{ij},56-d_{ij})\) is **1**. In fact every surviving triple has a cross block equal to a signed rotation of \(g\).

The receipt is **511** symmetric nonconstant blocks, **2269** bucket pairs, **1687** compatible affine triples, and **zero** allowable sign patterns. No numerical search tolerance is used. This proves the full order-19 exclusion, not just exclusion of a 57-cycle or of order 19 with fixed branches.

Together with the analytic cases above, every prime \(p\ge17\) is excluded in the gauge-permutation stabilizer (primes larger than 57 cannot occur). A transitive group on 57 branches would have order divisible by 19 by orbit-stabilizer, contradiction to Cauchy's theorem. None of this requires an automorphism of a completed Moore graph. Reproduce the finite case with:

    c++ -O3 -std=c++20 order19-signs.cpp -o order19-signs
    ./order19-signs
