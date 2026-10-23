# Session 35 independent zeta/theta audit

Verdict: PASS. No mathematical correction to session35-zeta/RESULTS.md is required. Its correction of THM-1926 is valid and narrowly scoped. Analytic results below apply to every positive integer dimension n, not just the finite test universe. No Git operations or repository edits were performed.

## 1. Finite Gamma obstruction

For a real skew M, Q=I+M^T M is real symmetric positive definite. Orthogonal diagonalization defines Q^{-s} unambiguously using real log(q_j). A finite sum of exponential functions is entire and Z_Q(-m)=sum q_j^m>0. Multiplication by Gamma gives the claimed nonremovable simple poles. In the length variable the residues, if needed, are 2(-1)^m pi^m tr(Q^m)/m! at s=-2m. The factor 2 is from Gamma(s/2); RESULTS does not give an incorrect length-variable residue.

The finite heat inversion obstructions follow directly from t tending to zero. With constant terms the right side diverges, without them it decays exponentially despite its polynomial prefactor. Positivity and nonempty finite spectrum suffice. The rational adjacency zeta has constant numerator 1, hence no finite zeros; zero eigenvalues do not create reciprocal-eigenvalue poles.

## 2. Infinite tower: convergence and remainder details

The block multiplication gives D(M)^2=diag(2M^2-I,2M^2-I). Thus the tower has squared-energy eigenvalues 2^k q_j with multiplicity 2^k. Its Dirichlet sum converges exactly on Re(s)>1, and the geometric sum gives the stated meromorphic continuation. The denominator derivative at every candidate s_k is log 2, so the residue is Z_Q(s_k)/log 2 unless canceled. Since Z_Q(1)>0, the real pole is genuine.

For fixed t>0, the k>=0 heat terms are bounded by n 2^k exp(-q_min 2^k t), and the k<0 terms of P are bounded by n t 2^k. These establish absolute convergence; on each compact subinterval of t>0 the bounds are uniform. Reindexing then proves P(2t)=P(t).

Write S(t)=sum_{r>=1} 2^{-r} h_Q(2^{-r}t). Absolute interchange of its exponential Taylor series is justified by
sum_j sum_{r>=1} 2^{-r} exp(|t|q_j 2^{-r}) <= sum_j exp(|t|q_j/2).
Therefore S(t) is entire, its coefficients are exactly tr(Q^m)/(m!(2^{m+1}-1)), and H=P/t-S is an exact identity, not a merely formal asymptotic expansion.

A useful stronger bound for real t>0 is
0 <= H(t)-P(t)/t+n <= t tr(Q)/3.
Moreover,
|H(t)-P(t)/t+n-t tr(Q)/3| <= t^2 tr(Q^2)/14.
These follow from 0<=1-exp(-x)<=x and |1-exp(-x)-x|<=x^2/2. They directly certify the O(t) remainder in RESULTS. No assertion that every candidate complex pole is genuine is made. The scalar Q=[1] example does make all of them genuine.

## 3. Theta lift and Epstein continuation in every dimension n>=1

The Fourier transform of exp(-pi t x^T Q x) is det(Q)^(-1/2)t^(-n/2) exp(-pi y^T Q^{-1}y/t). Gaussian Poisson summation proves the matrix theta identity. Because Q(DM)=diag(2Q,2Q) in the original integer coordinates, independent summation over the two Z^n factors proves exact theta doubling. An orthogonal spectral diagonalization alone would not justify that integer-lattice factorization, but the actual block identity does.

The lattice Dirichlet series converges absolutely for Re(s)>n by comparison with sum_{v!=0}|v|^{-Re(s)}. Termwise Mellin integration is valid there. Splitting at 1 yields exactly the two rational terms -2/s and 2 det(Q)^(-1/2)/(s-n) appearing in RESULTS. Exponential theta decay makes both remaining integrals entire by compact-uniform domination, including arbitrary s derivatives. Thus Lambda_Q has genuine simple poles with residues -2 at 0 and 2 det(Q)^(-1/2) at n. They are distinct for every allowed n>=1.

Applying this split formula to Q^{-1}, n-s verifies the functional equation with the determinant multiplier. At s=-2m (m>=1), Lambda_Q(s)=det(Q)^(-1/2)Lambda_{Q^{-1}}(n+2m)>0: the latter is in the positive, absolutely convergent Dirichlet region. Consequently the zero of E_Q is exactly simple, not merely at least simple. Explicitly,
E_Q'(-2m)=(-1)^m m! Lambda_Q(-2m)/(2 pi^m),
which is nonzero. At zero, the residue ratio gives E_Q(0)=-1. E_Q itself has only a simple pole at n, with residue 2 pi^(n/2)/(Gamma(n/2)sqrt(det Q)). No parity or even-order restriction on n is needed.

For B=Q/(det Q)^(1/n), det B=1. Hence Lambda_B(s)=Lambda_{B^{-1}}(n-s), and their sum has the claimed reflection. Multiplication by s(s-n) cancels both poles and is invariant under s->n-s. This is a deliberately symmetrized completion, not proof that the original Q is isodual. At n=1,Q=[1], E_Q=2 zeta. The polynomial hostile example has four explicit off-axis roots and is positive on the real axis; it correctly refutes only the insufficient inference, not RH.

## 4. THM-1926 correction and scope

For Paley tournaments of prime order p congruent to 3 modulo 4, write M=A-A^T and J for the all-ones matrix. One has A=(J-I+M)/2 and M^2=J-pI. On the orthogonal complement of the all-ones vector, M has eigenvalues +/-i sqrt(p), hence A has eigenvalues (-1+/-i sqrt(p))/2 and modulus sqrt(p+1)/2. The Perron adjacency eigenvalue is (p-1)/2, not part of that modulus assertion. The skew operator has an additional zero eigenvalue.

Exactly two adjacency-modulus statements in THM-1926 require this correction: the YAML status assertion and the Paley assertion in section 4. The THM-1880 related-item comment concerns char_S and its sqrt(p) modulus, so it should NOT be changed. Determinant identity, primitive closed-orbit Euler product, strong-component factorization, transpose/complement invariance, and low-degree closed-walk counts are unaffected. The poles associated to non-Perron adjacency eigenvalues have modulus 2/sqrt(p+1), not 1/sqrt(p). This is not a correction of the skew spectrum.

Independent exact matrix polynomial identities were checked at p=3,7,11,19,23,31,43. The derivation, rather than this finite list, explains the correction generally.

## 5. Sources independently inspected on 2026-10-23

- NIST DLMF https://dlmf.nist.gov/5.2 : confirms Gamma is zero-free, with simple poles and reciprocal-Gamma simple zeros at nonpositive integers.
- NIST DLMF https://dlmf.nist.gov/25.4 : equations 25.4.2 through 25.4.4 confirm the Riemann reflection and xi normalization cited.
- NIST DLMF https://dlmf.nist.gov/25.10 : distinguishes the negative-even trivial zeros from the unproved localization statement RH.
- NIST DLMF https://dlmf.nist.gov/20.7#viii : equation 20.7.32 at z=0,tau=it gives the scalar theta inversion. This is not by itself a citation for arbitrary matrix Q; that extension follows from the multivariate Gaussian transform and Poisson argument above.
- Elkies https://people.math.harvard.edu/~elkies/M229.18/zeta1.pdf : pages 3-4 supply the theta Mellin split and Poisson proof. Important notation caveat: these notes use xi for the meromorphic completed zeta without the polynomial s(s-1)/2; DLMF uses the entire xi. RESULTS is consistent with both and does not confuse their normalizations.

## 6. Reproduction

- Original script replay: original-replay.out, reproducing the 1099 labeled tournaments, Paley characteristic polynomials, hostile polynomial, scalar tower checks, and theta coefficients.
- Independent script: independent_check.py and independent_check.out. Added exact Paley adjacency/skew matrix identities through p=43; non-scalar three-to-six-dimensional theta coefficient enumeration; non-scalar heat tower and remainder tests; non-diagonal numerical Poisson inversion; numerical Epstein reflection and determinant-normalized symmetrization in dimensions 1 through 5.

Numerical checks are labeled as such and are not presented as rigorous numerical certificates or proofs of the general analytic statements. No blockers.
