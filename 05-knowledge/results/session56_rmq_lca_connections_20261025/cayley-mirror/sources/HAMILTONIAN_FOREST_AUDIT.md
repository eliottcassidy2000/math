# Independent audit of the Hamiltonian-size-biased forest identity

**Verdict: PASS.** The all-order proof in `session31-root/HAMILTONIAN_FOREST_SIZE_BIAS.md` is valid for labelled uniform tournaments, n>=2. No novelty verdict is asserted. Independent exact code uses parent-function forests and permutation Hamiltonian paths, not determinants or the original subset-DP implementation; every labelled tournament at n=2,3,4,5 passes. Output: `audit_hamiltonian_forests.json`.

## Proof audit

1. **Expectation is legal here.** A directed rooted forest contains at most one of the two possible arcs on each unordered pair. Its indicator monomial therefore uses independent edge variables. Summing forest monomials justifies replacement by fair or conditional edge probabilities. A generic determinant-expectation rule would be false; the proof does not use one.
2. **Size bias equals path conditioning.** Writing H as the sum over n! vertex permutations gives E[kappa F]=E[F | 1->2->...->n], because each specified path has probability 2^{-(n-1)}. Relabelling invariance is needed and holds for F.
3. **The conditional Laplacian is correct.** Each fixed adjacent edge changes the baseline weights 1/2,1/2 to 1,0. The endpoint diagonal changes are +1/2,-1/2, interior changes cancel, and the two off-diagonals are -1/2,+1/2. Thus the claimed B is exact and B*1=0.
4. **The quotient argument is correct.** J maps into the constant line, so its induced map on the quotient is zero even though B need not preserve the usual orthogonal complement. The conditional Laplacian induces nI/2+Bbar. Both full operators preserve the constant line. Therefore det(aI+tB)/a is exactly the forest determinant, with polynomial continuation through a=0. No unjustified orthogonal diagonalization is used.
5. **Tridiagonal characteristic polynomial is correct.** The recurrence det(xI-2B)=x f_{n-1}(x), f_0=1,f_1=x,f_j=xf_{j-1}+f_{j-2}, holds including n=2. The Chebyshev rescaling gives the claimed imaginary cosine eigenvalues. At even n the two zero eigenvalues need not yield two independent eigenvectors; the determinant argument remains valid.
6. **Positive covariance is coefficientwise, not pointwise tournament ordering.** Pairing roots proves nonnegative coefficients of the correction, as made explicit below. This does not imply F orders every pair of tournaments by H or that the forest count reconstructs H.
7. **Asymptotic regime is correctly stated.** For fixed t>0, the ratio is 1+1/n+O_t(n^-2). The t=0 equality is a distinct boundary; uniformity down to t=0 is not proved.

## Real polynomial form and strict coefficient covariance

The recurrence also has the elementary closed form

\[
f_m(x)=\sum_{s=0}^{\lfloor m/2\rfloor}\binom{m-s}{s}x^{m-2s}.
\]

With \(a=1+nt/2\), the audited identity is equivalently

\[
\mathbb E[\kappa F_T(t)]
 =\sum_{s=0}^{\lfloor(n-1)/2\rfloor}
 \binom{n-1-s}{s}\left(\frac{t^2}{4}\right)^s
 a^{n-1-2s}.
\]

The s=0 term is E[F]=a^{n-1}. Every remaining coefficient is nonnegative. The s=1 term is

\[
\frac{n-2}{4}\,t^2(1+nt/2)^{n-3}.
\]

Hence for n>=3 the covariance of kappa with the coefficient [t^j]F is strictly positive for every j=2,...,n-1; it is zero for j=0,1. This proves the proposed coefficient refinement with no additional hypothesis.

## Independent score-energy proof

Let \(S=\sum_v(d_v-(n-1)/2)^2\). Without conditioning, every d_v is Binomial(n-1,1/2), so E[S]=n(n-1)/4.

After conditioning on a specified Hamiltonian path, its two endpoints have means (n-1)/2 ±1/2 and variances (n-2)/4. Each of its n-2 interior vertices has mean (n-1)/2 and variance (n-3)/4. Therefore

\[
\mathbb E[\kappa S]
 =2\left(\frac{n-2}{4}+\frac14\right)
 +(n-2)\frac{n-3}{4}
 =\frac{n^2-3n+4}{4}.
\]

Subtracting E[S] proves exactly

\[
\operatorname{Cov}(\kappa,S)=-\frac{n-2}{2}.
\]

This independent calculation agrees with [t^2]F=(N^2-sum d_v^2)/2, N=n(n-1)/2, and the forest covariance (n-2)/4. Degree correlations between vertices are irrelevant because only the sum of individual squared deviations is used.

## Independent finite evidence

The independent checker enumerates all 2,8,64,1024 labelled tournaments at n=2,3,4,5. It counts H by vertex permutations and F by choosing, for each vertex, either a parent along an outgoing arc or root status, rejecting cycles. All expected/biased coefficient identities and the exact score-energy covariance agree. This is additional verification, not the reason the all-order proof is accepted.

## Consequence for the mirror note

The relationship between the tree coefficients and tournaments is not merely a shared integer: E[F_T(2)]=(n+1)^{n-1}=n!a_n. The H-size-biased correction supplies an exact Hamiltonian interaction. This is a genuine ensemble bridge, but it retains neither individual tournament identity nor Hamiltonian-spectrum realization. The directed triangle still has F(1)=7 and H=3, and H=63 is attained. No gap7/gap21 explanation or spectrum-completeness theorem follows.
