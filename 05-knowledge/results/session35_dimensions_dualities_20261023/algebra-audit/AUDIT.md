# Independent algebra audit: dimension horizon and marked tournament recursions

2026-10-23. Verdict: **PASS for the stated dimension theorem and E/D marked-resolvent formulas**, with the verification-hardening request below. No Git operations. This audit distinguishes mathematical proof from finite controls and does not certify historical repository inheritance claims or originality.

## 1. Dimension-positivity theorem

Audited `session35-root/DIMENSION_POSITIVITY_HORIZON.md`, its exact check script, and its output.

The stated fixed-order theorem is correct for **t >= 0**. For real t, simultaneous positivity at every order holds exactly at nonnegative integers. Negative t fails already at order one; this does not mean negative t fails every individual order (for example G_2(-1) is PSD). The note correctly retains the t >= 0 hypothesis for its fixed-order classification.

Independent derivation of the obstruction: put r=floor(t)+2 for positive nonintegral t. Restrict the S_k-indexed Gram matrix to permutations of r letters. Every unused letter contributes t, so the restricted matrix is t^(k-r)G_r(t). The sign vector has eigenvalue t(t-1)...(t-r+1), exactly one factor of which is negative. The prefactor is strictly positive. This establishes failure at every k >= r, not just at r.

The central convolution scalar is correctly normalized: expansion of the class function into irreducible characters gives scalar k! times character multiplicity divided by the Specht degree. Schur-Weyl, hook-content, and hook-length yield product_(i,j in lambda)(t+j-i). Interpolation is legitimate because equality holds at infinitely many integer dimensions. If t >= k-1, every individual factor is nonnegative; actual tensor representations handle all smaller nonnegative integers. At t=0 and k>=1 the whole matrix is zero.

Signed Burnside recurrence and the division by k! are correct. For t=7/2 and k=5 the alternating **trace** is -7/256; the corresponding convolution **eigenvalue** is 5! times this, namely -105/32. Distinguishing these avoids a normalization error. The first failure for 57/2 is order30. The proof does not depend on rational versus irrational parameter type.

The warnings about braided quantum dimension, retention of the original nonpositive object after doubling, absence of a tournament functor, and finite-depth positivity versus all-order positivity are appropriate. No category-equivalence, Lie-dimension, Collatz, or RH theorem follows.

## 2. E/D identities and independent structural checks

Audited `session35-recursions/RESULT.md` and its initial script/output.

For any real skew-symmetric M, its orthogonal real normal form consists of zero blocks and two-dimensional rotation blocks. Consequently

- p(-z)=(-1)^n p(z);
- F(z)=a_0/z + sum_j a_j z/(z^2+lambda_j^2), with a_j >= 0;
- F is odd and F(z)/z is rational in z^2.

This independently explains the branch independence of both D formulas. From D^2=diag(2M^2-I,2M^2-I), its squared characteristic polynomial equals the square of det((z^2+1)I-2M^2). Both candidates are monic, so equality has positive sign. The skew spectral blocks give det((z^2+1)I-2M^2)=2^n p(w)^2, w^2=(z^2+1)/2, for odd and even n alike. In particular the odd-n zero eigenvalue contributes the necessary factor z^2+1, not an erroneous square-root branch.

In the D marked solve, the off-diagonal contributions reduce to 2 u^T M R u=0 because M R is skew-symmetric. The remaining term is 2z u^T R u. The real spectral expansion above then proves F_D=(z/w)F(w) without choosing a branch or assuming w != 0 at the level of rational identities.

For E, writing the endpoint unknowns as x,y and a=x+y,b=x-y gives z a+b=2 and a=(z+2F)b-2F. The interior solution has marked sum F(1-b). Hence total a+F(1-b) is precisely the stated Möbius transform. Eliminating the interior block gives Schur determinant z^2+1+2zF. Thus both p_E and F_E are correct.

These arguments apply to arbitrary real skew matrices, not merely +/-1 tournament matrices; the independent controls exploit that extra freedom rather than repeating the producer's tournament enumeration.

## 3. Exceptional roots, cancellation, and minimality

“Continuation at apparent poles” must not mean that actual marked-resolvent poles disappear. The correct interpretation is equality in the rational function field followed by cancellation of common factors where possible. At poles that remain, no finite value is claimed.

Concrete exact controls:

- Directed three-cycle C has p_C=z(z^2+3) but F_C=3/z: the uniform mark misses its two nonzero eigenvalues.
- Its extension has F_E(C)=(5z^2+3)/(z(z^2+7)); the old roots of z^2+3 are absent from the marked resolvent even though p_E retains that factor.
- F_D(C)=6z/(z^2+1). At z=+/-i, w=0 produces genuine poles.
- For the two-vertex skew block J, F_D(J)=4z/(z^2+3). At z=+/-i, the apparent w=0 singularity is removable, and F_D(J)(i)=2i.

Thus branch points of the presentation are not themselves a reason to omit cases; cancellation and genuine poles must be distinguished. The producer was asked to make this distinction explicit.

The same-spectrum/different-extension example proves that unmarked spectrum alone is insufficient. It does **not** prove universal minimality or completeness of (p,F). The note expressly disclaims both and restricts E to the canonical source/interior/sink section. That is the correct scope.

## 4. Zeta and duality scope

The matrix identity Q_D=2 I_2 tensor Q gives Z_D(s)=2^(1-s) Z(s), an everywhere-nonzero multiplicative factor, so no new zeros are created. The tower geometric series has possible poles at s=1+2 pi i k/log2, with possible numerator cancellation. The note correctly says possible poles, not mandatory poles. None of this implies Riemann's functional equation or RH.

The positive-integer odd-core decomposition, theta heat scaling by4, and Mellin prefactor1/2 are consistent. The marked tournament Q-energy scaling is2 rather than4. Completion removes the trivial zeta zeros by gamma-pole cancellation for m>=1. If “m” is left unquantified in zeta(-2m), the conventional positive-integer restriction should be understood or added; m=0 is not a trivial zero.

## 5. Independent exact evidence and limits

`independent_check.py` does not import producer code, enumerate contents, call adjugates, or form matrix inverses. It uses:

1. Direct permutation-indexed rational Gram matrices and exact Schur-complement PSD elimination: 64 tests, orders1..4, quarter-integer parameters (negative cases tested only at order1).
2. A direct signed S_5 sum for the t=7/2 negative eigenvalue.
3. Indexed-entry constructions of E,D on 24 weighted integer skew matrices of sizes1..6, exact LU linear solves and determinants at three Pell-compatible (z,w) pairs: 288 checks.
4. Five symbolic pole/cancellation controls above.

These are independent finite error detectors, not substitutes for the all-order proofs. All checks use explicit exceptions and remain active under python -O. Normal and optimized outputs are compared byte for byte.

The producer's original recursion script used bare asserts; a hardening request was sent to replace these with explicit exceptions and compare optimized output. This is a verification robustness issue, not a counterexample to the formulas.

## 6. Producer revision and additive Cayley coordinate

The producer has replaced bare assertions with explicit RuntimeError checks and reported byte-identical normal/optimized outputs. The note now explicitly distinguishes removable presentation singularities from genuine meromorphic poles. These close the two verification/wording requests.

The added Cayley identity also passes: writing a=z^2+1,b=2z, the sum and difference of numerator and denominator of F_E give (a+b)(1+F) and (a-b)(1-F), respectively. Their ratio is ((z+1)/(z-1))^2 times C. This is a rational identity, including the usual canceled-value interpretation. The transitive p,F formulas hold at n=1 and n=2 and propagate by E, proving them for both parities. The all-depth additive law is correct; it is not a proof that D squares the Cayley coordinate.
