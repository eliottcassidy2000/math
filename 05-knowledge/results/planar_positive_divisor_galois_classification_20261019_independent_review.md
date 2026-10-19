# Independent review: positive-divisor Galois classification

Date: 2026-10-19. Verdict: **accepted; no mathematical gap found** in the reviewed statement and proof.

Reviewed producer artifact: `packet135-stage/05-knowledge/results/planar_positive_divisor_galois_classification_20261019.md`.

Reviewed SHA-256: `4c6a8a93c3e07c49a83aa4d172ff539cb10206960c9e1dbe7d58f44ad30e4db3`.

## Specific checks

1. **Critical values and indecomposability.** Substituting a critical point rho*zeta factors every primitive term by rho^(a+1)*zeta; the remaining alternating sum is the stated nonzero beta-integral value kappa. Thus the critical points have distinct values. In a hypothetical polynomial composition, every root over a critical point of the outer polynomial must be the same point, forcing the inner polynomial to be a translated pure power. Two outer critical points are then impossible; a single outer critical point forces the entire composition to have only one critical point. This proves the required indecomposability and hence primitivity.

2. **Group theorem.** Independently inspected Jones's primary PDF, Corollary 1.3 (page 160), and its stated dependence on the classification of finite simple groups. The hypothesis is a single cycle with at least three fixed points, without a prime-length requirement. The local cycles here are single cycles because no other critical point has their branch value. For m >= 3 the (a+1)-cycle has precisely m fixed points. The two m=2 cases are correctly handled by transpositions and primitivity. Parity distinguishes the symmetric cases, while the square discriminant gives the alternating upper bound.

3. **Actual fixed-time source group.** All nonzero differences of finite branch values are a nonzero algebraic-over-K factor times E^(-(a+1)/a); none can be a nonzero element of K. Hence B and B+lambda are disjoint. The intersection of the two finite Galois splitting fields can ramify only at infinity. Tame Riemann-Hurwitz over the algebraically closed constant field forces this intersection to be trivial. Galois linear disjointness then survives passage to the source root subfield C(s), preserving the full permutation group and irreducibility. This is a proof at every specified nonzero constant time, rather than an unsupported specialization of a generic-time claim.

4. **Arithmetic square class.** For a and r even, the monic derivative is (n/E^r) h^2. The resultant discriminant formula leaves square class (-1)^((n-1)/2)*n = (-1)^(a/2)*n. Its square factor remains nonzero after source substitution by separability. A constant square class cannot become a square in the pure transcendental field K(E,s) unless it was already a square in K. Since the arithmetic group contains the geometric A_n and lies in S_n, the discriminant gives exactly the stated alternatives.

5. **Individual coordinates and radicals.** The rational-map degrees m+1 and 2m+a+2 are correctly computed and neither is divisible by n. Therefore neither coordinate function belongs to the target field. Linear disjointness prevents either from entering the enlarged source field. Primitivity of the natural A_n/S_n action excludes proper intermediate root fields, so each coordinate generates the full degree-n extension. This also descends to K(E,s) by degree bounds. The natural actions are faithful, so the normal closure of each coordinate field has the claimed full group. The radical criterion yields exactly the S4 quartic exception; all remaining degrees are at least five.

## Scope and evidence

This review is of the universal proof, not a replacement of it by numerical or finite symbolic checks. No redundant checker run was required for this review. The source exponential and birational parameter identities are inherited from the previously reviewed packet 128. The theorem remains restricted to the positive-divisor locus and nonzero times in K; it does not assert numerical specialization, invariant-dependent times, or a general Jacobian-conjecture conclusion.

Primary group-theorem source independently inspected: Gareth A. Jones, *Primitive permutation groups containing a cycle*, Corollary 1.3, https://www.cambridge.org/core/services/aop-cambridge-core/content/view/BBF3021D8434AD2C0BB0D3D6CE637008/S000497271300049Xa.pdf/primitive-permutation-groups-containing-a-cycle.pdf .

No correction is requested. This review is staged separately from the producer packet; it may be copied into the final packet without modifying its mathematical assessment.
