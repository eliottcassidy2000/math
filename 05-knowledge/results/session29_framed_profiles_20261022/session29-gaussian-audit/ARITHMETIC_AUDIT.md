# Independent arithmetic nondegeneracy audit

Verdict: ACCEPTED as a proof using a cited algebraic-independence theorem. This review is separate from the interval certificate and leaves its frozen source intact.

Reviewed addendum SHA256: 805c2db2886598057c293a98eb993fd23eaf50389adc179d2e0d521ec2234a8a. Its reviewed contents are preserved in `arithmetic_frozen.md`.

## Rational operator

For an integral positive-definite Gram matrix G, the simultaneous linear equations H^T G=GH and tr H=0 have rational coefficients. Their rational kernel has a rational basis whose realification is the full Euclidean traceless self-adjoint space. The trace inner product is positive definite there: conjugation to a real symmetric matrix preserves trace of the square. Thus its rational matrix S is invertible. Each minimal-vector functional z^T GH z is rational on this basis. The spherical-average coefficient is rational. Consequently S^-1 B has rational entries and algebraic eigenvalues, independently of whether its shell is perfect. An irrational orthonormal change of coordinates cannot invalidate this similarity-invariant argument.

## General weights and reduction

For n=24m+16, k=n/2=12m+8, the extremal minimum squared norm is2m+2, so positive-degree harmonic theta vanishes at the cusp to order at least m+1. Dividing the degree2 and degree4 series by Delta^(m+1) yields holomorphic modular forms of weights−2 and0 respectively. Holomorphy away from the cusp follows because Delta has no zeros in the upper half-plane. The first quotient is zero; the second is constant. The resulting fixed angular operator has exactly the rational construction above.

A rational weight-k modular form has an E4,E6 polynomial with only even E6 exponents, since k is divisible by4. Put a=k/4=3m+2. Only c0 E4^a and c1 E4^(a−3)E6² survive two differentiations followed by E6=0. All higher E6 powers still vanish. For m=0, a=2 and the second monomial does not exist; setting c1=0 is correct. Since the actual theta value F(i)>0 and E4(i)>0, c0=F(i)/E4(i)^a>0. In particular this argument does not assume that a merely formal extremal modular form corresponds to an existing lattice.

Independently differentiating with Ramanujan's identities gives

D²(E4^a)|E6=0 = a(4a+1)E2² E4^a/36 +5a E4^(a+1)/36,
D²(E4^(a−3)E6²)|E6=0 = E4^(a+1)/2.

Substitution E2(i)=3/pi gives the producer D²F expression. The scalar Hessian normalization independently follows from

2t² F_tt/[n(n+2)] +t F_t/n,

using F_t=−2pi DF and F_tt=4pi² D²F. Dividing by b=(pi t)² Delta^(m+1) gives exactly the stated scalar formula. At1 the cancellation leaves

r(1)=2C·1728^(m+1)/[k(k+1)] − c0·1728^(m+1)/[8pi²E4(i)].

The second coefficient cannot vanish. An exact rational check at m=3 recovers both coefficients of the inherited rank88 formula, including the factor4 from y=4pi²E4(i).

## Algebraic independence and source checking

The cited Nesterenko Theorem1 applies to arbitrary complex q with0<|q|<1, not just algebraic q. At q=e^(−2pi), E6(i)=0, so the remaining triple must be algebraically independent. Then pi²E4(i)=9E4(i)/E2(i)² is transcendental: algebraicity would give a nonzero polynomial relation after clearing denominators. A nonconstant rational affine function of its reciprocal remains transcendental. Such a scalar cannot cancel an algebraic angular eigenvalue.

Primary source: Yu. V. Nesterenko, Modular functions and transcendence questions, Sbornik: Mathematics187:9(1996),1319–1348, DOI10.1070/SM1996v187n09ABEH000158. PDF: https://www.mathnet.ru/links/f3406fc7df8e07b1490855a647fef28b/sm158_eng.pdf .

Audit access detail: direct primary PDF retrieval timed out, and the independent indexed searches returned printed p1321–1322 rather than p1320. The parent/producer separately verified Theorem1 on p1320. Importantly, the independently retrieved primary p1322 Theorem2 also suffices for precisely this application: choose its three generators as q,E2(i),E4(i). All four Ramanujan values are algebraic over their generated field because E6(i)=0. The theorem's positive lower bound for nonzero integer polynomials in these generators rules out any algebraic relation. This verifies the needed special case from primary text without relying on a secondary paraphrase. Nesterenko's theorem itself remains CITED, not reproved or formally verified.

## Accepted conclusion and exclusions

Every actual extremal even unimodular lattice of dimension24m+16, m>=0, has a nonsingular determinant-one Gaussian Hessian at t=1. The argument requires neither perfection nor an angular eigenvalue calculation. It proves nonsingularity, not definiteness. It does not extend the rank88 three-critical-point shape theorem to other ranks and does not imply transcendence of the off-center minimum value.
