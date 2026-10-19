# Independent review: negative-divisor Galois classification

Date: 2026-10-19. Verdict: **accepted; no mathematical gap found and no correction requested**.

Reviewed artifact: `packet137-stage/05-knowledge/results/planar_negative_divisor_galois_classification_20261019.md`.

Reviewed SHA-256: `1a660b81bfc24ece05e36de0baf742543aa725f21d6b04fe14c11596c0cf1d23`.

## Focused findings

- **Kummer case r=2.** The equation is exactly Z^n=z^n−n lambda, n=b−1. At a simple zero of its right side, the geometric polynomial is Eisenstein. This gives the full degree after adjoining roots of unity as well. The splitting field presentation permits every cyclotomic automorphism to fix Z, E, and z; hence the displayed semidirect product has an actual section, rather than merely being an extension with an unspecified splitting. The cyclotomic action is the stated usual action.

- **Individual coordinates when n is composite.** The proof does not incorrectly invoke primitive cyclic actions. Since Z is transcendental over the enlarged constant field, equality of coordinate functions under Z -> zeta Z implies equality as rational functions. The finite pole set of chi or finite zero set of theta is precisely the set EV^b=1. Its invariance forces zeta^b=1, while zeta^n=1 and gcd(b,n)=1 force zeta=1. Thus both stabilizers are trivial, and each individual coordinate has full degree. This includes all composite n and descends to the original field by degree bounds. The n=1 exception is harmless.

- **Critical values, signs, and small cases for r>=3.** Substitution of rho*zeta into each primitive term gives rho^(b−1)*zeta^(-1). The remaining beta-integral factor has sign (-1)^(r−1) and magnitude exactly the stated kappa. Zero is critical only when b>=3, and its value is distinct from the nonzero critical values. Critical injectivity gives indecomposability as in packet 135. The cycle at zero fixes (r−2)b points; for b=2, the other cycles fix r−2 points. These are enough for Jones except for the explicitly treated cubic and quintic. The cubic has a transposition; the quintic has a single 3-cycle in a primitive group and therefore contains A5. I independently checked the cited small-cycle assertion in Kreso–Tichy Remark 3.11, pages 9–10 of https://arxiv.org/pdf/1601.07316 .

- **Fixed-time transfer and arithmetic group.** Nonzero critical-value differences scale as E^(-(b−1)/b), so no such difference is a fixed nonzero element of K. The disjoint finite branch loci, Galois intersection, and tame Riemann-Hurwitz argument from packet 135 therefore apply without a sign or exponent change affecting the conclusion. The actual source group and degree survive base change. When b and r are even, n is odd and the monic derivative is (n/E^(r−2)) times a polynomial square. The resultant discriminant leaves precisely square class (-1)^((n−1)/2)n = (-1)^(b/2−1)n. Separability preserves a nonzero square factor after source substitution; a constant square class is unchanged by the pure transcendental source field extension.

- **Coordinate nonfactorization for r>=3.** At each nonzero critical point, chi has pole order r, which is not divisible by ramification index r−1. Thus chi cannot factor through the primitive. Theta has its only finite pole at zero, but the primitive's zero fiber contains a nonzero point because its zero multiplicity b−1 is less than n. A factorization through the primitive would force another finite pole there. Both arguments are valid for the small cases too. Linear disjointness then prevents either function from becoming rational in the source field, and the primitive A_n/S_n action gives full individual generation. The normal closures and radical conclusions follow, with exactly the cyclic family and the stated cubic exception on this negative locus.

The review concerns the universal mathematical proof. No broad computational replay was performed. Source-flow existence and birational identities are inherited from the reviewed packet 128. The claims are restricted to the stated negative-divisor locus, nonzero scalar times in K, and unspecialized source variables; no general Jacobian-conjecture conclusion is involved.

This review is staged outside the producer packet and may be copied into it without changing the assessment.
