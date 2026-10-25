# Independent audit: squarefree projective-support supplement

Date: 2026-10-25. Verdict: PASS, with no required mathematical corrections. This is a support/divisibility correspondence, not descent of integer hourglass equations. Both arithmetic and geometric quantifiers are correctly separated. No files in the older absolute-geometry packet or any repository were changed; no primary source paper is included in this audit bundle.

## 1. Arithmetic proof audit

Fix a squarefree positive center e=product p_i with distinct p_i congruent to 1 modulo 4. The inherited Gaussian theorem uses fixed prime order and fixed conjugate choices. Every norm-e² Gaussian integer is, up to a unit, a product of the conjugate prime factors with exponents 1+j_i and 1-j_i, j_i in {-1,0,1}. The zero index alone has zero coordinate product. For a nonzero index, absolute coordinate product and norm recover the unordered absolute coordinate pair; its symmetries are exactly Gaussian units and conjugation. Thus no further gap identifications occur beyond j equivalent to -j. This establishes the actual bijection, rather than merely a count comparison. Reduction modulo three identifies exactly these equivalence classes with P^(k-1)(F_3), since F_3 has precisely two nonzero scalars.

At active p_i, only one of its two conjugate Gaussian primes divides the normalized Gaussian integer. If the rational prime divided either its real or imaginary coordinate, conjugation would force the other Gaussian prime as well, a contradiction. Therefore active p_i does not divide d. At inactive p_i, z contains the rational factor p_i, so d contains p_i². Since e is squarefree, truncating valuations at those of e² gives exactly gcd(e²,d)=(e/e_S)². Higher balanced valuation is harmless and genuinely occurs at e=65,d=3000. The count for exact nonempty support S is 2^|S| sign assignments divided by the free two-element global sign action, hence 2^(|S|-1). It follows in particular that full support means coprime to e.

Deleting inactive coordinates T removes the rational multiplier u=product_(i in T) p_i from z and u² from d, without changing the remaining Gaussian indices. Conversely u² dividing d forces every deleted coordinate inactive by the active valuation argument. Consequently the entire divisibility slice is exactly u² D_(e/u), not only a subset or a count equality. Independently, this slice equality follows directly from endpoint squares: u² divides both e²±d, so u divides their integer square roots. Division preserves every already-existing Schur relation, positivity, and distinctness. It does not construct a relation.

The all-zero index and k=0 are properly excluded from ordinary projective space and represented by an empty gap set. Deleting all coordinates is consistent with D_1 being empty. Squarefreeness is essential to the {-1,0,1} index model and to the displayed support/gcd formula; the note does not assert them for arbitrary prime powers.

## 2. Geometric proof audit

For any finite field, exact projective support S consists of nonzero coordinates on S modulo a common nonzero scalar. Fixing one coordinate to one gives (F_q^*)^(|S|-1). Scheme-theoretically this is a split torus; its closure in projective space is the coordinate projective subspace with supports contained in S. The sum of the stratum counts equals ((1+(q-1))^k-1)/(q-1)=1+q+...+q^(k-1), with the polynomial interpretation also valid at q=1. The reference to closure is to scheme strata, not to the discrete finite set of rational points with its induced topology.

At q=3 these toric strata agree with the arithmetic support partition. Coordinate-subspace containment and deletion have genuine divisibility meanings via the proved integer maps. General projective hyperplanes do not acquire such a meaning. Prime labels p_i and the geometric field-size variable q are different parameters; varying q is not an interpolation of integer gaps. Standard projective space over Z and its toric model exist, but no scheme of hourglass equations has been shown isomorphic to that model.

## 3. Point convention: primary-source and direct verification

Primary source checked directly: Oliver Lorscheid, Lectures on F1, https://lorscheid.org/notes/Lectures_on_F1.pdf , Lecture 1 section 3 (pointed monoids and local morphisms) and Lecture 3 section 6, printed page 16. The latter explicitly explains the discrepancy between expected closed-point sets and literal Hom(Spec F1,-), then motivates rank spaces and Weyl extension. This is a convention-dependent warning, not a universal assertion about every F1 framework.

In the specified free pointed-monoid projective model, a homogeneous prime is determined by which variables vanish. Removing the irrelevant ideal leaves exactly the nonempty supports. The residue group at S is free abelian of rank |S|-1. A local morphism to {0,1} sends all units to 1 and all nonunits to 0; this uniquely constructs the morphism from Spec F1 with that image point. Thus literal Hom has 2^k-1 elements in precisely this model.

There is also a direct chart proof that does not appeal to terminology. In chart i, x_i=1 and each other coordinate has image 0 or 1. On an overlap, the newly inverted coordinate must equal 1. Gluing therefore identifies two chart points exactly when their nonempty supports agree. Every support S occurs in precisely |S| charts and becomes one morphism. Singleton supports are the k closed points and have residue rank zero. At q=1 only singleton torus strata survive in the counting polynomial. Thus the three k-count notions coincide here, while raw Hom differs for k>=2. These claims do not mix unpointed and pointed monoid conventions, or replace raw Hom by the Weyl extension without notice.

## 4. Decisive hostile and loss boundary

The explicit e=65 triple indexed by (1,0),(0,1),(1,1) has gaps 4056,3000,3696. Those index vectors are projectively dependent, but none of those integers is the sum of the other two. In fact direct independent enumeration gives D_65={2016,3000,3696,4056}, with no distinct Schur triple anywhere in the set. This decisively refutes dependence-implies-Schur; it does not address the reverse implication globally. The projective indexing and support poset preserve their specified data, but not integer gap addition or Gaussian phase multiplication. No F1 hourglass descent, global hourglass exclusion, or RH claim follows.

## 5. Independent exact verification

The producer verifier and all three negative mutations were rerun from copies. Normal and optimized outputs are byte-identical. Wrong support count, wrong truncated gcd, and wrong torus count each raise RuntimeError in both modes.

A separately written implementation, independent_check.py, enumerates positive Pythagorean legs u<v satisfying u²+v²=e² and forms 2uv, rather than using the producer's Gaussian allocations or endpoint-square scan. It reproduces counts 0,1,1,4,4,13,40 at e=1,5,13,65,85,1105,32045. It verifies every support stratum and every square-divisibility deletion slice, including all sixteen slices at e=32045.

A second independent path constructs full scalar orbits on nonzero vectors rather than normalizing their leading coordinate. For k=1,...,5 and q=2,3,5,7 it checks projective totals and every support count. This is finite-prime-field testing; the general proof covers all prime powers. Direct gluing of 0/1 affine charts checks literal pointed-monoid Hom counts and overlap multiplicities through k=5.

Independent normal/-O outputs agree; an explicit negative exception is detected in both. All validation uses exceptions, not Python assertions. Producer SHA256SUMS verified. REVIEWED_SOURCE_HASHES.json pins the reviewed supplement and inherited arithmetic proof; the audit SHA256SUMS pins all local deliverables.

Reproduce: python replay_verify.py; python -O replay_verify.py; python replay_negative.py; python independent_check.py; python -O independent_check.py. The command python independent_check.py --negative must fail, as must its -O counterpart. No remote changes were made.
