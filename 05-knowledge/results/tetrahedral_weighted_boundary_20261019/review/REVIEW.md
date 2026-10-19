# Independent audit: tetrahedral weighted-boundary packet

Date: 2026-10-19. Scope: the paper, its finite verifier, its dependency ledger, and the linked primary sources.

## Verdict

**ACCEPTED, with the stated scope.** The central construction, finite normalization, divisor-class obstruction, finite-factorial-domination corollary, and global-standard-cusp classification are valid over C. The packet supplies a non-Keller example and an intrinsic obstruction for the specified field extension. It does not prove the general A4 case or the planar Jacobian conjecture.

This verdict applies to `A4_WEIGHTED_BOUNDARY.md` with SHA256

`4ac39451579cd75ef0026c4e45b198a7ce5c48dcf8bc7705d8b2ebae86afab97`.

Two wording corrections were requested and incorporated by the packet author:

1. Theorem 4 now says that the isomorphism over A² is **in the chosen target coordinates**. Without that phrase, a translate of the standard cusp would be a literal counterexample to the original phrasing, although not to its intended assertion.
2. Section 9 now places the general standard plane in H¹_et(Y_reg, μ₂), explicitly retaining the unit-squareclass versus Cl(Y)[2] alternatives. The class-group alternative alone is not supplied by THM-2655 in general.

No other substantive correction is required.

## Verification results

- Replayed a copy of the submitted verifier: 130/130 checks pass.
- Replayed that copy under Python optimization: 130/130 checks pass, with byte-identical stdout.
- Wrote and ran `independent_checks.py`, without importing the submitted verifier: 32/32 checks pass in normal and optimized modes, with byte-identical stdout.
- Verified all seven corpus dependency hashes against the current local source files.
- Verified the packet checksum manifest.

Independent checks include a companion-matrix trace calculation rather than the submitted Newton-identity calculation, local-coordinate valuation formulas, the monogenic étaleness criterion, the classical pair cubic resolvent, all order-three stabilizers in A4, and all Q8-to-V4 epimorphism kernels. These finite checks supplement the geometric audit below; they do not certify the geometric arguments by themselves.

## 1. Polynomial source rings and Zariski Main

The stable-polynomial-ring criterion is correct. Since u and v lie in A=C[r,s], their Jacobian matrix lies in M₂(A). The extended coordinate derivations supply its inverse if they preserve A. The determinant is therefore a unit of A, hence a nonzero complex constant. Conversely a constant nonzero determinant gives polynomial coefficients for both derivations.

R⊂A follows from integrality: an element integral over B is integral over A, and A is integrally closed in L. The Keller morphism Spec A→Spec B is étale and consequently quasi-finite. Zariski Main then identifies Spec A with an open of the normalization in L. The supplied reference gives precisely the affine normalization form needed here: [Stacks, Tag 03GT](https://stacks.math.columbia.edu/tag/03GT).

No unmentioned completeness, local-finiteness, analytic extension, or local-nilpotence condition is being used.

## 2. Boundary classes and codimension two

Proposition 2 is valid for a normal integral affine surface X with a dense open U≅A². Restriction injects O(X)* into O(U)*=C*, while constants give the reverse inclusion. Thus the unit boundary map is zero. Since Cl(U)=0, the Weil-divisor localization sequence identifies the free group on the omitted prime divisors with Cl(X).

There are no principal relations among those boundary primes: a rational function whose divisor is supported there becomes a regular invertible function on the normal open U, hence is constant. This proves independence as well as generation.

The treatment of isolated omitted points is correct. If no divisors are omitted, every height-one point of X lies in U. Normality then forces O(X)=O(U). Because both spaces are affine, the given open immersion is an isomorphism. The normal-domain intersection statement used in this argument is a special case M=R of [Stacks, Tag 0AVB](https://stacks.math.columbia.edu/tag/0AVB).

Therefore nonzero torsion in Cl(X) really excludes every affine-plane open, independently of coordinates. Conversely, torsion-freeness is only necessary; the packet does not promote it to sufficiency.

## 3. Finite factorial domination

Corollary 3 is valid. For a finite dominant map of normal integral varieties p:V→X of degree n, divisorial pullback is defined at height one by ramification indices, even if D is not Cartier and p is not globally flat. At a prime divisor of X, the equality Σeᵢfᵢ=n follows from the finite normalization over its DVR. Since Cl(V)=0, p*D is principal.

Pushforward of a principal divisor is the divisor of its field norm; the cited statement indeed concerns dominant proper maps and applies to this finite map: [Stacks, Tag 02RT](https://stacks.math.columbia.edu/tag/02RT). It follows that nD is principal and nCl(X)=0.

Combining bounded torsion with the free boundary basis forces no divisorial boundary for a Keller affine-plane open. The codimension-two argument then forces that open to be all of X. The resulting finite étale map A²→A² has degree one because C² is simply connected. The hypotheses require a genuinely finite dominant map from a factorial variety. An arbitrary rational parametrization, dominant map, or local factoriality statement would not suffice.

## 4. Normalization and exact ramification

The monic quartic presentation and the cone presentation are inverse ring presentations. Hence R is a normal domain, finite free of rank four over C[u,v], and is the integral closure in its fraction field. Normality can be obtained either from the finite-group invariant presentation or from the hypersurface with its unique codimension-two singular point.

An especially direct independent certificate is

∂f/∂a = 108a³−36ua−8v = 64(a³−2ab+c) = 64e in R.

The finite monogenic map is therefore étale exactly on D(e). This establishes the full ramification support without relying only on the discriminant.

On the chart a≠0, write t=b/a, so c=at². Then

- e=a(a−t)²;
- a²−b=a(a−t);
- h=1/(a−t);
- Δ=64a²t(a−t)³.

The generic points of E and C₀ lie in this chart. These identities independently give div(e)=2E, div(Δ)=C₀+3E, and ord_E(h)=−1. The support of e=0 is exactly E, including at the cone vertex. The displayed inverse chart formulas identify D(e) with A¹×Gm as schemes.

The source polynomial map has determinant 64e and contracts e=0 to the origin. Thus it is genuinely polynomial but is not Keller. Its restriction to A¹×Gm is genuinely étale. The finite completion has exactly the stated cusp branch curve, and deleting E produces the stated nonproper-value set. No conflation between these two notions is needed in the argument.

## 5. Monodromy, traces, and the two representations

The exact matrix calculation and the projective action prove that G has order 24 with scalar kernel {±I}, projective image A4, and H of order six. Its core is {±I}, so the quartic Galois closure has group A4. This is stronger than merely observing that the discriminant is a square.

The determinant-one kernel is Q8, while the full group is the reflection representation G4. The distinction is essential and is made correctly. The linked October 31 material explicitly contains the invariant polynomials and the character twist yielding Shephard–Todd group 4: [Elkies, Math 263x](https://abel.math.harvard.edu/~elkies/M263.12/index.html).

The inverse h equation, primitive-element recovery, and all four trace moments check. The companion matrix independently gives T₁=T₂=0, T₃=−64v/(9Δ), and T₄=256u/(27Δ). The three escaping roots have pole order 1/3 relative to the target parameter; the fourth root is finite at the generic cusp point. The numerator of the T₃ pole is not identically zero on the cusp.

A second algebraic route to the cyclic resolvent is useful. For the monic a-quartic, the classical pair resolvent is

(t+2u/9)³ + 64(u³−v²)/729.

After a linear rescaling, its field is exactly K(Δ^(1/3)). This agrees with the invariant-theoretic resolvent.

## 6. Class groups and boundary weights

The cone class group computation is correct. Localization at a gives a UFD. Its only missing prime is L₀, with div(a)=2L₀, so Cl(X) is generated by an element of order dividing two. Under the double-cover presentation C[s,t] over C[s²,st,t²], generic inertia along s=0 is trivial; therefore the pullback of L₀ is the single reduced line s=0. A principal L₀ would pull back to the divisor of a constant multiple of s, contradicting invariance under simultaneous sign change. Thus its order is exactly two.

The local c≠0 chart also verifies div(b−a²)=E+L₀. Consequently E is the nonzero class. The localization unit map for D(e) is multiplication by two, exactly as claimed. The minimal-resolution transition functions give O_P¹(−2), and the single weighted exceptional curve has matrix [−2]. The algebraic class-group proof is independent of the graph determinant, so no unsupported passage from an unweighted graph to torsion occurs.

For the actual cubic resolvent, Q8 acts freely off the origin. The quotient is the smooth locus of the stated D4 hypersurface. Its identification and resolution agree with [Mikosz–Weber, Section 7.1](https://link.springer.com/article/10.1007/s10801-020-00938-3). Line-bundle descent is legitimate because Pic(A²\{0})=0 and its units are constants. It gives Cl(Y)=Hom(Q8,C*)≅(Z/2)², with the cyclic action on the three nonzero characters. This class group is correctly kept distinct from Cl(X).

## 7. Global-standard-cusp classification

Theorem 4 is valid in the repaired coordinate phrasing. The following hypotheses and steps are load-bearing:

1. The whole affine branch support lies in the global standard cusp in the selected polynomial coordinates. A local cusp germ does not suffice.
2. The cubic quotient of the A4 closure is cyclic and unramified at every affine prime other than Δ. Factoring its Kummer representative in the UFD C[u,v] reduces it modulo cubes to Δ or Δ². Their fields are the same. Infinity creates no additional Kummer freedom, since every rational function is a quotient of affine polynomials.
3. Along Δ, characteristic-zero inertia is cyclic. Its nontrivial image in C3 forces order three in A4; hence its intersection with V4 is trivial. Outside Δ the closure has trivial inertia. Normal closure of an étale cover introduces no new ramification there.
4. The V4 cover is therefore unramified in codimension one over Y. Purity is applied only over the regular locus, where its hypotheses hold. The packet correctly avoids applying purity across the singular origin: [Stacks, Tag 0BMB](https://stacks.math.columbia.edu/tag/0BMB).
5. C²\{0} is simply connected, so the smooth quotient has fundamental group Q8. Complex finite-cover comparison identifies its finite étale covers with the finite topological covers; the cited theorem is exactly [Milne, Theorem 21.3](https://www.jmilne.org/math/CourseNotes/LEC210.pdf).
6. Every epimorphism Q8→V4 has kernel {±1}. The independent enumeration gives six marked epimorphisms and one kernel. Thus the connected V4 cover is unique as an unmarked cover over Y_reg. Its function field is unique over E, which also makes the isomorphism a K-isomorphism. There is no missing descent/twisting assumption here.
7. All four order-three subgroups of A4 are conjugate. Their fixed quartic fields are therefore K-isomorphic, and uniqueness of normalization finishes the classification.

Unramified omitted divisors of a quasi-finite source need not coincide with finite branch support. The final packet explicitly preserves this distinction. No replacement of branch support by the full nonproper-value set occurs in the classification proof.

## 8. Chau attribution and scope

The primary paper is correctly identified: Nguyen Van Chau, *A Remark on Vitushkin’s Covering*, Acta Mathematica Vietnamica 24 (1999), 109–115. Its Theorem 1.1 states that a polynomial complex-plane map with nonzero constant Jacobian is bijective if its restriction away from a curve homeomorphic to C is an unbranched covering. The cusp parametrization is a homeomorphism from C onto the standard cusp, so the packet’s all-degree special-case comparison is valid: [primary paper, Theorem 1.1](https://math.ac.vn/uploads/files/9901109.pdf).

The packet appropriately makes no historical novelty claim for the tetrahedral construction or classical tools, no general A4 exclusion, and no Jacobian-conjecture conclusion. Its substantive deliverable is the exact hostile example, the coordinate-independent weighted-class obstruction, and the explicitly delimited global cover classification.

## Reproduction

From the package directory in a repository checkout:

    python review/independent_checks.py
    python -O review/independent_checks.py
    python review/reproduction/verify_a4_weighted_boundary.py
    python -O review/reproduction/verify_a4_weighted_boundary.py

The review directory contains both independent outputs, the machine-readable independent result, the untouched copied companion and its outputs, and this audit. The final source identities are recorded in `independent-results.json`.
