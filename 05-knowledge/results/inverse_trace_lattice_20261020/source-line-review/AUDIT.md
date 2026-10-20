# Independent audit: source-line specializations and irreducible boundary

20 October 2026. **PASS under the retained global hypotheses and the explicitly cited external theorems.**

Reviewed addendum: `session11-jacobian/SOURCE_LINE_AND_IRREDUCIBLE_BOUNDARY.md`.

SHA-256: `4af75e8618e937dd345b94478ceb6e93d80cc8e41ca26e7cf5b0837a59ddad7e`.

This audit is separate from the frozen finite-selector/Euler-gate audit. The latter's files and receipts were not changed.

## Findings

1. **All scalar specializations have irreducible support.** Finite base change of the actual graph map G to Z=a is finite surjective from the source line onto V(R_a). The source line has nonempty irreducible curve image, because P(a,y) remains monic of positive degree in y. Thus R_a is neither zero nor a constant and has one irreducible support. At this stage its multiplicity is not yet known; the note explicitly preserves that distinction.

2. **Normalization is A¹, and the parametrization is birational.** Factoring the finite line map through normalization gives a finite dominant morphism A¹→the smooth normalization. Its projective completion is dominated by P¹, so has genus zero. Its affine complement is nonempty, and the inverse image of that complement is contained in the single point at infinity of P¹, so that complement is exactly one point. The normalized curve is therefore A¹. The lifted map is a polynomial q. Because the original line parametrization has nonzero derivative everywhere by the Keller determinant, the chain rule rules out every degree greater than one for q. Hence the original parametrization is the normalization, not a nontrivial multiple cover.

3. **No unibranch singularities.** Nonzero derivative gives a smooth local image for each normalization branch. If a target point had only one normalization preimage, finiteness would exclude an additional local branch, and the local inverse theorem on a nonzero coordinate derivative would make that point smooth. Thus singularities require distinct parameters. This is differential immersion, not an assertion that the entire normalization map is globally injective or an algebraic embedding.

4. **Multiplicity one is justified.** At a general smooth point p of D_a, there is a unique source-line preimage s, and consequently a unique preimage of w=(p,a) under G. Finiteness puts the full inverse image of a sufficiently small neighborhood of w inside an inverse-chart neighborhood of s. The Keller local inverse makes the image a smooth graph Z=x(F⁻¹(u,v)), and surjectivity makes this the entire germ of W. Because R defines a reduced hypersurface, it is locally a unit times this graph equation. Slicing Z=a cuts by a source coordinate, with nonzero differential, so the section has multiplicity one. This proves that R_a itself is irreducible, not merely that its support is irreducible. The argument uses the actual finite graph and étaleness essentially.

5. **Boundary location and smoothness.** On H=0, R_a=aA+B. The B-unit result at O excludes D_a∩O, and elsewhere on C the equation is exactly −B/A=a. The frozen gate bounds every ambient F-fiber on C by one point. Consequently the line normalization has at most one point over D_a∩C, and the preceding immersion argument makes D_a smooth at each such point. Therefore all singularities of D_a lie outside C.

6. **Every line must collide for d>1.** The cited one-line injectivity theorem applies to the affine source line x=a. Its contrapositive gives noninjectivity of the normalization parametrization for every a when F is not an automorphism. Since its normalization is A¹, the Euler formula is χ_c(D_a)=1−Σ(b_p−1), and a collision forces χ_c(D_a)≤0. Those collisions occur outside C by the preceding item. This does not imply generic smoothness of the specializations; a finite normalization can identify points in every scalar slice.

7. **Irreducible-boundary lower bound.** When C itself is irreducible, its polynomial parametrization gives A¹ normalization. The cited self-intersection theorem forces a normalization fiber of size at least two, so χ_c(C)≤0. Combining with the already reviewed gate gives #O≥d−1 and hence the stated degree-six lower bound of five distinct intersections. Reducible boundaries are not silently treated as irreducible, and branch counts are not replaced with delta invariants.

## External theorem verification

The versioned primary text [Z. Jelonek, A note on the Jacobian Conjecture, arXiv:2011.03472v3](https://arxiv.org/html/2011.03472v3), dated 1 February 2021, was independently opened and checked during this review.

Theorem 1.1 supplies the strong self-intersection exclusion. Its proof treats affine-line normalizations and topological affine lines, including the unibranch possibility; the manuscript is not relying merely on exclusion of an algebraically smooth A¹. Theorem 2.7 records Gwoździewicz's injectivity-on-one-line theorem, with the 1993 original paper listed as reference 4. The addendum accurately labels its use of the versioned preprint rather than claiming a new publisher-version comparison.

## Small exact controls

The separate independent verifier passes 29 checks in normal and optimized Python, with identical outputs and JSON. It checks an immersed nodal normalization, ramification of nontrivial polynomial-cover compositions, tangent nonreduced sections, the source-coordinate nonreduced hostile F(x,y)=(y²,xy), persistence of nodal identification in a product family, the constant-term-one Laurent surface, and the distinction between node and cusp Euler contributions.

The producer's separate verifier was replayed from a copy: 192 checks pass in normal and optimized Python, with byte-identical receipts matching the producer's final JSON. Its finite birational hostile has graph uZ²−v², special section −v², and Jacobian −2y²; it correctly exhibits why finite graph structure alone does not prove reduced source-line specializations.

These are finite algebra controls, not actual noninvertible Keller realizations or substitutes for the universal arguments and cited theorems. No heavy historical census was run.

Only this review staging directory was written. No producer proof, earlier audit, repository checkout, Git state, or remote was changed.
