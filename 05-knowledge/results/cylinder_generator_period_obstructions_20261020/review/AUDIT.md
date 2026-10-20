# Independent audit: fixed and moving cylinder obstructions

Date: 20 October 2026.

## Verdict

**PASS, under the stated hypotheses and cited classical inputs.** The fixed-cylinder finite-period theorem, moving-algebra Theorems A and B, quadratic-coordinate recognition Proposition C, and conditional actual-source Corollary D are sound. No counterexample or missing mathematical hypothesis was found.

The review independently confirmed a useful strengthening before freeze: the moving proof needs uniqueness of plane embedding only on a nonempty Zariski-open set of base values. Whole-algebra generation, rather than an unavailable exceptional-fiber embedding theorem, supplies the nonzero cotangent determinants on every remaining collision fiber. The producer incorporated that strengthening. One final wording clarification explicitly places the analytic period calculation inside the uniqueness open set.

This is a theorem about specified entire source algebras, not a proof of the complex Jacobian conjecture. In particular, Corollary D still requires a monic quadratic element of the actual source algebra. The existence of that element for arbitrary source algebras of conductor degree four or six has not been proved.

Reviewed source hashes are recorded in `source_hashes.json`. This reviewer wrote only in this review directory and did not edit the producer notes or repository.

## 1. Entire algebra and specialization

Let R=C[x], B=R[y], and E_g=R[y²,yg(x,y²)], with g(x,0)=c in C*. Its exact normal form is

    E_g = R[W] ⊕ y g(x,W) R[W],   W=y².

The relation Z²−Wg² is monic in Z. Reducing a representative to a(W)+Zb(W) and substituting proves that the claimed presentation has no extra kernel: its even and odd parts cannot cancel, and g is nonzero.

This also justifies specialization, including leading-degree drops and constant-g fibers. The specialized normal form maps to a(y²)+yg(ξ,y²)b(y²). Because g(ξ,0)=c≠0, the specialized polynomial is never identically zero. The same parity argument makes this map injective. Thus E_g/(x−ξ)E_g is exactly the image algebra C[y²,yg(ξ,y²)], not a larger nonreduced quotient.

More directly, R[P,Q]=E_g provides polynomial expressions for the original generators in P,Q and conversely. Those expressions specialize, giving equality of the specialized image algebras at every ξ. This is substantially stronger than C(x)[P,Q]=C(x)[y²,yg]. It supplies genuine closed embeddings of the same abstract fiber curve.

The fixed-cylinder tensor product has the same specialization property. For its node ideal I=R m, subtracting the endpoint evaluation of a polynomial in P,Q proves I=(P−p,Q−q) over R itself. Consequently the two generating cotangent classes define a surjective 2 by 2 matrix over R. Its determinant is a unit of R, hence a nonzero complex constant.

## 2. Fixed-cylinder finite periods

The normalization is C[y]. The conductor is a nonzero ideal of that ring. If r≠s map to the same point, every conductor element h satisfies h(r)=h(s) and r h(r)=s h(s), because both h and yh belong to the curve ring. Thus r and s lie in the finite conductor support.

An intrinsic automorphism extends to an affine transformation of the normalization and permutes that support. Two distinct points determine an affine transformation, so the intrinsic automorphism group is finite. This argument does not require the automorphisms to extend to the plane.

At the chosen singular image, the classes of the two plane coordinates form a two-dimensional cotangent basis. Each intrinsic automorphism identifies this space with the cotangent space at its image singular point by one fixed invertible complex matrix. Combining that matrix with the globally constant determinant above leaves only finitely many possible ambient Jacobian scalars. The periods transported from the corresponding finitely many intrinsic images likewise form a finite set. A complex polynomial in x with finite image is constant.

The finite-number-of-embedding-classes variant is also valid: choose one fixed representative for each class, together with its parametrization. There are then finitely many representative periods and cotangent identifications, each modified by only the finite intrinsic automorphism group. No assertion that all rational plane curves have finitely many embedding classes is being used.

The polynomial closed-form step is exact. The exterior derivative of Φ*(U dV)−δ W dZ vanishes when det DΦ=δ; polynomial one-forms on the affine plane admit polynomial primitives when closed. The primitive evaluates equally at the coincident image endpoints. The fixed-endpoint boundary term also vanishes because equality of the endpoint polynomials for Q implies equality of their x-derivatives.

No individual period needs to be nonzero. An independent control below has a genuine collision and zero chosen period; its zero derivative still excludes a nonzero base-only Jacobian.

## 3. Scope of the embedding input and parametrized correction

A fresh primary-source check verified the prime-exponent uniqueness theorem and the separate one-line injectivity theorem. The exact sources and access method are in `PRIMARY_SOURCE_CHECK.md`. The first is used only for equivalence of embedded images/defining ideals, not equality of arbitrarily rescaled equations or a relative ambient automorphism over R.

For a fixed nonconstant g, the conductor is g(y²)C[y]. Its finite support is symmetric and has at least two points. Every intrinsic affine transformation y↦αy+β preserving that support has β=0 by its centroid. Preservation of the complete conductor divisor then gives g(α²W)=g(W): the scalar multiple is determined to be one at W=0. Hence the intrinsic automorphism extends by (W,Z)↦(α²W,αZ).

This is enough to correct an ambient equivalence of unparametrized images into the exact required parametrization. In particular, one does not silently identify an arbitrary abstract automorphism with the identity, and one does not require the chosen ambient equivalences to vary algebraically or continuously with x.

## 4. Generic uniqueness and full-fiber scalar descent

Assume the generic polynomial g is squarefree and put

    D(x,W) = (F_W G_Z − F_Z G_W)(x,W,0).

On uniqueness fibers, representatives F,G and the corrected ambient automorphism agree modulo the defining equation. The equation and its gradient vanish at every singular point, so their derivative matrices agree there. Therefore D has the same nonzero value at all singular points of each such fiber.

Divide D by g over C(x). Away from the finitely many denominator, leading-coefficient, and discriminant exceptions, the remainder has degree below d and takes a common value at d distinct roots. It must be constant in W. Since this holds on an infinite open subset, every positive-degree coefficient of the rational remainder vanishes identically. Thus D≡λ modulo g for λ in C(x).

Write λ=U/V in coprime form. The unit constant coefficient makes g primitive over C[x], so Gauss's lemma upgrades divisibility of VD−U to C[x,W]. A zero ξ of V would make nonconstant g(ξ,W) divide −U(ξ)≠0. The hypothesis that every fiber retains at least one collision is exactly what excludes that possibility. Hence λ is polynomial.

At any exceptional fiber singular point, whole-algebra generation alone shows that P−P(q), Q−Q(q) generate the maximal ideal. The original classes W−r,Z are a basis of its two-dimensional cotangent space, even for a tangency or a multiple root of g. The derivative matrix of F,G is the change-of-generators matrix and has nonzero determinant. Consequently D(ξ,r)≠0 without any exceptional-fiber ambient automorphism. Since g divides D−λ over C[x,W], λ(ξ)≠0 on every base fiber. A complex polynomial without zeros is a nonzero constant.

This verifies the amended generic-only hypothesis. A generic exponent 2d+1 that is prime is sufficient; the leading coefficient may be nonconstant, and an exceptional fiber may have a composite leading exponent. The exact control g=xW⁵+W⁴+2 checks such an 11-to-9 drop. This is not a claim of uniqueness for the degree-four exceptional curve.

The distinction between a root-cover unit and a common scalar is essential. On W²−xW+1=0, W is a unit with inverse x−W, but it is not scalar in the quotient. In contrast, the proof first makes D a common value across all generic roots, then descends that value.

## 5. Repeated roots, moving endpoints, and coefficient contradiction

If g has a repeated generic root r, differentiation of g(x,r)=0 gives g_x(x,r)=0 because g_W(x,r)=0. Characteristic zero permits extending the base derivation to the algebraic root field. With t²=r≠0, Z, Z_x, and Z_y all vanish at ±t. The resulting Jacobian values are opposite, while a base-only Jacobian takes equal values there. This already forces that Jacobian to vanish and needs only membership, not generation.

For the squarefree case, the period identity is taken on a small analytic open set contained in the uniqueness, full-degree, and simple-root loci. Equality is pointwise for every base value in this set; no differentiability of the ambient automorphism choices is required. Both sides of the identity are analytic through their explicit polynomial integrands and chosen analytic endpoint functions.

For moving endpoints u,v with equal target values, the complete boundary correction after integration by parts is

    P(v) d[Q(x,v(x))]/dx − P(u) d[Q(x,u(x))]/dx.

This is zero because the two P values agree and the two total Q derivatives agree. Dropping the endpoint-velocity terms without this argument would be an error; the note keeps them through the total-derivative formulation. Independent tests use nonsymmetric moving endpoints with unequal speeds and nonzero common P and Q values.

The standard pair's Jacobian is −2 sum a_i' y^(2i+2). Integrating and dividing by the endpoint separation 2t gives the factors 1/(2i+3) in H. From H(r)=κ at every generic root, squarefreeness gives g|(H−κ). The quotient has degree at most one. Constant and linear coefficients give B=−κ/c and A=κa₁/c². Its leading coefficient gives

    a_d' = (2d+3) κ a₁ a_d / c².

Here κ is polynomial because the common scalar was proved constant. If a₁ is nonzero, the right side is nonzero of degree at least deg a_d, incompatible with a derivative, including the constant-a_d case. If a₁=0, A=0 and B≠0; successive coefficients force a₂,a₃,... to vanish. The d=1 endpoint case is separately covered. There is no division by a₁ or by a derivative, so no hidden degeneracy has been discarded.

## 6. Recognition and the actual-source consequence

Suppose S contains A₀=R[y²]. The even part of every element already lies in A₀⊂S; subtracting it shows S=A₀⊕yI for an ideal I of A₀. Because B is generated as an A₀-module by 1,y, testing fB⊂S reduces exactly to f,yf∈S. For f=e+yo, these conditions say e,o∈I. Thus the conductor is IB=I⊕yI.

If IB=CB is nonzero principal, parity preserves that principal ideal. Its generator satisfies C(x,−y)=±C(x,y), since the only units of B are nonzero constants. An odd generator C=yh would have h∈I, hence h∈yhB, impossible in this domain. Thus C is even. Contracting (g)B to A₀ gives (g)A₀ by parity decomposition, and therefore I=(g). This proves the entire-algebra recognition, not just a match of conductor ideals.

The final Keller condition at y=0 is valid: every element's y-derivative is divisible by g(x,0); a zero of that polynomial would make both P_y and Q_y zero at one source point. Thus g(x,0) must be a nonzero constant. The monic quadratic reduction uses a global triangular source shear of determinant one, and preserves the generic y-degree of the conductor. A nonunit quadratic leading coefficient would not permit this reduction.

For Corollary D, the cited frozen conductor packet supplies a nonzero principal conductor in its explicit monic, primitive actual-source gauge. After the permitted shear, Proposition C recognizes E_g. Degrees two, four, and six correspond to d=1,2,3 and are excluded by Theorem B. Degree zero gives S=B and an injective source line. Thus eight is a valid conditional lower bound. The added arithmetic consequence is also immediate: if the positive even conductor degree is m=2d, then m+1 must be composite, since a prime value would invoke Theorem B. This still assumes the extra monic quadratic element. The audit checked this application and its hypothesis preservation; it is not a replacement independent reproof of every theorem in the earlier conductor packet.

The same-conductor hostile is useful here. Take C(y)=(y²−1)(y²−4). The equalizer that glues −1 to 1 and −2 to 2 contains y². The equalizer that glues 1 to 2 and −1 to −2 has the same conductor C(y)C[y], but a monic quadratic y²+by+d would require simultaneously b=−3 and b=3. It has no such quadratic. These abstract equalizer rings are not asserted to be two-generated Keller source algebras. They demonstrate exactly why conductor support and multiplicity alone cannot establish the missing recognition hypothesis.

## 7. Exact controls and reproducibility

The independent `independent_checks.py` imports no producer code and uses explicit exceptions rather than Python assert. Its 396 exact checks over 43 declared controls cover:

- polynomial-ideal membership and even/odd conductor tests over Q[x,W], rather than only fraction-field divisibility;
- global inverse witnesses that survive specialization, including common-zero and degree-drop fibers;
- nonprincipal-conductor and same-conductor/different-pairing boundaries;
- Gauss divisibility with a nonunit leading coefficient, a failed nonprimitive upgrade, and the empty-fiber rational-pole counter-control;
- common scalar versus a nonconstant root-cover unit;
- representative changes by the defining equation on ordinary and exceptional fibers;
- nonsymmetric moving-endpoint cancellation and a deliberately nonzero missing-boundary trap;
- coefficient recurrences for d=1 through 8;
- nonlinear fixed-cylinder target automorphisms with signs and nonunit-integer scaling, a zero individual period, and localization-only generation;
- moving repeated roots and the global monic-quadratic source shear.

Normal and optimized independent JSON outputs are byte-identical. The producer checker was separately replayed in both modes: 205 exact checks over 26 declared families, reproducing its stored JSON byte-for-byte.

Reproduce from the workspace root:

    python3 session13-cylinder-review/independent_checks.py --output /tmp/cylinder-independent.json
    python3 -O session13-cylinder-review/independent_checks.py --output /tmp/cylinder-independent-O.json
    cmp /tmp/cylinder-independent.json /tmp/cylinder-independent-O.json
    python3 session13-cylinder-continuation/verify_cylinder_obstruction.py --output /tmp/cylinder-producer.json

These are FINITE-EXACT controls of identities and failure boundaries. They do not prove the cited all-embedding theorem, classify arbitrary actual source algebras, or independently establish the universal statements by testing examples.

## 8. Review boundary

The information transfer is:

    entire algebra generation
      + generic single-curve embedding uniqueness
      → common generic cotangent scalar
      + all-fiber cotangent invertibility and retained collisions
      → one constant scalar
      → moving-period derivative identity
      → polynomial coefficient contradiction.

The principal-conductor recognition theorem supplies a route into this chain only when a monic quadratic normalization coordinate is already present. Losing entire generation, replacing the scalar by an arbitrary unit on the root cover, ignoring exceptional fibers, or keeping only conductor multiplicities would break the argument. The revised note preserves each of these distinctions.

Recommended publication status: **PROVED from the explicitly CITED classical inputs; independently reviewed**, retaining every stated algebra, coordinate, and conditional-scope hypothesis. No priority claim or general Jacobian-conjecture resolution is warranted.

Final freeze: the two mathematical notes changed only their status lines after acceptance. The README additionally gained a review link and the verified independent-check counts. `promotion_check.json` verifies that undoing exactly these presentation changes recovers every reviewed source content hash after status-line normalization. `source_hashes.json` records the final promoted files.
