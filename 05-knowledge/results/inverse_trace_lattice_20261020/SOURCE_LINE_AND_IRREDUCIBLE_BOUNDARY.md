# Source-line specialization and a sharper irreducible-boundary budget

20 October 2026. **Separate addendum; independent review pending.** The finite-selector/Euler gate and local trace packet remain frozen. This note records two further constraints and labels the additional literature input explicitly.

## 1. Setup and conclusions

Retain the hypotheses and notation of `FINITE_SELECTOR_EULER_GATE.md`: an actual complex planar Keller map F=(P,Q), generic degree d≥3, primitive source coordinate x, P monic in its companion y, global regularity of Tr(x^k) for 1≤k≤d−2, and

    R(u,v,Z)=H Z^d+Hc₁Z^(d−1)+…+Hc_(d−2)Z²+A Z+B.

The graph map G=(P,Q,x) is finite and surjective onto W=V(R). For every a∈C put

    R_a(u,v)=R(u,v,a),    D_a=F({x=a}).

Then:

1. R_a is nonconstant and, up to a nonzero scalar, an irreducible polynomial.
2. The map y↦(P(a,y),Q(a,y)) is the finite birational normalization A¹→D_a and is an immersion everywhere. Thus every singularity of D_a has at least two branches; unibranch cusps are excluded.
3. D_a avoids the omitted set O=V(H,A), and on the nonproperness curve C=V(H),

       D_a∩C = {p∈C\O : −B(p)/A(p)=a}.

4. Using the one-line injectivity theorem cited below, d>1 forces every D_a to have a self-intersection. Therefore χ_c(D_a)≤0; all of these singular points lie outside C.
5. Separately, if C itself is irreducible, Jelonek's versioned theorem below gives χ_c(C)≤0. The frozen Euler gate therefore strengthens to

       d−1 ≤ #O = #V(h,A) ≤ deg(h) deg(A),

   where h=rad(H). In the degree-six sector, at least five distinct denominator/numerator intersections are required.

Statements 1–3 follow from the actual polynomial source and étaleness. Statement 4 and the exclusion of χ_c(C)=1 in statement 5 use the explicit external theorems in Section 5. Neither conclusion is inferred from a positive Hermitian trace form.

## 2. Finite parametrization, normalization, and immersion

Finite base change of G by the plane Z=a gives a finite surjective morphism

    {x=a} ≅ A¹ → V(R_a).

Set-theoretically, this proves V(R_a)=D_a. The image is nonempty, closed, and irreducible. It is a curve: P(a,y) is a nonconstant monic polynomial in y. Hence R_a is neither zero nor a unit, and its reduced support is irreducible. At this stage alone one has only R_a=c b_a^e for an irreducible b_a and e≥1.

Write f_a(y)=(P(a,y),Q(a,y)). Its derivative is never zero: if both P_y(a,y) and Q_y(a,y) vanished, the Jacobian of F would vanish. The finite map f_a factors through the normalization ν:˜D_a→D_a. The normalization is A¹. Here is a direct curve argument, avoiding an unstated injectivity assumption:

- The finite dominant map A¹→˜D_a extends to a nonconstant map P¹→the smooth projective completion of ˜D_a. Riemann–Hurwitz forces that completion to have genus zero.
- Because A¹→˜D_a is finite, the points over the complement of ˜D_a all lie at ∞ in P¹. Surjectivity on projective completions then makes that complement exactly one point.
- Thus ˜D_a≅A¹, and the lifted map is a polynomial q:A¹→A¹. Its coordinate factorization is f_a=ν∘q. If deg(q)>1, q′ has a complex zero; the chain rule would make both components of f_a′ zero there. Hence deg(q)=1.

Consequently f_a is finite birational and gives the normalization itself. Its nonzero derivative makes each individual normalization branch nonsingular. If a point had only one normalization preimage, the finite normalization map would be locally one branch, and the analytic inverse-function theorem on a nonzero component of f_a′ would make that point smooth. Thus a singular point has at least two normalization preimages. Distinct branches may be tangent; no transverse-node assumption is made.

## 3. Why the specialized polynomial has multiplicity one

It remains to justify e=1, which does not follow merely from irreducibility of the image.

Choose a general smooth point p of D_a where the finite birational parametrization has exactly one preimage s. In W over p, select w=(p,a). A source preimage of w must lie on {x=a}, so G^(-1)(w)={s}.

Since G is finite, one may choose sufficiently small analytic neighborhoods of w and s such that the full inverse image of the former lies in the latter. Since F is locally biholomorphic at s, the image of G there is the smooth graph

    Z=x(F^(-1)(u,v)).

Surjectivity and the single-point fiber imply this is all of the germ W_w. The polynomial R is irreducible over C[u,v,Z], hence defines W with its reduced structure; its local defining equation is therefore a unit times the equation of that smooth graph.

Restrict to Z=a. In the source coordinates transported by F, the resulting equation is a unit times x−a, which has a simple zero along the source line. Therefore R_a vanishes to order one at the general point of its only irreducible component. Thus e=1, proving statement 1.

This argument needs both finite G and the actual Keller local inverse. Finite surjectivity alone permits nonreduced hyperplane sections, while a formal root polynomial alone need not be the finite image of the source.

## 4. Selector fibers and strong coefficient controls

On H=0 the specialized polynomial is

    R_a = aA+B.

At O, A=0 but B≠0 by the frozen gate. Hence D_a∩O=∅. Away from O on C, A≠0, so aA+B=0 exactly when a=−B/A. This proves statement 3. Moreover every target point of C has at most one source preimage. Therefore the normalization of D_a has at most one point over each point in D_a∩C. The immersion argument makes D_a smooth there. Consequently all singularities of every D_a lie outside C.

Already a single specialization gives a useful actual-source test: B=R_0 must be nonconstant and irreducible, with a finite immersed A¹ normalization. An arbitrary coefficient choice B=1 cannot come from the specified source, regardless of its formal Newton polygon.

For example, the surface

    u Z^d+v Z+1=0

has invertible Z, since Z(uZ^(d−1)+v)=−1. Indeed it is isomorphic to A¹×G_m via (u,Z), with v=−uZ^(d−1)−Z^(-1). It cannot be the finite graph of a polynomial source A² with Z=x: the line x=0 would have a nonempty image, while R_0=1 has empty zero set. This is an exact hostile control for treating a locally admissible inverse polynomial as a polynomial Keller realization.

The uniform condition for every a is stronger: each R_a must be irreducible, and its singularities can only come from genuine collisions of distinct parameters. A candidate atlas with a reducible specialization, a constant specialization, or a unibranch singularity in a specialization is excluded before any global properness conclusion is attempted.

## 5. Versioned external inputs and the boundary inequality

The primary source checked for this note is Z. Jelonek, *A note on the Jacobian Conjecture*, arXiv:2011.03472v3 (1 February 2021), [versioned full text](https://arxiv.org/html/2011.03472v3). This citation is to that precise primary preprint version; no claim is made that a fresh publisher-PDF comparison was performed.

Its Theorem 1.1 states that the nonproperness set of a planar complex Keller map “cannot be a curve without self-intersections.” The proof explicitly uses polynomial parametrizations of its irreducible components and their A¹ normalizations. This excludes a unibranch cuspidal curve just as it excludes a smooth A¹: merely knowing a curve is not algebraically isomorphic to A¹ would not suffice here.

For irreducible C with A¹ normalization, let b_p be its number of normalization preimages above each singular point. Then

    χ_c(C)=1−Σ_p(b_p−1).

The stated theorem forces at least one b_p≥2, so χ_c(C)≤0. Combining this with the frozen gate yields

    #V(h,A)=(d−1)Σ_p(b_p−1)≥d−1.

The same primary source records Gwoździewicz's one-line injectivity theorem as Theorem 2.7: a planar Keller map injective on one affine source line is an automorphism. Applying it to {x=a}, d>1 forces the normalization map f_a above to identify two parameters. Therefore χ_c(D_a)≤0 as well. All these forced collisions occur outside C by Section 4. These cited statements are not claims of new classification results.

## 6. Next concrete target and scope

For a proposed even-degree-six low-trace atlas, the combined tests now require:

- h and A meet in at least five distinct points if the nonproperness curve is irreducible;
- B is nonconstant irreducible and has no zero at any of those points;
- every R_a is an irreducible curve with immersed A¹ normalization and a parameter collision outside C;
- each boundary component retains the local rank-four trace lattice's mandatory repair index 2a_D.

The new target is to derive an intersection upper bound from the all-a source-line specialization condition or from the source compactification. An upper bound below five, a forbidden specialization, or failure of B to be a unit at a required intersection would exclude that concrete degree-six sector. The all-a collision condition does not imply that a generic R_a is smooth: identifications in the finite normalization of W can persist across every parameter a, and their geometry at infinity has not been bounded. No such universal upper bound has been proved here, and no claim is made that degree six is the first globally open function-field degree.
