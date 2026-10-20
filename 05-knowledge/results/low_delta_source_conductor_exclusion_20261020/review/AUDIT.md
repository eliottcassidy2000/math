# Independent audit: low delta, degree gaps, and relative quadratic descent

20 October 2026.

## Verdict

**PASS under the exact stated hypotheses and the cited one-line theorem.** No mathematical gap was found in the frozen low-delta packet. The absolute degree-semigroup argument, quadratic classification, both sharp immersed delta-three examples, and actual-source globalization are sound.

Combined with the separately reviewed moving quadratic-source exclusion, the correct consequence is:

> For a hypothetical noninvertible complex Keller map in the specified monic primitive source gauge, the generic conductor degree is at least six. If it is exactly six, the generic source algebra cannot be quadratic; it has a generating presentation of degrees three and four over C(x).

This does **not** assert a lower bound of six on every special fiber, globalize the cubic/quartic presentation, or exclude the remaining generic cubic/quartic sector. It does not establish the Jacobian conjecture.

The reviewer wrote only this review directory. The frozen producer files and shared checkout were not edited. The source manifest verified before acceptance. The computational controls below are finite exact controls, not a substitute for the universal proof.

## 1. Trace duality and polynomiality of the conductor

Put D=k[u], u=f(y), with f monic of degree n. Since k(f,g)=k(y), the minimal polynomial R(u,v) of g over k(u) has degree n in v. Integrality of g over D and normality of D put the coefficients in D. Thus both claimed power bases really are D-bases, including when the displayed degrees of f and g have a common divisor. One must not replace this field-primitivity argument by a degree-gcd heuristic in the general theorem.

For a monic degree-n polynomial H, the functional Tr(h(theta)/H'(theta)) is the coefficient of degree n−1 of the remainder of h. The pairing matrix on 1,theta,...,theta^(n−1) has zero entries on one side of its anti-diagonal and ones on that anti-diagonal. Its determinant is a unit. This proves the two trace-dual formulas over D itself, not only over its fraction field, and establishes double duality for the free lattices in use.

The step proving polynomiality is particularly important and is valid: B-dagger is contained in A-dagger, so the particular element 1/f' belongs to (1/R_v)A. Therefore c=R_v/f' belongs to A. No immersion or prior divisibility hypothesis is needed.

For b in B, the displayed conductor equivalence follows by testing Tr against both factors. Explicitly, bB is contained in A if and only if Tr(b B A-dagger) is contained in D, if and only if b A-dagger is contained in B-dagger. After the monogenic substitutions this is (b/c)A contained in B, equivalent to b/c in B because 1 belongs to A and A is contained in B. The zero element b presents no exception.

The determinant comparison gives Norm(c)=I² with no residual sign: both degree-n monogenic discriminants contribute the same sign. Constant factors in a nonmonic choice of defining equation would matter, but the packet fixes R monic in v. Finally, factorization after extension of constants shows deg_u Norm(h)=deg_y h for every nonzero polynomial h. Smith normal form therefore gives deg c=2δ. The degree-one case has c a nonzero constant and δ=0 as required.

## 2. The trace-to-residue step has the correct factor

For a in A, Tr(a/R_v) lies in k[u]. The identity

    (a/c) dy = (a/R_v) du

is exact, since du=f' dy and c=R_v/f'. There is only one place above u=infinity, namely y=infinity. Thus the residue of the traced differential is the residue of the original differential, with no extra ramification factor.

The factor can be checked directly. In local coordinates v=z^n, a monomial z^m dz becomes (1/n) z^(m−n+1) dv. Trace kills it unless n divides m−n+1; if it survives, the trace contributes n and cancels 1/n. The downstairs differential has a v^(−1) term exactly when m=−1, with the same coefficient. This establishes the required identity term by term. Characteristic zero allows the formal parameter change and division by n. Constant extension is harmless because the relevant point at infinity is rational; vanishing descends.

A polynomial in u times du has residue zero at infinity. Consequently every a in A satisfies Res_infinity(a dy/c)=0. If deg a=deg c−1, that residue is minus the ratio of leading coefficients and is nonzero. This proves the missing degree 2δ−1 for δ>0.

There is no hidden choice of independent local dualizing generators in this argument. Its strength is the specific global differential dy/c obtained from the plane presentation and its polynomial trace. The existence of a principal conductor alone would not justify its residue property; Section 6 below gives an exact hostile control.

## 3. Codimension, gaps, and the delta-three generation argument

For each d, the dimension increment of A intersect k[y]_(degree at most d) is either zero or one, and equals one exactly when d occurs as the degree of an element of A. The filtered quotients inject into B/A and eventually exhaust that finite-dimensional vector space. Thus the number of missing degrees is exactly δ. Equivalently, leading-term reduction gives a basis of B/A indexed by those degrees. Finite codimension supplies cofiniteness; multiplicativity supplies closure under addition.

For a numerical semigroup with largest gap F and δ gaps, the injection s ↦ F−s from the semigroup elements in [0,F] into its gaps gives F+1−δ at most δ. Hence F at most 2δ−1. The residue gap realizes that bound. Equality also proves complementarity of each pair s,F−s, so the claimed symmetry follows without importing the stronger full plane-semigroup theorem.

The low-genus conclusions are now forced:

- δ=0 gives B itself.
- δ=1 has gaps {1}.
- δ=2 has gaps {1,3}; otherwise 1,2,3 would already be three gaps.
- δ=3 has gaps either {1,3,5} or {1,2,5}.

In the last nonquadratic case, choose h and k of degrees three and four. The extension degree [k(y):k(h,k)] divides both three and four, so it is one. The basis 1,k,k² over k[h] has leading degrees respectively 0,4,8, distinct modulo three. Therefore no leading cancellation between its summands is possible, and k[h,k] has exactly the semigroup generated by three and four. Its codimension is three. Inclusion into A with equal codimension proves equality of the entire algebras.

The further normal form follows by translating the cubic to remove its quadratic term, subtracting constants, and then subtracting a multiple of the cubic from the quartic to remove its cubic term. It is a statement about a selected generating presentation. It does not imply an ambient polynomial change carrying an arbitrary original generating pair to this pair, nor a C[x]-defined change when the ground field is C(x).

## 4. Quadratic algebras and the geometric sharp examples

Once t² belongs to A, the even part of any element already belongs to A. Its odd part therefore also belongs to A. The set I of odd coefficients is a nonzero ideal in k[t²], hence principal. This gives the whole algebra A=k[t²,t p(t²)], not only an inclusion or conductor match. The quotient dimension is deg p.

For b=e(t²)+t o(t²), testing b and tb is enough to test bB contained in A, because B is generated by 1,t over k[t²]. It gives precisely divisibility of e and o by p, proving conductor p(t²)B.

Over C the derivative criterion for immersion is exactly p(0) nonzero. A distinct nonzero root of p glues the two parameters ±sqrt(a); its multiplicity is their intersection multiplicity. Thus degree-two p gives precisely two ordinary nodes or one ordinary tacnode in the immersed case. After undoing the translation, both nodal pairs have the same sum. This conclusion does not hold for arbitrary nonplanar equalizers.

Both delta-three examples pass independently:

- For (y³,y⁴+y), the relation is v³−u(u+1)³ and c=(y³+1)². The three distinct roots of y³+1 form one fiber; the three branch slopes are distinct. The point is an ordinary triple point. The semigroup is generated by three and four.
- For (y³−y,y⁴), the relation and conductor in the packet are exact. The conductor (y²−1)(y⁴+1) is squarefree. The three pairs have sums 0,sqrt(2),−sqrt(2), are disjoint, and the determinant of their tangent vectors is 8(y−t), nonzero. These are exactly three ordinary nodes. Different pair sums independently preclude a monic quadratic.

The two sharp examples establish only realizability as immersed polynomial plane curves, not as slices of an actual noninvertible Keller map.

## 5. Specialization, saturation, and globalization

The relevant inherited source hypothesis is the actual complex Keller map, with P monic in y and C(x,P,Q)=C(x,y). The monicity ensures that every P_a still has the same positive degree n.

For clarity, the source proof can be expanded as follows. Each slice is finite onto its image. Its normalization is A¹: its function field is a rational subfield of C(y), and the unique point above infinity leaves a polynomial normalization parameter. Write the normalization-factor map as h:A¹→A¹. Both P_a' and Q_a' are divisible by h'. The Keller identity prevents their common vanishing anywhere, so h' is a nonzero constant and h is linear. Thus every slice is birational and Q_a has minimal-polynomial degree n over C(P_a).

Globally S is C[x,u,v]/(R), with R monic in v. The characteristic polynomial of multiplication by Q on C[x,y] over C[x,P] is R, because generic primitivity makes its minimal degree n. Specialization commutes with the characteristic determinant. The slice minimal degree n just proved forces the specialized characteristic polynomial to be the slice minimal polynomial. Consequently S/(x−a)S embeds in C[y] for every a. This explicitly justifies the quotient-to-image identification; no arbitrary base-change theorem for conductors is needed.

If (x−a)b lies in S, its zero image after specialization lies in (x−a)S by that injection. Cancellation in the ambient polynomial domain gives b in S. Factoring any nonzero polynomial in C[x] proves full base saturation.

Given a generic monic quadratic, its constant term can first be subtracted in C(x). Clearing a common denominator then gives α₀y²+b₀y in S. Dividing the gcd of α₀,b₀ is permitted by saturation and gives αy²+by in S with gcd(α,b)=1. No inversion of the remaining leading coefficient has yet occurred.

If α(a)=0, coprimality gives b(a) nonzero. The specialized image algebra therefore contains y and is C[y]. Hence the original Keller map is injective on the affine source line x=a. The original Gwoździewicz theorem, freshly checked, gives a polynomial automorphism. In the alternative noninvertible case, α has no zero in C and is constant. This is exactly the missing justification for global monicity. An automorphism has C[P,Q]=C[x,y] and supplies a quadratic trivially.

Nothing here assumes constant slice δ, a nonvanishing conductor leading coefficient, flatness of the conductor quotient, or constancy of a rationally chosen shift. The field being C and the one-line theorem are essential at the final step.

## 6. Three sharp failure boundaries

### 6.1 Principal degree-four conductor does not force the residue identity

Let E be the subalgebra gluing 1 to 2 and −1 to −2. Its conductor is

    (y²−1)(y²−4) C[y]

and its codimension is two. To see the conductor necessity, if bB lies in E then b(a)=b(b) and a b(a)=b b(b) at each distinct glued pair; both values must be zero. The converse is immediate.

Nevertheless E has no nonconstant element of degree at most two: the two evaluation-difference conditions are independent on the span of y,y². It contains y³−7y. The residue of

    (y³−7y) dy / ((y²−1)(y²−4))

at infinity is −1, not zero. Thus the plane trace condition is a real additional ingredient. This E cannot be two-generated over the field by the theorem under review. Local nodality and a principal conductor of degree 2δ do not supply the global plane trace condition.

### 6.2 Generic primitivity does not imply saturation

Take S=C[x,xy²,y³]. The generic algebra contains y², but y² is not in S even though xy² is. The presentation has relation u³−x³v². At x=0 the abstract quotient still has a nonzero nilpotent u, while its polynomial image sends u to zero. This is precisely the failure of specialization injectivity. Monicity and slice primitivity in the actual source proof cannot be discarded.

### 6.3 Even saturation plus immersed primitive slices does not suffice without Keller

The sharper example is

    P=y³,  Q=xy²+y,  S=C[x,P,Q].

P is monic. Every slice is immersed, since P_y vanishes only at y=0, where Q_y=1. Every slice is primitive: for x nonzero the degrees are three and two, and at x=0 the second coordinate is y. The characteristic-polynomial argument therefore makes S base-saturated. Its monic relation and conductor are

    R=v³−3xuv−x³u²−u,
    c=x²y²+xy+1.

The generic algebra contains y²+y/x and has δ=1. At x=0 it is C[y] and δ=0. On the generic fiber the conductor roots form the single collision pair, whose sum is −1/x. Every monic quadratic in the algebra must therefore have linear coefficient 1/x, so no global monic quadratic exists.

Here J=−3y⁴. The map has an injective x=0 source line but is not Keller, so the external theorem cannot apply. This isolates the final logical dependency with all preceding hypotheses still present.

## 7. Exact combination with the cylinder exclusion

Use the hypotheses of the audited moving-source packet, not merely a conductor-support analogy. The inherited actual-source theorem supplies the nonzero principal conductor. A global monic quadratic permits the determinant-one triangular source shear that makes it y². The recognition proposition then identifies the **entire** source algebra as

    C[x,y²,y g(x,y²)],

with g(x,0) a nonzero constant. The Keller condition gives this last property by evaluating the y-derivatives at y=0. For a noninvertible map every slice retains a collision by the one-line theorem.

The already reviewed exclusion applies when the generic exponent 2d+1 is prime, where d is the degree in y² of g. Its argument only needs embedding uniqueness on a nonempty generic open set; exceptional fibers are covered by whole-algebra cotangent generation. Thus the nonconstant leading coefficients and potential degree drops allowed in the present globalization do not violate its hypotheses.

Consequently generic conductor degrees two and four give δ at most two, hence a generic quadratic, hence a global quadratic, hence the excluded d=1 or d=2 source sectors. Degree zero also gives an automorphism. Since the conductor degree is even, the remaining noninvertible generic degree is at least six.

At degree six, δ=3. If the generic algebra were quadratic, globalization would invoke the excluded d=3 sector (2d+1=7). Therefore it must have semigroup generated by three and four and be generated by a cubic and quartic over C(x). The normalization translation and chosen generators may have rational x-coefficients. This audit does not globalize them or extend the prime-exponent quadratic uniqueness input to arbitrary cubic/quartic families.

The bound is **generic**. The present argument does not rule out a degree-two or degree-four conductor specialization inside a higher-degree nonquadratic generic family.

## 8. Independent controls, software trap, and reproducibility

The independent checker imports no producer code and uses explicit exceptions, so Python optimization cannot disable the checks. It computes multiplication matrices, characteristic polynomials, actual field traces by modular inversion, norm determinants, conductor quotients, and their degree filtrations for twelve independently specified pairs. They include cusps, tacnodes, mixed singularities, translated quadratics, both sharp delta-three examples, a hidden quadratic in an original degree-(3,6) generating pair, and a degree-(4,6) cancellation example with semigroup generated by 4,6,7.

The finite algebra modulo the conductor is generated to closure and row-reduced in descending degree. Thus the gap checks do not incorrectly infer a degree semigroup merely from the original generator degrees. Additional controls cover local trace/residue coefficients for ramification indices one through eight, the complete low-genus numerical alternatives, all three failure boundaries above, and tangent determinants.

A broader resultant control exposed an implementation sign pitfall in the installed SymPy 1.14.0: its direct resultant(t³,t⁵+1,t) gives −1, although the defining Sylvester determinant and the root product give +1. The tested family (y⁴,y⁶+y) encounters this ordering. The checker therefore first reduces the second divided difference modulo the first, which is monic. This preserves the mathematical resultant exactly and avoids the library's degree-swap sign issue. It is not a defect in the trace/conductor proof, and the frozen producer examples do not meet this odd-degree swap case. A separate exact receipt records the discrepancy.

Normal and optimized outputs must agree byte-for-byte; the producer is also replayed independently in both modes and compared with its frozen receipts. Counts, hash provenance, and comparison results are recorded in REVIEW_SUMMARY.json.

Reproduce from the workspace root:

    python3 session13-low-delta-review/independent_checks.py --output /tmp/low-delta-review.json
    python3 -O session13-low-delta-review/independent_checks.py --output /tmp/low-delta-review-O.json
    cmp /tmp/low-delta-review.json /tmp/low-delta-review-O.json
    python3 session13-low-delta/verify_low_delta.py

Recommended status: **proved within the stated scope, with the one-line theorem cited, and independently reviewed.** No proof repair or scope contraction is required. An optional exposition improvement is to write the monic hypersurface presentation explicitly before the specialization-injectivity assertion.

## 9. Final promotion verification

The producer subsequently added the explicit combined corollary as Section 6.1, pinned to cylinder dependency commit `9605e62f324416ddb1a448c2e923d12e7806166f`, and promoted the status and README. That final wording is accepted. Removing the new outcome item, restoring the original status, and replacing Section 6.1 by the original two dependency-pending paragraphs recovers the originally reviewed main-file hash exactly. Thus the original proof body did not otherwise change. The new corollary has also been checked against Theorem B, Proposition C, and Corollary D read directly from the pinned commit. It preserves the generic, normalization-degree, whole-generation, and non-globalization limitations identified above. The four frozen control/source files remain byte-identical. `promotion_check.json` records the final main and README hashes; `prepromotion_source_hashes.json` preserves the original input provenance, and `source_hashes.json` records the accepted final state.
