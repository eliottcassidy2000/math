# All quadratic-semigroup Keller sources are excluded

20 October 2026. **PROVED; independently audited PASS in the stated scope.** The proposal author originated the all-m mechanism. An independent reviewer reconstructed and audited the proof, supplied the shorter filtration induction and hostile controls, and checked the primary sources. The proposal author then independently reviewed the complete written proof and freshly checked the two principal cited sources. No full planar Jacobian-conjecture claim or historical priority claim is made.

## Exact theorem

Let R=C[x], K=C(x), and P,Q in R[y]. Assume P has positive y-degree and nonzero constant leading coefficient (in particular P may be monic), and J(P,Q)=P_x Q_y-P_y Q_x is a nonzero complex constant. Then

    Gamma(K[P,Q]) != <2,m> for every odd integer m>=3.

No separate base saturation, special-fiber constancy, global quadratic generator, conductor monicity, prime-m assumption, or unique-embedding hypothesis is needed.

The useful interpretation is that a proper primitive generic source algebra containing a quadratic is impossible in this monic Keller setting. This is a restricted source-algebra result, not an assertion that every possible Keller source contains a quadratic.

## 1. Constant-leading reduction is global

Suppose p has y-degree n>=1 and leading coefficient a in C*, while q has positive y-degree d and leading coefficient b(x). Then

    [y^(n+d-1)] J(p,q) = -n a b'(x).

Since n+d-1>=1, a constant Jacobian forces b'=0 and b in C*. The degree of p_x is at most n-1, so the other Jacobian product cannot contribute at the displayed degree.

Whenever one positive y-degree divides the other, subtract the unique constant multiple of the appropriate power of the lower-degree coordinate from the higher-degree coordinate. This is a polynomial target automorphism over C, preserves R[p,q] and K[p,q], and preserves the nonzero constant Jacobian up to a constant after swaps or scalings. It strictly reduces the sum of degrees. The unchanged coordinate has constant leading coefficient; applying the same coefficient calculation forces the new positive-degree coordinate to do so too.

A zero-degree coordinate cannot be an escape: if r=r(x), then J(r,q)=r'(x)q_y is a unit. Hence r is affine with nonzero slope and q is linear in y with constant nonzero coefficient. This gives K[r,q]=K[y], whose semigroup is N, not <2,m>. A constant coordinate would make the Jacobian zero. Thus reduction ends at a pair p,q of positive degrees neither of which divides the other. Scale both by constants to make them monic. All these operations are over C, not merely K.

## 2. Exact classical input and applicability

For a characteristic-zero plane curve with one place at infinity, the actual two coordinate pole orders start a characteristic degree sequence generating its semigroup. For a nonprincipal pair with semigroup conductor c,

    delta_0,delta_1 <= c+1;
    delta_i <= c-1 for i>=2.

Every primitive semigroup generator occurs in this sequence. These are **CITED** properties, not a newly proved classification. Direct primary-source checks are recorded in PRIMARY_SOURCES.md.

They apply over an algebraic closure of K. Indeed Gamma(A)=<2,m> supplies h,r in A=K[p,q] of degrees 2,m. The intermediate-field degree [K(y):Frac(A)] divides both 2 and m by the tower law, and is one. Since y is integral over K[h], K[y] is the finite normalization of A. The only point of the normalization where a polynomial in y can have a pole is infinity, and h has a pole there. Hence the curve has exactly one place at infinity.

Constant extension does not change the degree semigroup. For each degree bound, the defining finite linear systems for polynomial-coefficient cancellation are over K; extending scalars does not change their ranks. Equivalently, the monic degree-filtration basis of A and its missing monomial complement remain a basis after scalar extension. The same coprime-degree argument proves geometric primitivity. No source-coordinate coefficient is descended from an algebraic-closure normal form.

## 3. Terminal degrees: one is m, the other is even and smaller

The conductor of <2,m> is c=m-1: the gaps are exactly 1,3,...,m-2. The element m is primitive, since all positive semigroup elements below it are even. It must occur in the characteristic sequence. It cannot occur at index i>=2, because the later-entry bound is m-2. Therefore m is one of the actual terminal coordinate degrees.

The other degree is at most m. It cannot equal m, because equal degrees are principal and would have been reduced. Every positive semigroup element below m is even. Thus, after a possible swap,

    deg_y p = n=2k < m = deg_y q.

This argument includes composite m and terminal pairs with nontrivial degree gcd. It does not assert that p,q themselves have coprime degrees or that the characteristic sequence has only two entries.

**Important hostile:** two even coordinate degrees alone do not prevent primitive field generation. The pair t^4,t^6+t has degrees 4,6, while its cancellation (t^6+t)^2-(t^4)^3=2t^7+t^2 has degree 7; the field degree divides both 4 and 7 and is one. It is the primitive-m occurrence plus the later-entry bound, not parity of the initial degrees alone, that excludes two-even terminals here.

## 4. The rational quadratic center is polynomial

Choose a monic quadratic h=y^2+u y+v in A with u,v in K. Put s=u/2 and t=y+s. Subtracting the K-constant v-u^2/4 gives t^2 in A.

Every element a of A of degree less than m is in K[t^2]. To see this without any full normal-form classification, its positive degree is even; subtract its leading coefficient times the corresponding power of t^2. Repeat at strictly smaller degree. All remaining degrees below m are again even, so the process ends at a constant.

In particular the actual terminal monic p is a monic degree-k polynomial in t^2. Its y^(n-1) coefficient is exactly n s. All other powers of t^2 have degree at most n-2 and cannot affect that coefficient. Since p belongs to R[y],

    s = [y^(n-1)]p / n belongs to R.

Thus y=t-s(x) is a genuine polynomial source shear, of Jacobian one. In the new coordinates pbar(x,t)=p(x,t-s(x)) belongs to R[t] and is even in t: its equality to a polynomial in t^2 over K forces all its odd coefficients to vanish. Its even coefficients are automatically in R. This uses no assertion that the selected generic quadratic itself already lies in R[p,q].

## 5. A genuine affine line after shear

Write C(x)=pbar(x,0) and qbar(x,t)=q(x,t-s(x)). The source shear preserves the Jacobian, so

    j=J(pbar,qbar)|_(t=0)=C'(x) qbar_t(x,0).

Both factors belong to C[x]. Their product is a nonzero constant, so each is a nonzero constant. In particular C(x)=a x+b with a in C*. The restriction of (pbar,qbar) to the actual affine line t=0 is therefore injective, because its first component is ax+b.

**CITED:** Gwozdziewicz, *Injectivity on one line*, Theorem 1.1, states that a polynomial map of the plane over an algebraically closed characteristic-zero field, with nonzero constant Jacobian, is a polynomial automorphism if its restriction to one affine line is injective. Apply this theorem to (pbar,qbar), not directly to the possibly nonlinear graph y=-s(x) in the original coordinates.

It follows that C[pbar,qbar]=C[x,t], and hence K[pbar,qbar]=K[t]. Translation preserves degrees in the normalization variable, while all previous target reductions preserved A. Its semigroup would therefore be N, contradicting <2,m>. This proves the audited statement.

## 6. Hostile examples and scope checks

- The saturated family R[y^3,xy^2+y] has generic semigroup <2,3> and a monic first generator, but J=-3y^4, not a nonzero constant. Its generic quadratic center is y+1/(2x), so monicity of the initial first coordinate alone cannot globalize the center. The constant-Jacobian leading-coefficient gate is load-bearing: the second leading coefficient here is x.
- The stronger saturated hostile R[y^4,xy^3+y^2+y] has generic semigroup <3,4>, no global monic cubic, immersed primitive noninjective slices, and J=-4y^6. It is outside both the present semigroup scope and Keller hypothesis. The new argument does not retract it.
- The moving quadratic family p=y^2, q=y(y^4-s(x)y^2+c), c!=0, is immersed on every slice but has J=2s'(x)y^4. Its even-coordinate critical line has C'(x)=0, precisely violating the unit identity needed in Section 5.
- The same first-order nodal lift used by the constant-leading packet has J=1+x(4-16y). It agrees with a Keller identity on x=0 but is not a Keller map globally. Formal all-orders lifts with unbounded y-degree do not supply polynomial counterexamples, and the terminal degree reduction is not a process on such unbounded series.
- The proof makes no statement about arbitrary higher-multiplicity semigroups, arbitrary orders without an actual Keller pair, or nonconstant leading coefficients. It imports no unique-embedding theorem and does not require m prime.

## 7. Review, computation, and scope

The proof above was reconstructed from the proposal and independently checked against the original one-line theorem, original characteristic-coordinate paragraph, and author-hosted characteristic-degree bounds. The normal-form descent was replaced by the direct degree-filtration induction of Section 4, reducing dependence on predecessor classifications.

`verify_quadratic_semigroup.py` verifies the universal coefficient mechanisms in finite symbolic universes, the two-even primitive hostile, the saturated non-Keller hostiles, the moving-family Jacobian, and a nonlinear-shear positive control. These computations corroborate the displayed identities; the all-m conclusion is the proof in Sections 1-5, not a finite test. All 111 controls pass normally and with Python -O, producing identical output.

The independent review verdict is **PASS**, with two essential safeguards: use primitive m plus the later-degree bound to eliminate even-even terminals, and apply the one-line theorem only after the polynomial source shear has converted the graph into an actual affine line.

## 8. Equivalent quadratic-source formulation

Suppose the generic source algebra A is primitive and proper in K[y], and contains a quadratic. Its degree semigroup is numerical: finite normalization K[y] gives finite codimension, or this follows from a common nonzero conductor element. Since 2 occurs, every nonnegative even degree occurs. Let m be the least odd degree. Properness excludes m=1, since a degree-one element generates K[y]. Then m>=3 and the degree semigroup is exactly <2,m>. Thus the theorem excludes every such proper quadratic source algebra, not just prime m or bounded conductor.

This reformulation does not claim that an arbitrary possible noninvertible Keller source contains a quadratic. The remaining higher-multiplicity semigroup sectors require separate arguments.
