# Constant-leading Keller reduction repairs cubic/quartic globalization

20 October 2026. **PROVED; independently reviewed PASS in the stated monic/semigroup scope.**

See [the independent audit](review/AUDIT.md) for the proof, source checks and a separate direct terminal-coefficient exclusion. This packet does not alter the frozen predecessor proofs.

## Exact theorem and consequence

Put R=C[x] and K=C(x). Let P,Q belong to R[y], with P monic in y of positive degree, and suppose

    J(P,Q)=P_x Q_y-P_y Q_x ∈ C*.

If the degree semigroup of A=K[P,Q] is <3,4>, then there exist monic polynomials H,L in R[y], of respective y-degrees three and four, such that

    R[P,Q]=R[H,L].                                      (T)

In fact the proof does not need a separate base-saturation hypothesis. It uses the frozen generic reduced-pair classification, the constant leading coefficient forced by the Jacobian identity, and an explicit polynomial recovery formula in the exceptional (4,6) class.

Combined with the independently audited moving-monic obstruction, (T) implies that no such P,Q exist. Combined further with the audited low-delta actual-source theorem, a noninvertible primitive monic Keller source has generic source-conductor degree at least eight, not merely six. The independent direct coefficient proof in Section 7 also excludes both terminal branches without invoking the moving-period obstruction. The conductor degree is not the field degree of the Keller map. The planar Jacobian conjecture remains open. No historical priority claim is made for these reductions or the scoped conductor bound.

## Inheritance and loss ledger

Closest proved mechanism: the moving monic cubic/quartic obstruction, repository commit 129154b04024f1246570592ac8df4db249712bc4, together with its generic reduced-pair classification.

Canonical hostile: R[y^4,xy^3+y^2+y] is saturated, has a principal conductor, immersed primitive noninjective slices, and generic semigroup <3,4>, but has no global monic cubic. Its displayed Jacobian is -4y^6, and its cubic generator has nonconstant leading coefficient. The new gate excludes precisely that lost leading-coordinate phenomenon for an actual Keller generating pair.

Corrected near miss: a first-jet Jacobian identity at one fiber gives no contradiction. Section 8 retains an explicit first-jet lift and an all-orders formal lift for every immersed polynomial curve. The formal construction is not required to keep a uniform finite y-degree bound; the explicit nodal completion below does not.

Least-used sidecar: the leading y-coefficient of the *entire generating pair*, not the leading coefficient of a generic chosen cubic.

Concept board: top Jacobian coefficient; constant-coefficient target reduction; generic degree filtration; exceptional embedding label; polynomial source translation; unit recovery coefficient; local formal lifting versus whole-polynomial identity.

Anchor: eliminate the degree-drop mechanism under the actual Keller identity. Niche: distinguish a generic normal form from its integral coordinates. Wildcard: all-orders formal Keller lifting explains why arbitrary finite jets are too weak without a degree bound.

## 1. Top-coefficient lemma

Let p,q∈R[y], with p of degree n≥1 and leading coefficient a∈C*, and q of positive degree m and leading coefficient b(x). Then

    [y^(n+m-1)]J(p,q)=-n a b'(x).                      (1)

Indeed deg_y p_x≤n-1, so deg_y(p_x q_y)≤n+m-2; the only term at the displayed degree is the leading term of -p_y q_x. If J(p,q)∈C*, its positive-degree coefficient in (1) vanishes. Characteristic zero gives b∈C*, since it is a nonzero leading coefficient.

Thus an entire Keller generating pair with one constant leading coefficient has two constant leading coefficients.

## 2. Degree reduction takes place over C, not merely over K

Start with the given P,Q. Their y-degrees are positive: if Q has degree zero in y, then A=K[P], whose degree semigroup is a cyclic semigroup nN, not <3,4>. By Section 1 both leading coefficients are nonzero complex constants.

Whenever one positive degree divides the other, subtract the appropriate constant multiple of a power of the smaller-degree coordinate from the larger-degree coordinate. This is a target triangular automorphism over C, and preserves the Jacobian up to a nonzero constant after possible swaps and scalings. It strictly reduces the sum of y-degrees. The unchanged coordinate still has constant leading coefficient, so Section 1 forces the new nonconstant coordinate's leading coefficient to be constant too.

A coordinate of degree zero cannot occur at any stage, since the K-algebra and its noncyclic semigroup have not changed. More concretely, if a terminal coordinate is r(x), then J(r,q)=r'(x)q_y∈C* forces both factors to be nonzero constants. Thus r is linear in x and q is linear in y with constant nonzero leading coefficient, so R[r,q]=R[y], again incompatible with semigroup <3,4>. The process therefore terminates with a nonprincipal generating pair p,q, meaning neither y-degree divides the other. Rescale both by nonzero constants to be monic. All reductions preserve the entire R-algebra exactly and preserve nonzero constant Jacobian. No x-dependent target coefficient has entered this reduction.

## 3. The only terminal degrees are (3,4) and (4,6)

This is the precise frozen classical input. For a characteristic-zero plane polynomial curve with one place at infinity, a nonprincipal characteristic degree sequence of semigroup conductor C has its first two degrees at most C+1. The characteristic gcd sequence stops once its gcd is one, and its entries include every primitive generator of the degree semigroup.

Here the semigroup is <3,4>, with conductor six. The possible positive degrees at most seven are 3,4,6,7. The nondividing pairs are (3,4), (4,6), (3,7), (4,7), (6,7). Any pair involving seven has gcd one and therefore no later characteristic entry; the resulting semigroup fails to contain both three and four. The only terminal possibilities are consequently

    (deg_y p,deg_y q)=(3,4) or (4,6),                  (2)

up to swapping.

The bound and characteristic-sequence properties are CITED, not newly proved here. They are the identical input already used and independently audited in the fixed and moving cubic/quartic packets. The remaining calculations below are explicit. Passage from K to its algebraic closure preserves the degree filtration: each filtered membership/rank calculation is a linear system over K. No source coordinate or normal-form coefficient is assumed to descend merely from an algebraic-closure classification.

In the (3,4) case, the entire reduced pair itself proves (T).

## 4. Exact exceptional-class extraction over K

Suppose the terminal monic pair has degrees (4,6). The degree semigroup supplies monic elements of degrees three and four in A. Their generated subalgebra has K-basis 1,k,k² over the cubic polynomial ring, with leading degrees 0,4,8 in distinct residues modulo three. Hence its degree semigroup is <3,4> and its codimension in K[y] is three, the same as A. The inclusion is consequently equality, so these elements generate A over K. Complete the cubic's degree-three translation over K and adjust the quartic to obtain

    t=y+s,   h=t^3+a t,   k=t^4+b t^2+c t,
    A=K[h,k],   s,a,b,c∈K.                            (3)

At this point s and the other coefficients are genuinely allowed to be rational. Polynomiality will be proved, not assumed.

The filtered spaces are

    A_≤4=K⊕K h⊕K k,
    A_≤6=K⊕K h⊕K k⊕K h^2.

Because p and q are monic of degrees four and six, respectively, they have unique forms

    p=k+βh+d,
    q=h^2+ηk+γ₀h+e₀.                                (4)

Replacing q by q-ηp-(e₀-ηd) gives h²+γh, where γ=γ₀-ηβ. This is used initially only over K.

For this generic comparison, also replace p by p-d. The normalized pair is then exactly k+βh, h²+γh; both changes are affine over K and preserve the K-algebra. This subtraction is not asserted to preserve the (x,y)-Jacobian. Retain (4) for the original terminal pair. Applied to the normalized pair, the frozen direct degree-separated-basis calculation gives

    β=0,   γ=3c/2,   b=2a/3,   c≠0.                  (5)

For auditability, its logic is recalled. The K[p]-basis 1,q,q²,q³ can be replaced unimodularly by 1,q,r,qr whenever r=q²-A(p)q-B(p). If β≠0, r=q²-p³ has degrees 0,6,11,17 in distinct residues modulo four, giving delta seven. With β=0, subtracting (4a-3b)pq leaves a degree-nine term 2γ-3c; if nonzero the corresponding basis gives delta six. When γ=3c/2, a further k² subtraction leaves a degree-seven term -c(2a-3b)/2; if nonzero the basis gives delta five. If c=0, both candidate coordinates are even in t and cannot have fraction field K(t). These contradictions to delta three prove (5). The frozen packet contains the complete polynomial coefficients and exact identity for the survivor.

Since β=0, equations (4)-(5) now read

    p=k+d,
    q=w+ηk+e,   w=h²+(3c/2)h.                        (6)

No intrinsic automorphism of A, additional source reparametrization, or unspecified ambient equivalence is inserted in (6); these are identities for the same terminal p,q and the same t.

## 5. Every rational normal-form coefficient in (6) is polynomial

Write p=y⁴+p₃(x)y³+⋯. Since k=t⁴+b t²+c t and t=y+s, the y³ coefficient of p=k+d is exactly 4s. Therefore

    s=p₃/4∈R.                                       (7)

The translation t=y+s(x) is a polynomial source shear of Jacobian one, not an arbitrary rational source change. Write p in this t-coordinate. Its t², t and constant coefficients are respectively b,c,d, so all belong to R. Equation (5) gives

    a=3b/2∈R.                                       (8)

Consequently h,k,w belong to R[y]. In the same coordinate, the t⁴ coefficient of w is 2a, while k is monic degree four. Equation (6) therefore gives

    η=[t⁴]q-2a∈R.

Because h,k,w vanish at t=0, the constant coefficient of q in t is e, so e∈R. Thus every operation in (6) now lies over R and

    R[p,q]=R[k,w].                                   (9)

This is the coefficient-level repair of the putative pole at a quadratic degree-drop fiber.

## 6. Immersion makes the recovery coefficient a unit

If c(ξ)=0 for some ξ∈C, evaluate in the source at x=ξ,t=0. From

    k_t=4t³+(4a/3)t+c,
    w_t=(2h+3c/2)(3t²+a)

both derivatives vanish. Equation (6), now known over R, gives p_y=q_y=0 there, contradicting J(p,q)∈C*. Hence c has no complex zero, and

    c∈C*.                                          (10)

The frozen recovery identity is

    w²-k³-2akw+(2a²/3)k²
      +(2a³/27+3c²/4)w-(a⁴/9+ac²)k=(c³/8)h.        (11)

Since c³/8 is a nonzero constant, (11) recovers h in R[k,w]. Conversely w∈R[h,k], so

    R[p,q]=R[k,w]=R[h,k].                            (12)

Both h and k are globally monic in y, of degrees three and four. This proves (T) in the remaining terminal case.

An alternative weaker finishing route is saturation: after (7)-(8), h belongs to R[y]∩S_K and hence to S; a saturated monic cubic/quartic subalgebra then equals S. The direct proof (9)-(12) is preferable because it eliminates that extra dependency and shows exactly why the recovery denominator cannot vanish.

## 7. Independent direct exclusion of both terminal branches

The independent reviewer supplied a second endpoint. It uses the same constant-leading reduction and generic terminal classification but no moving-period argument. The universal Jacobian identities are checked independently in `review/independent_checks.py`.

### 7.1 Terminal degrees (4,6)

Use the polynomial coordinate t and exact pair from (6), where c∈C*. Put A=a'. The coefficient of t⁷ in J(p,q) gives η'=-A. Its t⁵ coefficient then gives d'=2aA/9, and the t³ coefficient gives

    e'=a²A/3+2aηA/9.

After these three substitutions, the entire Jacobian, including its constant coefficient, is

    J(p,q)=-(a'c²/2)t.                              (13)

It cannot be a nonzero constant. This excludes the exceptional terminal case independently of cubic recovery, saturation, conductor monicity, or collision periods.

### 7.2 Terminal degrees (3,4)

Center the monic cubic by the polynomial source shear obtained from its quadratic coefficient. The entire terminal pair has the form

    p=t³+a t+d,
    q=t⁴+ηt³+b₀t²+c₀t+e₀,    a,d,η,b₀,c₀,e₀∈R.

The t⁵ coefficient of J is -3η', so η is constant. The constant target shear q↦q-ηp removes this cubic term. Rename the remaining coefficients to obtain

    p=t³+a t+d,   q=t⁴+b t²+c t+e.                 (14)

The four positive-degree coefficient equations are

    4a'-3b'=0,
    -3c'+4d'=0,
    2ba'-ab'-3e'=0,
    ca'-ac'+2bd'=0.

Thus b=4a/3+B for B∈C, and

    d'=3c'/4,
    e'=(4a/9+2B/3)a',
    ((2a+3B)c)'=0.                                 (15)

The constant coefficient reads

    j=(3c²/8-4a³/27-Ba²/3)'.                       (16)

Write (2a+3B)c=D∈C. If D≠0, both polynomial factors are units and a,c are constant, contrary to j≠0. If D=0 and c is not identically zero, then a=-3B/2 is constant, so (16) would make c² have degree one, impossible for a polynomial. If c=0 identically, (16) would make 4a³/27+Ba²/3 have degree one. For constant a its derivative is zero; for nonconstant a that polynomial has degree 3 deg(a)≥3. This too is impossible.

Hence neither terminal pair is Keller. This independently proves the nonexistence consequence of (T), and combines with the prior low-delta actual-source theorem to give the scoped generic source-conductor lower bound eight.

## 8. Why a single fiber, and even formal jets, do not replace the whole identity

For p(y)=y⁴ and q(y)=y²+y, put

    a(y)=4y²-2y+1,   b(y)=2.

Then a q'-p'b=1. The explicit first-order lift

    P=y⁴+x(4y²-2y+1),
    Q=y²+y+2x

satisfies

    J(P,Q)=1+x(4-16y),

and therefore J(P,Q)|_(x=0)=1, although the special image is the immersed nodal curve C[y⁴,y²+y]. This is not an entire Keller map and is not a conductor-six example.

More generally, if p,q∈C[y] have gcd(p',q')=1, choose A,B∈C[y] with A q'-p'B=1. There exist formal series

    P=Σ_(n≥0) p_n(y)x^n,  Q=Σ_(n≥0)q_n(y)x^n ∈ C[y][[x]],
    p₀=p, q₀=q,

with J(P,Q)=1 exactly. Inductively, after coefficients below n have been fixed, the coefficient of x^(n-1) is

    n(p_n q'-p' q_n)+E_n(y),

where E_n depends only on earlier coefficients. Put R_n=(1_(n=1)-E_n)/n and take p_n=A R_n, q_n=B R_n. This solves every order. Differentiation and products are well-defined x-adically.

There is also a closed-form universal construction. Put C(y)=A B'-A'B and let F∈x C[y][[x]] be the unique formal solution of F+C F²/2=x. Equivalently F is the formal expansion of (sqrt(1+2Cx)-1)/C, understood by its polynomial coefficients also when C=0. Then P=p+AF and Q=q+BF satisfy

    J(P,Q)=F_x(1+CF)=1.

The F_y terms cancel identically. The implicit equation has a unit derivative at F=x=0, so this series exists with polynomial y-coefficients. For the explicit nodal control C=4-16y, its coefficient of x^n has exact y-degree n-1; the corresponding P coefficient has degree n+1. Thus this particular completion has unbounded y-degrees.

This completion is in C[y][[x]], not necessarily C[[x]][y] and certainly not necessarily C[x,y]: y-degrees can grow unboundedly. It need not retain the generic semigroup or any finite monic coordinate degree. Thus the correct obstruction is the whole polynomial generating pair with controlled finite degrees; arbitrary local jets alone cannot produce it.

## 9. Scope and next checks

The degree-drop alternative of the frozen globalization packet remains a real phenomenon for arbitrary immersed saturated orders, but cannot occur in the generic <3,4> sector of an actual monic Keller generating pair. The hostile is not retracted; the missing hypothesis is now supplied.

The first residual generic conductor degree becomes eight. The symmetric numerical-semigroup alternatives there, their plane realizability, and their entire embedding behavior need separate analysis. No claim about them is made here.

The proof imports only the exact bounded characteristic-sequence input and the explicit (4,6) generation calculation already audited in the predecessor. The moving entire-generation obstruction is a separate audited alternative endpoint used only after (T). It does not circularly assume (T), and it is not a dependency of the direct exclusion in Section 7.
