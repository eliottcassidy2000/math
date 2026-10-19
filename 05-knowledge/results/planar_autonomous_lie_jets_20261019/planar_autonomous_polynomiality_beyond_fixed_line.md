# Polynomial autonomous Hamiltonian flows beyond a fixed source line

19 October 2026. **Proof assembled from explicitly cited classical inputs; independent review pending.** No historical-novelty or Jacobian-conjecture resolution is claimed. The contribution to the present research program is the exact bridge from a point-local autonomous model to a global coordinate-fibration obstruction, together with explicit examples separating it from arbitrary polynomial symplectic maps.

## Main theorem

Let K be a field of characteristic zero, m=(x,y), and H in K+m^3. Use

    {A,B}=A_x B_y-A_y B_x,       D_H={-,H}.

For every nonzero scalar lambda in K the original m-adic exponential

    F_lambda=(exp(lambda D_H)x, exp(lambda D_H)y)

is defined in K[[x,y]]^2. The following statements are equivalent:

1. F_lambda has polynomial coordinates for one nonzero lambda.
2. D_H is locally nilpotent on K[x,y].
3. Every scalar-time map is a polynomial automorphism, and the full flow is polynomial in (s,x,y).
4. H=f(q), where q is a polynomial coordinate of K[x,y] and f is a one-variable polynomial over K.

Constant H is included: the flow is the identity and any coordinate q works. No assumption that the map fixes a specified line is made.

A polynomial coordinate means that some p in K[x,y] satisfies K[p,q]=K[x,y]. If a={p,q} is its nonzero constant Jacobian, the full flow in these coordinates is exactly

    (p,q) -> (p+s a f'(q),q).

Thus the polynomial autonomous locus is the conjugate-shear locus. Polynomiality at a single nonzero time cannot be an accidental cancellation confined to that time.

## 1. Formal convergence and the global invariant

D_H(m^j) is contained in m^(j+1). Consequently every fixed total-degree coefficient in exp(sD_H)x and exp(sD_H)y receives only finitely many contributions and is a polynomial in s. The inverse at scalar time lambda is the formal map at -lambda. Hamiltonian derivations preserve the Poisson bracket, so

    det JF_lambda=1,     F_lambda(0)=0,     JF_lambda(0)=I.

Also D_H H=0 gives H(F_lambda)=H as a formal identity. If F_lambda is polynomial, that is an identity of polynomials on the entire plane.

This last identity is the global information missing from arbitrary finite jets. By Shpilrain--Yu, Corollary 1.7, a characteristic-zero planar polynomial map with invertible Jacobian matrix which fixes a nonconstant polynomial is an automorphism. For nonconstant H this makes F_lambda a polynomial automorphism. The identity case is harmless: if D_H is nonzero, the lowest-degree nonzero term of F_lambda-id is lambda times the corresponding lowest-degree term of (D_Hx,D_Hy), so F_lambda cannot be the identity.

## 2. Necessity over C: invariant polynomial forces an elementary map

We use the classical plane-automorphism dichotomy: a polynomial automorphism of C^2 is polynomially conjugate either to an elementary triangular map, or to a composition of generalized Hénon maps of degrees at least two.

The Hénon case has no nonconstant invariant polynomial. The exact input needed here is the polynomial degree statement in Bedford--Smillie, Proposition 4.2: for a nonconstant polynomial h and a generalized-Hénon composition g, some h composed with g^n has a single pure y-power as its highest total-degree term. Thereafter its degree multiplies by deg(g)>1 under further iteration. Thus h composed with g=h is impossible. This uses the algebraic conclusion of their proposition, not a claim that all invariant curves or arbitrary maps have already been classified by this argument.

Therefore F_lambda is conjugate to

    T(u,v)=(a u+P(v), b v+c).

The conjugacy sends the fixed origin to a fixed point p. The derivative at p is exactly I because JF_lambda(0)=I. Hence a=b=1, and the fixed-point equation gives c=0. Translating p to the origin if desired gives P(0)=P'(0)=0. Since the flow is nonidentity, P is nonzero.

Write h=H composed with the inverse conjugacy. Then

    h(u+P(v),v)=h(u,v).

If h has positive u-degree n, with leading coefficient A_n(v), the coefficient of u^(n-1) in the difference is n P(v) A_n(v), which is nonzero in characteristic zero. Contradiction. Thus h=f(v). Pulling v back by the polynomial conjugacy gives a polynomial coordinate q with H=f(q).

This also proves local nilpotence directly: in companion coordinates (p,q), D_H is a f'(q) partial_p. There is no need to assume that a formal logarithm commutes with a change of expansion center. The identity and constant-H cases were already separated.

## 3. Exact ground-field statement

The argument above first proves (1) implies (2) over C. For arbitrary characteristic-zero K, take the subfield L generated over Q by the finitely many coefficients of H, lambda, and the asserted polynomial pair F_lambda. It embeds in C. Every coefficient identity defining F_lambda remains true after this embedding, because each formal coefficient uses finitely many arithmetic operations. The complex argument gives integers a,b with

    D_H^a x=0,       D_H^b y=0

after scalar extension. Injectivity of L into C gives the same polynomial identities over L and hence over K. A derivation nilpotent on the two algebra generators is locally nilpotent on the whole polynomial ring, by the iterated Leibniz rule. Thus no descent of the complex conjugating automorphism is assumed.

To pass from (2) to (4) over K itself, use Rentschler's theorem in its characteristic-zero field form: a nonzero locally nilpotent derivation of K[x,y] has kernel K[q] for a polynomial coordinate q. Since D_H H=0, H lies in K[q]. The zero derivation means H is constant. This is an explicit classical descent input, not an inference from the complex conjugacy.

Conversely (4) gives (2) and (3) from the displayed shear formula, and (3) gives (1). This completes the equivalence.

## 4. A singularity obstruction with no fixed-line hypothesis

Suppose H is nonconstant and the origin is an isolated point of the reduced critical locus

    V(H_x,H_y).

Then no nonzero scalar-time formal flow is polynomial.

Indeed, if H=f(q), the coordinate q has nowhere-vanishing gradient: the polynomial identity {p,q}=a is a nonzero constant. Therefore the critical locus is precisely V(f'(q)). As H is in K+m^3, the origin is critical, so f'(q(0))=0. After extending scalars to an algebraic closure, the entire coordinate fiber q=q(0), an affine line in polynomial coordinates, is a critical curve through the origin. This contradicts isolation. Multiplicity of a root of f' changes the critical scheme's multiplicity, not the existence of this reduced curve.

More generally, the reduced critical locus of an allowed polynomial flow is a union of disjoint coordinate fibers. This gives a global geometric certificate to test before calculating long time-series coefficients.

## 5. Two exact examples

### A polynomial flow fixing a parabola

Set q=y+x^2 and H=q^3. Its flow is

    X=x+3s q^2,
    Y=y-6s x q^2-9s^2 q^4.

The invariant is Y+X^2=q. The map has Jacobian one and fixes the parabola q=0 pointwise. It is a polynomial autonomous flow in a nonlinear coordinate; the earlier restriction to a specified line was not intrinsic.

### Arbitrarily deep agreement with a polynomial flow can conceal failure

For M>=3 put

    H_M=y^3+x^M,      H_0=y^3.

H_M has isolated critical point (0,0), since its two derivatives are M x^(M-1) and 3y^2. Therefore its flow is nonpolynomial at every nonzero scalar time. Nevertheless D_(H_M)-D_(H_0)=-M x^(M-1) partial_y, so the two derivations induce the same operator modulo m^(M-1). Their exponentials agree modulo m^(M-1): they have identical coordinate terms through total degree M-2. The H_0 flow is the polynomial shear

    (x+3s y^2,y).

The first discrepancy is exact:

    F_(H_M),s - F_(H_0),s
       = (0,-s M x^(M-1)) modulo m^M.

Higher occurrences of the perturbation in the exponential raise total degree by at least one more. Thus the agreement depth is exactly M-2 for s nonzero. This is an explicit finite-jet obstruction within polynomial Hamiltonians, independent of any unresolved Keller-map example.

## 6. Polynomial autonomous logarithms are a strict sufficient subclass

Every formal symplectic map tangent to the identity has a unique filtration-raising formal logarithm, which is a Hamiltonian derivation in two variables with a Hamiltonian in m^3, unique after fixing its constant term. This follows coefficientwise from the formal exp/log correspondence; the closed formal one-form integrates coefficientwise in characteristic zero.

If a normalized polynomial Keller map has a polynomial formal Hamiltonian logarithm, the theorem proves that it is a polynomial automorphism of conjugate-shear type. There is no justification for assuming this polynomial-log hypothesis for a general Keller map.

In fact it already fails for the explicit polynomial automorphism

    F(x,y)=(x+y^2, y+(x+y^2)^2).

Let h(x,y)=(y,x+y^2). Then F=h^2, and deg(F^n)=4^n by the leading-term recurrence. It is tangent to the identity and symplectic, but is of Hénon type. It has no nonconstant polynomial invariant, so its unique formal Hamiltonian logarithm cannot be polynomial. It is a composition of two polynomial shear flows; autonomy with a polynomial generator is the extra, genuinely restrictive condition.

The route now has a precise boundary: polynomial invariant or polynomial autonomous logarithm supplies global structure; arbitrary finite point-jets do not. Products with different invariants may leave this class immediately.

## 7. Primary sources and the exact external inputs

- Vladimir Shpilrain and Jie-Tai Yu, *Polynomial Retracts and the Jacobian Conjecture*, Corollary 1.7, https://arxiv.org/pdf/math/9701210 . Read directly: characteristic-zero K, invertible polynomial Jacobian matrix, nonconstant fixed polynomial imply automorphism.
- Shmuel Friedland and John Milnor, *Dynamical properties of plane polynomial automorphisms*, Ergodic Theory Dynam. Systems 9 (1989), 67--99, https://doi.org/10.1017/S014338570000482X . Elementary/generalized-Hénon conjugacy dichotomy. Also restated at the opening of the Bedford--Smillie paper below.
- Eric Bedford and John Smillie, *Polynomial diffeomorphisms of C^2: currents, equilibrium measure and hyperbolicity*, Invent. Math. 103 (1991), 69--99, Proposition 4.2 and its proof on pp. 79--81, https://www.math.stonybrook.edu/~ebedford/PapersForM655/BS1.pdf . Read directly; only the polynomial degree-growth conclusion is needed.
- Rudolf Rentschler, *Opérations du groupe additif sur le plan affine*, C. R. Acad. Sci. Paris Sér. A-B 267 (1968), A384--A387. The characteristic-zero-field triangularization statement is explicitly recalled in the introduction of Masakazu Suzuki, *Sur les opérations holomorphes du groupe additif complexe sur l'espace de deux variables complexes*, https://www.numdam.org/item/10.24033/asens.1335.pdf . The original 1968 note was not opened in this session; the stated field scope was checked in Suzuki's primary research paper.

## Scope and next step

This classifies polynomiality of single autonomous formal Hamiltonian flows with H in K+m^3. It does not classify arbitrary Hamiltonians with linear/quadratic parts, arbitrary compositions, rational/algebraic time maps, or arbitrary planar Keller maps. It is not a proof of the planar Jacobian conjecture.

A next useful test for any proposed local-to-global route is to record which global object survives: an exact polynomial invariant, a coordinate fibration, a bound on pole orders at infinity, or a common compactification. Replacing that object by only its finite germ loses the decisive hypothesis, as H_M demonstrates.

## 8. Common invariants for polynomial maps and for finite words

The argument actually gives a broader statement. Let F be a polynomial Keller map over K with F(0)=0 and JF(0)=I. Then the following are equivalent:

- F fixes a nonconstant polynomial;
- F is a polynomial automorphism whose unique formal logarithm is a polynomial locally nilpotent derivation;
- F is polynomially conjugate over K to a shear.

Over C, the first implication follows directly from Sections 1--2, applied to the fixed polynomial rather than an assumed Hamiltonian. In centered conjugate coordinates the shear has form (u+P(v),v) with P(0)=P'(0)=0. Its polynomial logarithm is P(v) partial_u. Conjugation is an isomorphism of completed local rings at the corresponding fixed points, so formal logarithm uniqueness identifies the transported polynomial derivation with log(F). To descend over K, note that each coefficient of log(F) lies in K; its vanishing above the complex polynomial degree therefore descends. Its local nilpotence descends on x and y as before, and Rentschler gives a shear conjugacy over K. The other implications follow by preserving the shear coordinate. The identity is included.

Consequently any finite product of formal symplectic maps tangent to the identity, whose complete product preserves a nonconstant polynomial P, has a conjugate-shear result whenever that result is polynomial. The individual factors need not have polynomial logarithms. The invariant must be an exact identity for the whole word; a truncated invariant does not suffice.

There is an especially transparent Hamiltonian subcase. Let H_i be polynomial elements of K+m^3 and suppose one nonconstant polynomial P satisfies {P,H_i}=0 for every i. Over K(x,y), dP is nonzero and every dH_i is proportional to dP. Hence {H_i,H_j}=0 for every i,j. Their Hamiltonian derivations commute, and any finite scalar-time word is exactly

    exp(D_B),       B=sum_i lambda_i H_i.

Thus its polynomiality is equivalent to B being a polynomial in a polynomial coordinate. This does not require the H_i to be polynomial functions of one preselected H. A common polynomial invariant is sufficient; arbitrary different-invariant words such as the Hénon example do not satisfy it.
