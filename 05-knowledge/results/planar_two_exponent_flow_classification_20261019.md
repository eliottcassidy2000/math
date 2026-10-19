# Two-exponent algebraicity and rational-map classification

2026-10-19. **PROVED:** let K have characteristic zero, m≥2, k≥0, and a=k−2m. For H=x^k t^m(1+x²t), every nonzero scalar-time source Hamiltonian flow has algebraic coordinate images exactly in the following cases:

1. a=0;
2. a>0 and a divides m;
3. a=−b, where b≥2 is a proper positive divisor of m+1.

Outside this locus each coordinate is individually transcendental. Within it the flow is a rational map exactly when a=0 or (m,k)=(3,4). All other algebraic cases are nonrational as maps; separate nonrationality and exact degree of each individual coordinate are not asserted here. Every nonzero-time rational map in the stated locus is nonpolynomial. This classifies one two-exponent family, not arbitrary polynomial carriers, compositions, or planar JC.

## 1. The original source flow and its invariants

Work in original polynomial coordinates x,t with {A,B}=A_xB_t−A_tB_x. Put

    u=x²t,  q=1+u,  E=H=x^a u^m(1+u),  δ={−,H}.

The last expression for E is a function-field identity; H itself is polynomial even when a<0. The source derivation satisfies δx∈t^(m−1)K[x,t] and δt∈t^m K[x,t]. It therefore raises t-order by at least m−1≥1 on K[x][[t]]. For every scalar λ∈K the exponential C_λ=exp(λδ) is well-defined and continuous, with inverse C_(−λ). Let X=C_λ(x), T=C_λ(t), U=C_λ(u)=X²T.

Direct differentiation gives

    δE=0,   δu=−a xE,
    E=X^k T^m(1+X²T).                              (1)

If either X or T is algebraic over K(x,t), the invariant equation makes the other algebraic too: it is a nonzero positive-degree polynomial in the remaining coordinate, including k=0. Thus excluding simultaneous algebraicity proves individual transcendence.

## 2. The scaling and trace bridge

We state the details needed to extend the [one-exponent classification](planar_monomial_carrier_flow_trichotomy_20261019.md). The dilation (x,t)↦(zx,z^(−2)t) multiplies H by z^a and the Hamiltonian vector field by weight a+1. Consequently

    X_λ(zx,z^(−2)t)=z X_(λz^(a+1))(x,t),
    T_λ(zx,z^(−2)t)=z^(−2) T_(λz^(a+1))(x,t).       (2)

If a≠−1 and a nonzero scalar-time coordinate pair is algebraic, (2) makes the full formal-time pair algebraic over K(x,t,τ): substitute τ=λz^(a+1), eliminate the finite algebraic extension in z, and clear coefficients to a nonzero relation in K[x,t,τ,Y]. The topology is controlled in

    A=K[x,τ][[t]] ∩ K[x,t][[τ]] ⊂ K[x][[t,τ]].

Both flow coordinates lie in A. Substitution τ↦λz^(a+1) is injective coefficientwise in t, for either sign of a+1, because each coefficient is a polynomial in τ. Regrouping the vanishing relation by τ gives the same relation for the formal-time exponential. Thus no unjustified scalar-time/formal-time interchange is used.

Assume also a≠0. Let C be an algebraic closure of the initial-value constant field K(x,t), choose η∈C* with η^a=E, and in the target flow put v=U_τ and w=η/X_τ. Then

    w^a=v^m(1+v),   dτ=−w dv/(aηE).                (3)

The function v is nonconstant because its initial derivative is −a xE≠0. The target curve has function field M=C(v,w) of degree |a| over C(v): the valuation of v^m(1+v) at v=−1 is one, giving full Kummer degree, with w inverted when a<0. Algebraicity of the formal-time coordinates makes L=C(τ,X_τ,T_τ) finite over C(τ). Both M and L have transcendence degree one over C, and L is finitely generated, so L/M is finite separable.

Tracing (3) gives

    d Tr_(L/M)(τ)=[L:M]·[−w dv/(aηE)].

Trace commutes with differentiation by the sum-over-embeddings formula in characteristic zero. Hence w dv must have a rational primitive on the target curve. No descent of exactness to the original constant field is needed.

## 3. Exactness is equivalent to a rational equation

The cyclic automorphisms w↦ζw allow projection of a primitive to the character of w. That eigenspace is wC(v), so exactness is equivalent to the existence of R∈C(v) satisfying

    a v(1+v)R′+[m+(m+1)v]R=a v(1+v).               (4)

This follows by differentiating wR and using a w′/w=m/v+1/(1+v).

A finite pole of R away from 0,−1 is impossible, since the derivative term has strictly higher pole order. A pole of order r>0 at 0 requires ar=m by the leading coefficient, and a pole at −1 requires ar=1.

For positive a not dividing m, neither possibility occurs: a=1 already divides m, so it is excluded from this case. Thus R must be polynomial. Its leading degree d≥2 would give a nonzero term of degree d+1, impossible on the quadratic right side. For R=A v+B the coefficient equations are

    mB=0,
    (a+m)A+(m+1)B=a,
    (a+m+1)A=a,

which contradict a≠0. Therefore positive a requires a|m.

For a=−b<0 there are no finite poles either. If a polynomial solution has degree d≥2, cancellation of its highest term requires

    bd=m+1.

If d≤1 the same linear coefficient equations are contradictory. Hence b must divide m+1 with quotient at least two. The remaining b=1 case is excluded by a residue: then w=1/[v^m(1+v)], and w dv has the nonzero simple-pole coefficient (−1)^m at v=−1. A derivative of a rational function has zero coefficient of (v+1)^(−1), so no rational primitive exists.

This proves the necessary divisor conditions for a≠0,−1. Sufficiency will be proved constructively in Sections 5–6. The resonant scalar-time case a=−1 requires the next argument because scaling cannot vary time there.

## 4. The exceptional scaling weight a=−1

Now E=u^m(1+u)/x and δu=u^m(1+u), so

    U=exp(λu^m(1+u)∂_u)(u)∈K[[u]],
    U=u+λu^m+O(u^(m+1)).                           (5)

If the source coordinates were algebraic, U would be algebraic over K(u,E). Substitution t=u/x² embeds the source completion into K(x)[[u]]. Here x is transcendental over K((u)), by coefficientwise comparison of a putative polynomial relation in x, and therefore E=u^m(1+u)/x is transcendental over K((u)). Expand a cleared algebraic relation for U as a polynomial in E. At least one nonzero coefficient polynomial in u,Y annihilates U, proving algebraicity over K(u).

Write the elementary partial fraction decomposition as

    1/[u^m(1+u)] = B′(u) + c[1/(1+u)−1/u],
    c=(−1)^m ≠0,   B∈K(u).                         (6)

For example B integrates the finite sum Σ_(j=0)^(m−2) (−1)^j u^(j−m), whose exponents are at most −2. This explicitly gives (6).

The flow preserves its vector field, so U′/[U^m(1+U)]=1/[u^m(1+u)]. Define

    r=(1+U)u/[U(1+u)],
    h=−[B(U)−B(u)]/c.

Equations (6) and vector-field preservation give dh=dr/r. Extend constants to an algebraic closure. If h were nonconstant in the algebraic function field containing u,U, that field would be finite over the rational field in h. The norm n of r to that field would satisfy n′/n=d for a positive integer d, because dr/r=dh. This is impossible: a rational logarithmic derivative is O(1/h) at infinity. Thus h is constant, and so is r. Equation (5) makes r have constant term one; hence r=1. But r=1 gives (1+U)u=U(1+u), or U=u, contradicting the nonzero λ coefficient in (5).

Therefore both coordinates are individually transcendental at a=−1 for every nonzero scalar time.

## 5. Positive divisors: explicit polynomial time primitives

Suppose a>0 divides m, and put r=m/a. Define the source rational parameters

    z=xu^r,   s=1/z.

Then E=z^a(1+u), and

    u=E/z^a−1,
    x=z^(m+1)/(E−z^a)^r,
    t=(E−z^a)^(2r+1)/z^(2m+a+2).                  (7)

Thus K(x,t)=K(E,z)=K(E,s), with E,s algebraically independent. Direct differentiation gives

    δz=z^(m+a+2)/(E−z^a)^r,
    δs=−1/[s^a(Es^a−1)^r].                        (8)

The polynomial

    P_E(s)=−∫s^a(Es^a−1)^r ds                     (9)

is defined by termwise polynomial integration over K(E), has degree d=m+a+1, and satisfies δP_E(s)=1. Additive integration constants are irrelevant.

The source exponential extends to K(x)((t)): the derivation raises t-valuation even on negative Laurent orders and rational x-coefficients. Consequently it agrees with its rational-function-field extension. Applying it to the rational source function P_E(s), whose second derivative under δ is zero, gives directly

    P_E(C_λ(s))=P_E(s)+λ.                           (10)

This is a nonzero polynomial relation at every scalar λ. Hence C_λ(s), and therefore X,T through (7), are algebraic. No root is chosen separately from the existing source exponential.

## 6. Negative proper divisors: a second primitive family

Suppose a=−b, where b≥2 divides m+1 with r=(m+1)/b≥2. Set

    z=x/u^r.

Then

    u=1/(Ez^b−1),
    x=z/(Ez^b−1)^r,
    t=(Ez^b−1)^(2r−1)/z²,                         (11)

so K(x,t)=K(E,z). Direct differentiation gives

    δz=−1/[z^(b−2)(Ez^b−1)^(r−2)].                (12)

The polynomial

    P_E(z)=−∫z^(b−2)(Ez^b−1)^(r−2) dz             (13)

has degree d=m−b and satisfies δP_E(z)=1. The same source-Laurent argument gives P_E(C_λ(z))=P_E(z)+λ for every scalar λ, proving algebraicity of both coordinate images via (11).

For completeness, a general specialization safeguard also applies to both constructions. If algebraic formal-time coordinates are first described by Q(x,t,τ,Y), remove the greatest common power of τ−λ from all coefficient polynomials before specializing. Division of the vanishing relation is valid in the domain K[x,τ][[t]], and evaluation τ=λ is defined coefficientwise there. The resulting specialized polynomial is nonzero and still annihilates the scalar-time coordinate. Thus even a formal-time relation that initially degenerates at a selected λ cannot create an unproved exceptional scalar time. The explicit identities (10) and (13) already avoid this degeneracy directly.

## 7. Rational translations of a polynomial primitive

**Lemma.** Let F be any characteristic-zero field, P∈F[v] have degree d≥2, and λ∈F* be nonzero. There is no R∈F(v) satisfying P(R(v))=P(v)+λ.

**Proof.** A finite pole of R over an algebraic closure of F would give a finite pole of P(R), impossible. Thus R is polynomial. Degree comparison makes R=αv+β with α≠0, and leading coefficients give α^d=1. If α=1, the coefficient of v^(d−1) in P(v+β)−P(v) is d p_d β, so β=0, contrary to λ≠0. If α≠1, the affine transformation R has finite order N because α is a root of unity and its geometric translation sum vanishes. Iteration gives P(R^N(v))=P(v)+Nλ, impossible in characteristic zero. ∎

Apply the lemma over F=K(E) to (9) or (13). If the source map were rational, its transformed primitive parameter would lie in K(E)(s) or K(E)(z), respectively, contradicting the lemma whenever d≥2. Conversely a rational parameter image reconstructs a rational source map using the displayed birational inverse formulas. Positive-divisor cases always have d=m+a+1≥4 and hence give nonrational algebraic maps.

For negative proper divisors, d=m−b. Degree one occurs precisely when m=b+1 and b|(m+1)=b+2. Since b≥2, this forces b=2,m=3, and then k=2m−b=4. Every other negative proper-divisor case has d≥2 and is nonrational as a map.

The preceding argument concerns rationality of the pair. It does not by itself show that each separate algebraic coordinate is nonrational or compute its algebraic degree.

## 8. The two rational loci and nonpolynomiality

When a=0, H=f(u) with f(u)=u^m(1+u) and f′(0)=0. The exact source formula is

    D=1−λx f′(u),
    X=x/D,   T=tD².                                (14)

For nonzero λ, D is nonconstant and coprime to x, so X is nonpolynomial. This is the locus k=2m.

In the unique additional rational case (m,k)=(3,4), formula (12) gives δz=−1 for z=1/(x³t²). Thus C_λ(z)=z−λ. In source coordinates put

    D=1−λx³t²,   B=(1+u)D²−u.

Equations (11) become

    X=xD/B²,   T=tB³/D².                           (15)

These are rational source images, with source-unit denominators D,B. They have a rational inverse at time −λ and preserve H. For λ≠0, D is nonconstant; gcd(D,t)=gcd(D,u)=1 and B≡−u mod D, so gcd(D,tB)=1. Hence T is nonpolynomial. This establishes nonpolynomiality of the rational map without an appeal to a global polynomial inverse theorem.

All source flows have formal Jacobian one. The rational examples are only rational-domain maps, and the nonrational algebraic examples lie on algebraic covers. None is a polynomial Keller counterexample.

## 9. Interpretation, verification and limits

The divisor conditions come from actual valuations of the invariant curve and exactness of its time differential. The algebraic cases have explicit rational invariant parameters above. Rational invariant fibres alone are insufficient: nonzero logarithmic residues or other nonexact time differentials still obstruct algebraic scalar time. The theorem uses exactness and the scaling/trace bridge rather than a genus-only inference.

The [companion checker](../../04-computation/experiments/planar_two_exponent_flow_classification_20261019.py) passes 372 exact controls in normal and optimized Python. [Canonical output](planar_two_exponent_flow_classification_20261019.json), [optimized output](planar_two_exponent_flow_classification_20261019_optimized.json), [run receipt](planar_two_exponent_flow_classification_20261019_run_receipt.json), [full independent review](planar_two_exponent_flow_classification_20261019_independent_review.md), and [independent universal review](planar_two_exponent_flow_classification_20261019_bridge_review.md) accompany the proof. The checker verifies the general Hamiltonian/scaling identities, exactness coefficient conditions, the exceptional partial fractions, parameter inverses and time primitives for representative admissible pairs, and the extra rational example. Finite controls do not replace the all-integer pole and primitive proofs.

The classification is restricted to H=x^k t^m(1+x²t), m≥2,k≥0, over characteristic-zero fields and nonzero scalar times. It leaves individual degree classifications in most algebraic cases, arbitrary polynomial carriers, different-invariant compositions, source-chart entry, and planar JC open. No historical novelty claim or external integration-classification theorem is used.
