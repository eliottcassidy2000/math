# Scalar-time trichotomy for monomial carriers x^k t²(1+x²t)

2026-10-19. **PROVED, characteristic zero:** for H_k=x^k t²(1+x²t), k≥0, every nonzero scalar-time source Hamiltonian flow has algebraic coordinate images exactly when k∈{4,5,6}. At k=4 it is rational but nonpolynomial. At k=5 and k=6 both coordinates are algebraic and individually nonrational, with joint extension degrees four and five respectively. At every other k both coordinates are individually transcendental. This classifies this one monomial family, not arbitrary carriers, compositions, or planar JC.

## 1. Source meaning and basic identities

Let K be any characteristic-zero field. All coordinates x,t are original affine polynomial coordinates; {A,B}=A_xB_t−A_tB_x. Put

    q=1+x²t,  u=x²t,  a=k−4,  E=H_k=x^a u²(1+u).

The expression for E in x,u is a rational-field identity when a<0; H_k itself is always polynomial. Its Hamiltonian derivation δ={−,H_k} has δx∈tK[x,t] and δt∈t²K[x,t]. Hence exp(λδ) exists on K[x][[t]] for each scalar λ∈K, and is an automorphism with inverse exp(−λδ). Write X_λ,T_λ for its coordinate images and U_λ=X_λ²T_λ. The theorem concerns precisely these scalar-time source exponentials.

Invariance and direct differentiation give

    E=X_λ^k T_λ²(1+X_λ²T_λ),
    δu=−a xE.                                      (1)

If either X_λ or T_λ is algebraic over K(x,t), the invariant equation makes the other algebraic as well. Indeed, with one nonzero coordinate adjoined it is a nonzero positive-degree polynomial in the other; for k=0 the equation in X is E−T²−X²T³=0. Thus it suffices to exclude simultaneous algebraicity to prove individual transcendence.

## 2. From one scalar time to algebraic formal time

**Scaling lemma.** Suppose k≠3 and λ≠0. If both scalar-time coordinates are algebraic over K(x,t), then the full formal-time coordinates X_τ,T_τ are algebraic over K(x,t,τ).

**Proof.** The dilation ρ_z(x,t)=(zx,z^(−2)t) multiplies H_k by z^(k−4). The bracket has dilation weight one, so with b=k−3,

    X_λ(zx,z^(−2)t)=z X_(λz^b)(x,t),
    T_λ(zx,z^(−2)t)=z^(−2) T_(λz^b)(x,t).           (2)

These identities follow termwise on the polynomial iterates δ^n(x),δ^n(t). Applying the invertible dilation to a nonzero polynomial algebraic relation for a scalar-time coordinate gives a nonzero relation for the corresponding coordinate at time λz^b over K(x,t,z).

Since b≠0, the field extension K(x,t,z)/K(x,t,τ), with τ=λz^b, is finite algebraic. Transitivity of algebraicity, or elimination of z, yields a nonzero relation Q(x,t,τ,Y)∈K[x,t,τ,Y] for that coordinate. For the elimination step one may first use an abstract τ and a chosen root z of z^b=τ/λ; negative b is handled by clearing powers of z. The relation evaluated at time τ is valid as a source t-adic identity.

Here is the needed topology comparison. Work in the restricted double-series ring

    A=K[x,τ][[t]] ∩ K[x,t][[τ]] ⊂ K[x][[t,τ]].

The substitution τ↦λz^b is injective on K[x,τ][[t]] into K[x,z,z^(−1)][[t]]: coefficientwise in t it is an injective substitution in a polynomial in τ, for either sign of b. Thus the eliminated relation, known after this substitution, already vanishes before substitution. Each t-coefficient of the source exponential is a polynomial in τ, because δ raises t-order. Each τ-coefficient is a polynomial in x,t, because δ is a polynomial derivation. These two descriptions are compatible coefficient arrays, and the same remains true after substitution in the finite polynomial Q. Equality in the source t-adic arrangement says every coefficient of every monomial x^i t^j τ^n vanishes. Regrouping by τ therefore gives the identical relation for the formal-time exponential in K[x,t][[τ]]. This proves the lemma without assuming a map between unrelated completions. ∎

## 3. Algebraic time forces an exact differential on the invariant curve

Assume a≠0 and that the formal-time flow is algebraic. Set F=K(x,t), now regarded as a field of initial constants, and extend these constants to an algebraic closure C. The initial E is a nonzero element of C. Choose η∈C* with η^a=E. In the formal-time flow set

    v=U_τ,   w=η/X_τ.

Then

    w^a=v²(1+v),
    dτ=−w dv/(aηE).                                (3)

Since δu=−a xE≠0, v is nonconstant as a formal function of τ. The function field M=C(v,w) is the invariant-curve function field. It has degree |a| over C(v): for a>0 the polynomial W^a−v²(1+v) is Eisenstein at v+1; for a<0 apply the same argument to w^(−1), using the valuation at v+1 or its reciprocal. Equivalently the valuation of v²(1+v) at v+1 is one, giving full Kummer degree in either sign.

Let L=C(τ,X_τ,T_τ). Algebraicity of the formal-time flow makes L finite over C(τ), so it has transcendence degree one over C. As M⊂L also has transcendence degree one and L is finitely generated, L/M is finite. Characteristic zero makes it separable. Equation (3) is an equality of differentials over C in L.

Trace that equality down to M. Differentiating the sum over the embeddings shows that trace commutes with differentials. Since the right side belongs to M,

    d Tr_(L/M)(τ)=[L:M]·[−w dv/(aηE)].

Thus w dv has a primitive in M. No descent to the original constant field is needed; the obstruction below works already over the enlarged algebraically closed field C.

## 4. The exactness test is a rational first-order equation

Let n=|a|. The cyclic automorphisms w↦ζw, ζ^n=1, act on M and send w dv to ζw dv. If B∈M satisfies dB=w dv, project B to that character:

    B_1=(1/n)Σ_(j=0)^(n−1) ζ^(−j) σ^j(B).

Then dB_1=w dv, and the character eigenspace in M is w C(v), including n=1. Consequently B_1=wR(v) for some R∈C(v). Logarithmic differentiation of w^a=v²(1+v) gives the necessary and sufficient equation

    a v(1+v)R′+(2+3v)R=a v(1+v).                  (4)

We now show (4) has no rational solution for positive integers a≥3 or negative integers a.

A pole away from 0 and −1 is impossible: the derivative term has strictly higher pole order than the other terms. A pole of order r>0 at 0 requires ar=2 by its leading coefficient; at −1 it requires ar=1. Neither is possible for a≥3 or a<0. Hence R must be polynomial.

If R has degree d≥2, its leading term on the left has degree d+1 with coefficient proportional to ad+3. Since the right side has degree two, ad+3=0. For integer a in the stated range this leaves only a=−1,d=3. In every other case R=A v+B has degree at most one. Comparing its constant, linear and quadratic coefficients yields

    2B=0,  (a+2)A+3B=a,  (a+3)A=a.

These force A=0 and then a=0, a contradiction. For the exceptional a=−1,d=3, the operator in (4) sends v^d to (2−d)v^d+(3−d)v^(d+1). A general cubic has zero constant coefficient, and its v coefficient must be −1, forcing its v² coefficient on the left to be −2, whereas the right side requires −1. This is again impossible.

Therefore algebraic formal time is impossible unless a∈{1,2}, apart from the excluded a=0 case. Combined with the scaling lemma, this proves individual transcendence for every k≥0 outside {3,4,5,6}.

For context, at a=1 the primitive is P(v)=v³/3+v⁴/4. At a=2 one may write w=v√(1+v), with primitive (2/5)(1+v)^(5/2)−(2/3)(1+v)^(3/2). The constructive cases below do not rely on an unproved converse to the trace argument.

## 5. The dilation-resonant case k=3

Here a=−1, E=u²(1+u)/x, and δu=u²(1+u). Thus

    U=exp(λu²(1+u)∂_u)(u)∈K[[u]],
    U=u+λu²+O(u³).                                 (5)

If the source coordinates were algebraic, U would be algebraic over K(x,t)=K(u,E). To descend this assertion to K(u), substitute t=u/x² in the source completion, obtaining an injection into K(x)[[u]]⊂K(x)((u)). The indeterminate x remains transcendental over K((u)): a polynomial relation in x with Laurent-series coefficients has all coefficients zero by coefficientwise comparison. Since E=u²(1+u)/x, E is likewise transcendental over K((u)). After clearing denominators, a nonzero relation Q(u,E,U)=0 can therefore be expanded as a polynomial in E; at least one nonzero coefficient polynomial Q_j(u,Y) annihilates U. Hence U is algebraic over K(u).

The one-variable flow preserves its generating vector field, giving

    U′ u²(1+u)=U²(1+U).

For h=1/U−1/u and r=(1+U)u/[U(1+u)], partial fractions give

    dh=dr/r.                                        (6)

We justify the obstruction algebraically. Extend constants to an algebraic closure of K. If h were nonconstant in the algebraic one-variable function field K(u,U), that field would be finite over K(h). The derivation d/dh extends uniquely in characteristic zero. Equation (6) would give (dr/dh)/r=1. The norm n of r to K(h) would then satisfy n′/n=d for a positive integer d. A rational logarithmic derivative is O(1/h) at infinity and cannot be a nonzero constant. Thus h is constant. Equation (6) then makes r constant as well; constants in an algebraic one-variable function field over an algebraically closed constant field are precisely that field.

Writing h=c gives U=u/(1+cu), hence

    r=[1+(c+1)u]/(1+u).

Constancy forces c=0, so U=u. This contradicts the nonzero λ coefficient in (5). Therefore U, and consequently each source coordinate, is transcendental at every nonzero time when k=3.

## 6. The algebraic cases k=4 and k=5

At k=4, E=u²(1+u), and the [common-u rational formula](planar_common_u_rational_flows_20261019.md) gives

    D=1−λxu(3u+2),   X=x/D,   T=tD².

The flow is rational. For λ≠0, D is nonconstant and coprime to x, so X is nonpolynomial. T itself is polynomial; no individual nonpolynomiality claim is made for T.

At k=5, [the quartic source-flow theorem](planar_quartic_algebraic_source_flow_20261019.md) proves both coordinates algebraic and individually nonrational, with joint extension degree four, using P(U)=P(u)−λE² and irreducibility. Its branch is the same original source exponential used here.

## 7. The quintic algebraic case k=6

Put z=x³t and E=x⁶t²q. Then

    u=E/z²−1,
    x=z³/(E−z²),   t=(E−z²)³/z⁸,
    δz=z⁶/(E−z²).                                  (7)

These are function-field identities, giving K(x,t)=K(z,E). The rational function

    A_E(z)=1/(3z³)−E/(5z⁵)

satisfies δA_E(z)=1. The source exponential extends to its fraction field and to K(x)((t)): the derivation raises t-valuation by at least one even for negative Laurent orders, and its action on rational x-coefficients has the same estimate. Thus for Z=exp(λδ)(z), one has A_E(Z)=A_E(z)+λ. This can also be obtained by applying the exponential to a cleared rational identity: differentiation terminates after its first derivative. Consequently

    Q(Z)=(3E−5z²−15λz⁵)Z⁵+5z⁵Z²−3Ez⁵=0.          (8)

For λ≠0 this quintic is irreducible over K(z,E). Indeed regard it first as linear in E:

    3(Z⁵−z⁵)E + 5z⁵Z²−5z²Z⁵−15λz⁵Z⁵.

The two coefficients are coprime in K(z)[Z]. At any root Z=ζz of Z⁵−z⁵, with ζ⁵=1 in an algebraic constant extension, the second coefficient is 5z⁷(ζ²−1−3λz³), which is nonzero since z is transcendental and λ≠0. A primitive linear polynomial in E is irreducible over K(z)[Z,E] by Gauss's lemma. Viewed as a polynomial in Z it is primitive over K(z)[E], since its Z² coefficient 5z⁵ is a unit there. A second Gauss argument proves irreducibility over K(z,E), of degree five.

Hence Z has degree five over the source field. Its coordinate images satisfy

    X=Z³/(E−Z²),   T=(E−Z²)³/Z⁸.                    (9)

These are algebraic. Conversely Z=X³T, so their joint field has degree five. If X were rational, Z³+XZ²−XE=0 would contradict degree five. If T were rational, divide (E−V²)³−TV⁸ by Q(V). The remainder has degree less than five and its V⁴ coefficient is exactly 3E≠0, independently of T. It would be a nonzero polynomial of degree below five annihilating Z, again impossible. Thus each coordinate is nonrational.

Every step concerns the existing source branch Z; quintic irreducibility is not used to choose an arbitrary root. The formal inverse exponential does not become a rational inverse map.

## 8. Scope and verification

The theorem gives a full algebraic/rational/transcendental trichotomy for the individual monomial family H_k=x^k t²q. It does not classify arbitrary polynomials R in H=t²qR, different-invariant words, or filtration limits outside previous hypotheses. It gives no polynomial Keller counterexample and no resolution of planar JC.

The [companion checker](../../04-computation/experiments/planar_monomial_carrier_flow_trichotomy_20261019.py) passes 118 exact controls in normal and optimized Python. [Canonical output](planar_monomial_carrier_flow_trichotomy_20261019.json), [optimized output](planar_monomial_carrier_flow_trichotomy_20261019_optimized.json), [run receipt](planar_monomial_carrier_flow_trichotomy_20261019_run_receipt.json), [full independent review](planar_monomial_carrier_flow_trichotomy_20261019_independent_review.md), and [independent bridge review](planar_monomial_carrier_flow_trichotomy_20261019_bridge_review.md) accompany the proof. The checker tests scaling, the rational exactness operator and its exceptional resonance, the k=3 differential identity, and the k=6 birational, primitive, quintic and remainder identities. Universal exclusions use the scaling/trace and rational-pole proofs, not a finite range of tested exponents. The proof is self-contained; no binomial-integration classification or analytic transcendence theorem is invoked.
