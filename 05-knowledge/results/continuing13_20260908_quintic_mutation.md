# A polynomial shear reveals moving punctures and critical-value collisions

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.**
This constructs globally submersive first functions on fixed W2 with rational
mates, including an explicit genuine source quintic. Their original affine
fibres form a non-isotrivial family of punctured rational curves. Their smooth
projective completions all have genus zero. The full unit response detects
distinct target values while forgetting how many components share a value.
There is no polynomial mate or general Jacobian-conjecture conclusion.

## 1. Inheritance and the operation

The supplier is the [complete univariate family](continuing13_20260908_univariate_fixed_unit.md),
independently recovered from the cubic extension of the previous quadratic
classification. The [previous hidden-arm construction](continuing12_20260908_fixed_unit_hidden_arms.md)
and [pair-degree corollary](continuing12_20260908_pair_degree.md) showed that
every operation on one fixed unit can miss ambient geometry. The present
operation changes the first function using an actual polynomial repair; it
is not another Weyl operation on the old unit.

The complete original-source component theorem and its canonical connection
are [vertical component jets](planar_jc48_sep06_torsion.md), using THM-3412,
`01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md`,
and THM-3770,
`01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md`.
The [coefficient-span lemma](continuing11_20260908_quadratic_family.md) supplies
the Weyl consumer after all fibres and rational constants are paid below.

The source is a polynomial identity J(A,Q)=A^3. The map sends it to
T=lambda*A+Q and the rational primitive 1/(2A^2). It preserves the actual
source Jacobian predicate and, under the stated boundary conditions, global
submersion. It changes the first function, target fibre supports, and the
generic puncture configuration. The sidecars are the full chart, the actual
Eu jet, and the critical values of Q. The cheap hostile is lambda=-f(-2h),
which leaves a rational mate but makes the whole Eu component critical.
Another is confusing a non-isotrivial punctured curve with varying completed
genus. The live board is original fibres, primitive degree, target collisions,
complete unit modules, both boundary charts, and configuration moduli.

## 2. The complete family and global submersion

Fix a squarefree polynomial f in C[z] of degree q>=2 with f(0)=0. Choose h
such that a0=f(-2h)!=0, and lambda such that lambda*(lambda+a0)!=0. Set

    u=x-h, z=u-2h+u^3 t, A=u f(z),
    Q(z)=integral_0^z f(v)^2 dv,
    T=lambda*A+Q(z), G=1/(2A^2).                    (1)

These are literal formulas in the original source C[x,t]. Since f is
squarefree, its zero at 0 is simple. Write b=f'(0)!=0. The source t-degree
of T is exactly 2q+1, with leading coefficient

    lead(f)^2*u^(3(2q+1))/(2q+1).                  (2)

In the actual (u,z) coordinates J_(x,t)(u,z)=u^3. Hence

    J_(x,t)(A,Q)=u^3*f(z)*Q'(z)=A^3,
    J_(x,t)(T,1/(2A^2))=1.                        (3)

On A!=0, dA and dQ are independent, so dT cannot vanish. On a component
z=rho with f(rho)=0, Q'=0 and dT=lambda*dA is nonzero, since f is
squarefree and u cannot vanish there. On Eu={u=0}, the actual source
relation gives z=-2h and z_x=1; therefore

    T_x=a0*(lambda+a0), T_t=0.                     (4)

This is nonzero exactly under the stated additional condition. In
particular u and z must not be treated as independent on Eu.

For the full second chart use x=1/R, t=-R^2-R^4 B and

    M=-h^2(3-hR)-B(1-hR)^3,
    z=R M, u=(1-hR)/R,
    A_inf=(1-hR)M [f(RM)/(RM)].                   (5)

The bracketed expression is a polynomial. Since f(0)=0, Q has a zero of
order three at 0. Thus T_inf=lambda*A_inf+Q(RM) is a polynomial, and on
the ENTIRE added divisor R=0 it has tangential derivative -lambda*b!=0.
The two actual charts prove that T is a global submersion on fixed W2.
For this fixed admissible f,h, the two excluded lambda values are exact:
lambda=0 is critical on every root component of f, and lambda=-a0 is
critical on Eu. They may still admit the rational identity (3).

## 3. The whole field, every fibre, and the moving punctures

The full field is C(A,z), since u=A/f(z) and
t=(z-u+2h)/u^3 recover the original coordinates. Also A=(T-Q(z))/lambda,
so C(x,t)=C(T)(z). The derivation D_T is nonzero on z. Its rational
constant field is exactly C(T).

Let Z be the q distinct roots of f, z0=-2h, and

    S={Q(z0)} union {Q(rho):rho in Z}.              (6)

For c outside this FINITE SET, the complete source fibre has coordinate ring

    C[z, f(z)^-1, (c-Q(z))^-1],
    u=(c-Q(z))/(lambda*f(z)),
    t=(z-u+2h)/u^3.                               (7)

This is the WHOLE source fibre: a point with A=0 lies on Eu or one of the
root components, whose target values are exactly in S. Conversely the
inverses in (7) recover every original-source point of this fibre. The
identities do not merely parametrize a chosen dense subset.

The polynomial Q has degree 2q+1 and Q'=f^2. Thus for c outside S, c-Q
has 2q+1 distinct roots, disjoint from Z. The fibre is P1 minus exactly

    r=q+(2q+1)+1=3q+2                             (8)

points: q fixed finite points Z, one fixed point at infinity, and the
2q+1 moving roots of Q(z)=c. The projective genus is always zero.

These PUNCTURED curves are non-isotrivial, even without labels. Fix one
reference puncture set P of size r. Any isomorphism between two smooth
punctured rational curves extends to a Mobius transformation of their
unique smooth projective completions. Choose three distinct fixed points
from Z union {infinity}; this is possible because q>=2. Their images in P
have at most r(r-1)(r-2) possibilities and determine the transformation.
Consequently at most that many different source puncture sets in this
family can be isomorphic to P. The sets for distinct c are different:
after removing the fixed Z and infinity, their remaining sets are the
disjoint root sets of the different equations Q(z)=c. Therefore the
family has infinitely many isomorphism classes; indeed each class
contains at most r(r-1)(r-2) parameters outside S. No chosen labeling or
unpaid assertion about a moduli space is needed for this proof.

For ANY c, the locus A!=0 in T=c is still exactly (7), and is nonempty
and irreducible. Its closure is one irreducible component E_c^reg. Any
other component lies in A=0, hence is exactly Eu if c=Q(z0), or a root
component E_rho if c=Q(rho). Each root component z-rho is primitive
linear in t, because rho!=z0; all are irreducible. The source submersion
already proved makes the fibre reduced and its components disjoint.
This pays ALL fibres, including every possible collision among S.

## 4. Exact unit supports, hidden multiplicities, and cyclicity

The rational primitive G is regular on E_c^reg. At a root component
E_rho, Q(z)-Q(rho)=O(A^3), since Q'=f^2 has a double zero there. With
g_c=T-c, one has g_c=lambda*A+O(A^3), so the COMPLETE scalar principal
part is

    pp_Erho(G)=lambda^2/(2g_c^2).                  (9)

At Eu the actual relation z=z0+u+u^3t gives

    A=a0*u+f'(z0)*u^2+O(u^3),
    Q(z)-Q(z0)=a0^2*u+a0*f'(z0)*u^2+O(u^3)
              =a0*A+O(A^3).

Thus g_c=(lambda+a0)*A+O(A^3), and

    pp_Eu(G)=(lambda+a0)^2/(2g_c^2).               (10)

There is NO simple-pole coefficient in (9) or (10). This follows from the
absence of an A^2 term in the target reparametrization, rather than from
checking only a leading valuation. Every E_c^reg retains principal part
zero before quotienting by the common diagonal.

Write e_c for the number of A=0 components assigned to c in S. The
complete original-source response torsion at c has e_c full arms. Its
unit has one nonzero principal-coefficient vector at order two, so it
generates EXACTLY ONE of those arms under the canonical Weyl algebra.
Across distinct supports the unit parts separate by polynomial CRT, and
the exact conclusions are

    total ambient torsion arms = q+1,
    total unit-generated arms = |S|,
    unit generates all torsion iff |S|=q+1,
    Ann_C[T]([1])=( product_(c in S)(T-c)^2 ).     (11)

The product has distinct supports but exponent TWO at each support; its
degree is 2|S|, not an unqualified global unit order two. It actually
annihilates: P(T)=product_(c in S)(T-c) vanishes on each reduced component
of A=0, so A divides P(T) in C[x,t] and P(T)^2 G is polynomial. Exactness
follows from the nonzero order-two tuple and the regular component in
each fibre, with the already paid rational constant field.

The pointed unit Weyl module depends only on S: at each support it is a
single principal-part tower with distinguished element g_c^-2, after
rescaling its nonzero coefficient vector. It does not record e_c.
Critical-value collisions and collision with Q(z0) are therefore the
exact lost coordinate for ambient cyclicity. There is no assumption that
these target values are generically or universally distinct.

## 5. Intrinsic rational pair degree and an explicit quintic

Over C(T), A is transcendental and the polynomial equation
Q(z)=T-lambda*A has degree 2q+1 in z. The usual irreducibility of
Q(Z)-Y gives

    [C(A,z):C(T,A)]=2q+1.

Since G=1/(2A^2), one also has [C(T,A):C(T,G)]=2. Therefore

    [C(x,t):C(T,G)]=2(2q+1).                      (12)

Any other rational mate with a nonzero constant Jacobian is a nonzero
constant multiple of G plus a rational function of T, since ker D_T=C(T).
It gives the same embedded pair field. Thus (12) cannot be lowered by
choosing another rational mate.

For the concrete choice h=1, lambda=1, f(z)=z(z-1), take

    u=x-1, z=u-2+u^3t,
    Q=z^5/5-z^4/2+z^3/3,
    T=u z(z-1)+Q(z), G=1/[2u^2z^2(z-1)^2].       (13)

Here a0=6, S={0,1/30,-256/15}, and all three values are distinct.
The source t-degree is FIVE, leading coefficient u^15/5. The full
torsion and the unit-generated torsion both have three arms; the exact
scalar annihilator is

    [T(T-1/30)(T+256/15)]^2.                       (14)

The rational pair degree is TEN. The ordinary source fibres outside S
are non-isotrivial eight-punctured rational curves. This gives an actual
first function at the [live quintic frontier](planar_jc48_sep06_board.md),
with a paid rational mate and integral obstruction. It has NO polynomial
mate, since (14) is nontrivial and the unit is nonzero. The statement does
not weaken the preceding global quartic-pair exclusion or claim this is
the least source degree for non-isotrivial affine fibres.

For a collision control keep f=z(z-1), lambda=1 and choose z0 to be any
root of 6z0^2-15z0+10=0, with h=-z0/2. Then Q(z0)=Q(0)=0, while
f(z0) and 1+f(z0) are nonzero. Thus the global submersion remains valid,
ambient torsion still has three arms, and the unit generates only two.
This is a genuine admissible wall, not a singular parameter boundary.

## 6. Exact controls and scope

The paired engine imports no producer and uses explicit raising gates.
It checks the symbolic Jacobian, complete second chart and boundary,
literal source Eu jet, inverse fields and fibres, pure order-two transfer,
the complete explicit quintic and its target values, an admissible complex
collision modulo its exact quadratic equation, and the two forbidden
lambda controls. Additional finite f choices and pole-module controls
support the analytic all-q and all-operator conclusions; they do not
replace those proofs. Normal and optimized outputs/certificates are
compared as raw LF bytes. No parameter sampling is called a classification.

The strongest paid progress is a real global quintic with variable affine
fibre geometry and a complete target-collision response law. Arbitrary
global cubics/quintics, source unit order one on W2, and JC remain OPEN.

## Independent acceptance

The [independent audit](continuing13_20260908_quintic_mutation_audit.md) accepts
this scoped theorem without mathematical repair. Its 341 exact
gates pass in normal and optimized modes with identical raw LF transcripts
and certificates where applicable. The parent repeated those checks after
filing the frozen sources under the repository reproduction paths.
