# A fixed cubic-order unit with arbitrarily many hidden source arms

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
This classifies an entire univariate ansatz of global functions on fixed W2.
For every original source degree q>=2 it has exact unit order three and
the same pointed two-arm Weyl module, while its full source torsion has q
arms and its rational-pair field degree is 2q+1. Cubic source degree already
gives three ambient arms. The result does not classify arbitrary cubic
first functions or settle the search for unit order one outside this ansatz.

## 1. Inheritance and the whole family

The closest mechanisms are the source component-jet theorem in
[Vertical component jets](planar_jc48_sep06_torsion.md), the coefficient-span
lemma in [the fifth-order family](continuing11_20260908_quadratic_family.md),
and the fixed-unit construction in
[Continuing12](continuing12_20260908_fixed_unit_hidden_arms.md).
The recovered cubic sidecar is
[THM-3975 / danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate](../../01-canon/theorems/THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate.md).
That result concerns a different completion, with source-form divisor D
rather than W2's 2D; its chosen first function has no rational mate. It
cannot supply a W2 cubic by forgetting the completion. Likewise the
[prime-leading filter](planar_jc48_sep08_prime_leading.md) explicitly needs
a polynomial mate and is not a filter for arbitrary rational mates here.

The present family grew from a cubic probe and the root agent's univariate
extension. The live concepts are boundary normal order, labelled critical
values, full versus generated torsion, exact rational fields, and pair
degree. The hostile is a source submersion with an added critical point.
The corrected near miss is discarding a regular component or requiring
every critical value to be nonzero instead of retaining their full vector.
The least-used operation is integrating the square of one univariate
polynomial before forming the global first function.

Use the actual fixed charts

    x=1/R, t=-R^2-R^4B, omega=dx wedge dt=R^2 dR wedge dB.

For an arbitrary complex h, put

    u=x-h, z=u-2h+u^3t, F=u f(z),                      (1)

where f is a nonconstant complex polynomial. The exact whole-family
globality condition is f(0)=0. Within that global family, F is everywhere
submersive if and only if f is squarefree and f(-2h)!=0.

For the main module theorem assume additionally q=deg f>=2. These
conditions imply h!=0 and f'(0)!=0; they are not omitted extra hypotheses.
The original source t-degree of F is exactly q.

All response statements use the original ring
C[x,t]/D_F C[x,t], where D_F G=F_xG_t-F_tG_x. They do not replace it by
O(W2) or by a rationally substituted affine chart. No external priority
claim or general Jacobian-conjecture consequence is asserted.

## 2. Exact globality and both-chart submersion

Define J=h^2(3-hR)+B(1-hR)^3. The complete coordinate identity is

    z_inf=-R J, u=(1-hR)/R.                            (2)

Thus z itself is global. The only possible boundary pole of u f(z) is
f(0)/R. This proves globality exactly when f(0)=0. If f(z)=sum a_i z^i,
the complete extension is

    F_inf=(1-hR) sum_(i>=1) a_i(-1)^i R^(i-1)J^i,
    F_inf|D=-f'(0)(3h^2+B).                            (3)

In the actual source, u=0 implies z=-2h and F_x=f(-2h), F_t=0. Hence
submersion requires f(-2h)!=0. On u!=0, (u,z) are valid coordinates
with J_(x,t)(u,z)=u^3. The two derivatives of u f(z) in these coordinates
are f(z) and u f'(z), so their common vanishing is equivalent to a repeated
root of f. Every such root is realized by actual source points with u!=0.
This proves necessity and sufficiency of squarefreeness in the open chart.
Finally f(0)=0 and squarefreeness give f'(0)!=0, and (3) makes the added
boundary tangential derivative the nonzero constant -f'(0). All points
have now been paid, proving the exact iff in Section 1.

## 3. Rational primitive, complete fibres, and pair degree

Let Q be the unique polynomial satisfying

    Q'=f^2, Q(0)=0, g=F,
    G=Q(z)/g^3.                                       (4)

In (u,z) coordinates,

    J_(x,t)(F,Q(z))=u^3 f(z) Q'(z)=u^3 f(z)^3=F^3,

so J_(x,t)(F,G)=1 exactly. The third scalar repair g^3G=Q(z) is a
polynomial in the original source and is also global on W2.

The whole function field is C(F)(z): u=F/f(z), then
t=(z-u+2h)/u^3. The derivation on z is F^3/f(z)^2!=0, so its full rational
constant field is exactly C(F), with no implicit field extension.

For every c!=0, the complete source fibre F=c has coordinate ring

    C[z,f(z)^-1], u=c/f(z), t=(z-u+2h)/u^3.            (5)

Indeed u and f(z) are units since their product is c, and these inverse
formulas recover the original generators. This pays the entire fibre,
not only a dense open. Every nonzero fibre is irreducible, with rational
smooth projective completion. Its punctures are the q roots of f and
infinity, independently of c.

Since Q has degree 2q+1, (4) also gives the exact function-field degree

    [C(x,t):C(F,G)]=2q+1.                              (6)

This is the degree of the nonconstant polynomial map z -> Q(z) over
C(F), since F^3G=Q(z). It is not an assertion that a rational pair is
defined or finite at every point of W2.

## 4. Every source component and every scalar principal part

The zero fibre has exactly q+1 reduced, irreducible, pairwise disjoint
source components:

    E_u: u=0;  E_rho: z=rho for each root rho of f.     (7)

For each root, z-rho=u^3t+u-2h-rho is primitive linear in t, since its
constant coefficient at u=0 is -2h-rho!=0 by f(-2h)!=0. This proves
irreducibility and separation from E_u. Distinct roots give disjoint
components. Together with (5), the inherited component-jet theorem gives
exactly q full torsion arms in the original source module, all at g=0.

The complete principal parts of the primitive (4) are

    E_u:   Q(-2h)/g^3+f(-2h)/g^2,
    E_rho: Q(rho)/g^3.                                (8)

At a simple root rho, Q'=f^2 has a double zero, so Q(z)-Q(rho) is
divisible by (z-rho)^3. Since g=u f(z) is a local uniformizer there,
the remaining quotient is regular. Thus there are no hidden g^-2 or
g^-1 terms on those components. The root rho=0 has Q(0)=0 and is regular;
this component remains present before quotienting by the common diagonal.

For E_u put alpha=-2h. The actual source jet is
z=alpha+u+u^3t, and it must be retained. Taylor expansion gives

    Q(z)=Q(alpha)+f(alpha)^2u+f(alpha)f'(alpha)u^2+O(u^3),
    g=u f(z)=f(alpha)u+f'(alpha)u^2+O(u^3).

Therefore Q(z)-Q(alpha)-f(alpha)g is divisible by u^3. Since
f(alpha)!=0, its quotient by g^3 is regular at E_u. This proves its
complete entry in (8), including the absence of a simple pole.

## 5. Exact fixed pointed module in all degrees q>=2

Use the regular component E_0 as reference after retaining all labels.
Then (8) represents the distinguished unit as

    theta=A/g^3+B/g^2,
    B=(f(-2h),0,...,0),
    A=(Q(-2h),Q(rho_1),...,Q(rho_(q-1))).              (9)

B is nonzero. At least one nonzero root of f has Q(rho)!=0. Otherwise
Q would vanish at every root of f, and Q'=f^2 would make every such
zero have multiplicity at least three. Hence f^3 would divide Q. But
3q>2q+1=deg Q for q>=2, a contradiction. This argument does not require
all critical values to be distinct or nonzero.

Consequently A and B are independent for every admitted parameter, and
theta has exact scalar annihilator (g^3). By the coefficient-span lemma
the unit generates exactly two complete torsion arms, regardless of q.
The quotient of the full source torsion by this generated module has
q-2 additional full arms.

Let D=C<g,partial>/(partial*g-g*partial-1), where partial is the canonical
connection on principal parts. The full exact left annihilator is

    I=D*g^3+D*(1+partial*g+(1/2)partial^2*g^2).         (10)

This gives the same pointed module (D/I,[1]) for every q>=2 and every
admitted f,h. For completeness, in the free C[partial]-basis A/g,B/g,
the two elements g theta and g^2 theta have columns

    [[-partial,1],[1,0]],                             (11)

whose determinant is -1. Thus they form another free basis. Both
operators in (10) annihilate theta by direct differentiation. Conversely
normal-order any Weyl operator modulo D*g^3 as
a_0(partial)+a_1(partial)g+a_2(partial)g^2. The second relation in (10)
eliminates the constant term. Independence of the two columns in (11)
forces both remaining coefficients to vanish if the operator kills theta.
This proves (10) for arbitrary operators, beyond finite controls.

Every j-th canonical derivative of theta has exact scalar order j+3.
All finite sequences of multiplication and canonical differentiation
therefore have identical relations across this family, even as both
ambient arm count q and rational-pair degree 2q+1 grow without bound.
Additional source response classes or fibre data are needed to see them.

## 6. Cubic entry, sharp boundaries, and the order-one stopping reason

The first new ambient gain occurs at q=3. For example take

    f=z(c+d z^2), c*d*h*(c+4dh^2)!=0,
    F=u z(c+d z^2),
    Q=d^2 z^7/7+2dc z^5/5+c^2 z^3/3.                 (12)

The zero fibre has four components and three full source arms. At each
nonzero root rho, Q(rho)=8c^2rho^3/105!=0. The unit still generates the
two arms and exact ideal (10). At q=2 the construction recovers the a=0
part of the already classified seventh-order quadratic pencil, so the
new theorem agrees with the complete quadratic two-arm ceiling.

The boundary q=1 is different. For f(z)=a z, Q=a^2z^3/3 and the only
root is the regular zero root. The component-coefficient space has only
one dimension, and the unit generates one arm. The strict degree argument
in Section 5 becomes equality and must not be used there. Repeated roots
of f create actual open-source critical points; a root at -2h creates a
source critical line; f(0)!=0 creates an actual boundary pole.

Critical-value collisions alone do not cause these failures. For
f=z(z-1)(z-a), a=(7+sqrt(-7))/14, the exact control gives Q(1)=0 but
Q(a)!=0. This globally submersive example retains both coefficient
directions despite a nonzero root sharing the regular root's Q-value.

The initial cubic unit-order-one search has a precise stopping reason in
this whole ansatz: (8)-(9) force a nonzero cubic scalar principal part
for every admitted q>=2. No parameter adjustment or further derivative
can turn its distinguished unit into order one. This does not exclude
order one for a different original-source cubic family. A narrower cubic
probe with a fourth-order coefficient was overtaken by this stronger
all-degree structure and is not used as a dependency.

## 7. Reproduction and finite universe

The runner imports no mathematical producer. It checks the full source
and added-chart identities, actual source coordinate jets, polynomial
and rational Jacobians, complete nonzero-fibre inverses, all root principal
parts, and formal Weyl actions. Its declared degree universe is q=2,...,8
with f=product_(j=0)^(q-1)(z-j), h=1. Additional exact controls pay the
symmetric cubic, the critical-value collision above, and every excluded
boundary. Derivatives j=0,...,8 are checked. The all-degree, all-parameter,
full-operator proofs are analytic, not conclusions from this finite bank.

Reproduce after relocation:

    python 04-computation/continuing13_20260908_univariate_fixed_unit.py
    python -O 04-computation/continuing13_20260908_univariate_fixed_unit.py

All 222 gates remain active under optimization. The runner writes its
deterministic certificate into 05-knowledge/results after relocation,
or beside itself in the external packet. Normal and optimized raw LF
outputs and certificate bytes are compared before freeze. No maintained
repository file or frozen supplier was edited to produce this packet.

## Independent acceptance

The [independent audit](continuing13_20260908_quintic_mutation_audit.md) accepts
this scoped theorem without mathematical repair. Its 341 exact
gates pass in normal and optimized modes with identical raw LF transcripts
and certificates where applicable. The parent repeated those checks after
filing the frozen sources under the repository reproduction paths.
