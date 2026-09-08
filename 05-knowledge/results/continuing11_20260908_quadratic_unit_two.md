# A globally submersive quadratic on W2 realizes unit order two and both torsion arms

**Status: PROVED + FINITE-EXACT; independently audited.**
This supplies a positive answer to the preceding source-linear classification's
explicit higher-degree question. The surface remains fixed at W2, and the first
function has genuine source-t degree two. No polynomial mate or solution of
JC(2) is claimed. The response ring is always the original C[x,t].

## 1. Inheritance, the actual object, and the theorem

Use the actual surface and charts of
`planar_jc48_sep08_dg_quadratic.md`:

    W2=(P1_x x P1_Z) minus{Z=x^2},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^2-r^4b, D=W2 minus U0={r=0},
    omega=dx wedge dt=r^2 dr wedge db.

The complete quadratic rational-mate criterion is inherited from
`planar_jc48_sep08_exact_pencils.md`; its fifth-order pencil
u^5 span{1,u} is the closest mechanism. The source-linear gap is
`continuing10_20260908_dg_linear_all_m.md`: a globally submersive rationally
integrable source-linear function on any W_m cannot have distinguished unit
order two. That was a theorem about one source layer, not an exclusion of
higher-degree functions. The canonical response and component-jet theorem is
`planar_jc48_sep06_torsion.md`, with **THM-3412**,
`01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md`,
and **THM-3770**,
`01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md`.

The canonical hostile is the older order-two source-linear v=x+x^3t, smooth
on U0 but critical at the added boundary. Its h(b) thickenings still have a
boundary critical point. The corrected near miss is importing affine
submersivity into the global surface. The least-used operation here couples a
translated repeated root in an exact quadratic pencil to the normal derivative
on the actual boundary; no source translation is called a global automorphism.

Set

    u=x-1, w=1+u^2t, y=u w, L=(u+1)w-4,
    F=u(u+1)w^2-4u w=u w L,
    G=((1-2u)w+2)/(3u^2w^2).

Equivalently the complete original polynomial is

    F=x(x-1)^5 t^2+2(x-1)^3(x-2)t+(x-1)(x-4).

**Theorem.** F is global and has no critical point anywhere on W2. It has
the rational mate G with J_(x,t)(F,G)=1. In

    C_F=C[x,t]/D_F(C[x,t]), D_F=F_x partial_t-F_t partial_x,

the distinguished unit theta=[1] has exact annihilator(F^2). The complete
torsion submodule has two component principal-part arms, and theta generates
both under multiplication by F and the canonical transverse connection nabla.
More precisely, with g acting as F and

    D_alg=C<g,nabla>/(nabla g-g nabla-1),

there is an isomorphism of left modules

    tors_(C[F])(C_F) = D_alg theta ~= D_alg/(D_alg g^2).

The proof below retains every original-source component. It does not identify
C_F with a quotient of O(W2), on which this Hamiltonian derivation need not be
regular because omega vanishes at D.

The concept board is: exact quadratic pencils; both physical charts; fixed
rational constants; all three special components; unequal scalar principal
parts; and the canonical Weyl action.

## 2. Global regularity and all-point submersivity

In the full boundary chart put

    zeta=r[1+b(1-r)^2].

Literal substitution gives

    F_inf=(1-r)(zeta^2-4), y_inf=(1-r)(2-zeta).

Both are polynomial on all Uinf. Along the entire added D one has

    F=-4, F_r=4, F_b=0.

Thus there is no boundary critical point, including any special value of b.
On the omitted line u=0 in U0, the original derivative is F_x=-3. On u!=0,
(u,y) are valid coordinates because J_(x,t)(u,y)=u^3, and

    F=(1+1/u)y^2-4y,
    partial_u F=-y^2/u^2,
    partial_y F=2(1+1/u)y-4.

The first derivative can vanish only at y=0, where the second is-4. These
charts and the u=0 line cover the source and the added boundary. This proves
global submersivity without a finite point search or a generic-parameter claim.

The source discriminant at fibre value c is

    P^2-4N(Q-c)=4u^5[(c+4)u+c],

where F=N t^2+P t+Q. Hence this actual global function lies on the inherited
exact pencil u^5 span{1,u}; the rational primitive is also checked directly.

## 3. Rational primitive, constants, and the three-component fibre

Solving the displayed expression for F gives

    u=y^2/(F-y^2+4y),
    C(x,t)=C(u,y)=C(F)(y).

The derivation satisfies

    D_F y=-u y^2=-y^4/(F-y^2+4y).

Its rational constants are therefore exactly C(F). With F held fixed,

    G=F/(3y^3)+2/y^2-1/y

has D_F derivative1. Substituting y=u w gives the original rational formula
above, so no parameterwise primitive is being mistaken for a rational source
function. Its only affine poles are on u=0 or w=0.

The fibre F=0 has exactly the three factors u,w,L. They are pairwise
comaximal: w at u=0 is1; L at u=0 is-3; L at w=0 is-4. They are reduced and
irreducible. The factor w=1+u^2t is primitive and linear in t. So is

    L=(u+1)u^2t+u-3,

since u^2(u+1) and u-3 are coprime. Call their components E_u,E_w,E_L.

There is no other reducible fibre on the ORIGINAL U0. For c!=0, F-c is
primitive in C[u][t]: the only common root of its leading and linear
coefficients could be u=0, but its constant coefficient there is-c. Its
discriminant4u^5[(c+4)u+c] is nonsquare in C(u). For c=-4 it is a nonzero
multiple of u^5, still nonsquare; for other nonzero c it has the odd-order
root u=0 and a distinct moving root. The quadratic is therefore irreducible.
This is not a statement about all W2 fibres: D is an additional component
at the boundary value F=-4.

## 4. Exact unit order and the complete scalar principal parts

The polynomial

    P2=F^2 G=L^2[(1-2u)w+2]/3

satisfies D_F P2=F^2. It is an actual original-source polynomial, so
F^2 theta=0. On E_w, where u is a unit, G has a pole of order two with nonzero
coefficient; on E_L it is regular. The fibre F is a uniformizer on each
component. Thus F G has a genuine simple pole on E_w and is regular on E_L.
Every rational primitive of F differs from F G by H(F), by the constant-field
identity. Cancelling the E_w pole with H creates an uncancelled pole on E_L.
Hence no polynomial primitive of F exists and F theta!=0. This proves the
complete C[F] annihilator(F^2), with no degree bound on a possible primitive.

The full scalar principal parts in the SAME variable g=F are

    E_u: 9/g^2+4/g,
    E_w: (32/3)/g^2+4/g,
    E_L: 0.

For E_u one must first retain w=1+u^2t. Then g=-3u+u^2+O(u^3) and
G=u^-2-(2/3)u^-1+O(1), giving the stated scalar part. For E_w, the valid
(u,w) chart gives g=-4u w+u(u+1)w^2 and the rational formula for G gives the
second part. The remaining denominators u,w are units along E_L.
The source separately checks both leading and simple coefficients by exact
limits with the correct local coordinates. Neither component is discarded.

## 5. A single actual unit generates the two complete torsion arms

The polynomial gradient is a unit ideal by the affine submersion proof, and
the rational constants are C(F). The hypotheses of the inherited component-jet
theorem are therefore satisfied. Since only F=0 is reducible on U0 and it has
three components, its complete torsion is

    (C^3/C(1,1,1)) tensor g^-1 C[g^-1].

The canonical connection differentiates scalar principal parts. Represent the
component quotient by triples with final coordinate zero, and put

    A=(9,32/3,0), B=(4,4,0).

Their determinant in this two-dimensional quotient is-20/3, so they are
independent. The actual unit is theta=A/g^2+B/g, and

    g theta=A/g, (g nabla+2)theta=B/g.

Repeated derivatives supply every negative power in both independent arms.
Thus theta generates the whole torsion under the actual Weyl action. In
particular its j-th canonical derivative has exact order j+2, but the result
here retains more than that single order statistic.

To identify the full left annihilator, normally order a Weyl operator with
powers of g on the right. Modulo the left ideal D_alg g^2 it has a unique
representative a(nabla)+b(nabla)g. If this kills theta, its B-component is
a(nabla)(1/g). The functions nabla^j(1/g)=(-1)^j j! g^(-j-1) are independent,
so a=0. The A-component then similarly forces b=0. Conversely g^2 kills
theta. This proves the stated isomorphism D_alg/D_alg g^2 with the full
torsion module, not only a surjection or a finite jet comparison.

## 6. Scope, reproduction, and the connection that succeeded

The example resolves the explicit question whether the source-linear unit
order-two gap persists into actual global DG quadratics: it does not. It
realizes the order on fixed W2 with a first function smooth everywhere. It
does not contradict the complete source-linear classification, produce a
polynomial/global mate, or establish a Keller entry mechanism.

The successful transfer retains the full quadratic pencil and both charts,
then passes through the paid constant field to actual component principal
parts. The scalar repair order alone would lose the two independent component
directions. Restoring them reveals the cyclic Weyl presentation of the entire
torsion. This is an explicit realization of inherited mechanisms, with no
external-priority claim for torsion or Weyl module theory.

The standalone source checks the original and full second-chart identities,
the all-point gradient mechanism, both rational field inverses, the exact
mate and polynomial annihilator brackets, all three component controls, the
entire discriminant pencil, every scalar principal coefficient, and literal
finite controls for the Weyl commutator and normal forms through derivative6.
These validate the formulas; the unbounded statements use the proofs above.
Run normally and with `-O` after filing under04-computation. It imports no
other producer and writes a standard certificate under results.

Both modes pass56 always-active exact gates with identical raw LF stdout;
certificate regeneration is unchanged. Frozen SHA256:

    Source 09e286549c66c0944c7e7dc98d16aefd66fd37a4092224dbdf7c2075177b4742
    Output 70167bb62e53ec201dc9050404423a41033269e80e4bfe2b95c6f44ec9cf0642
    Certificate ea3124d9334fb1fe4e4c80892f34ee353e948d3a8c2ae72c16ae4d1afee2f5ef

The independent audits below are accepted; their all-degree proofs pay the
separate source, principal-part and Weyl-module gates.

## Accepted independent audits

- [quadratic unit two audit](continuing11_20260908_quadratic_unit_two_audit.md).
- [quadratic weyl audit](continuing11_20260908_quadratic_weyl_audit.md).

The parent accepted these reviews without mathematical repair. Source, output
and certificate bytes remain frozen; this filing changes status and routing only.
