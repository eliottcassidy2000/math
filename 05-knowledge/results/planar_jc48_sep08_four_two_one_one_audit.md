# Independent audit of the complete binary 4+2+1+1 boundary exclusion

**Status: INDEPENDENT ANALYTIC / SOURCE / NORMAL-AND-OPTIMIZED REPLAY PASS.**

This audits [the primary proof](planar_jc48_sep08_four_two_one_one.md) and
its complete frozen source/output. I independently reconstructed the
original shifted section constraints, both infinity charts and their
normalized branches, the generic-integrality argument, the full divisor
of the relative differential, the complete primitive space, and the two
coefficient obstructions. Both fresh source modes reproduce the frozen
output exactly. A separate reconstruction without importing the producer
checked thirteen further symbolic identities, including the full local
faces and a bracket reduction obtained without polynomial division.

The producer applied both harmless prose clarifications: in the
explicit infinity coordinate used below the lower numerator equals
-lambda at the origin, and the first-jet failures give nonzero residues,
whereas the earlier regular alternatives give the holomorphic
contradiction. Neither changes an equation used in the proof or any
frozen source/output byte. The final primary also explicitly allows
the primitive coefficients in the algebraic generic constant field,
as justified below. I read and accepted those final clarifications.
The audit accepts the mathematical argument;
status promotion and Git remain parent-owned.

## 1. Accepted scope and recovered dependencies

The object is the fixed DG surface

    W=(P1_x x P1_z) minus {z=x^2},
    t=1/(z-x^2), omega=dx wedge dt,

with genuine global H in L2 and L in L1, and F=H^2+L.
Its complete binary-octic leading section has four distinct points with
multiplicities 4,2,1,1. The accepted rational conclusion is necessary:
a rational mate forces the double point to be the original infinity,
the finite fourfold point to be the midpoint of the two simple points,
and L to be constant. It does not assert that every H satisfying
these necessary conditions has a rational mate. An explicit family
attains them at every finite fourfold position.

For all placements, no polynomial mate of any degree exists. This
includes a mate not global on W. The rational exclusion is restricted
to the branches actually excluded in the proof; the stated rational
positive prevents conflating the two conclusions. No general planar
Jacobian conclusion or automatic envelope representation is claimed.

The exact leading bridge and the complete differential classification
are the proved suppliers
[leading_exactness](planar_jc48_sep08_leading_exactness.md) and
[boundary_exactness](planar_jc48_sep08_boundary_exactness.md).
The full rows are from
[dg_quadratic](planar_jc48_sep08_dg_quadratic.md).
The local value/jet arguments use
[shared_roots](planar_jc48_sep08_shared_roots.md) and
[quartic_common_root](planar_jc48_sep08_quartic_common_root.md).
The closest complete-genus/primitive-space predecessor is
[four_three_one](planar_jc48_sep08_four_three_one.md), with its
[independent audit](planar_jc48_sep08_four_three_one_audit.md).
No reserved supplier is needed.

If the double point is finite, dx/sqrt(N) has a genuine simple pole
with nonzero residue there. That immediately covers every placement
except double infinity. In that remaining placement deg N=6, and the
proved 4+1+1 row gives exactly the midpoint condition. Scaling the
original x and z coordinates together normalizes the nonzero
half-separation; target scaling removes the leading scalar. These maps
have nonzero constant source Jacobian. The surviving finite position
p stays arbitrary. Writing u=x-p does not assert an automorphism of W.

## 2. Complete all-p rows and the actual entry conditions

After these legitimate scalings, N=u^4(u^2-1). The original numerator
description gives, with x=u+p,

    H=[Q z^2+(P-2Qx^2)z+(N-Px^2+Qx^4)]/(z-x^2)^2.

Each numerator row has degree at most four in x. The full conditions
are equivalent to

    deg P<=4,
    Q=(A4-1)u^2+(A3-2pA4+4p)u+C0,
    P=sum A_i u^i.

The coefficients of u^6 and u^5 in the last numerator row give a
triangular system with unit pivots for Q2 and Q1. All remaining
coefficients are unconstrained. This both proves completeness and
prevents a generic-p pivot assumption. For L=Mt+R the two numerator
rows have degree at most two and force

    M=sum B_i u^i,
    R=B4 u^2+(B3-2pB4)u+E0.

The frozen source checks every complete row, not just selected
functions in its span.

At each simple boundary zero the local differential is regular on
every normalized branch, for arbitrary lower coefficients. At the
double infinity point every possible regime is regular after
multiplication by the actual canonical factor r^2. The M-unit and
normal-unit cases are already regular. When both values vanish and
the normal derivative also vanishes, the generic tangent quartic has
simple nonzero roots because its tangential quadratic coefficient is
nonzero. The unweighted local logarithms become regular after the
r^2 factor. This exhausts the local degree; it is not a statement
based on one favorable branch.

Thus possible primitive poles occur only at the finite fourfold
point. If M(0) is a unit, its multiplicity four is not any balanced
value 3j, so the complete M-unit analysis is regular. If M(0)=0 but
P(0) is a unit, the shared normal-unit analysis is regular. In either
case eta is holomorphic on each compact normalized generic component
and cannot have the required nonconstant rational primitive.
This step is valid before proving geometric integrality: the relative
form is nonzero on every component meeting W away from D, and no
generic component equals the fixed divisor D or S. The finitely many
exceptional denominator fibres of a rational G can be avoided.

After both values vanish, the finite first-jet tangent calculation
gives actual nonzero logarithmic residues unless

    P'(0)=M'(0)=0.

Here the leading tangential quadratic coefficient of the numerator
is zero because its boundary multiplicity is four. Factoring the
tangent polynomial by Z^2 leaves a quadratic with nonzero constant
if P'(0)!=0; when P'(0)=0 and M'(0)!=0 it leaves a nonzero linear
factor. Generic fibre values give simple nonzero roots. These
arguments are residue obstructions to a rational primitive, not
polynomial-source criticality arguments.

Finally T2=-M/(4N) is universally exact in the radical field under a
rational mate. Normalized trace makes it rational-exact. Its two
simple-root residues force M(1)=M(-1)=0. The value/jet conditions
and the full degree-four bound now give exactly

    M=lambda*u^2*(u^2-1),
    R=lambda*u^2-2p*lambda*u+e.

This is not already a zero-count contradiction:
M/N=lambda/u^2 really is exact. If lambda=0, the global row forces
L constant. The remaining analysis retains every parameter with
lambda!=0.

## 3. Generic geometric integrality is paid

The fibre calculation is on the genuine generic curve. For
lambda!=0, L is nonconstant. A nontrivial polynomial composition
F=h(K) must have outer t-degree two or four.

For an outer quadratic, completing its square produces a polynomial
K1 with F=K1^2+constant. Choose the leading sign so K1+H has
t-degree two. Then

    (K1-H)(K1+H)=L-constant

has right side of degree at most one in t. Therefore K1=H and
L is constant, a contradiction. For an outer quartic the inner
polynomial has t-degree one. Comparing leading coefficients forces
N to be a polynomial square up to a nonzero scalar, contrary to
its two simple roots. No assertion that nonconstant L alone excludes
composition is used.

I re-read the primary
[Arzhantsev--Petravchuk, Theorem 1 and Lemma 3, PDF pages 2 and 4](https://arxiv.org/pdf/math/0608157v2).
They identify noncomposition with closedness and relative algebraic
closure of the generated rational-function field. In characteristic
zero the resulting regular function-field extension gives a
geometrically integral generic fibre. This is the same cited route
already recovered by the all-finite 4+3+1 supplier. It does not assert
irreducibility of every special fibre.

## 4. The full local divisor, independently reconstructed

Let eta=omega/dF on the compact geometric generic curve with
F=zeta. On the finite chart put s=z-x^2. The cleared equation is

    E=(N+sP+s^2Q)^2+s^3(M+sR)-zeta*s^4,
    eta=s^2 du/E_s,

up to the consistent choice of orientation. At u=0 set s=u^2 Z.
The complete face is

    P_zeta(Z)=(-1+cZ+dZ^2)^2-lambda Z^3-zeta Z^4.

Its constant is one and its leading coefficient is d^2-zeta.
A repeated nonzero root must annihilate the zeta-independent
polynomial ZP0'-4P0, whose constant is minus four. It follows
that there are four simple nonzero roots over the generic constant
field. They give four unramified branches with parameter u and
exhaust local Weierstrass degree four. The relative differential
has exact order -2 at each: s^2 has order four and E_s order six.
Thus a hypothetical primitive has at most one simple pole at each
of four points. Denote their sum by E0.

At u=1 or -1, use the full numerator v=N+sP+s^2Q as a tangential
coordinate. Its tangential derivative is a unit. Since M vanishes
there, the generic leading equation is

    v^2+(constant-zeta)s^4+terms of greater weight=0.

Both branches have v of order two in s. E_v also has order two,
so eta is a unit. H=v/s^2 is bounded. This remains true if the
normal coefficient of the numerator vanishes: u-u0 then has order
two rather than one, with a nonzero generic coefficient. No simple
boundary point contributes a hidden zero or pole.

For the deleted infinity point I independently used the actual
coordinate

    r=1/x, s=z^{-1}-r^2,

and recovered the complete numerators, with N,P,Q,M,R evaluated at
u=1/r-p,

    N_inf=r^8N+s(2r^6N-r^4P)+s^2(r^4N-r^2P+Q),
    M_inf=-r^4M+s(R-r^2M),
    E_inf=N_inf^2+s^3M_inf-zeta*s^4.

Consequently N_inf(r,0) starts with r^2,
N_inf,s(0,0)=2-a, and M_inf(0,0)=-lambda. The latter sign has no
effect on the nonzero-unit argument.

When a!=2, solve N_inf=0 for its moving analytic centre s=psi(r).
Its leading term is -r^2/(2-a). To check that higher units have not
been dropped, set

    s=r^2[z0+r(z1+Y)],
    z0=-1/(2-a),
    z1=-(2ap+b)/(a-2)^2.

The independent expansion gives

    [r^3]N_inf=(2-a)Y,
    [r^6]E_inf=(2-a)^2Y^2-lambda*z0^3,
    [r^3](E_inf)_s=2(2-a)^2Y.

The two roots of the middle equation are simple and nonzero.
They produce two unramified branches, exact relative-form order
2+4-3=3 at each, and H has order -1 in r. The local degree in s
is two, so no third local branch has been omitted.

When a=2, the normal term has higher Newton weight. The complete
face under r=tau^3, s=tau^4 Z is

    [tau^12]E_inf=1-lambda Z^3,
    [tau^8](E_inf)_s=-3lambda Z^2.

Its primitive Newton edge has coprime exponents 3 and 4: the three
Puiseux determinations form one normalized branch. The displayed
nonzero coefficients exclude hidden splitting or cancellation.
The actual differential has order 6+8+2-8=8 and H has order -2
in tau. This includes the full a=2 stratum, not merely its generic
coefficients.

On the retained divisor D, the independently computed full values are

    H_D=(2-a)b_D-3ap^2+2bp-c+d+10p^2-1,
    L_D=-lambda*(b_D+3p^2-1).

For a!=2, F_D has degree two and two generic transverse points.
For a=2, its degree is one with slope -lambda and it has one
generic transverse point. At each such point the exact volume
r^2 dr wedge db_D gives eta order two. The deleted infinity
branches are distinct from these finite-b_D points.

At every other point of the smooth generic fibre in W, the source
volume and dF are nonzero units in local smooth coordinates; eta
is a unit. All boundary points have now been included. Its total
canonical degree is therefore

    a!=2: 2*3+2*2-4*2=2,
    a=2:  8+2-4*2=2.

Geometric integrality then gives genus exactly two in both regimes
for every retained parameter with lambda!=0.

## 5. The complete primitive space, including descent

Every rational primitive can have poles only in E0, with multiplicity
at most one. A pole elsewhere would give a pole of its derivative,
contrary to the complete divisor just computed. Since deg E0=4
and g=2, Riemann--Roch gives dimension exactly three.

The functions 1, 1/u and H/u all belong to L(E0):

* At the active points H is bounded and u is a parameter.
* The affine line u=0 is contained in the fixed fibre F=d^2, so
  no affine denominator remains on the generic fibre.
* At the simple boundary points u is nonzero and H is bounded.
* At every retained D point H is global and 1/u has a zero.
* At deleted infinity, H/u has orders zero and one in the two
  regimes respectively, by the actual H orders computed above.

These are all points, including all affine denominators. The
three functions are independent over the geometric generic
constant field: after multiplying a relation by u, its t-degree
is at most two, below the irreducible generic equation of degree
four. The t^2 coefficient first eliminates H; then 1 and u are
independent. Membership and exact dimension prove completeness,
not only the absence of a primitive in a chosen finite ansatz.

The basis is defined over C(zeta). For a rational source mate,
uniqueness of its expansion makes its coefficients Galois-fixed,
so they belong to C(zeta). Alternatively the coefficient
contradiction below works over the larger algebraic constant field:
the fibre derivation kills every element algebraic over C(zeta)
in characteristic zero. Thus no unproved constant-field descent is
required.

## 6. Independent bracket derivation

This part was reconstructed directly, without importing the
producer or using polynomial-remainder computation. Write M=L_t
and N=[t^2]H. Then

    J(F,1/u)=F_t/u^2,

whose t^3 coefficient is 4N^2/u^2=4u^6(u^2-1)^2. Also the exact
identity is

    J(F,H/u)
      =J(L,H)/u+2(zeta-L)H_t/u^2+HM/u^2
        +2(F-zeta)H_t/u^2.

The last term vanishes on the actual generic fibre. The remaining
polynomial has degree at most two in t. Its t^2 coefficient is

    (2NM'-MN')/u-3NM/u^2
      =lambda*u^4*(u^2-1)*(3-u^2).

Thus for a primitive C+A/u+B H/u, the t^3 coefficient forces A=0
and the t^2 coefficient forces B=0. The remaining constant has
zero derivative and cannot be a mate. The two factors are
nonzero polynomials for every permitted parameter, including all
p,a,b,c,d, so no specialization has been lost.

When lambda=0, L is constant. The identity
J(H^2+L,G)=2H J(H,G) excludes polynomial G in the original ring.
It does not exclude all rational G, as the literal positive
control demonstrates.

## 7. Sharp control and frozen replay

The family

    H=(u^2-1)(1+u^2t)^2+d,
    G_H=-1/[2u(1+u^2t)]

satisfies every original shifted global row for every finite p.
Its original source Jacobian is one; H^2+e has mate G_H/(2H).
This is in the exact binary 4+2+1+1 class, with the double point
at original infinity. It pays the polynomial/rational distinction
within the theorem's own class.

The source has an explicit symbolic universe: all free shifted
coefficients before entry, every parameter after entry, both
a!=2/a=2 regimes in the analytic proof, and the two universal
primitive coefficients. The numeric-looking genus/order gates
are arithmetic controls of the independently audited geometry;
they do not themselves prove branch exhaustion or generic genus.
No ambient-degree or mate-degree census is used. Every gate is
an always-active exception check.

Independent reproduction:

    python3 04-computation/planar_jc48_sep08_four_two_one_one.py
    python3 -O 04-computation/planar_jc48_sep08_four_two_one_one.py

Both fresh outputs equal the frozen output byte for byte:
**41 gates, 283 bytes**.

* Source: 4,581 bytes, SHA256
  be559b1eb0b15a8a55179d72e9defb01fb62cec8d4341f2b1a9c0f6ba1318805.
* Output SHA256:
  a6a4dcf34b277956677a86916a1d224becff92b514fc7837ca1ca0f896b95fc1.
* Semantic trace:
  6410484364b2dc0ab1c3412baed101976b764e4c4410567f7dc2957a0ee94818.

The separate fresh reconstruction checked thirteen identities in
the full local numerators, moving centre, Newton faces, actual
D restrictions, and direct bracket formulas. Its analytic content
and exact formulas are recorded above, so the acceptance does not
depend on retaining a temporary script.

**Final mathematical acceptance: PASS.** There is no remaining
scope or mathematical correction. The original source/output
remain untouched; the parent owns primary status and checkpointing.
