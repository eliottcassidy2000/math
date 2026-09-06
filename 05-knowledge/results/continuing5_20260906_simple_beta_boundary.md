# The full two-interlacer model has simple beta roots

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.** In the full weak anchored model, every beta root is simple.
Consequently its fixed-original-phase extrema occur only at a zero beta
root or a shared beta/interlacer root. Compactness gives a uniform positive
gap between distinct beta roots, without an explicit numerical constant.
The zero-root boundary admits an exact fixed-atom reduction.

This does not prove simplicity in the B-only model, strict interlacing, a
strict reference for every coefficient fibre, or the remaining response
sign. General anchored sign and Laurent noncancellation remain OPEN.

## 1. Inheritance and precise model

The incoming proved boundary-extremum theorem is
`05-knowledge/results/long_frontier2_sep06_boundary.md`, with its independent
audit. It reduces fixed-phase extrema to zero/repeated beta roots and shared
interlacer roots. The closest earlier mechanism is the weak residue-measure
representation and its exact Gram determinants in
`05-knowledge/results/continuing3_20260906_residue_floor.md`, Sections 1--2,
with `continuing3_20260906_residue_floor_audit.md`. That result also proves
the stronger floor y>161875/888583>9/50. The proof below uses only the exact
moment identities, so it does not depend on that numerical floor or its
tail calculation.

The canonical hostile is the C-only repeated configuration
(0,0,3,5,5). The corrected near miss is to regard every discriminant zero
as a possible full-model boundary. The underused sidecar is the ordinary
three-by-three Gram determinant of **D/B residue moments**, not Newton
moments of D. The ordered-resultant selector in
`continuing5_20260906_pencil_selector.md` retains the admissible shared-root
endpoints once a strict reference is supplied.

Let B have five nonnegative real roots, counted with multiplicity, whose
sum is13 and square sum59. Define x=e3, y=e4, z=e5 and

    B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
    C(v)=v^4-12v^3+45v^2-(2/3)xv+(3/7)y,
    D(v)=v^4-11v^3+36v^2-(5/12)xv+(1/7)y.

Assume both monic quartics C and D weakly interlace B in the usual ordered
root sense, allowing equalities. Denote this full model by K. All
coefficients x,y,z are nonnegative because they are elementary symmetric
functions of nonnegative roots.

The concept board is: repeated beta nodes; forced common factors; separate
positive residue measures; a low-order Gram obstruction; compact root
geometry; and zero-atom peeling. The map below keeps the actual repeated
root and both interlacers before eliminating coefficients. It preserves
the full predicate; its scalar determinant is used as a necessary condition,
never as an equivalent replacement for the geometry.

## 2. Three necessary moment determinants

For either A=C,D, weak interlacing forces at least multiplicity m-1 in A
at a beta root of multiplicity m. After cancellation, A/B therefore has
only simple poles. The ordered root signs give nonnegative residues, and
the monic degree difference gives total mass one. Thus

    C(v)/B(v)=sum_i w_i/(v-beta_i),
    D(v)/B(v)=sum_i h_i/(v-beta_i),
    w_i,h_i>=0, sum_i w_i=sum_i h_i=1.

Canceled nodes may have weight zero; the measures need not have the same
weights. All their ordinary moment matrices are PSD. Their shifted moment
matrices are PSD because every beta node is nonnegative. This argument
includes repeated roots directly and does not assume simplicity in order
to prove it.

Division at infinity gives

    mu0,...,mu5 =
      1, 1, 3, x/3-16, 16x/3-373-4y/7,
      54x-59y/7+z-3969,

    nu0,...,nu4 =
      1, 2, 7, 7x/12-19, 115x/12-632-6y/7.

Here mu are moments of C/B and nu are moments of D/B. The required exact
identities are

    det(mu_(i+j+1))_(i,j=0)^1 = (x-75)/3,                 (1)

    det(mu_(i+j+1))_(i,j=0)^2
      = -x(x-75)^2/27 +(15/7)(x-75)y
        -(16/49)y^2 +(x-75)z/3,                         (2)

    Delta_D := det(nu_(i+j))_(i,j=0)^2
      = -49x^2/144 +(269/4)x -(18/7)y-3132.             (3)

Each determinant is nonnegative in K. The source reconstructs the formal
series and all determinants rather than inserting them as its input.

## 3. Simplicity theorem

**Theorem.** Every B in K has five distinct roots.

First suppose zero is repeated. The constant and linear coefficients of B
give z=y=0. Equation (1) gives x>=75. Equation (2) then becomes
-x(x-75)^2/27>=0. Since x>=75>0, this forces x=75. But (3) equals
-37/16 at (x,y)=(75,0), a contradiction. This is a self-contained exclusion
of a repeated zero; the inherited strict y floor is a consistent stronger
consequence of the same measures.

Now suppose r>0 is repeated. Weak interlacing forces C(r)=D(r)=0.
These two equations are linear in x,y with coefficient determinant r/12,
so they imply exactly

    x=(24/7)r^3-36r^2+108r,
    y=3r^4-28r^3+63r^2.                                (4)

Substitution in the necessary equation B'(r)=0 gives

    B'(r)=(4/7)r^2(2r^2-14r+21)=0.

Therefore

    r=(7-sqrt(7))/2 or r=(7+sqrt(7))/2.                 (5)

Modulo 2r^2-14r+21=0, the coefficients and determinant reduce to

    x=126-12r,
    y=735/4-49r,
    Delta_D=5r-75/4.                                   (6)

At the upper root in (5), y=49/4-(49/2)sqrt(7)<0, contradicting
nonnegative beta geometry. At the lower root,

    Delta_D=-5/4-(5/2)sqrt(7)<0,

contradicting the ordinary D/B Gram condition. These exhaust the possible
positive repeated roots. Together with the zero case, this proves the
theorem, regardless of any proposed higher multiplicity.

The C-only hostile is retained exactly:

    B(v)=v^2(v-3)(v-5)^2,      (x,y,z)=(75,0,0),
    C(v)/B(v)=2/(3v)+1/(3(v-3)).

C is a lawful weak interlacer with positive residue measure, but
D(5)=-25/4 and Delta_D=-37/16. Thus dropping the second interlacer really
does destroy the simplicity conclusion; it is not merely an omitted proof
hypothesis.

## 4. Uniform mutual separation and the reduced extremum set

The ordered beta vectors lie in the compact set

    0<=beta1<=...<=beta5<=13,
    sum beta_i=13, sum beta_i^2=59.

The constraints that the associated monic quartics weakly interlace are
closed: their coefficients depend continuously on the beta vector, and
the ordered real root inequalities persist under limits. Hence K, viewed
as a set of ordered beta vectors, is compact. It is nonempty; exact controls
below verify both (84,35,0) and (84,35,1).

The continuous function

    gap(beta)=min_(1<=i<5)(beta_(i+1)-beta_i)

is strictly positive everywhere in K by the theorem. Its minimum is
therefore some eta>0. **There exists a uniform positive mutual beta-root
gap on the whole full model.** No explicit value or effective certificate
for eta is supplied. This does not assert distance from zero or separation
between beta and interlacer roots. In particular, (84,35,0) is an admissible
simple-zero example.

Now fix s>0 and keep the original equation P(-s)=0. For a nonempty full
model slice K_CD(s), the incoming proved theorem places every maximum and
minimum of Q(-s) at a zero/repeated beta root or a shared C/D root. The
simplicity theorem deletes only the repeated-root alternative. Consequently,
with

    E_0(s)={ (x,y,z) in K_CD(s): z=0 or Res(B,C)=0
                                     or Res(B,D)=0 },

E_0(s) is nonempty and

    max_(K_CD(s)) Q(-s)=max_(E_0(s)) Q(-s),
    min_(K_CD(s)) Q(-s)=min_(E_0(s)) Q(-s).              (7)

Every maximizer and minimizer lies in E_0(s). The incoming constructive
fixed-phase move likewise terminates on this reduced set. Because B is
simple, a shared interlacer root is exactly disappearance of that node's
weight in its residue measure; it is not a collision of beta nodes.

Neither resultant equality alone nor its sign imposes the full predicate.
In a fibre supplied with a simultaneous strict reference, the ordered
middle-root selector intersects both PSD intervals and z>=0, rejecting
strictly exterior determinant zeros. Fibres without such a reference,
including potentially collapsed fibres, still require the full degree-eight
PSD decoder. Simplicity does not supply a strict reference in every fibre.

Equation (7) is a boundary reduction, not the missing sign on that boundary.
For the B-only model the repeated-root alternatives remain necessary.

## 5. The zero-root boundary peels off two fixed positive atoms

In K with z=0, simplicity makes zero a simple root. Thus y>0 and

    B(v)=v G(v),
    G(v)=v^4-13v^3+55v^2-xv+y

has four strictly positive simple roots. Set

    Ctilde(v)=v^3-(45/4)v^2+(75/2)v-(5/12)x,
    Dtilde(v)=v^3-(32/3)v^2+(197/6)v-(23/72)x.

Direct coefficient identities give

    C(v)/B(v)=3/(7v)+(4/7) Ctilde(v)/G(v),
    D(v)/B(v)=1/(7v)+(6/7) Dtilde(v)/G(v).              (8)

The remaining measures are positive and normalized: the zero weights are
exactly3/7 and1/7, leaving masses4/7 and6/7. Thus the zero-root full model
reduces faithfully to the quartic G and its two cubic weak interlacers.
Conversely, if G has four positive roots and both displayed monic cubics
weakly interlace it, their positive residue representations combined with
the fixed zero atoms prove that C,D weakly interlace vG. The anchors remain
the sum13 and square sum59 of those four positive roots. Multiplicities may
be allowed when formulating this converse; the proved full-model theorem
then excludes them.

This reduction removes a node in the residue representation, not a row
from the original Laurent normalization. Keep exactly

    beta=t^-1(1+13t+55t^2+xt^3+yt^4+zt^5),
    C_raw=t^-1(1+12t+45t^2+(2x/3)t^3+(3y/7)t^4),
    D_raw=t^-1(1+11t+36t^2+(5x/12)t^3+(y/7)t^4),
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).

In particular q_(-1)=28. At z=0 the original phase equation becomes

    P(-s)=182-20020s+2002x s^2-3432y s^3=0.             (9)

One must not replace O14,E14 by lower-height binomial carriers, omit the
mixed response, or treat every quartic model point as an actual factorial
row. Formula (8) is a lower-dimensional geometry packet; a uniform sign
proof on (9) remains a next test.

## 6. Exact controls, reproduction and stopping point

The source uses no repository or producer imports. It reconstructs both
formal residue sequences through moment8, checks their products with the
original denominator, all three Gram identities, the complete forced
repeated-root equations, and both exact quadratic candidates. It retains
the C-only hostile, proves both peeling identities by polynomial equality,
and independently rebuilds the original P,Q carriers and inverse carry.

For each of (x,y,z)=(84,35,0),(84,35,1), exact Sturm counts and square-free
tests verify the beta roots. Sylvester's criterion on the full ordinary
five-by-five residue matrices gives the leading principal minors

| Coefficients | C/B minors | D/B minors |
|---|---|---|
| (84,35,0) | 1,2,11,188,7125 | 1,3,26,705,30600 |
| (84,35,1) | 1,2,11,170,6069 | 1,3,26,646,24617 |

For these simple real beta roots, Vandermonde congruence identifies the
strictly positive residue weights, verifying both full weak-model controls
without a numerical root test. In particular the zero-root boundary is
nonvacuous and a positive distance of all beta roots from zero would be false.

    python continuing5_20260906_simple_beta_boundary.py
    python -O continuing5_20260906_simple_beta_boundary.py

Both runs pass **89 always-active exact gates**. The source configures
stdout LF explicitly. Normal and optimized subprocess stdout were captured
as raw bytes, are identical and contain no carriage returns. No text
normalization was used to claim byte identity.

- Source SHA256: `979a4540c112b049e032bd1320ed78e7991a5996b0a5a4991cd88e185bbcc185`.
- Output SHA256: `8f65fb22877e43128560c0dcc3f4fb4c346111adb0d2175324a72cff4a73b7d6`.

The all-model simplicity and compactness conclusions are analytic, not
inferred from these two controls. The precise next objects are the fixed
zero-atom boundary and admissible shared-root endpoints; no general sign,
strict-reference or LRC consequence is claimed. All artifacts were created
outside the repository, with no shared file or Git changes.

Independent [proof and exact referee](continuing5_20260906_simple_beta_boundary_audit.md) passes.

The [raw-byte checkpoint manifest](continuing5_20260906_manifest.json) pins
the filed source, report and identical normal/optimized transcript. Any
candidate-report hashes above identify the pre-promotion bytes.
