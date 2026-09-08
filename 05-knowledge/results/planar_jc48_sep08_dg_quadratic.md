# The DG source filtration and recovered exclusion through cubic degree

**Status: PROVED FILTRATION + RECOVERED COROLLARY + INDEPENDENT AUDIT PASS.**
The [independent audit](planar_jc48_sep08_dg_quadratic_audit.md) accepts the
filtration, full low-degree dependencies, and their extension obstruction.
The filtration calculation concerns the specified explicit surface. The
quadratic and cubic Keller rigidity used below is recovered proved canon, not a new
planar Jacobian theorem. JC(2) remains OPEN.

## 1. Inheritance and the actual object

Work over C on the audited surface

    W=(P1_x x P1_z) minus {z=x^2},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^2-r^4 b, D={r=0},
    omega=dx wedge dt=r^2 dr wedge db.

The two source charts and their transition are from the
[explicit DG surface](planar_jc48_sep06_dg_surface.md). The incoming
[complete source-linear carrier exclusion](continuing10_20260907_dg_linear_carrier.md)
classifies all global functions of source-t degree at most one and excludes
every regular constant-Jacobian mate, with no degree bound on that mate.
Its critical-free rational-primitive example is the retained hostile:
neither boundary separation nor a generic A1 fibre pays global regularity.

Before attempting a quadratic hyperelliptic calculation, the inheritance
search recovered [THM-2071, quadratic-fibre rigidity](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md).
That theorem already proves every planar Keller pair with one coordinate
quadratic in one linear source variable is a polynomial automorphism.
Its centered parity argument and central-binomial noncancellation close
all degrees of the other coordinate. The companion
[THM-2063, one-fibre-linear pairs](../../01-canon/theorems/THM-2063-one-fiber-linear-planar-keller-pairs.md)
covers lower degree. A second inheritance search then recovered
[THM-2118, all-degree cubic Faber closure](../../01-canon/theorems/THM-2118-all-degree-cubic-faber-boundary-flux-coprimality.md).
Its all-degree boundary/flux coprimality removes the centering poles;
the remaining coefficient poles and polynomial leading coefficient are
then controlled before applying
[THM-2102, power-free weighted faces](../../01-canon/theorems/THM-2102-power-free-weight-face-and-first-defect-descent.md).
Thus the cubic layer is also already closed in the plane. These are the closest proved mechanisms, with their
full slugs retained to avoid numeric-ID collisions.

The corrected near miss would be to start a new search in a carrier already
excluded by that canon. The least-used sidecar is the pole of the original
coordinate x along D, which turns the recovered plane automorphism into
an obstruction to a pair of global functions on W.

## 2. The complete filtration in every degree

Define L_n={F in O(W): deg_t(F|U0)<=n}, for each integer n>=0.
Then a basis is

    E_(a,j)=x^a t^(n-j)(1+x^2 t)^j,
    0<=a<=2n, 0<=j<=n.                                  (1)

In particular dim L_n=(2n+1)(n+1), and the actual second-chart formula is

    E_(a,j)=(-1)^n r^(2n-a) b^j(1+r^2 b)^(n-j).          (2)

Formula (2) proves regularity everywhere on the boundary chart, including
the locus 1+r^2 b=0; no rational denominator has been discarded.

For completeness of the basis, use the original second projective
coordinate z=x^2+1/t on the dense overlap. Every polynomial F of t-degree
at most n has a unique representation

    F=Q(x,z)/(z-x^2)^n,    Q in C[x,z], deg_z Q<=n.       (3)

Indeed Q is obtained by replacing t by 1/(z-x^2) and clearing that fixed
denominator. Near the generic point of D, r=1/x and z are regular
coordinates and

    F=r^(2n)Q(1/r,z)/(r^2 z-1)^n.                       (4)

The denominator is a unit along D. Thus regularity forces deg_x Q<=2n:
the nonzero leading coefficient, a polynomial in z, cannot vanish at the
generic point of D. Conversely all monomials with these two degree bounds
give precisely (1), whose complete regularity was already checked by (2).
Independence follows from independence of x^a z^j in C[x,z]. This proves
the full space, not a selected span of native generators.

The filtration respects multiplication, L_m L_n subset L_(m+n).
Here L_m L_n denotes the linear span of pairwise products; in fact equality
holds, since the bidegree boxes in (3) add and every monomial in the larger
box splits into one monomial in each smaller box. The denominator powers
also add. No statement about associated-graded multiplication is needed.
The union of the L_n is O(W), because restriction of every global function
to the affine chart U0 is an ordinary polynomial in x,t.

## 3. The fifteen-dimensional quadratic layer explicitly

Write F=A(x)t^2+B(x)t+C(x), and write A_i,B_i,C_i for coefficients of
x^i. The complete globality criterion is

    deg A<=8,
    B=2A_8 x^6+2A_7 x^5+sum_(i=0)^4 B_i x^i,
    C=A_8 x^4+A_7 x^3+(B_4-A_6)x^2+(B_3-A_5)x+C_0.   (5)

All nine coefficients A_0,...,A_8, the five coefficients B_0,...,B_4,
and C_0 are free. One direct proof is to retain separately the b^2, b,
and constant coefficients after substituting t=-r^2-r^4b. The first
gives deg A<=8, the second the two forced top coefficients of B, and the
third the four displayed coefficients of C. The converse is literal
Laurent cancellation. This is also the n=2 case of (1)-(4).

The full boundary restriction, a useful coordinate retained for later
searches, is

    F|D=A_8 b^2+(2A_6-B_4)b+(A_4-B_2+C_0).             (6)

Thus the quadratic layer contains boundary separators and is larger than
the old collapsed carrier. Its exclusion below does not follow from
boundary collapse or absence of critical-free members.

## 4. Recovered exclusion through cubic degree with an unrestricted mate

**Corollary.** If F,G in O(W) satisfy dF wedge dG=lambda omega with
lambda!=0, no nonconstant member alpha F+beta G of their constant output
pencil belongs to L_3. More generally this conclusion holds in every
linear source direction after an affine source change on U0.

Suppose otherwise. Complete the nonzero output direction (alpha,beta) to
an invertible constant target matrix. On U0 the resulting polynomial
pair has nonzero constant Jacobian and one coordinate of t-degree at most
three. THM-2063, THM-2071 or THM-2118, in the exact forms linked above, makes that pair
a polynomial automorphism of A2. Undoing the target matrix preserves that
conclusion. Hence its polynomial inverse expresses x as a polynomial in
F and G. Since F and G are global on W, that expression is global on W.
On the dense U0 it equals x, so it is the rational function x everywhere.
But x=1/r has a genuine simple pole along the nonempty divisor D. This
contradiction proves the corollary. An affine source change leaves the
original x a polynomial in its two new coordinates, so the identical
inverse argument proves the final statement.

The same proof applies to any variety containing this affine plane as a
dense open if some original polynomial coordinate fails to extend
regularly: a plane automorphism whose two components extend globally
would extend every original coordinate. This elementary extension
observation is separate from any classification of Keller maps.

The scope is exact. Neither F nor G is required to have bounded degree
except for the one specified pencil member. Functions such as the global
t itself do have polynomial mates on U0, for example (t,-x); the mate's
pole on D is essential. Products and powers of the linear generators
populate L_2 but do not create a globally regular pair. Quartic and higher
source degree, and the identification of W with an actual finite Keller
envelope, remain OPEN here. This result does not strengthen THM-2071's
already proved planar scope or claim literature priority.

## 5. Concept comparison and reproduction

| Concept | Map and preserved predicate | Loss, sidecar, and next test |
|---|---|---|
| Moving source | Actual restriction O(W) to C[x,t] retains polynomiality and bracket | A completed carrier flow need not land in this global ring |
| Changed global carrier | F maps to Q in (3), retaining the exact two degree bounds | Chart denominators restored by the polynomial expression (2) |
| Infinity and meridians | A new coordinate pair would give actual covering data | This ring filtration alone gives no whole-support passport |
| Collision and torsion | Polynomial inversion would restore the missing x coordinate | Generic rational primitives still require every special-fibre principal part |
| Next search | L_3 has dimension28 and is excluded; L_4 has dimension45 by (1) | Any candidate needs two global coordinates and the unrestricted mate equation |

The verifier tests the actual two-chart identities, complete independent
Laurent-kernel dimensions through n=3, the fifteen-parameter formula and
its boundary restriction, and the positive affine-plane/negative boundary
pole control. These finite checks support the displayed identities, not
the unbounded theorem of THM-2071. Its proof is an inherited dependency.

    python3 -B 04-computation/planar_jc48_sep08_dg_quadratic.py
    python3 -B -O 04-computation/planar_jc48_sep08_dg_quadratic.py

The independent audit checks the full all-degree inherited proofs, including
THM-2118's cubic pole regimes and THM-2102's direct weighted descent.
It also replays this source normally and under optimization:122 always-active
gates with the same409-byte output. The source's mention of the quadratic
dependency records the earlier recovered layer; the cubic consequence is
analytic and is not inferred from those finite gates. Frozen SHA256:

    source b170cb02793e5eea5423aa1a838ee296406e993c18b26c97a033d15a01f9d568
    output 98dbafefbba323c2401b5af303de6ca4c5b2fa8ac37cb28b04acb7baca0761a8

A targeted quartic inheritance search found monic and finite-pole subcases,
but no whole one-coordinate, unrestricted-mate quartic theorem. The results
with both coordinates of degree at most four or five impose a second
coordinate bound and are not substituted for that missing theorem. Thus
L_4 is the first unresolved layer in this specific recovered route; some
of its subfamilies may already be excluded by the named prior work.
