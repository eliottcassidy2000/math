# Sharp genus and rational mates on every graph-complement surface W_m

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The proof and exact controls below concern the specified surfaces and their
genuine quadratic source functions.  The polynomial-mate conclusion is a
recovered corollary of proved planar quadratic rigidity.  It is not a new
plane Jacobian theorem, and JC(2) remains OPEN.  Source/output are frozen
with the final replay pins below.

## 1. Statements, inheritance, and scope

Work over C.  For every integer m>=1 put

    W_m=(P1_x x P1_z) minus {z=x^m},
    t=1/(z-x^m),       omega=dx wedge dt.

Let H be global on W_m and have **degree exactly two** in t on the
source chart A2_(x,t).  Write

    H=N(x)t^2+P(x)t+Q(x),       N!=0.

A rational mate means G in C(x,t) with J_(x,t)(H,G)=1.  Any nonzero
constant Jacobian is reduced to this normalization by scaling G.  No
regularity of G on W_m or degree bound on G is assumed.

**Genus bound.** If H has a rational mate, every geometric generic
component of H has genus at most m-1.  This bound is attained for every
m, by actual global H with an explicit rational mate.

There is also a complete equality criterion when m>=2.  Express H in
its unique global numerator form

    H=[A(x)z^2+B(x)z+C(x)]/(z-x^m)^2,
    deg A, deg B, deg C <=2m,
    D0=B^2-4AC.

Then the following are equivalent:

1. H has a rational mate and its geometric generic component has genus
   m-1.
2. For one fixed p in C, with u=x-p, the **whole pencil** satisfies

       D0+4lambda N
          =u^(2m+1)[a(lambda)u^(2m-1)+b(lambda)],     (1)

   where a,b are affine polynomials in lambda, neither identically zero.

Here and below generic means outside a finite set of fibre values.  In
(1) both a(lambda) and b(lambda) are nonzero generically.  Proportional
D0 and N are allowed.  Condition (1) also supplies a rational mate:

    G=2(2Nt+P)/[(2m-1)b(H)(x-p)^(2m)].                (2)

The equality criterion is not asserted for m=1, where square genus-zero
fibres introduce additional cases.  Nor is a root translation asserted
to extend to an automorphism of W_m: u=x-p is a coordinate in the
radical-curve calculation.  Criterion (1) is tested on the actual global
H, and its sufficiency follows from the literal source formula (2).

**Recovered polynomial exclusion.** No such genuine quadratic H has a
polynomial mate in C[x,t], even if the proposed mate need not be global
on W_m.  This is stronger than excluding a pair of global functions,
but follows directly from
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and the explicit global coefficient space below.  No new rigidity
theorem for arbitrary plane quadratic functions is claimed.

The closest geometric mechanism is the complete m=2 filtration in
[dg_quadratic](planar_jc48_sep08_dg_quadratic.md); its proof extends to
the actual degree-m graph.  The closest genus mechanism is the now
proved and independently audited
[genus_capacity theorem](planar_jc48_sep08_genus_capacity.md), itself
connected to
[the exact degree-eight classification](planar_jc48_sep08_boundary_exactness.md).
The same-partition rational hostile in
[five_one_one_one, Section 6](planar_jc48_sep08_five_one_one_one.md)
is the m=2 member of the sharp family below.  The corrected near miss
would be to import the m=2 elliptic bound without retaining m, or to
exclude rational mates from the polynomial theorem.  The least-used
sidecars are the **whole discriminant pencil**, the degree of the
actual field map, and the mate's poles.

The five live concepts are global section spaces, primitive pole degree,
branch count, fixed repeated roots in a pencil, and rational versus
polynomial mates.  The map from H to its discriminant retains the generic
fibre and the relative form; just retaining its degree loses all position
and exactness data.  The sparse equality calculation restores those data
and, unusually, supplies the actual rational mate.  The direct source
Jacobian of (2) is the cheapest decisive test of this extra connection.

## 2. The actual charts and the complete filtration

There are two affine charts covering W_m:

    U0=A2_(x,t),
    Uinf=A2_(r,b),
    x=1/r,       t=-r^m-r^(2m)b,
    b=-x^m-x^(2m)t,
    D=W_m minus U0={r=0}.

On the projective-coordinate overlap,

    b=z/(1-r^m z),       z=b/(1+r^m b).

These formulas include the loci t=0 and 1+r^m b=0 in their respective
affine charts; their regularity is not inferred only on a smaller torus.
The original source two-form is

    omega=r^(2m-2) dr wedge db.                       (3)

In particular it is nonvanishing at D when m=1.  One must not import the
m=2 vanishing order two into this case.

For every integer n>=0 define

    L_n={F in O(W_m): deg_t(F|U0)<=n}.

The complete basis is

    E_(i,j)=x^i t^(n-j)(1+x^m t)^j,
    0<=i<=mn,       0<=j<=n.                          (4)

Thus dim L_n=(mn+1)(n+1).  Its full second-chart expression is

    E_(i,j)=(-1)^n r^(mn-i)b^j(1+r^m b)^(n-j),         (5)

which proves global regularity, including at every boundary point.

For completeness, any polynomial of t-degree at most n has a unique
representation

    F=R(x,z)/(z-x^m)^n,       deg_z R<=n.

Near the generic point of D this becomes

    F=r^(mn)R(1/r,z)/(r^m z-1)^n.

The denominator is a unit there.  The leading x-coefficient of R is a
nonzero polynomial in z, so generic regularity forces deg_x R<=mn.
Conversely the monomials in this full degree box give (4) and satisfy
(5).  Their independence is that of the numerator monomials x^i z^j.
This proves completeness, not only an exhibited family of sections.
Also the union of the L_n is O(W_m), since every global function
restricts to an ordinary polynomial on U0.

For n=2 the numerator notation in Section 1 consequently gives

    N=A x^(2m)+B x^m+C,
    P=2A x^m+B,
    Q=A,
    P^2-4NQ=D0=B^2-4AC.                              (6)

All three numerator polynomials have degree at most 2m, so

    deg N<=4m,       deg D0<=4m.                      (7)

The three coefficient boxes are independent; their total dimension is
3(2m+1).  They include lower-degree functions too, which is why the
genuine-degree hypothesis N!=0 is kept separate.

## 3. Generic fibres and the genus bound

Put V=2Nt+P.  Over the x-line, the generic fibre H=lambda has the
actual quadratic presentation

    V^2=Delta_lambda(x):=D0(x)+4lambda N(x).            (8)

Away from N=0, the inverse is t=(V-P)/(2N); hence this is a birational
description of each relevant generic component.  Values at which a
vertical line x=constant is a whole fibre component are finite and are
discarded.  Indeed such a line requires N=P=0 and Q=lambda at one of
the finitely many zeros of N.  Thus (8) accounts for every geometric
generic component, not only one selected branch.

On a fibre, the relative differential of omega is

    eta=-dx/V.                                       (9)

The equality J(H,G)=1 says dG=eta on every generic component where G
is defined.  There are only finitely many components on which a
denominator of a fixed rational G can vanish identically, since that
denominator has finitely many irreducible factors.  These cause only
finitely many exceptional fibre values.  There is no hypothesis that G
is regular on either affine chart.

If Delta_lambda is square in C[x], each of its two geometric components
is rational.  Their union is not treated as a connected double cover.
Otherwise its normalized curve is connected and has the radical field
C(x)(sqrt(Delta_lambda)).  In this field dx/V is exact by (9).

Here is the short capacity argument from the cited genus supplier.  Let
F be a nonsquare polynomial of degree n and suppose dx/sqrt(F) is exact.
A finite double root gives nonzero residues and is forbidden.  A root of
multiplicity j>=3 allows total primitive pole degree j-2: for odd j
there is one normalized point, and for even j there are two points each
allowing degree j/2-1.  Simple roots are regular.  Thus the finite
primitive pole budget is

    P_fin=sum_(j>=3)(j-2).                            (10)

For odd n>=3 there is one infinity point, with differential order n-3.
The primitive's local degree there is n-2, so P_fin>=n-2.  A pattern
with two or more roots has strictly smaller budget; the polynomial must
be a pure power and the genus is zero.

For even n>=4 the two infinity points each force local degree n/2-1.
One uses either local degree, not their sum, since the primitive may
have different values there.  Write s for the number of simple roots,
h_o for the number of odd high roots, and h_e for the number of even
high roots.  The connected branch count is B=s+h_o=2g+2 and

    P_fin=n-s-2h_o-2h_e>=n/2-1,
    B+h_o+2h_e<=n/2+1.

Some high root is required; hence B<=n/2 and

    g<=floor((n-4)/4).                               (11)

Degrees zero and one have genus zero.  Degree two is never exact,
because of a finite double root or the two nonzero residues at infinity.
Square fields have genus zero separately.

By (7)-(8), n<=4m.  The odd-degree and small-degree cases have genus
zero, while (11) gives g<=m-1 in every remaining case.  This proves
the genus bound, including m=1, and every component qualification above.

## 4. Equality, fixed roots, and the rational mate

Assume m>=2 and equality g=m-1.  The generic discriminant must have
degree exactly 4m: any smaller even degree gives at most m-2 in (11),
and odd degree gives genus zero.  Equality in the branch-budget argument
then forces

    h_o=1,       h_e=0,       s=2m-1.

Thus Delta has a unique high root p of multiplicity 2m+1, and all its
other 2m-1 roots are simple.  The full extremal theorem in
[genus_capacity, Section 3](planar_jc48_sep08_genus_capacity.md)
adds the position equations.  For clarity, its complete primitive-space
argument is as follows.  Write u=x-p and

    Delta=u^(2m+1)D(u),       deg D=2m-1,
    D squarefree,       D(0)!=0.

With X=1/u and Y=V/u^(2m), the actual normalized model is

    Y^2=T(X)=X^(2m-1)D(1/X),
    dx/V=-X^(2m-2)dX/Y.

It has a unique infinity point, and the differential has its only pole
there, of order 2m.  Its primitive is regular on the smooth affine
curve and has pole order at most 2m-1.  The complete function space is

    L((2m-1)infinity)=span{1,X,...,X^(m-1),Y}.

Indeed its normal affine ring consists of A(X)+YB(X), and X,Y have
pole orders 2,2m-1.  Even and odd leading pole orders cannot cancel.
Taking the odd part under Y->-Y therefore forces a primitive proportional
to Y.  Its derivative forces T'=constant*X^(2m-2), so

    Delta=u^(2m+1)(a u^(2m-1)+b),       ab!=0.          (12)

This derives the coefficient condition as well as the multiplicities.

It remains to show that the root p is fixed across the whole pencil.
That statement is not obtained merely by naming the high root separately
on each fibre.  Set G0=gcd(N,D0), N=G0 N1, D0=G0 D1.  If N,D0 are
not proportional, N1,D1 are coprime and

    W=D1' N1-D1 N1'!=0.

For a repeated root of D1+4lambda N1, coprimality implies N1!=0 there,
and the root must be a zero of the fixed polynomial W.  Each such root
determines at most one lambda.  For the finitely many roots of G0,
the residual factor is nonzero generically, again by coprimality.
Consequently all generic repeated roots of Delta_lambda are exactly
the fixed roots of G0, with their fixed multiplicities.  The unique
root of multiplicity 2m+1 is therefore one fixed p.  If N,D0 are
proportional, Delta_lambda is a scalar multiple of N and the same
conclusion is immediate.  This case includes D0=0 and is not removed
by a division by W.

For this p, (12) holds generically.  Two distinct generic fibre values
now show that N and D0 themselves belong to the vector space

    span{(x-p)^(4m),(x-p)^(2m+1)}.

This proves the entire-pencil condition (1); a(lambda),b(lambda) are
affine and both generically nonzero.  The repeated-root argument also
shows that generic branch counts are stable outside a finite set.

Conversely assume (1).  Its 2m-1 nonzero residual roots are simple,
and its high root has odd multiplicity 2m+1.  Thus the connected
double cover has 2m branch points and genus m-1.  More than genus is
recovered: on H=lambda,

    d[2V/((2m-1)b(lambda)u^(2m))]=-dx/V.

This follows by differentiating (1); equivalently, for
R=2/((2m-1)b u^(2m)), one has

    Delta R'+Delta' R/2=-1.

Substitute lambda=H to obtain (2).  The derivatives of b(H) in a tangent
direction of H vanish, so (2) has dG=eta on generic fibres.  Therefore
J(H,G)=1 as a rational identity on the source.  Its denominator is not
identically zero: b is a nonzero affine polynomial, H is nonconstant,
and u is a nonzero rational function.  This proves the full iff without
assuming global regularity of G, or assuming that a scalar fibrewise
primitive automatically extends without a formula.

An independent source control uses coordinates (x,V): here
H=(V^2-D0)/(4N), and the source Jacobian is 2N times the (x,V)
Jacobian.  It directly verifies (2) with all four coefficients of
N,D0 in the displayed two-dimensional space left free.

## 5. Actual sharp families and the field-degree sidecar

For every m>=1, delta!=0, and arbitrary beta,q in C, put

    h=x^m+x^(2m)t,
    H=(x^(2m)+delta*x)(1+x^m t)^2+beta*h+q.            (13)

These are actual global functions.  In the full second chart,

    h=-b,
    H=(1+delta*r^(2m-1))b^2-beta*b+q.                 (14)

No point in that chart has been deleted to make the formulas regular.
The numerator coefficients are

    A=x^(2m)+delta*x+beta*x^m+q,
    B=-beta*x^(2m)-2q*x^m,
    C=q*x^(2m),

all of degree at most 2m.  The source leading coefficient is

    N=x^(4m)+delta*x^(2m+1)!=0.

Let e=2m-1 and u=x^(-e).  Then

    J_(x,t)(u,h)=-e,
    omega=-du wedge dh/e,
    H=(1+delta*u)h^2+beta*h+q.

It follows directly that

    G_sharp=1/(e*delta*h),       J(H,G_sharp)=1.       (15)

The field map is not silently treated as birational.  One has
C(x,t)=C(x,h), and, with v=1/x,

    v^e=u.

Over C(h,u), this has degree e: the polynomial v^e-u is Eisenstein
at u in C(h)[u].  This includes composite e, not only prime exponents.
For m=1 the map is birational.  For m>=2 it has degree 2m-1 and is
not a source automorphism.  The rational mate (15) has the declared
pole along h=0, consistent with the polynomial exclusion below.

The complete discriminant pencil is

    Delta_lambda=(4lambda+beta^2-4q)x^(4m)
                    +4delta(lambda-q)x^(2m+1).       (16)

Both coefficients are nonzero generically.  Hence (16) has one odd root
of multiplicity 2m+1 and 2m-1 simple residual roots, giving genus m-1.
This includes the proportional case beta=0.  It proves sharpness by
actual global functions and actual rational mates, for every integer m.
For example m=3 already supplies genus two, refuting a uniform elliptic
bound without the surface parameter.  The m=2 case is the inherited
five-one-one-one rational control, now placed in its full family.

## 6. Recovered exclusion of every polynomial mate

Suppose instead that a polynomial G in C[x,t] has nonzero constant
Jacobian with the genuine quadratic H.  The exact unrestricted-mate
statement of
[THM-2071, Sections 3 and 4](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
forces

    N=N0 in C*,
    Q-P^2/(4N0) affine in x with nonzero slope.        (17)

No bound on deg_t G or deg G enters that theorem.  We only need to
check what (17) means in the complete global coefficient box (6).

From

    A x^(2m)+B x^m+C=N0,
    deg A,deg B,deg C<=2m,

one first obtains deg A<=m: if deg A>m, its leading term in
A x^(2m) has degree greater than 3m and cannot cancel either other
term.  Put E=A x^m+B.  Since C=N0-x^m E has degree at most 2m,
one has deg E<=m.  Thus the full remaining form is

    B=-A x^m+E,
    C=N0-x^m E,
    P=A x^m+E,
    Q=A,
    deg A,deg E<=m.                                 (18)

If deg A=k>=1, then deg P=m+k>m, with no possible cancellation by E.
Consequently Q-P^2/(4N0) has degree 2(m+k)>1.  This contradicts (17).
If A is constant, Q is constant.  A nonconstant P then gives centered
degree 2 deg P>=2; a constant P gives centered degree zero.  Neither
case has nonzero affine slope.  This proves the polynomial exclusion.
The reasoning includes m=1: a nonconstant A would give deg P=2,
and the constant-A alternatives remain the same.

This is a recovered corollary of the plane theorem, not a replacement
proof of it.  The globality of H provides the incompatible coefficient
space, while G is entirely unrestricted except for polynomiality and
its constant Jacobian.  The genuine-degree and polynomiality boundaries
are both sharp:

* The global degree-one function t has polynomial mate -x.
* The genuine global quadratic t^2 has rational mate -x/(2t), although
  it has no polynomial mate.  Its geometric generic components are
  rational, so it also controls the split case in Section 3.
* Family (13) has rational mates and maximal genus for every m.

No exclusion for arbitrary rational mates is claimed; their complete
extremal family is part of the positive result above.  No assertion is
made that a generic low-genus discriminant already supplies a mate.
Exactness and the whole-pencil equations are essential.

## 7. Exact controls and audit boundary

The standalone source imports no inherited mathematical implementation.
Its complete declared finite universes are:

* Every basis monomial in every box m=1..6, n=0..3, with literal full
  chart substitution.
* Independently, all source monomials x^i t^j with 0<=i<=4m and
  0<=j<=2 for m=1..4, retaining every negative Laurent coefficient.
  Kernel dimension, basis rank, and kernel containment are all checked.
* Actual sharp families for m=1..8 with free delta,beta,q, including
  both source coefficients and boundary charts, the cleared rational
  Jacobian, full discriminant pencil, and the proportional beta=0 case.
* The rational mate for the entire four-coefficient sparse pencil,
  m=2..6, checked directly in (x,V), together with the radical primitive
  and residual-binomial root-separation identities.
* Every post-cancellation degree pattern in (18) for m=1..12, together
  with its full symbolic coefficient cancellation.

The finite tests do not prove uniformity in m or all-polynomial
exactness.  Those are paid by the arguments in Sections 2-6.  No generic
smoothness or exactness is inferred from a numerical sample.  Positive
and hostile controls retain lower-degree polynomial mates, genuine
quadratic rational mates, a special repeated-root pencil fibre that is
not generic, the proportional-pencil boundary, and a genus-two exact
curve.  Missing intermediate coefficients in (12) are paid by the
complete primitive space, not by a multiplicity-only inference.

Reproduce from the repository root:

    python3 04-computation/planar_jc48_sep08_dg_genus.py
    python3 -O 04-computation/planar_jc48_sep08_dg_genus.py

Every gate raises an explicit exception on failure; none uses Python
assert.  The producer's normal, optimized, and frozen output are
byte-identical: **1,550 gates**, 856 output bytes.  Frozen SHA256 pins:

    source (9,148 bytes):
      45f25cbe2dea8486f8988136a23d11d386672ac9a6445d9110c929711fd8cff8
    output (856 bytes):
      386c8905a444b7184dffdff3dfc1694f841a8e7b27f3504026dc02875fd2d9f9
    semantic trace:
      c0e665ed29f7969341556b17978f101f4db88223d0ddc9bb5c4ba3691a0d5f9b

The finite controls retain 480 literal basis identities, four complete
quadratic Laurent kernels of dimensions 9,15,21,27, and 806 complete
nonconstant-A degree patterns.  The source, output, and proof are now
frozen for independent audit.  The [independent complete audit](planar_jc48_sep08_dg_genus_audit.md)
accepts the entire all-m argument, literal shifted source identity,
complete finite universes and both1550-gate replays. Root has promoted
the primary while preserving every frozen source/output byte.
