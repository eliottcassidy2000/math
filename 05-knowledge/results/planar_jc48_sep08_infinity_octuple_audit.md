# Independent audit of the infinity-octuple exclusion

**Status: full independent analytic/source audit PASS.**
Root read the complete [primary proof](planar_jc48_sep08_infinity_octuple.md)
and [source](../../04-computation/planar_jc48_sep08_infinity_octuple.py),
and recovered the exact M-unit and shared-root local dependencies.
The new theorem excludes polynomial mates of unrestricted degree in
the ORIGINAL source plane. Rational mates survive in its composite
case; they are included as an explicit hostile.

## 1. Surface, section, and volume transformations

Product inversion on P1 times P1 preserves the graph z=x² and is an
actual automorphism of its complement W. Its source-chart expression
u=1/x,T=-x²-x⁴t is rational, and its inverse gives the exact Jacobian
dx wedge dt=u² du wedge dT. This weight is essential. The proof never
imports the unweighted finite-octuple mate exclusion through this map.

The old graph equation becomes -(zeta-u²)/(u²*zeta). Squaring it gives
the stated positive numerator transformation for H; its first power
gives the negative sign for the L numerator. These transformations
retain the full section boxes. On the new boundary, an original
constant octic becomes alpha*u⁸. The complete affine fibre of that
octic restriction is alpha*h² plus the full six-dimensional L1 space.

The source pays this dimension statement by a rank-nine restriction
matrix on all15 section monomials, and by an invertible six-column
L1 coefficient matrix. Multiplication by the nonzero graph equation
injects L1 into the restriction kernel; its dimension is six, so this
is the entire kernel. No first normal coefficient is omitted.

The auditor independently transformed the OLD six-element global basis.
In the order 1,t,xt,x²t,x+x³t,x²+x⁴t, its images are respectively
1,-h,-v,-1-w,-uT,-T. This checks the actual coefficient-space bijection,
including the constant shift, by a different route. It also verifies
h=-old t and v=-old x*old t, which are load-bearing in the final
polynomial step.

On the new second chart the weighted volume is dr wedge db, with no
pole. On u=0 it vanishes to order two. These chart facts make the
compact-component holomorphic test legitimate: the only possible poles
of the generic relative form lie at the one boundary point. A generic
component cannot be contained in the fixed zero divisor of the volume.
Thus the relative form is nonzero on every generic component, even if
the generic fibre is reducible. No generic irreducibility is assumed in
the early local steps.

## 2. Unit cases and exact local degree

For M a unit and octic multiplicity8, the inherited full local
classification has no balanced integer case8=3j. A normal unit gives
exactly two local determinations. For positive normal order, the local
degree is three: either one low branch and two cancellation determinations,
or the dominant simple cubic. The earlier local proof makes every
unweighted form regular in these unbalanced cases. Multiplication by
u² improves the orders. This is a use of the proved LOCAL classification,
not its unweighted global conclusion.

When M is not a unit but N_s is a unit, the inherited shared-root
normal-unit lemma again makes every unweighted local form regular.
Only after removing both unit cases does the exact equation become
E(0,s)=(k0²+l0-c)s⁴. Its generic nonzero coefficient pays the local
degree four used in all subsequent ledgers. The source checks degrees
two, three and four separately; it never counts a branch away from
the boundary as an additional local branch.

## 3. First normal order one

If kxt is nonzero, the low s=uZ equation has a quadratic factor with
constant kxt² and discriminant of nonzero c-slope4kxt². Its two
nonzero simple roots give unweighted logarithmic branches, hence
weighted regular forms of order one. They are actual branches, with
local coordinate u, not just formal roots counted at an endpoint.

The two other determinations cancel alpha*u⁸+kxt*u*s at s~u⁷.
The normal derivative has order one; M-cs has boundary order n=1,2,3,4
because L is nonconstant and its normal terms start later here.
The weighted form has exponent (9-n)/2 in u. Every value is positive,
and when the split ramifies, u=tau² adds the derivative order needed
for regularity. These two determinations and the low pair exhaust the
degree-four restriction. Thus this entire case is holomorphic on the
compact normalized fibre and cannot have a rational primitive.

## 4. The linear M jet and its only balanced contact

With kxt=0 and lxt nonzero, the possible remaining normal orders are
2,3,4,6; the last coefficient2alpha is forced when the lower ones
vanish. There is always one simple low s~u branch with weighted
regular form.

At normal order2, a middle s~u³ branch is a simple nonzero root of
Z²(k2²+lxt Z), and the last two determinations cancel at s~u⁶.
The middle unweighted logarithm becomes regular. The final weighted
exponent5/2 becomes order six after quadratic normalization. All four
determinations are accounted for.

At normal order3, the remaining cubic at s=u⁵Z is
(alpha+k3 Z)²+lxt Z³. Its constant is nonzero. Coefficient comparison
rules out a triple root: the supposed cube identities require9=12
after canceling nonzero factors. At a simple root the weighted form
has order one. At a double root, the actual scaled germ Q has a unit
second derivative in Z and Q_c=-u⁴Z⁴. Its analytic critical centre
therefore has critical value with exact derivative -u⁴z(u,c)⁴.
Since its initial centre is nonzero, generic critical order lambda
is at most four. This is a generic-parameter bound, not a finite
truncation assumption.

Parametric Morse normalization gives a unit times u du/sqrt(-R).
For lambda1,2,3 the normalized orders are respectively2,0,0, all
regular. For lambda4 there are two smooth branches of order-1 with
nonzero leading coefficients. Those coefficients are their nonzero
residues, so exactness fails. Higher jets cannot cancel a simple-pole
residue supplied by a nonzero leading unit.

At normal orders4 and6 the remaining cubic is alpha²+lxt Z³, with
three simple nonzero roots and weighted order one. Together these
cases exhaust all possible lower normal jets. They either make the
entire relative form holomorphic or exhibit an unavoidable nonzero
residue; neither outcome permits a rational primitive.

## 5. The remaining order-two coefficient

After both first jets vanish, k2 nonzero gives the complete two low
double-pole branches of the fixed-zero proof. The weight u² makes
them regular of order zero. The high cancellation centre at s~u⁶
has n=2,3,4 and weighted exponent3-n/2. On the actual normalizations
these are orders2,4,1. The source checks the leading suppliers and
these normalized values, while the degree-four restriction proves
exhaustion. Consequently k2 must vanish in a rational-mate candidate.

The original constant-L case was already excluded as a polynomial
case by J(H²+L,G)=2H*J(H,G). Thus no missing infinite-order boundary
case is concealed by the n=2,3,4 list.

## 6. The entire quadratic field and genus-two form

The remaining h,v functions generate the original rational field.
The correct volume becomes (1/h) dh wedge dv. The full F, including
every remaining coefficient, is a(h)v²+2b(h)v+d(h). With Y=h(av+b),
the original generic fibre and relative form are exactly

    Y²=h²[b²+a(c-d)],       eta=-dh/(2Y).

The field map is birational when a is nonzero; zeroes or poles of a
at individual points do not invalidate its function-field inverse.
If l2 is nonzero, the polynomial on the right has degree five and
leading coefficient -l2*alpha². Its coefficient of c is
B=h(k3²h+l2), with only simple roots. Outside B=0, repeated roots
require critical points of the nonconstant rational function Q0/B.
There are only finitely many. At a persistent zero shared by Q0 and B,
the derivative has nonzero c coefficient B', so the root is simple
generically. The source's special persistent-second-root control is
valuable: one must not incorrectly assume Q0 and B are coprime.

This pays generic squarefreeness and a geometrically integral genus-two
normalization. At finite branch points dh and Y both have order one.
At the unique infinity point they have orders-3 and-5, so dh/Y has
order two. The form is everywhere holomorphic and nonzero, excluding
a rational primitive. The argument covers all k3 and l3 values when
l2 is nonzero.

The auditor also reconstructed the quadratic directly from its first
and second derivatives in v, then independently checked its c coefficient,
quintic leading coefficient and original-source t=0 derivatives. These,
together with the old-basis transformations, give ten independent
symbolic identities. The explicit source critical point for l2 nonzero
is a valid weaker polynomial check; the genus-two argument proves the
stronger rational exclusion.

## 7. Logarithmic and composite boundaries

If l2=0 and k3 is nonzero, the exact form is -dh/(2h*y). The quadratic
right side at h=0 has nonzero c coefficient k3², so it is a unit for
generic c. Each of the two unramified points there has a nonzero
logarithmic residue. This remains true after extending the generic
constant field, including constant quadratic extensions; geometric
integrality is not required in this residual case.

If l2=k3=0 but l3 is nonzero, the exact linear-field form is
-dh/(l3*h), again with a nonzero residue. If l3 also vanishes, F=P(h)
has quartic P. The proof returns to h=-old t in the original polynomial
ring. Its Jacobian has the nonconstant polynomial factor P'(-t), so
cannot be a nonzero constant. No polynomiality is transported through
the rational inverse map in this last step.

The actual hostile F=t⁴,G=-x/(4t³) has Jacobian one and the required
infinity octic. It prevents strengthening the whole theorem to a
rational-mate exclusion. The combined finite/infinite conclusion is
therefore precisely: every single-point boundary octic is excluded
for a polynomial mate, with no mate-degree bound.

## 8. Verification and acceptance

The entire source uses exact symbolic algebra and always-active checks;
the complete analytic branch arguments supply its unbounded claims.
Normal and optimized replays each pass141 gates and match the399-byte
frozen output. Source:12,949 bytes, SHA256
`8fc3fcf4c719deb653970090d8ca63589a2a3649074639d59cf96ef764bde333`.
Output SHA256:
`96805cf8e6722abd373477ed3ab6f1af62b83c50b7c18d1a0c7aba476978a33e`.

No mathematical or source correction was needed. Root accepts primary
promotion with its original polynomial-ring scope and every retained
unit, branch and constant-extension boundary. Other octic multiplicity
partitions and the general Jacobian conjecture remain open.
