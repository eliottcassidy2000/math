# Independent audit: exact radical pencils and quadratic rational mates

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
This audit accepts the complete six-space theorem, its all-member strengthening,
parameter-rational primitive formulas, and the necessary-and-sufficient rational
mate criterion in the producer's stated quadratic domain. It also accepts the
complete DG quadratic consumer, its sharp genus-one boundary, and the
constant-L quartic equivalence. It makes no claim about globally regular mates,
quartics with nonconstant L, or the planar Jacobian conjecture.

Auditor: `three_ray_geometry`, independently of the root producer. The producer
is [exact_pencils](planar_jc48_sep08_exact_pencils.md); its source and frozen
output have the same stem under `04-computation` and this directory. The source
was read in full, then independently run in normal and optimized modes. Both
outputs equal the frozen output byte for byte: **64 always-active exact gates**.
No source or mathematical correction was required. One pre-promotion prose
clarification was requested and checked: the generic-root argument now says
“Every residual root is simple and moves,” preserving its next sentence about
possible fixed simple gcd roots.

## 1. Domain, inheritance and generic-root exhaustion

The pencil is a two-dimensional complex linear subspace of the polynomials of
fixed affine x-degree at most eight. Exactness refers to the differential
`dx/sqrt(N)` in the actual field `C(x)(sqrt(N))`, including the square case where
this is just `C(x)`. The zero polynomial is excluded. Generic exactness means a
nonempty Zariski-open subset of the projective parameter line. These hypotheses
are used; this is not a classification of a nonlinear family, two isolated exact
polynomials, or a differential after an unweighted projective coordinate change.

I checked the use of the proved
[boundary_exactness classification](planar_jc48_sep08_boundary_exactness.md),
which covers all 67 multiplicity partitions of affine degree zero through eight.
Its 17 exact types include four position-dependent cases. The producer does not
replace those conditions by their multiplicity patterns. Its pole-capacity
antecedents remain the named
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and
[THM-2723, split-exact-square-prefix-rational-primitive-pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
The new pencil argument is proved directly; no external-priority claim is
needed or accepted by this audit.

For independent A,B write `A=R A0`, `B=R B0` with coprime residuals. The rational
function A0/B0 is nonconstant and its Wronskian is a nonzero polynomial. At a
repeated root of A0+cB0, B0 cannot vanish; the root must be a zero of that fixed
Wronskian and determines at most one c. Thus only finitely many c have a
repeated residual root. Avoidance of the finitely many roots of R and stability
of the positive generic residual degree each exclude at most finitely many
additional values. This remains correct when either endpoint has smaller degree
or one residual is constant. It does not require individual moving roots to be
rational functions of c.

Consequently all generic repeated roots are fixed gcd roots, whereas every
residual root is simple and nonconstant. Fixed simple gcd factors are retained
at this stage. Of the inherited exact types, precisely

    1; 3+1; 5+1; 7+1; 4+1+1; 6+1+1; 4+3+1; 5+1+1+1

can have a positive number of moving simple roots. Counting every possible
positive allocation of their simple roots gives twelve allocations. No other
exact type can support an independent pencil: all its root locations would be
fixed and its polynomial would vary only by scalar.

I checked each elimination, including the allocations with extra fixed simple
roots. Type 1 gives the whole affine-linear space. A sole high root of order
3,5,7 and one simple root gives `u^m span{1,u}` by dimension, for a fixed
`u=x-p`. For type 4+1+1 the midpoint equation kills the residual quadratic's
linear coefficient identically on the space; dimension then gives
`u^4 span{1,u^2}`. That residual space has gcd one, so no extra fixed simple
root remains. For 5+1+1+1 the two elliptic coefficient conditions similarly
give `u^5 span{1,u^3}`, again with residual gcd one.

The 6+1+1 condition is the genuine quadratic equation `3b^2=4ad` on the
residual coefficients. Its rank-three projective conic contains no projective
line. The producer's direct proof pays both coefficient charts: with a0,a1
nonzero its mixed equation is

    -3(a0*b1-a1*b0)^2/(a0*a1)=0,

forcing proportionality; when a0=0 the endpoint equation forces b0=0 and
d0 nonzero, then the polar equation forces a1=0 and its endpoint equation
forces b1=0. These are all possibilities, including zero b-coordinates. For
4+3+1 the two repeated roots p,q are fixed, and
`3/(p-q)+1/(p-r)=0` fixes the simple root as `(4p-q)/3`; this is proportional
rather than an independent pencil. There is no missing degree-drop or infinity
chart: all conclusions concern the stated affine differential and its complete
normalized curve, and the generic degree was fixed before using the table.

## 2. Universal primitives and parameter descent

I differentiated all six universal formulas independently at the level of a
quadratic differential field. For `y^2=D`, differentiation of `yR` gives
`(D R'+D' R/2) dx/y`; thus this check uses no analytic choice of square root.
For the three spaces `u^(2k+1)(a u+b)`, set `v=y/u^k`, so
`v^2=a u^2+b u`. The coefficient identity

    D*(u^-j)' + D'*(u^-j)/2
      = (1-j)*a*u^(1-j) + (1/2-j)*b*u^-j

shows that the stated recurrence cancels all coefficients except u^-k. The
terminal coefficient `Ck=-2/((2k-1)b)` has the correct sign and normalization.
The linear primitive `2y/a` and the two sparse formulas

    -y/(b*u^3),       -2y/(3*b*u^4)

also differentiate to `du/y` in their respective fields. These expressions are
rational in a,b,u,y. They adjoin neither sqrt(a) nor sqrt(b), which is essential
for the subsequent sufficiency direction over the parameter field.

For exceptional a=0 or b=0 members, the residual degree drops or roots collide.
Every nonzero resulting member is a constant or a pure power of exponent
1,3,4,5,6,7,8. All are exact and exponent two never arises. This pays the
all-member strengthening without claiming that a single displayed formula is
regular at every exceptional parameter. The pole-bearing formulas remain
legitimate rational functions at generic parameter.

## 3. The complete quadratic iff

The general domain is `H=N(x)t^2+P(x)t+Q(x)` with N nonzero and both N and
`D0=P^2-4NQ` of degree at most eight. Genuine t-degree two is required. In the
generic field, `c=H` and `z=2Nt+P` give

    z^2=D0+4cN,       t=(z-P)/(2N).

These are inverse function-field identities. Zeros of N do not remove places
from the normalization or erase residue obligations. With c held fixed,
`dt/dx=-H_x/H_t`, hence a Jacobian-one mate satisfies

    dG=-dx/z.

The sign agrees with `J(H,G)=H_x G_t-H_t G_x`.

Necessity survives specialization. A nonzero rational denominator has only
finitely many irreducible polynomial factors. A factor contained in a fibre
component of H=c determines that one constant c, so only finitely many c can
make the denominator vanish identically on a component. At the remaining
ordinary complex fibres the displayed differential is exact. When D0,N are
independent this gives generic exactness of precisely their projective pencil.
When `D0=kappa N`, choosing `kappa+4c` nonzero and scaling by its complex square
root gives exactness of `dx/sqrt(N)`. If N is square, the generic geometric
fibre may split; the argument on each component uses its chosen sign. There is
no hidden connectedness assumption or unjustified descent from a geometric
primitive to the original parameter field.

For independent D0,N, the pencil classification expresses them in one fixed
translated sparse basis. Their coefficient functions a(c),b(c) are affine in c.
Independence prevents any coefficient required in a generic denominator from
vanishing identically. The universal primitives therefore yield a rational
`B(c,x,z)` with `d_x B=dx/z` at fixed c. Substituting

    G=-B(H,x,2Nt+P)

is legitimate in `C(x,t)`: H is nonconstant, the generic denominator factors
remain nonzero field elements, and poles on special fibres are allowed. The
chain rule with H fixed gives the claimed Jacobian. This explicit construction,
not pointwise existence alone, pays rational parameter descent.

In the proportional case, every primitive can be chosen as `sqrt(N) R(x)`.
For nonsquare N, take its odd part under the quadratic involution; for square
N divide a rational primitive by a chosen rational square root. In either case

    N R' + N'R/2=1.

The formula `G=-NR/(2Nt+P)` is entirely rational and has the required
Jacobian. At fixed H, use `z^2=(kappa+4H)N` and the preceding identity to
obtain `dG=-dx/z`. This includes kappa=0, geometrically split fibres, and a
nontrivial constant field over `C(c)` without adjoining sqrt(c). These checks
establish both directions of the stated iff and their disjoint case split.

## 4. Actual DG scope, genus and constant-L consumer

I checked the claimed application against the proved
[dg_quadratic filtration](planar_jc48_sep08_dg_quadratic.md). Every actual
L2 section is

    [a(x)z^2+b(x)z+c(x)]/(z-x^2)^2,  deg(a),deg(b),deg(c)<=4.

Putting `t=1/(z-x^2)` gives

    N=a*x^4+b*x^2+c,   P=2a*x^2+b,   Q=a,
    D0=b^2-4ac.

Thus the required degree-eight bounds cover all fifteen global coefficients.
The cancellation is an identity, not a truncation or an imposed sparse chart.
Sections whose t-degree is below two are explicitly outside this theorem.
Rational mates need not be globally regular; earlier global regular-pair
obstructions are not contradicted.

The genus statement follows from the complete inherited table and the six
spaces. Their only positive-genus possibility is the sparse elliptic curve
`Y^2=a+bX^3`; squares give rational components. In the explicit actual global
family

    h=x^2+x^4*t,
    H=(x^4+delta*x)*(1+x^2*t)^2+beta*h+q,
    G=1/(3*delta*h),     delta!=0,

I checked the Jacobian, the other-chart expression
`H=(1+delta*r^3)b^2-beta*b+q`, and the complete discriminant

    (4c+beta^2-4q)*x^8 + 4delta*(c-q)*x^5.

The rank determinant is nonzero exactly when beta is nonzero (delta is already
nonzero). Thus both the independent elliptic pencil and its proportional
boundary occur in the actual global ring. The mate's pole is retained; the
result does not provide a global regular Keller pair.

For `F=H^2+l0`, l0 constant, the equivalence of rational mates is exact:
if `J(F,G)=1`, then `J(H,2HG)=1`; if `J(H,G)=1`, then
`J(F,G/(2H))=1`. The nonconstant factor 2H also rules out a polynomial mate
of F. Nothing here treats a nonconstant L or rationality of an infinite formal
quartic expansion.

## 5. Source universe and decisive controls

The source uses always-active exceptions, exact SymPy rational cancellation,
and no imports from another producer. Its finite checks have the scope claimed:

* All 17 inherited exact types and all 12 positive moving-root allocations are
  retained, including allocations later ruled out analytically.
* The conic's full coefficient charts, nondegeneracy, and a Wronskian identity
  are checked. The generic-root theorem is analytic, not inferred from a sample.
* All six symbolic primitive identities and one genuine quadratic mate for each
  pencil are checked. These latter controls satisfy the general discriminant
  degree bounds; they are not all mislabelled global DG sections.
* The unrestricted fifteen-parameter DG formula cancels discriminant degrees
  12,11,10,9 and checks its degree and parameter derivative exactly.
* The proportional 5+3 mate, the actual global elliptic family, its rank drop,
  and its other-chart expression are checked.
* The two exact 6+1+1 octics whose sum is not exact prevent replacement of the
  exactness locus by a linear span. The generic discriminants of `t^2+x^2`
  and `t^2+x^3` prevent treating the leading exactness gate as sufficient.

The complete pencil exhaustion, generic specialization, and field descent are
proved analytically; 64 gates are reproducible supporting identities, not a
finite search standing in for these implications. No unresolved mathematical
obligation or failed hostile control remains within the stated theorem.

## 6. Independent replay and frozen acceptance

Run from the repository root:

    python3 04-computation/planar_jc48_sep08_exact_pencils.py
    python3 -O 04-computation/planar_jc48_sep08_exact_pencils.py

I ran both commands independently, saved separate outputs outside the repository,
and compared both bytes against the frozen output. Both pass 64 gates and
produce exactly 596 bytes. Source and output were not edited by this auditor.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/planar_jc48_sep08_exact_pencils.py` | 7430 | `040e6af0681f8b2d04704418d6e3fc179df94ae3ffc6454cbc14f3647da59b80` |
| `05-knowledge/results/planar_jc48_sep08_exact_pencils.out` | 596 | `be1ad4ecb256073ab2e92d2e1d8b48eadaa6ef669c06678d689dcc2a8b1469e1` |
| Primary proof at acceptance, before status promotion | 15245 | `2f22afe5f0ec3c66b956def6dac3b6bf2e4e893b24439a7f671b59d34eaf1193` |

The primary pin records the clarified text read here while its header is still
RESERVED. A subsequent status/link promotion is owned by root and does not
change this mathematical acceptance or the frozen source/output pins.
