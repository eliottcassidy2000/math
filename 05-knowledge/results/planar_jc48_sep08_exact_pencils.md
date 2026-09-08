# All exact radical pencils through degree eight and quadratic rational mates

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The first theorem classifies linear spaces of differentials. The second
gives a full rational-mate criterion for the stated quadratic class,
including every global quadratic on the fixed DG surface. Rational
mates may have poles; this is not a global regular pair or a solution
of the planar Jacobian conjecture.

## 1. The complete pencil theorem

Let V be a two-dimensional complex linear subspace of C[x] consisting
of polynomials of degree at most eight. Say a nonzero N is exact if
`dx/sqrt(N)` has a primitive in `C(x)(sqrt(N))`. Suppose a generic
member of V is exact, meaning every member in a nonempty Zariski-open
subset of the projective line P(V) is exact.

**Pencil theorem.** For some fixed p in C and `u=x-p`, V
is exactly one of the following six spaces:

    span{1,x};
    u^3 span{1,u};
    u^5 span{1,u};
    u^7 span{1,u};
    u^4 span{1,u^2};
    u^5 span{1,u^3}.                                  (1)

Conversely every nonzero member of each space in (1) is exact, including
the degree-drop and root-collision members. An overall nonzero scalar
does not change the space. Only affine translation of the one-variable
differential is used; no projective change of its weight is asserted.

The closest proved mechanism is the complete
[degree-eight radical differential classification](planar_jc48_sep08_boundary_exactness.md),
which pays all67 multiplicity partitions and the four position conditions.
Its primitive-degree antecedents are
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and [THM-2723, split-exact-square-prefix-rational-primitive-pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
The present operation asks which projective lines can remain inside that
proved locus. Targeted exact/synonym searches recovered no equivalent
pencil classification in current canon. This is not an external-priority
claim.

The live concepts are a fixed gcd; genuinely moving simple roots;
individual exactness versus a whole pencil; rational parameter descent;
and the complete quadratic discriminant. The source is a linear space
of polynomial differentials, the target is a sparse two-dimensional
space in (1), and the map retains the gcd and the position equations.
It loses the original basis of the pencil, which is restored for the
mate construction. The cheapest hostile is a line joining two exact
points on the6+1+1 conic: its intermediate members are not exact.

## 2. Generic repeated roots are fixed, and this exhausts the table

Choose independent A,B spanning V. Write `A=R A0`, `B=R B0`, with
R their polynomial gcd and `gcd(A0,B0)=1`. The rational function
`A0/B0` is nonconstant, since A,B are independent. Its derivative has
nonzero numerator

    W=A0' B0-A0 B0'.                                  (2)

A repeated finite root of `A0+c B0` must be a root of W. At such a
root B0 cannot vanish, by coprimality, and the root determines at most
one value c. Hence the residual pencil is squarefree for generic c.
It also avoids every root of R generically: at a fixed root of R the
two residual values are not both zero, so at most one c is excluded.
The residual generic degree is positive; degree zero would make A,B
dependent. Its leading coefficient drops for at most one c.

Consequently every repeated root of the generic full pencil is a fixed
root of R, with fixed multiplicity. Every residual root is simple and
moves. Some simple roots of R can initially be retained too; the
argument below eliminates such extra fixed roots when necessary.

The complete proved table has17 possible exact types. The only ones
with a simple root available to move are

    1; 3+1; 5+1; 7+1; 4+1+1; 6+1+1;
    4+3+1; 5+1+1+1.                                  (3)

There are12 possible allocations of at least one of their simple roots
to the moving residual polynomial. No type without a simple root can
occur in an independent pencil.

For type1, V has degree at most one and dimension two, hence is
`span{1,x}`. For `m+1`, m=3,5,7, the repeated root p is fixed and
`V subset u^m span{1,u}`; dimension two gives equality.

For `4+1+1`, the fourfold root is fixed. Write a generic member as
`u^4(a u^2+b u+d)`. The proved midpoint condition is b=0. Since the
coefficients depend linearly on the pencil parameter and vanish
generically, b vanishes identically on V. Thus
`V=u^4 span{1,u^2}`. Its residual space has gcd1, so the alternative
allocation with an additional fixed simple root is impossible.

For `5+1+1+1`, the same argument uses the exact elliptic condition:
after factoring the fixed u^5, the cubic's u^2 and u coefficients
both vanish. Hence `V=u^5 span{1,u^3}`, again leaving no additional
fixed simple root.

For `6+1+1`, the sixfold root is fixed and the residual quadratic
`a u^2+b u+d` must satisfy

    3b^2=4ad.                                        (4)

This nonsingular projective conic contains no projective line. Here is
an elementary check retaining all coefficient charts. If two endpoint
vectors have nonzero a-coordinates a0,a1, substitute
`di=3bi^2/(4ai)` into the mixed coefficient of (4). It becomes

    -3(a0 b1-a1 b0)^2/(a0 a1)=0,

so the vectors are proportional. If an endpoint has a=0, then b=0
and d!=0. Its mixed equation forces a=0 for the other endpoint,
whose own equation then forces b=0. They are again proportional.
If necessary choose a different pair of endpoints so the first chart
applies. Thus no two-dimensional residual space can satisfy (4).

Finally `4+3+1` has two fixed repeated roots p,q. Its sole simple
root r would have to satisfy `3/(p-q)+1/(p-r)=0`; this fixes r as
well. All members are then scalar multiples of one polynomial,
contradicting dimension two. This exhausts (3), including every
possible fixed-simple-root allocation, and proves necessity in (1).

## 3. Universal primitives over the parameter field

The following formulas prove sufficiency and also retain the coefficient
field needed for the quadratic consumer. They use rational functions
of a,b and the radical y itself; no square root of a or b is adjoined.
Take a,b generic in any characteristic-zero field containing C.

For `y^2=a x+b`, use `B=2y/a`.

For `y^2=u^(2k+1)(a u+b)`, k=1,2,3, put `v=y/u^k`, so
`v^2=a u^2+b u`. Define

    C_k=-2/((2k-1)b),
    C_j=-2j*a*C_(j+1)/((2j-1)b),   j=k-1,...,1.

Then a primitive is

    B=v sum_(j=1)^k C_j u^(-j).                       (5)

Indeed for `D=a u^2+b u`,

    d(v R)=(D R'+D'R/2) du/v,
    D (u^-j)'+D'u^-j/2
       =(1-j)a u^(1-j)+(1/2-j)b u^-j.

The recurrence cancels every term except u^-k, proving `dB=du/y`.

For the remaining two lines use

    y^2=u^4(a u^2+b):   B=-y/(b u^3);
    y^2=u^5(a u^3+b):   B=-2y/(3b u^4).                (6)

Direct differentiation proves both identities. In the elliptic case
the inversion `X=1/u`, `Y=y/u^4` gives `Y^2=a+bX^3` and
`du/y=-X^2 dX/Y`, explaining the second formula.

The displayed denominators are valid for generic pencil parameters.
For the exceptional a=0 or b=0 members, the polynomials become the
pure powers of exponents1,3,4,5,6,7,8, or a nonzero constant. The
proved pure-power primitive covers each one; exponent2 never occurs.
Thus every nonzero member of every line (1), not merely a generic one,
is exact.

## 4. Complete rational-mate criterion for a quadratic discriminant of degree eight

Let

    H=N(x)t^2+P(x)t+Q(x) in C[x,t],   N!=0,
    D0=P^2-4NQ,    deg N<=8,    deg D0<=8.             (7)

The last two degree bounds are hypotheses for this general statement.
They will be automatic for the actual global DG layer below.

**Quadratic criterion.** A rational G with `J(H,G)=1`
exists if and only if one of the following mutually exclusive cases
holds:

* D0 and N are independent, and their span is one of the six spaces (1).
* `D0=kappa N` for a complex constant kappa, and `dx/sqrt(N)` is
  exact in `C(x)(sqrt(N))`.

The second case includes D0=0. Its condition has the complete17-type
classification and position equations of the inherited theorem.

On a generic fibre H=c, the rational coordinate

    z=2Nt+P,        z^2=D0+4cN=:D_c                   (8)

identifies its function field with the quadratic discriminant field.
No points over zeros of N are thrown away: this is a function-field
identity, and the differential and primitive are subsequently understood
on the complete normalized curve. A rational mate would restrict to

    dG=-dx/z.                                        (9)

There are only finitely many bad fibre values at which specialization
of its denominator can vanish on an entire component. Thus every
generic discriminant D_c must be exact. If D0,N are independent,
the pencil theorem applies and gives the first case. If proportional,
choose any ordinary complex c with kappa+4c!=0 and scale (9) by its
nonzero square root; this gives exactness of `dx/sqrt(N)`. This also
covers geometrically split generic fibres when N is square, by using
each component with its chosen square-root sign.

For sufficiency in the independent case, use (5)--(6) or the linear
formula with a(c),b(c) the actual affine coefficient functions of D_c.
The result `B(c,x,z)` is rational in c,x,z and satisfies
`d_x B=dx/z` on `z^2=D_c`, with c held fixed. Substitute

    G=-B(H,x,2Nt+P).                                 (10)

Its denominators are nonzero rational functions: H is nonconstant,
the required generic coefficient functions are not identically zero,
and N,u,z are nonzero in C(x,t). At fixed H the derivative is exactly
`-1/(2Nt+P)`, so the chain rule gives `J(H,G)=1`. This proves
rational descent in the parameter; pointwise primitives alone would
not have paid it.

For sufficiency in the proportional case, choose a primitive in the
form `sqrt(N) R(x)` with R rational. If N is nonsquare, average any
primitive against the quadratic involution to make it odd, so it has
this form. If N is square, divide its rational primitive by the chosen
rational sqrt(N). In either case

    N R'+N'R/2=1.                                    (11)

Then the entirely rational formula

    G=-N R/(2Nt+P)                                   (12)

has Jacobian one. For example, at fixed H, differentiate (12) using
`z^2=(kappa+4H)N`; equation (11) gives `dG=-dx/z`.
No Galois or geometric-connectedness assumption is hidden in (12).

## 5. The entire global DG quadratic layer and its sharp elliptic boundary

For the actual surface

    W=(P1_x x P1_z) minus {z=x^2},   t=1/(z-x^2),

the [proved complete filtration](planar_jc48_sep08_dg_quadratic.md)
says every `H in L2` has the unique form

    H=(a(x)z^2+b(x)z+c(x))/(z-x^2)^2,
    deg a,deg b,deg c<=4.

Substituting `z=x^2+1/t` gives

    N=a x^4+b x^2+c,   P=2a x^2+b,   Q=a,
    P^2-4NQ=b^2-4ac.                                 (13)

Thus both N and D0 have degree at most eight. Formula (13) pays all
fifteen global coefficients at once; equivalently the full coefficient
formula cancels discriminant degrees12,11,10,9. Every genuine t-degree
two global H is therefore covered by the complete criterion in Section4.
Lower t-degree functions are outside that statement and already have
their separate filtration results.

Every generic geometric component of such a quadratic H with a rational
mate has genus at most one. In the independent case this follows from
the six lines: five give genus zero, and the final line gives the
elliptic curve `Y^2=a+bX^3`. In the proportional case it follows from
the complete inherited exactness table; its only elliptic type is
5+1+1+1 with that same sparse cubic equation. Square discriminants are
handled as rational components, not as connected double covers.

The genus-one possibility is realized inside the actual global ring.
Put `h=x^2+x^4t`, take delta!=0 and arbitrary beta,q, and define

    H=(x^4+delta*x)(1+x^2t)^2+beta*h+q,
    G=1/(3delta*h).                                  (14)

Then `J(H,G)=1`. In the actual other chart,

    H=(1+delta*r^3)b^2-beta*b+q,

so H is global. Its discriminant pencil is

    D_c=(4c+beta^2-4q)x^8+4delta(c-q)x^5.              (15)

For beta!=0 it spans the elliptic line; at beta=0 it is proportional
to `N=x^5(x^3+delta)`. Thus both ranks of the criterion meet in the
same exact rational family. The [five-plus-three-simple-root note](planar_jc48_sep08_five_one_one_one.md)
independently records this rational boundary of its polynomial exclusion.
The globality of G is not asserted; its poles are essential.

If `F=H^2+l0` with l0 constant, F has a rational mate if and only if
H does: multiply a mate of F by 2H in one direction and divide a mate
of H by 2H in the other. Hence Section4 also completely classifies
rational mates in this constant-L global quartic subfamily. Polynomial
mates there are already excluded by the nonconstant Jacobian factor2H.
This does not classify quartics with nonconstant L.

## 6. Hostiles, validity boundary and reproduction

The exact octics

    N_plus=x^6(x^2+2x+3),
    N_minus=x^6(x^2-2x+3)

each satisfy the6+1+1 residue equation. Their sum does not. Thus exact
endpoints alone do not define an exact pencil. The first failed
implication is replacement of a nonlinear exactness locus by its span.
The repaired statement requires generic exactness of the whole line;
the conic argument explains the obstruction.

The proportional boundary is also essential: `H=x^5(x-1)^3 t^2`
has the rational mate

    (8x^2-4x-1)/(3x^4(x-1)^2 t),

although its fixed5+3 pattern supports no independent exact pencil.
Finally `H=t^2+x^2` and `H=t^2+x^3` have generic discriminants with
two or three simple roots and are excluded; merely having a constant
t-leading coefficient is not sufficient.

The next connection is the quartic inverse's later coefficients.
This complete quadratic discriminant criterion pays a fixed-fibre
primitive classification. It does not assert that an infinite formal
primitive for a deformed quartic converges, is rational, or is global.
Those are separate obligations.

Reproduce from the repository root:

    python3 04-computation/planar_jc48_sep08_exact_pencils.py
    python3 -O 04-computation/planar_jc48_sep08_exact_pencils.py

The explicit finite universe consists of all17 inherited exact types and
all12 moving-root allocations, the six universal symbolic primitives,
one genuine quadratic mate for each line, the full fifteen-parameter DG
discriminant, proportional and elliptic-rank controls, and the named
nonlinear-span and generic-discriminant hostiles. No producer module is
imported. Always-active exact checks support the identities; the complete
line theorem and rational descent are analytic, not finite sampling.

All64 gates pass in byte-identical normal and optimized replays. Source
and output are frozen. The [independent complete audit](planar_jc48_sep08_exact_pencils_audit.md)
accepts the complete six-space theorem, all-member extension, rational
parameter descent, both ranks of the quadratic iff, full DG identity,
genus bound and constant-L consumer. The fixed-versus-residual-root
wording was clarified before freeze; no mathematical correction remains.

04-computation/planar_jc48_sep08_exact_pencils.py: 7430 bytes; SHA256 040e6af0681f8b2d04704418d6e3fc179df94ae3ffc6454cbc14f3647da59b80

05-knowledge/results/planar_jc48_sep08_exact_pencils.out: 596 bytes; SHA256 be1ad4ecb256073ab2e92d2e1d8b48eadaa6ef669c06678d689dcc2a8b1469e1
