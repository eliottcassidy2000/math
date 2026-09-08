# Independent audit: the all-degree radical genus bound and extremal locus

**Status: PROVED / INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.**
This audit accepts the sharp genus bound in every degree and the complete
extremal coefficient classification in degrees divisible by four, with the
producer's stated hypotheses. Exactness is a property of `dx/sqrt(N)` in the
actual normalized radical field. No arbitrary two-variable Jacobian mate,
global regular mate, or planar Jacobian consequence follows from this audit.

Auditor: `three_ray_geometry`, independently of the root producer. I read the
entire [genus_capacity proof](planar_jc48_sep08_genus_capacity.md) and source,
reconstructed the analytic argument, and ran independent normal and optimized
replays. Both reproduce the frozen **5,081-gate** output byte for byte. A
separate additive partition construction recovers all 7,334 declared partitions,
1,103 necessary capacity survivors and all six extremal rows. No mathematical,
source or prose correction was required.

## 1. Types and the all-degree local argument

The actual field is `K=C(x)(y)`, where `y^2=N` for a nonzero polynomial N. If N
is square over C, K is `C(x)`; the two formal components of a reducible equation
are not incorrectly counted as one connected double cover. Its genus is zero.
For nonsquare N, K is the connected degree-two function field and all local
orders refer to its smooth projective normalization. Affine x is fixed except
for a translation in the extremal classification. The hypotheses do not
silently replace the differential by an unweighted projective transform.

The local-order and primitive-degree mechanism is inherited from the proved
[boundary_exactness classification](planar_jc48_sep08_boundary_exactness.md),
with named antecedents
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and
[THM-2723, split-exact-square-prefix-rational-primitive-pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
Here it is proved again with no degree cutoff. The finite partition enumeration
is a consistency control and is not used to justify the unbounded implication.

At a root of odd multiplicity j, the normalized parameter has `x-p=z^2`.
The differential has order `1-j`; when j>=3 any primitive can have pole order
at most j-2. A simple root is regular. At an even root j the cover is unramified
with two points and differential pole order j/2 at each. Multiplicity two gives
a nonzero simple-pole residue, hence cannot occur in an exact polynomial. For
even j>=4, each primitive pole has order at most j/2-1. These two points
therefore contribute j-2 in total, matching the odd-root contribution. A
primitive cannot acquire a pole at a place where the differential is regular,
because in characteristic zero differentiating a pole raises its order by one.

Thus for n>=3 the only possible primitive poles are finite and their total
degree is bounded by

    P = sum over root multiplicities j>=3 of (j-2).

A primitive is nonconstant because `dx/y` is nonzero. Its map to P1 has degree
equal to its total pole degree, so this upper bound can be compared with any
one local mapping degree of that same primitive.

For odd n, the unique infinity point has differential order n-3. Integrating a
nonzero leading coefficient gives primitive local mapping degree exactly n-2.
Thus `P>=n-2`. If there is more than one distinct root, a simple root costs one
in `n-P` and each high root costs two; with multiplicity two already forbidden,
`P<=n-3`. Therefore exact odd degree n>=3 forces one distinct root, a pure odd
power, and genus zero. This argument also rules out all-simple odd patterns
without treating them as a separate numerical case.

For even n>=4 there are two unramified infinity points. Each has differential
order n/2-2 and primitive local degree n/2-1. The correct conclusion is
`P>=n/2-1`; these local degrees are not added, since the two primitive values
may differ. The proof retains this distinction. Degrees zero and one are
rational and exact; the nonexact degree-two cases are also covered by the
inherited elementary classification (a finite logarithmic pole for the square
case, or infinity logarithmic poles for the two-simple-root case).

## 2. Branch arithmetic and the equality pattern

Write s for the number of simple roots, h_o for odd roots of multiplicity at
least three, and h_e for even roots of multiplicity at least four. For even
n and nonsquare N the branch count has no infinity contribution and equals

    B=s+h_o=2g+2,       P=n-s-2h_o-2h_e.

The infinity degree bound is equivalent to

    B+h_o+2h_e <= n/2+1.

There must be a high root: otherwise P=0 contradicts the positive required
local degree. Consequently `h_o+2h_e>=1`, so `B<=n/2` and

    g <= floor((n-4)/4).

The square case is genus zero and separately satisfies this bound for all
stated even degrees. No branch-count formula for a connected double cover is
applied to a square polynomial.

At n=4m, m>=2, equality in genus gives B=2m. The preceding inequalities force
`h_o=1`, `h_e=0`, and `s=2m-1`. Degree then fixes the sole high multiplicity
as 2m+1. Thus the unique possible extremal pattern is

    (2m+1, 1 repeated 2m-1 times).

In particular no additional even high root or fixed simple factor can be
hidden in an extremal polynomial. This is only the necessary multiplicity
pattern; its sufficiency still requires the position argument below.

## 3. Complete primitive space, not merely a residue test

Move the high root to u=0 and write `N=u^(2m+1)D(u)`, where D is squarefree,
`deg D=2m-1`, and D(0) is nonzero. The change

    X=1/u,       Y=y/u^(2m)

retains the weighted differential exactly:

    Y^2=P(X)=X^(2m-1)D(1/X),
    dx/y=-X^(2m-2)dX/Y.

P is squarefree of odd degree 2m-1 with nonzero constant and leading terms.
This is a smooth affine hyperelliptic curve with one point at infinity. Its
infinity point corresponds to the original high root. The differential's only
pole is there, of order 2m, so any primitive is regular on the entire affine
curve and has pole order at most 2m-1 at infinity.

I checked the use of the full affine ring here. Smoothness implies normality,
and the global regular functions on this affine curve are precisely
`C[X,Y]/(Y^2-P(X))`; a rational function regular at all affine places lies in
this ring. Every such element is uniquely `A(X)+YB(X)` with polynomial A,B.
At infinity, X has pole order two and Y has pole order 2m-1. Even and odd pole
orders cannot cancel, nor can lower monomials cancel the leading pole within
either summand. The complete pole space is therefore

    L((2m-1) infinity)=span{1,X,...,X^(m-1),Y}.

The differential changes sign under the hyperelliptic involution. Taking the
odd part of any primitive keeps its derivative and leaves only a constant
multiple C*Y in this pole space. Consequently

    (C/2)P'(X)=-X^(2m-2).

C cannot vanish. Integration forces `P(X)=a+bX^(2m-1)`, a,b nonzero, and the
inverse coordinate change gives exactly

    N=(x-p)^(2m+1) [a(x-p)^(2m-1)+b].

Conversely `-2Y/((2m-1)b)` differentiates to the required form, proving the
full iff rather than just necessity. The restriction m>=2 correctly excludes
the additional square genus-zero extremal branch at n=4. The nonsquare pole
space itself also works at m=1, as stated.

No period assumption is omitted: the complete odd pole space reduces every
possible rational primitive to C*Y, so the coefficient condition pays more than
vanishing residues. This is the decisive mechanism excluding wrong-position
polynomials with the same genus and pole capacity.

## 4. Sharpness and hostile controls

For every r>=1, a,b nonzero, the polynomial

    N=u^(r+2)(a*u^r+b),       n=2r+2

has r simple nonzero residual roots. The high root is branched precisely when
r is odd. The resulting branch count gives genus `floor((r-1)/2)`, exactly
the claimed bound for that even n. Direct differentiation with

    R=-2/(r*b*u^(r+1))

gives `N R'+N'R/2=1`, so `sqrt(N)*R` is a primitive in the actual field.
This proves sharpness also in degrees congruent to two modulo four, where the
producer makes no complete extremal classification claim.

The literal `N=x^7(x^5+1)` has six branch points and genus two; its primitive
`-2sqrt(N)/(5x^6)` is exact. Thus the degree-eight genus-one restriction cannot
be transported to all degrees without retaining n. The source checks this
actual identity, not just its branch count.

The opposite control retains the extremal multiplicity pattern but takes
`P(X)=1+X+X^(2m-1)`, m>=2. Its derivative has a nonzero constant term and cannot
satisfy the necessary sparse equation. Smoothness is proved for every m, beyond
the finite source checks: putting d=2m-1, a repeated root would obey both
`X^d=-X/d` and `1+X+X^d=0`, hence `X=-d/(d-1)`. Its derivative equation would
also demand `X^(d-1)=-1/d`, impossible because this X is negative real and
d-1 is even. Constant and leading terms are nonzero. This gives a genuine
same-genus normalized-curve hostile, not a singular degeneration that escapes
the theorem's hypotheses.

## 5. Independent finite reconstruction and source audit

I read the full 4,102-byte source. All checks use an always-active `need`
function that raises on failure and exact symbolic cancellation. Its partition
generator enumerates each partition once in descending order. It excludes
multiplicity two and square polynomials before imposing the connected-curve
capacity bound; none of its 1,103 survivors is asserted exact merely because
it passes that bound.

As an independent enumeration, I constructed sets of ascending partitions by
adding each possible part size to previously obtained totals. This uses a
different traversal from the producer's recursive first-part generator. It
recovers exactly 7,334 partitions in degrees 3..24 and 1,103 necessary capacity
survivors. After the filters, every odd survivor is pure. Every even survivor
obeys the genus bound. The extremal rows, reordered as decreasing partitions,
are precisely

    n=4:  (3,1)
    n=8:  (5,1,1,1)
    n=12: (7,1,1,1,1,1)
    n=16: (9,1,1,1,1,1,1,1)
    n=20: (11,1,1,1,1,1,1,1,1,1)
    n=24: (13,1,1,1,1,1,1,1,1,1,1,1).

The source additionally verifies the rational primitive identity for r=1..14
(even degrees 4..30), the complete monomial pole-space lists for m=1..8, every
intermediate derivative coefficient in those lists, and the literal high-genus
and wrong-position controls. These finite checks corroborate the analytic
formulas and scopes; they do not replace the unbounded argument or claim a
complete exactness classification outside the stated 4m extremal locus.

## 6. Independent replays and frozen pins

Both independent commands completed successfully from the repository root:

    python3 04-computation/planar_jc48_sep08_genus_capacity.py
    python3 -O 04-computation/planar_jc48_sep08_genus_capacity.py

Their separately saved outputs are byte-identical to the frozen 528-byte
output, including all 5,081 gates and the gate-label digest
`073cf1ef4384cb17dff9b6a72e290f1ecdff5078ab350223e13c7003c445c1b1`.
No producer file was edited by this auditor.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/planar_jc48_sep08_genus_capacity.py` | 4102 | `9def2dae7f430cb5b7611124eec7422dda7d2861c566f922401a8cb542f21a1f` |
| `05-knowledge/results/planar_jc48_sep08_genus_capacity.out` | 528 | `d870a22095b9f7144780b20338deb45c61befda31018b46a33a4a79728ade07f` |
| Primary proof read at acceptance, before status promotion | 9400 | `d97309b9e2421015c433ac5e71938f0471b51fcca92088713329ab7c4101f5f7` |

Root owns any subsequent status/link promotion. This audit's mathematical
acceptance and source/output pins are final; no correction remains pending.
