# A sharp genus bound for exact radical differentials in every degree

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
These are theorems about `dx/sqrt(N)` on its actual normalized radical
curve. They do not assert Jacobian-mate existence for an arbitrary
two-variable polynomial. Exact differentials of arbitrarily high genus
are explicitly realized, so the degree bound is essential.

## 1. Statements and inheritance

Let N be a nonzero complex polynomial of degree n. Choose y with
`y^2=N`, let `K=C(x)(y)`, and let g be the genus of its smooth projective
curve. Say N is exact if `dx/y=dB` for some B in K. If N is square,
K is C(x) and g=0; one does not count both components of `y^2=N` as
one connected curve.

**Genus theorem.** If N is exact, then:

* For odd n>=3, N has only one distinct root and g=0.
* For even n>=4,

      g <= floor((n-4)/4).                            (1)

The second bound is sharp for every even n>=4. There is no exact
degree-two polynomial; degrees zero and one have genus zero and are
exact.

**Extremal theorem.** Let n=4m with m>=2. An exact N has
genus g=m-1 if and only if, for some p in C and nonzero a,b,

    N=(x-p)^(2m+1) [a(x-p)^(2m-1)+b].                 (2)

This is a complete coefficient and position classification, not only
a multiplicity pattern. It recovers the exceptional elliptic degree-eight
class when m=2 and supplies extremal examples of every higher genus.

The closest mechanisms are the all-degree local-order and primitive-map
argument already used in
[the complete degree-eight classification](planar_jc48_sep08_boundary_exactness.md),
and its antecedents
[THM-2071, quadratic-fiber-square-parity-gate](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and [THM-2723, split-exact-square-prefix-rational-primitive-pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
Here their necessary capacity inequality is retained in terms of branch
count, instead of enumerating a bounded list of partitions. The equality
case then supplies a complete primitive space. No external-priority
claim is made.

The live concepts are labelled local poles, primitive degree, branch
count/genus, the equality pattern, and the involution's odd pole space.
The source is an exact radical differential, the target is (1) or (2),
and the map retains the pole orders and branch labels. Passing the
inequality alone loses position and period information; the extremal
odd primitive space restores exactly those missing equations.

## 2. The genus inequality and its equality pattern

At a finite multiplicity-two root, `dx/y` has nonzero simple-pole
residue on each unramified normalization branch. Such roots are
forbidden by exactness, including when N is square.

First take N nonsquare. At an odd root of multiplicity j, its normalized
parameter z has `x-p=z^2`, so `ord(dx/y)=1-j`. If j>=3, a primitive
can have pole order at most j-2 there. At an even root j>=4 there are
two points, each with differential pole order j/2 and primitive pole
order at most j/2-1. A simple root is regular. Thus the total possible
primitive pole degree at the finite points is

    P=sum_(j>=3) (j-2).                               (3)

For n>=3, the differential has no pole at infinity, so a primitive
cannot have a pole there either. If n is odd, infinity is one branch
with differential order n-3, forcing primitive local mapping degree
n-2. Hence P>=n-2. Any multiplicity pattern with more than one root
has P<=n-3, and only a pure odd power survives. This argument has no
degree cutoff.

For even n>=4, there are two infinity points, each with differential
order n/2-2. Each forces local primitive degree n/2-1, whence

    P>=n/2-1.                                        (4)

The two local degrees are not added: the primitive need not have the
same value at the two points.

Let s be the number of simple roots, h_o the number of odd roots of
multiplicity at least three, and h_e the number of even roots of
multiplicity at least four. The connected double cover has branch count

    B=s+h_o=2g+2,
    P=n-s-2h_o-2h_e.

Together with (4) these identities give

    B+h_o+2h_e <= n/2+1.                              (5)

At least one high root is present, because otherwise P=0 contradicts
(4). Thus `h_o+2h_e>=1`, and `B<=n/2`. Since B=2g+2, this is
exactly (1). If N is square, g=0 already satisfies (1); no connected
double-cover argument is applied to that case.

Now suppose n=4m, m>=2, and g=m-1. Then B=2m, and (5) forces
`h_o+2h_e<=1`. It follows that `h_o=1`, `h_e=0`, and s=2m-1.
There is exactly one repeated root, it has odd multiplicity2m+1,
and every other root is simple. This proves the necessary extremal
partition

    (2m+1,1,...,1), with2m-1 simple roots.              (6)

The position information is still unpaid by (6); Section3 supplies it.

## 3. The complete extremal primitive space

Translate the unique high root to u=x-p=0, and write

    N=u^(2m+1) D(u),
    deg D=2m-1, D(0)!=0, D squarefree.

With `X=1/u` and `Y=y/u^(2m)`, the actual curve and differential are

    Y^2=P(X):=X^(2m-1)D(1/X),
    dx/y=-X^(2m-2) dX/Y.                              (7)

The polynomial P is squarefree, has degree2m-1, and has nonzero
constant and leading coefficients. The affine curve is nonsingular.
There is a unique infinity point on this odd-degree hyperelliptic
model, corresponding to the original high root. The only pole of
the differential in (7) is there, of order2m. A primitive must be
regular everywhere else and have pole order at most2m-1 there.

The affine coordinate ring is normal and equals
`C[X,Y]/(Y^2-P(X))`. Every function regular off this infinity point
therefore has the form `A(X)+Y B(X)` with polynomials A,B. The pole
orders of X and Y are2 and2m-1. Even and odd pole orders cannot
cancel; within either part its leading monomial has unique largest
pole order. Consequently the complete space is

    L((2m-1)infinity)=span{1,X,...,X^(m-1),Y}.           (8)

The differential (7) is odd under `Y->-Y`. Averaging a primitive
against this involution makes it odd without changing its derivative,
so (8) forces that primitive to be C*Y. Equation (7) now becomes

    (C/2) P'(X)=-X^(2m-2).                            (9)

It follows that `P(X)=a+bX^(2m-1)`, with a,b nonzero. Undoing the
coordinate change is exactly (2). Conversely this P makes
`-2Y/((2m-1)b)` a primitive, proving the full iff.
The same argument works for m=1 in the nonsquare genus-zero case;
the theorem starts at m=2 so square genus-zero polynomials do not
introduce an additional extremal branch.

## 4. Sharp examples in every even degree and the high-genus hostile

For any integer r>=1 and a,b nonzero, put

    N=u^(r+2)(a u^r+b),    n=2r+2.

Its residual r roots are simple and nonzero. It has r branch points
there and an additional branch at u=0 precisely when r is odd. Thus
its genus is `floor((r-1)/2)=floor((n-4)/4)`. A primitive is

    B=-2sqrt(N)/(r b u^(r+1)).                        (10)

For example, direct differentiation of yR with
`R=-2/(r b u^(r+1))` gives
`N R'+N'R/2=1`. This proves sharpness in every even degree, including
degrees congruent to two modulo four, where no complete extremal
classification is asserted here.

The explicit degree-twelve polynomial

    N=x^7(x^5+1)

has genus two and exact differential with primitive
`-2sqrt(N)/(5x^6)`. Hence an unqualified extension of the degree-eight
elliptic bound to all degrees is false. The first missing coordinate
would be n, and the repaired bound is exactly (1).

Conversely the partition (6) alone is not enough. If the transformed
polynomial is `P(X)=1+X+X^(2m-1)`, m>=2, then its derivative has a
nonzero constant coefficient and cannot satisfy (9). This is a smooth
same-genus hostile with the correct pole budget but wrong positions.
Indeed a repeated root, with d=2m-1 odd, would satisfy both
`X=-d/(d-1)` and `X^(d-1)=-1/d`, impossible since the first value is
negative real and d-1 is even.
The complete odd pole space, rather than the intermediate genus count,
is what finishes the extremal theorem.

## 5. Reproduction and scope

Run from the repository root:

    python3 04-computation/planar_jc48_sep08_genus_capacity.py
    python3 -O 04-computation/planar_jc48_sep08_genus_capacity.py

The finite universe is all7,334 partitions in degrees3..24. Explicit
filters remove multiplicity-two roots and square polynomials before
testing the connected-curve capacity inequality; all1,103 capacity
survivors are labelled necessary-only. The source checks every one
against the uniform arithmetic, including all six extremal partitions
in degrees4m through24. It also checks actual sharp families in every
even degree4..30, complete monomial pole spaces for m=1..8, and the
named high-genus and wrong-position hostiles. The unbounded proofs are
analytic; no finite partition enumeration substitutes for them.

All5,081 always-active gates pass in byte-identical normal and optimized
replays. The [independent complete audit](planar_jc48_sep08_genus_capacity_audit.md)
accepts the all-degree inequality, complete equality locus, odd primitive
space, sharp examples and both hostile boundaries. A separate additive
partition DP reproduces all7334 partitions,1103 capacity survivors and
every extremal row. Source and output are frozen and independently accepted.

04-computation/planar_jc48_sep08_genus_capacity.py: 4102 bytes; SHA256 9def2dae7f430cb5b7611124eec7422dda7d2861c566f922401a8cb542f21a1f

05-knowledge/results/planar_jc48_sep08_genus_capacity.out: 528 bytes; SHA256 d870a22095b9f7144780b20338deb45c61befda31018b46a33a4a79728ade07f
