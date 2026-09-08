# Exactness of dx/sqrt(N) for every polynomial of degree at most eight

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**

This classifies a differential, not Jacobian mates. No theorem about a
full two-variable polynomial is used or proved here. The separately [proved and audited leading-coefficient expansion](planar_jc48_sep08_leading_exactness.md)
pays the necessary Jacobian-mate application; this classification remains
a statement about the differential alone. All root positions below are in the fixed affine x-line;
an arbitrary projective change of x changes the differential's weight.

## 1. Complete statement and scope

Let `N in C[x]` be nonzero, with `0 <= n=deg N <= 8`, and choose a square
root `y` in an algebraic closure of `C(x)`. Put

    K=C(x)(y),       y^2=N,       eta=dx/y.

Call `eta` exact if `eta=dB` for some `B in K`. Replacing the chosen square
root by its negative changes the primitive's sign and not exactness.
If `N` is a square, **K is C(x)**: one does not apply a connected quadratic
cover argument to the two components of `y^2=N`.

Write a multiplicity partition using only distinct finite roots. The
following table is exhaustive. A displayed partition with an extra
condition is exact precisely when that condition holds. Every partition
not displayed is nonexact, for every choice of its distinct roots.

| Degree | Exact partitions and any necessary additional condition |
| --- | --- |
| 0 | Nonzero constants |
| 1 | `1` |
| 2 | None |
| 3 | `3` |
| 4 | `4`, `3+1` |
| 5 | `5` |
| 6 | `6`, `5+1`, `3+3`; `4+1+1` subject to (1) |
| 7 | `7` |
| 8 | `8`, `7+1`, `5+3`; `6+1+1` subject to (2); `4+3+1` subject to (3); `5+1+1+1` subject to (4) |

Here the precise labelled conditions are as follows.

* For `4+1+1`, let `p` be the fourfold root and `q,r` the two simple
  roots. Exactness is equivalent to

      2p=q+r.                                           (1)

* For `6+1+1`, let `p` be the sixfold root and `q,r` the simple roots.
  Put `d1=p-q`, `d2=p-r`. Exactness is equivalent to

      3(d1+d2)^2=4d1*d2.                                (2)

* For `4+3+1`, let `p,q,r` have multiplicities `4,3,1`, respectively.
  Exactness is equivalent to

      3/(p-q)+1/(p-r)=0,
      equivalently (p-q)+3(p-r)=0.                       (3)

* For `5+1+1+1`, translate its fivefold root to `u=x-p=0` and write

      N=u^5(a*u^3+b*u^2+c*u+d),
      a*d!=0,   a*u^3+b*u^2+c*u+d squarefree.

  These hypotheses include the full scalar in `N` and ensure the three
  simple roots are distinct and nonzero. Exactness is equivalent to

      b=c=0.                                            (4)

In particular the exact locus in (4) is nonempty: `N=x^5(x^3+1)` belongs
to it. It is an elliptic curve example with an exact differential of the
second kind, not a genus-zero reduction.

The closest recovered mechanism is the rational-map primitive-degree
lemma in [THM-2071, quadratic-fiber-square-parity-gate, Section 5](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md),
also retained by [THM-2723, split-exact-square-prefix-rational-primitive-
pole-capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
The present finite degree classification applies the same degree comparison
on the actual radical curve, then spends the remaining residues or the
complete elliptic primitive space. No literature-priority claim is made.

The five live concepts are labelled local residues, the global degree of a
primitive, genus-zero parametrization, the elliptic pole space, and the
separate asymptotic leading-coefficient map. The canonical hostile is
`N=x^4(x-1)^4`: its two residues cancel in their sum but each is nonzero.
The corrected near miss is to demand a primitive in `C(x)` when it only
belongs to `C(x)(sqrt N)`. The underused sidecar is the local degree of a
primitive at a zero of its derivative at infinity.

## 2. Local orders, connectedness, and the complete partition reduction

First suppose `N` is nonsquare, so its function field is a connected
quadratic extension. All statements are on its smooth projective
normalization. At a finite root `p` of multiplicity `m`, the local order
of `eta` and the maximum pole degree of a hypothetical primitive are:

| Multiplicity | Points over p | Order of eta at each point | Total primitive pole capacity |
| --- | --- | --- | --- |
| `m` odd | One, `x-p=z^2` | `1-m` | `max(m-2,0)` |
| `m` even | Two, parameter `x-p` | `-m/2` | `2 max(m/2-1,0)` |

The coefficients of these leading local terms are nonzero. The table
comes directly from `dx=2z dz` in the ramified case. A rational derivative
has no simple pole; if a primitive has pole order `e>=1`, its derivative
has pole order exactly `e+1`. Therefore any finite root of multiplicity
two immediately excludes exactness. At odd-multiplicity roots the local
series has only even powers of `z` multiplying `dz`; there is no residue.

At infinity the analogous calculation is

    n odd:   x=z^(-2),  ord_infinity eta=n-3;
    n even:  x=z^(-1),  ord_eta=n/2-2 at each of two points.  (5)

For `n>=3`, there is no pole of `eta` at infinity. Thus a hypothetical
nonconstant primitive has no pole there and has total map degree at most

    P = sum over finite multiplicities m>=3 of (m-2).       (6)

At the infinity point its local mapping degree is `n-2` for odd `n`,
or `n/2-1` at each of the two infinity points for even `n`. Each individual
local degree is bounded by the total degree of the same map. Hence

    P >= n-2             for odd n>=3,
    P >= n/2-1           for even n>=4.                     (7)

The two even-degree local requirements are **not added together**: the
two points need not have the same value under the primitive.

If odd `n>=3` has at least two distinct roots, then either all roots are
simple and `P=0`, or `P<=n-3`, contradicting (7). Thus only a pure odd
power can survive. For even `n=4,6,8`, removing roots of multiplicity two
and imposing (7) leaves, in the nonsquare case, precisely

    n=4:  3+1;
    n=6:  5+1, 4+1+1, 3+3;
    n=8:  7+1, 6+1+1, 5+3, 5+1+1+1, 4+3+1.                (8)

These lists can be checked directly by partitions; the companion source
enumerates every partition, including the rejected ones, not only (8).
There are `1,1,2,3,5,7,11,15,22` partitions in degrees zero through eight,
67 in all. All distinct complex root positions with any one partition
share the pole-order reduction. The later formulas retain their positions.

For `n=2`, a nonsquare has two simple poles at infinity by (5), so it is
nonexact. The square case has a finite multiplicity-two root and is also
nonexact. Degree one has a pole of order two at infinity and is exact;
it is not subject to (7). Nonzero constants are exact on `P1_x`.

If `N` is square, write `N=R^2` and work directly in `C(x)`. Multiplicity
two roots of `N` still give simple poles of `dx/R`. After their removal,
the only square patterns of degree at most eight are the pure powers
`4,6,8` and `4+4`. For the latter, writing
`N=kappa*(x-p)^4*(x-q)^4`, the residue at `p` is

    -2/[sqrt(kappa)*(p-q)^3] != 0.                        (9)

Thus `4+4` is excluded. This agrees with the full rational-primitive
lemma: if `1/R` has a rational primitive and `R` is nonconstant, then
`R` has exactly one root, of multiplicity at least two. The connected
quadratic capacity count is not used to prove this square case.

## 3. Explicit pure-power and genus-zero primitives

If `N=kappa*(x-p)^m`, with `m=1` or `3<=m<=8`, a primitive is

    2/[sqrt(kappa)*(2-m)] * (x-p)^(1-m/2).                (10)

For odd `m` it belongs to the radical field, as required; for even `m`
it is rational in `x`. If `N=kappa` is constant, use `x/sqrt(kappa)`.
The omitted exponent `m=2` is the logarithmic obstruction.

For every two-odd-root case in (8), write

    N=kappa*(x-p)^(2a+1)*(x-q)^(2b+1),
    k=a+b in {1,2,3},   d=p-q!=0,
    z^2=(x-q)/(x-p).

Then `x-p=d/(z^2-1)`, the function field is `C(z)`, and

    eta= -2/[sqrt(kappa)*d^k]
          * (z^2-1)^(k-1) z^(-2b) dz.                    (11)

Consequently an explicit primitive is

    -2/[sqrt(kappa)*d^k] *
    sum_{j=0}^{k-1} binom(k-1,j)(-1)^(k-1-j)
        z^(2j-2b+1)/(2j-2b+1).                           (12)

Every denominator in (12) is an odd, nonzero integer. This proves
exactness uniformly, including `3+3` and `5+3`; it is not a residue-only
guess. It also exhibits why the radical field cannot be discarded.

It remains to inspect the even root in the three exceptional genus-zero
partitions. Put

    u=x-p,   d1=p-q,   d2=p-r,
    D(u)=(u+d1)(u+d2)=u^2+s*u+t,
    s=d1+d2,  t=d1*d2,  v^2=D(u).                        (13)

Here `d1*d2*(d1-d2)!=0`. Each of these curves has exactly two odd branch
points and hence genus zero. At their odd roots there are no residues;
there are no poles at infinity. The two points over the even root have
opposite residues, so one exact scalar condition pays both. Alternatively
the explicit primitives below prove sufficiency without invoking the
residue criterion for `P1`.

For `4+1+1`, `eta=du/[sqrt(kappa)*u^2*v]`. Its residue vanishes exactly
when the coefficient of `u` in `D(u)^(-1/2)` vanishes, namely `s=0`.
This is (1). A primitive is

    -v/[sqrt(kappa)*t*u].                                (14)

For `6+1+1`, `eta=du/[sqrt(kappa)*u^3*v]`. The relevant coefficient is

    [u^2]D(u)^(-1/2)=(3s^2-4t)/(8*t^(5/2)).              (15)

It vanishes precisely under (2). In that case a primitive is

    v/sqrt(kappa) * [-1/(2t*u^2)+3s/(4t^2*u)].           (16)

For `4+3+1`, choose `q` to be the triple root. Now
`eta=du/[sqrt(kappa)*u^2*(u+d1)*v]`. The logarithmic derivative at zero
of the unit `(u+d1)^(-1)D(u)^(-1/2)` is

    -(d1+3d2)/(2d1*d2).                                  (17)

Thus the residue condition is (3). Under it, put

    B=-1/(d1*d2),   A=-2/[d1*d2*(d1-d2)].

A primitive is

    v*(A*u+B)/[sqrt(kappa)*u*(u+d1)].                     (18)

All denominators used here are nonzero under the stated distinct-root
hypotheses. Formulas (14), (16), and (18) are checked for all parameters
by the identity

    d(v*R) = [D*R'+D'*R/2] du/v.                         (19)

Rational-coefficient positive controls are respectively

    N=u^4(u^2+1),          primitive -sqrt(u^2+1)/u;
    N=u^6(u^2+2u+3),       primitive sqrt(u^2+2u+3)*(u-1)/(6u^2);
    N=u^4(u+3)^3(u-1),     primitive sqrt((u+3)(u-1))*(u+2)
                                    /[6u(u+3)].           (20)

The simple roots in the middle control are complex and distinct; no real
root-location assumption is present in the theorem.

## 4. Complete elliptic primitive space and its exact exceptional locus

Let `N=u^5(a*u^3+b*u^2+c*u+d)` have precisely the stated `5+1+1+1`
partition. Set

    X=1/u,       Y=y/u^4.

The exact curve and differential become

    Y^2=a+bX+cX^2+dX^3,       eta=-X^2 dX/Y.              (21)

The cubic is squarefree and `a*d!=0`: these follow from the original
distinct nonzero simple roots. Thus its smooth completion has one point
at infinity, four simple branch places of the degree-two map to `P1_X`,
and genus one. At finite points `dX/Y` is regular, including `Y=0` by
using `Y` as a local parameter. At infinity

    ord X=-2,       ord Y=-3,       ord eta=-4.

Therefore every rational primitive, if it exists, is regular on the
entire affine cubic and has pole order at most three at infinity.

This primitive space is exactly

    L(3 infinity)=span_C {1,X,Y}.                         (22)

To prove completeness directly, the affine cubic is smooth and its
coordinate ring is normal. A function without finite poles belongs to
`C[X,Y]/(Y^2-a-bX-cX^2-dX^3)` and has a unique expression
`A(X)+B(X)Y`. The possible pole orders of its two summands are even and
odd, respectively, so their top poles cannot cancel. The bound of three
forces `deg A<=1` and `deg B<=0`, proving (22). This pays all rational
primitives, not only a trial ansatz.

If a primitive is `k0+k1 X+k2 Y`, the involution `Y -> -Y` shows
`k1=0`, since `eta` is anti-invariant and `k1 dX` is invariant. The
remaining coefficient equation is

    (k2/2)*(b+2cX+3dX^2)=-X^2.                           (23)

It forces `k2=-2/(3d)` and `b=c=0`, and these conditions also suffice.
The primitive in the original field is

    -2Y/(3d) = -2*sqrt(N)/(3d*u^4).                      (24)

For example `N=u^5(u^3+1)` is exact. By contrast both
`u^5(u^3+u^2+1)` and `u^5(u^3+u+1)` are squarefree away from the
fivefold root and nonexact. Their differential has **zero residues
everywhere**: residue vanishing alone is insufficient in genus one.
The first failed implication is the attempted extension of the genus-zero
partial-fraction criterion without paying the primitive space.

## 5. Boundaries, connection map, and what is lost

This theorem is exhaustive for one algebraic differential with `deg N<=8`.
The source is the full labelled polynomial `N`, the target is the smooth
projective radical curve with its specified differential, and the map
adjoins the actual square root. The proof retains root positions,
multiplicities, residues, infinity orders, connected components, and the
field in which a primitive must lie. The partition filter alone loses
root positions; conditions (1)--(4) and the explicit formulas restore
exactly the missing data in the surviving strata.

Two scope hostiles are important.

* `F=x^3 t^2`, `G=1/(x^2 t)` has Jacobian one. Its leading radical
  differential `dx/x^(3/2)` has primitive `-2/sqrt(x)`. Requiring that
  primitive to lie in `C(x)` would incorrectly reject this actual rational
  pair. The extension field is load-bearing.
* `F=t^2+x^3` has constant leading coefficient, so the corresponding
  leading differential is exact. Nevertheless it has no rational mate:
  on the generic elliptic fibre `t^2+x^3=c`, the nonzero relative
  differential `-dx/(2t)` is holomorphic on the smooth compact curve.
  An exact rational differential with no poles would have a constant
  primitive. Thus a passed leading gate is not a supplier of a full mate.

The degree bound belongs to the affine polynomial `N`, not merely to an
unweighted divisor on `P1`. For instance `dx/x^2` is exact, whereas
`dx/[x^2(x-1)^2]` is not. Sending one finite root to infinity by a
projective substitution changes `dx` and must carry that weight. No
automorphism of the DG surface is claimed by an affine translation used
to express the formulas in `u=x-p`.

The repaired research direction is to combine this exact leading gate
with the actual formal coefficient map and then inspect higher
coefficients only in the surviving strata. Neither this note nor its
finite controls prove that any of those strata contains a full Jacobian
mate or a Keller counterexample. The general Jacobian conjecture remains
open.

## 6. Exact universe, reproduction, and audit boundary

The standalone source is
[planar_jc48_sep08_boundary_exactness.py](../../04-computation/planar_jc48_sep08_boundary_exactness.py).
It enumerates all 67 multiplicity partitions in degrees zero through
eight and compares each degree's complete survivor set with the theorem.
It verifies the local/infinity orders, every orientation of (12), all
three exceptional genus-zero primitives and residue factors, the
elliptic inversion and complete coefficient equations, and named
positive and negative controls. The formulas, not sampling, cover all
complex parameter values allowed in the theorem. No imported producer
or finite search over possible rational mates is used.

Reproduction, from the repository root:

```bash
python3 04-computation/planar_jc48_sep08_boundary_exactness.py
python3 -O 04-computation/planar_jc48_sep08_boundary_exactness.py
```

All gates use explicit exceptions rather than Python assertions, so the
optimized replay retains every check. The normal and optimized replays
are byte-identical: **195 gates**, 632 output bytes. The frozen output is
[planar_jc48_sep08_boundary_exactness.out](planar_jc48_sep08_boundary_exactness.out).

SHA256 pins:

* Source, 10,124 bytes:
  `bdc72806752cdfe31e2250b48cc6e26b19ead0bcadfc7fdaea887377ee579af0`.
* Output, 632 bytes:
  `c8ea8699d0d7f13eafa3961c017fea1bfe39d1c26a26cf231b609dc9d5105d66`.
* Semantic control manifest:
  `91100476475b4851378fa4b6a21cf69ecd002bdabc77cc6938a750397213564e`.

The [independent complete analytic/source audit](planar_jc48_sep08_boundary_exactness_audit.md)
accepts every partition, the square-field separation, all labelled residue
conditions and the complete elliptic primitive space. A separate
multiplicity-count enumeration reproduces all67 partitions, and independent
Fraction polynomial arithmetic verifies the exceptional primitives.
The source and output are frozen and independently accepted.
