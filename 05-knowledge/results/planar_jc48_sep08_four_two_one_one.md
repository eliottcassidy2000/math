# Every (4,2,1,1) quartic boundary placement

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The result concerns complete global square-prefix sections on the fixed
DG surface. It makes no claim that arbitrary Keller envelopes have this
form, and no conclusion about the full planar Jacobian conjecture.

## 1. Statement, inheritance and the leading entry

Keep the original surface and volume

    W=(P1_x x P1_z) minus {z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4b_D, D={r=0},
    omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L2`, `L in L1`, `F=H^2+L`. Suppose the leading
binary octic of `H` has four distinct zeros with multiplicities
four, two, one and one.

**Theorem.** A rational mate forces the double point
to be the original infinity, the finite quadruple point to be
the midpoint of the two simple points, and `L` to be constant.
No polynomial mate of any degree exists. Actual constant-`L`
rational mates attain the stated leading position condition for
every original finite quadruple-point position.

The closest mechanisms are [leading exactness](planar_jc48_sep08_leading_exactness.md),
the complete [degree-at-most-eight table](planar_jc48_sep08_boundary_exactness.md),
[shared-root jets](planar_jc48_sep08_shared_roots.md), and
[M-unit local forms](planar_jc48_sep08_quartic_common_root.md).
The complete global coefficient supplier is
[the DG quadratic filtration](planar_jc48_sep08_dg_quadratic.md).
The complete primitive-space method and nonsquare noncomposition
sidecar are inherited from
[the all-finite (4,3,1) proof](planar_jc48_sep08_four_three_one.md).
No pending result is a dependency.

The live concepts are leading positions, actual canonical weight,
two infinity regimes, complete Riemann--Roch spaces and global
regularity. The source-to-target map preserves the original generic
fibre and restricted volume. A list of local principal parts loses
affine points and the actual second chart; both are paid below.
The canonical hostile is the rational family in Section 6. The
least-used sidecar is the degree drop of `F|D`: it changes the
infinity branch at exactly the same parameter value, preserving
the genus. Treating only the generic infinity coefficient would
miss an entire residual stratum.

If the double point is finite, `du/sqrt(N)` has a nonzero
logarithmic residue there, so leading exactness excludes every
rational mate. With that point at infinity, `deg N=6`. The
proved degree-six `(4,1,1)` row of the exactness table requires
the quadruple point to be the midpoint of the other two roots.
No rational mate exists when this position condition fails.

Actual source/surface scalings then normalize the nonzero
half-separation and leading scalar. Put `u=x-p`, retaining
arbitrary `p`, and obtain

    N=u^4(u^2-1).

This is an ordinary polynomial translation for calculation, not
an arbitrary projective or surface automorphism. The leading
condition is sharp: `d(sqrt(u^2-1)/u)=du/[u^2 sqrt(u^2-1)]`.

## 2. Complete global rows and necessary active jets

The full original global coefficient box gives

    P=sum_(i=0)^4 A_i u^i,
    Q=(A4-1)u^2+(A3-2p A4+4p)u+C0,
    M=sum_(i=0)^4 B_i u^i,
    R=B4 u^2+(B3-2p B4)u+E0,
    H=N t^2+P t+Q, L=M t+R.                         (1)

Indeed the three numerator rows `Q,P-2Qx^2,N-Px^2+Qx^4`
have degree at most four in `x=u+p`. The last row's degree-six
and degree-five coefficients force precisely the two displayed
coefficients of `Q`; the other rows then have the required
degrees. Conversely the complete section box forces
`deg Q<=2,deg P<=4`. The linear numerator rows `R,M-Rx^2`
have degree at most two and force precisely `R` in (1).
All remaining parameters, including `p`, are free at entry.

The two finite simple points have regular relative forms by
the complete M-unit/shared-root lemmas. Infinity has leading
multiplicity two, and its extra canonical factor `r^2` makes
every branch regular for every coefficient in (1): the normal
and M-unit cases are regular already, and the shared tangent
quartic has only logarithms before weighting. Thus all possible
primitive poles lie over the finite quadruple point.

If that point has an M-unit, `m=4` is never its balanced value
`3j`, so all forms are regular. A shared normal unit is regular
as well. Otherwise the finite first-jet obstruction forces

    P(0)=P'(0)=M(0)=M'(0)=0.                         (2)

Failure of the value conditions makes a primitive holomorphic
and constant on each compact generic component, contrary to the
actual nonzero relative form. Once those values vanish, failure
of the first jets gives a nonzero logarithmic residue. Therefore
(2) is necessary under a rational mate.
This step does not assume geometric integrality.

The complete inverse coefficient `T2=-M/(4N)` and normalized
trace require `M du/N` rational-exact. Its residues at the two
simple roots force `M(1)=M(-1)=0`. Combined with (2) and the
full degree-four bound this gives

    M=lambda u^2(u^2-1),
    R=lambda u^2-2p lambda u+e.                     (3)

If `lambda=0`, the complete row makes `L` constant. Assume
`lambda!=0`, and absorb `e` into the formal fibre constant.
Unlike a mere zero-count contradiction, the survivor in (3)
really passes this exactness gate: `M/N=lambda/u^2` has rational
primitive `-lambda/u`.

The entire remaining family is consequently

    P=a u^4+b u^3+c u^2,
    Q=(a-1)u^2+(b-2pa+4p)u+d,
    H=N t^2+P t+Q,
    L=lambda[u^2(u^2-1)t+u^2-2pu].                  (4)

No assumption that `F|D` is constant has been made; in fact its
degree is positive throughout this family.

## 3. The full canonical divisor, including the infinity degeneration

As in the cited all-finite (4,3,1) proof, a nontrivial polynomial
composition of the quartic in `t` has outer degree two or four.
Completing a quadratic outer polynomial to a square and comparing
the unique square prefix forces `L` constant. An outer quartic
would make `N` a polynomial square, contradicted by the two
simple roots. Thus for `lambda!=0` the classical closed-polynomial
criterion gives geometric integrality of the generic fibre.
Use its compact normalization over the algebraic closure of
`C(zeta)`.

At the active finite quadruple point, put `s=u^2 Z`, where
`s=z-x^2`. The complete leading equation is

    (-1+c Z+d Z^2)^2-lambda Z^3-zeta Z^4.            (5)

It has four simple nonzero roots for generic `zeta`. The
repeated-root eliminant `ZP0'-4P0` has constant minus four,
and its leading coefficient is `d^2-zeta`. All four branches
are actual unramified branches and exhaust Weierstrass degree
four. Each has relative-form order `4-6=-2`, so there are four
possible primitive simple poles. Let their sum be `E`, degree four.

At either finite simple point, (3) makes the lower leading
section vanish. There are two simple generic branches after
centering the numerator of `H`. If its normal coefficient is
nonzero then `u-u0~s`; otherwise `u-u0~s^2`. In either case
the derivative of the full equation with respect to `u` has
order two in `s`, as does the numerator of `eta=s^2 ds/E_u`.
The form is a unit. Also `L` and therefore `H` are bounded
on those branches.

At the original infinity, the transformed lower section is a
unit: in the inverted boundary coordinate its value is
`-lambda!=0`. The normal coefficient
of `H` there is `2-a`, as follows from the actual second chart.
There are two distinct regimes, both of which are necessary:

* If `a!=2`, the normal-unit cancellation has two branches
  `s~r^2`, with splitting order three in `r`. The full equation
  derivative with respect to `s` has order three. The actual
  form `r^2 s^2 dr/E_s` has order `2+4-3=3` on each branch.
  The polynomial `F|D` has degree two in `b_D`, because its
  quadratic coefficient is `(2-a)^2`. Its two generic transverse
  points each contribute order two from the volume factor.
* If `a=2`, the lower-unit Newton edge has three determinations
  forming one actual normalized branch, with
  `r=tau^3,s=tau^4`. The normal term has strictly higher
  weight. The form has order `6+8+2-8=8`. Now `H|D` is
  constant and `F|D` has degree one, with coefficient
  `-lambda`; its single generic transverse point contributes
  order two.

These local faces exhaust the respective Weierstrass degrees.
The nonzero unit coefficients and the generic fibre parameter
exclude further degenerations. On the smooth generic source
part the volume and relative differential are units. Thus the
canonical degrees in the two cases are respectively

    2*3+2*2-4*2=2,    8+2-4*2=2.

The compact generic curve has genus two for every parameter
in (4) with `lambda!=0`.

## 4. The complete three-dimensional primitive space

Since `deg E=4>2g-2`, Riemann--Roch gives `dim L(E)=3`.
A complete basis is

    1, U=1/u, V=H/u.                               (6)

At each active branch `u` is a parameter and `H` is bounded,
so both nonconstant candidates have poles of order at most one.
The line `u=0` is the constant fibre `F=d^2`, hence it has no
ordinary affine points on the generic fibre. This checks all
affine denominators, not only their boundary principal parts.
At the two simple points `u` is nonzero and `H` is bounded.
On the retained divisor `D`, both `U` and `V` vanish, since
`1/u` has a simple zero and `H` is global there.

Finally, at the deleted infinity point the identity
`H^2=zeta-L` gives `H~s^-1/2`, since the transformed lower
section is a unit. In the first regime `H` has order minus one
in `r`, canceled by `1/u~r`. In the second regime it has order
minus two in `tau`, while `1/u` has order three. Thus `V`
is regular in both regimes, and so is `U`.

The three functions in (6) are linearly independent even over
the algebraic closure of `C(zeta)`: a relation multiplied by
`u` would be a polynomial of degree at most two in `t`, below
the irreducible fibre degree four. Its coefficient of `t^2`
forces the coefficient of `V` to vanish; then the independence
of `1,u` forces the other two. Membership plus the exact
dimension proves completeness, rather than a bounded search.

## 5. Two bracket coefficients exclude every rational mate

A hypothetical primitive is `C(zeta)+A(zeta)/u+B(zeta)H/u`,
even allowing coefficients in the algebraic closure of `C(zeta)`.
Reduce its Jacobian with `F` modulo the original equation
`F-zeta`, keeping the entire formal parameter. The coefficient
of `t^3` is

    4A u^6(u^2-1)^2.

A constant bracket forces `A=0`. With that substitution, the
coefficient of `t^2` is

    B lambda u^4(u^2-1)(3-u^2).

Since `lambda!=0`, this forces `B=0`. The resulting primitive
is constant on the fibre and has zero Jacobian, not one. This
contradiction excludes every rational mate for nonconstant `L`.

For constant `L`, `J(F,G)=2H J(H,G)` excludes polynomial
mates because `H` is nonconstant. The finite-double and
nonmidpoint leading obstructions in Section 1 cover every
remaining boundary location and position.

## 6. Rational sharpness and exact controls

For every original finite position `p`, put

    q=1+u^2t,
    H=(u^2-1)q^2+d,    G_H=-1/(2u q).

These satisfy the original global rows and `J(H,G_H)=1`.
The leading polynomial is exactly `u^4(u^2-1)`, so
`H^2+e` has rational mate `G_H/(2H)` in the claimed boundary
partition. Its denominator is essential.

The [exact companion](../../04-computation/planar_jc48_sep08_four_two_one_one.py)
controls the full all-p coefficient equations, the necessary
residues and their nonzero exact survivor, both actual infinity
regimes, the complete bracket coefficients and literal rational
examples. It does not enumerate a bounded degree of possible
mates. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_four_two_one_one.py
python3 -B -O 04-computation/planar_jc48_sep08_four_two_one_one.py
```

Both modes pass **41 gates** and reproduce the identical
[283-byte frozen output](planar_jc48_sep08_four_two_one_one.out).
The source is 4,581 bytes, SHA256
`be559b1eb0b15a8a55179d72e9defb01fb62cec8d4341f2b1a9c0f6ba1318805`.
Output SHA256:
`a6a4dcf34b277956677a86916a1d224becff92b514fc7837ca1ca0f896b95fc1`.
Semantic gate SHA256:
`6410484364b2dc0ab1c3412baed101976b764e4c4410567f7dc2957a0ee94818`.

The [independent full audit](planar_jc48_sep08_four_two_one_one_audit.md)
accepts the complete coefficient class, both infinity regimes, the
full genus-two primitive space and both bracket obstructions. Separate
original-chart reconstructions and both41-gate replays agree.
