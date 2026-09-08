# Independent audit of the fixed finite-six, constant-D quartic exclusion

**Status: PASS — independent analytic, exact-source, and normal/optimized/frozen-output audit.**
The theorem accepted here is the precise fixed-point theorem in
[the primary note](planar_jc48_sep08_const_d_quartic.md). Its parameters
are arbitrary complex constants subject to `a*m0!=0`; a proposed mate
may be any element of `C(x,t)` with nonzero constant source Jacobian.
This is not an exclusion of every common-root quartic or a change of
coordinates moving an arbitrary sixfold point to zero.

## 1. Complete section space and actual surface

I independently recovered the section-space reduction before reading
the final exact rank controls. Start with an arbitrary polynomial in
the full `O(4,2)` box, restrict it to `z=x^2`, and impose the boundary
section `a*x^6`. Writing `s=z-x^2` gives exactly

    N=a*x^6+s[B0+B1*x+B2*x^2+C1*x^3+(a+C2)*x^4]
                   +s^2[C0+C1*x+C2*x^2].

In the actual boundary chart `x=1/r`, `t=-r^2-r^4*b_D`, its
restriction is

    H|D=(a-C2)*b_D+C0-B2.

The general `O(2,1)` section similarly gives

    M=m0+m1*x+d*x^2+n1*x^3+n2*x^4
                                      +s(n0+n1*x+n2*x^2),
    L|D=-n2*b_D+n0-d.

Since `deg H|D<=2` and `deg L|D<=1`, a constant `H|D^2+L|D`
forces both restrictions to be constant. Thus `C2=a`, `n2=0`;
the resulting spaces have respectively six and five parameters.
They are precisely the spaces displayed in primary (1), including
the two lower normal coefficients `beta0,beta1`. The final source
also checks spanning, independence, and every linear constraint in
the full fifteen- and six-monomial section boxes. The proof does not
substitute a parameter sample for completeness.

The source and target of this reduction are sections and actual
global functions on the fixed surface `W=(P1 x P1)\{z=x^2}`. It
preserves the full functions, both surface charts, and the literal
volume form. In particular, no translation of `x` is used. Such a
translation with `t` fixed generally fails to extend across `D`.
The finite multiplicity-six point here really is `(x,z)=(0,0)`;
the other octic zero has multiplicity two at infinity.

## 2. Complete infinity supplier and lower local cases

I checked both complete infinity numerators, including the extra
terms contributed by `beta0,beta1`. They do not change the quadratic
part of `Ni`. The tangent quartic is exactly

    [a-c1*Z+(c0-b)*Z^2]^2-n1*Z^3+(n0-d-c)*Z^4.

Its nonzero constant term is `a^2`; its generic leading coefficient
is `F|D-c`. The repeated-root eliminant `Z*T0'-4*T0` has constant
term `-4a^2`, independent of `c`. Hence only finitely many original
fibre values can fail to have four simple nonzero tangent roots.
Each root lifts to an actual smooth branch `sigma=r*Z(r)`, and the
degree-four Weierstrass restriction at `r=0` proves exhaustion.
The actual factor `r^2` in the volume form gives order **one** for
the relative differential on all four branches.

The same-component use of this zero is valid. On every generic
component, `x` is nonconstant: for fixed nonzero `x`, the leading
coefficient of `F` in `t` is `a^2*x^12`; at `x=0`,

    F=(beta0*t+c0)^2+m0*t+n0

is also nonconstant. A compact component therefore has a pole of
`x`. It cannot meet the finite part of `D`, where `F` is constant,
so it meets one of the four infinity branches. A rational primitive
of the nonzero relative differential must consequently have local
degree two there and pole degree at least two on that same component.
This statement needs no generic irreducibility assumption.

The inherited local M-unit lemma is correctly used with `m=6`.
The `beta0!=0` case has `j=0` and local Weierstrass degree two;
`beta0=0,beta1!=0` has `j=1`. Both are regular local cases. After
both vanish, `b=0` gives `j>=3`, again regular. The balanced
`b!=0` case has degree three and cubic face

    (a+b*Z)^2+m0*Z^3.

Its discriminant and double-face parameterization in the primary
are correct. A simple face is regular; the remaining simple branch
at `Z=q/4` is also regular. At the double face, the critical-value
derivative `R_c=-x^4*z(x,c)^4` ensures generic order at most four.
Normalized orders `1,2,3,4` give total possible primitive pole
budgets `0,0,1,2`, respectively. In order two there are two actual
simple poles with nonzero residues, already inconsistent with
exactness. In order three the single ramified branch has a double
pole, with primitive budget one. Thus every earlier layer fails
the infinity degree-two gate. Neither Puiseux determinations nor
local degrees on different components are improperly added.

Two full-family prose errors were identified during review and
are repaired in the accepted version: the value of `F` at `x=0`
retains `beta0`, and the degree-three Weierstrass claim is made
only after the lower normal coefficients vanish. The final source
includes literal checks of both corrections. No remaining repair
was requested.

## 3. Independent critical-centre and residue reconstruction

In a separate scratch calculation, I formed `E=N^2+s^3*M-c*s^4`
directly from the literal sections and set `s=x^4*Z`. I then solved
the critical-value coefficient of order one for `m1`, the critical
centre coefficient of order one for `z1`, and repeated this process
at orders two and three. This recursion did not import the producer
or begin from its displayed tuning formulas. It recovered exactly

    m1=-4*b*c1/(3*q),
    d=4*(2*b^2*q-3*b*c0*q-3*c1^2)/(9*q),
    n1=4*c1*(b^2*q-3*b*c0*q-c1^2)/(9*b*q),

as well as the stated `z1,z2,z3,r4,r5`. Thus the equations are
successive necessary conditions for critical order four, not just
one convenient family of sufficient tunings.

The residue formula is also independently checked. If

    Q=R+h*(Z-z)^2+O((Z-z)^3),
    R=r4*x^4+r5*x^5+O(x^6),  h=h0+h1*x+O(x^2),
    Z=z+e*x^2+f*x^3+O(x^4),

then the actual equation gives

    e^2=-r4/h0,
    f/e=-h1/(2*h0)+r5/(2*r4).

In this chart the relative form is exactly `Z^2*dx/Q_Z`.
Its residue divided by its nonzero double-pole coefficient is

    2*z1/q-h1/(2*h0)-r5/(2*r4).

The cubic term in `Z-z` starts too late to change this coefficient.
My recursive calculation gives the primary's numerator and its
derivative `2*c1*q^4/b` with respect to the **original** generic
fibre value. Hence generic vanishing of residues forces `c1=0`,
then `m1=n1=0` and `d=8*b^2/9-4*b*c0/3`. The even family is
indeed the entire remaining layer. Its two actual branches are
even series in the parameter `x`, so their residues vanish.

## 4. Literal elliptic transport and complete primitive space

I independently reconstructed the two birational coordinate changes
and their volume-form coefficients. With

    u=1/x,  v=x+x^3*t,  w=v-u/q,  z=v*w,
    u=q*(z-w^2)/w,      v=z/w,

the exact remaining original fibre is

    A(z)*w^2=P2(z),
    A=4*a^2*z+C,        P2=3*a^2*z^2-B*z+c-K0,

with the same `B,C,K0` as the primary. The actual forms are

    omega=-q^2*(z-w^2)/w^2 * dz wedge dw,
    eta=q^2*(Az-P2)/(2*A^2*w^3) * dz.

These signs agree with `dF wedge eta=omega`. In particular the
transport preserves the differential, not just the equation of an
abstract genus-one curve.

Over the algebraic closure of the generic constants field, `P2/A`
has four simple odd valuations: its two distinct numerator zeros,
the distinct denominator zero, and infinity. The nonzero fibre-value
slopes in the exact source pay all generic distinctness conditions.
This proves geometric integrality and, by the degree-two map to
the `z` line, genus one. The coordinate change is birational on a
dense part of the actual generic curve; no component is lost in
inverting `x` or `w`.

The complete polar support agrees by direct local calculation.
At the two zeros of `P2`, use `w` as parameter: `dz` has order one
and the form has order minus two. At the root of `A`, a parameter
has `ord(z-z_A)=2`, `ord(w)=-1`, so the form is regular. At
infinity, `ord(z)=-2`, `ord(w)=-1`, again giving regularity.
There are no other poles. Thus an exact primitive can have only
simple poles at the two points `w=0`.

I checked the elementary involution proof of the complete allowed
function space. Its even part is a rational function of `z`;
branch-point poles would have even order and are forbidden, so
this part is constant. The odd part is `w*S(z)`. The only possible
poles of `S` are simple poles at the two zeros of `P2`; it must
vanish at the zero of `A` and at infinity. Consequently
`S=K*A/P2`, and the entire space is exactly

    G=K/w+constant.

Equivalently it is the two-dimensional space of functions with at
most those two simple poles on this genus-one curve. The argument
allows `K` in the enlarged generic constants field, so it does not
silently assume constants of a proposed primitive are complex
numbers independent of `c`.

Exact differentiation now requires

    q^2*(Az-P2)=-K*(P2'*A-P2*A').

The quadratic coefficient forces `K=-q^2/(12*a^2)`, independently
of `c`. With this value, the remaining constant coefficient has
nonzero `c` slope `-2*q^2/3`. This is a genuine generic identity
obstruction even in the enlarged constants field. Residue freedom
and the permitted pole budget therefore do not supply exactness.

## 5. Exact replay, controls, and limits

I read the full final source. Its checks are always active under
`-O`, it imports no inherited mathematical implementation, and its
polynomial identities retain symbolic parameters. Its finite
controls include all section boxes, both surface charts, every
critical layer, the asymmetric order-four residue obstruction,
the even second-kind case, actual birational and differential
identities, and named rational-mate hostiles outside the theorem.
In particular `h^4+h` with `h=x^2+x^4*t` has the wrong octic
multiplicity for this theorem, while an affine critical point alone
does not rule out a rational mate. Those boundaries are kept distinct.

Independent commands, from the repository root:

```sh
python3 -B 04-computation/planar_jc48_sep08_const_d_quartic.py
python3 -B -O 04-computation/planar_jc48_sep08_const_d_quartic.py
```

Both independent outputs and the frozen output are byte-identical:
**97 always-active exact gates, 462 bytes**. Reviewed pins:

- [Primary proof](planar_jc48_sep08_const_d_quartic.md), before audited-status
  promotion: 17,250 bytes, SHA256
  `9c8b08a077c0625d00c3fd4760f15b20d54c97d16e4f62daa4a0717d91d33827`.
- [Exact source](../../04-computation/planar_jc48_sep08_const_d_quartic.py):
  13,598 bytes, SHA256
  `cfc977adb5dbfb0b1eb3f0118669551bd40dde758d2bab112734b9bf4f5a2194`.
- [Frozen output](planar_jc48_sep08_const_d_quartic.out): SHA256
  `7b38b00e6ac5c6282122cfafb3d3fcfbc99c05cb3bb1fd52d299e0d85919e92f`.
- Semantic digest:
  `7602857f67214f96f1000d116cdf16c3d80036bb0fb83ea44daf64655de095d6`.

The accepted result closes the specified fixed finite-six,
constant-D, M-unit family. It neither normalizes other boundary
points nor asserts a general common-root exclusion. The earlier
shared-root stopping example is correctly identified as logarithmic;
the genuinely second-kind tuning is handled by the exact elliptic
primitive-space obstruction. No source or frozen predecessor was
edited in this audit.
