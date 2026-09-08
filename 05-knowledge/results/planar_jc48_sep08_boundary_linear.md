# Every global coordinate linear in the boundary-chart base has no Keller mate

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; controls FINITE-EXACT.**
The statement below has no bound on source fibre degree or on the degree
of a proposed mate. It concerns the specified surface and chart; JC(2)
and the general quartic stratum remain OPEN.

## 1. Statement and inheritance

Use the surface and full charts established in
[the DG filtration](planar_jc48_sep08_dg_quadratic.md):

    W=(P1_x x P1_z) minus {z=x^2},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^2-r^4 b, b=-x^2(1+x^2t),
    v=x(1+x^2t), rb=-v,
    omega=dx wedge dt=r^2 dr wedge db.

Here D={r=0} is the added boundary, whereas E={x=0} is a curve
inside the original affine source U0, outside the other chart.

**Claim.** If a global F on W has degree at most one in r in Uinf,
then no G in C[x,t] satisfies J_(x,t)(F,G)=lambda for any nonzero
lambda in C. In particular no global mate exists. This includes all

    F=r b C(b)+B(b),                  C,B in C[b],       (1)

without any bound on their degrees. Formula (1) is the entire stated
class, not a sufficient family chosen from it.

The closest proved mechanism is the rational-primitive pole-degree
argument in [THM-2071 / quadratic-fiber-square-parity-gate,
§5](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md),
combined with the special-fibre repair obstruction in
[the source-linear carrier](continuing10_20260907_dg_linear_carrier.md).
The operations here exchange the base and fibre in the second full
chart. The Jacobian becomes lambda r^2, so the numerator in the relative
primitive is quadratic and its denominator is A^3. The globality
sidecar b|A then identifies the only possible pole support. No theorem
requiring a constant Jacobian in that chart is silently applied.

The live concepts are full-chart regularity, rational primitive residues,
the degree of a rational map, and principal parts on two components of
one fibre. The hostile is F=b^4-v: it has a global quadratic approximate
root, no affine-source critical point, and an actual rational Jacobian
mate, but no polynomial mate. The corrected near miss is to infer a
regular mate from rational integrability or from a global square prefix.
The least-used sidecar is the other component of the special fibre.

## 2. Globality and the relative primitive

Write F=A(b)r+B(b) with A,B in C[b]. Substitution on U0 gives
A(-x^2(1+x^2t))/x+B(-x^2(1+x^2t)). Its only possible negative power
of x is A(0)/x. Thus globality is equivalent to A(0)=0, proving (1).

If A=0, both source derivatives of F=B(b) vanish on E, because those
of b=-x^2-x^4t do. This rules out any nonzero constant Jacobian.
Assume henceforth A!=0. The independent pair F,b generates C(r,b),
with r=(F-B)/A. If a rational mate G exists, the chain rule gives

    (partial G/partial b)|_F
        =lambda (F-B(b))^2/A(b)^3.                       (2)

In particular differentiating twice with respect to the independent
variable F shows that 1/A^3 has a rational primitive in C(F)(b).

Here is the complete pole-degree argument. If a polynomial A of degree
M>=1 has s distinct roots of multiplicities m_i, a rational primitive
of 1/A^3 has poles of orders 3m_i-1 at precisely those roots. It has no
polynomial part of positive degree, since its derivative vanishes at
infinity. Subtract its value at infinity. The resulting nonzero
rational map has degree 3M-s and a zero of order 3M-1 at infinity.
The latter cannot exceed its degree, so s<=1. This argument works over
an algebraic closure of C(F); it does not require G to be polynomial.
Since b divides A, necessarily

    A=a b^m,                   a!=0, m>=1.              (3)

If m>=2, the term r b^m has x-order 2m-1>=3 at E, while B(b)-B(0)
has x-order at least two. Both source derivatives of F therefore vanish
on E. Only A=a b can survive even this elementary source test.

## 3. Exact rational-integrability boundary

Put g=F-B(0) and h(b)=B(b)-B(0). With A=a b, equation (2) becomes

    G_b|_g =lambda (g-h(b))^2/(a^3 b^3).                (4)

Write h=h_1 b+h_2 b^2+.... The residue at b=0 is

    (lambda/a^3)(h_1^2-2g h_2).

Because g is transcendental, a rational primitive exists if and only if
h_1=h_2=0. Necessity is the zero-residue criterion. Sufficiency is given
by the explicit primitive below. Thus the remaining rationally
integrable class is precisely h in b^3 C[b].

For this class define polynomials with zero integration constants

    L(b)=integral h(b)/b^3 db,
    M(b)=integral h(b)^2/b^3 db.

Then

    G0=(lambda/a^3)[-g^2/(2b^2)-2g L(b)+M(b)]           (5)

satisfies the desired rational Jacobian equation. Substitute g=a rb+h
to obtain

    G0=-lambda r^2/(2a)+Q(b,rb),                       (6)

where Q is a polynomial: in the expansion of (5) the only potentially
problematic term besides r^2 is r h/b=(rb)(h/b^2), and h/b^2 is a
polynomial. In particular G0 is regular on Uinf and has an exact double
pole, with coefficient -lambda/(2a), along E in U0.

Every other rational mate is G0+H(F) for H in C(T). Indeed the kernel
of the nonzero rational derivation (a b/r^2) partial_b on C(F)(b) is
exactly C(F), in characteristic zero. This also specifies the complete
freedom available for a pole repair; arbitrary functions of the source
are not allowed.

## 4. The same fibre blocks every pole repair

The fibre F=B(0) contains the two distinct curves

    E={x=0},
    Gamma={1+x^2t=0} in U0, equivalently {b=0,r!=0}.

The order of F-B(0) along each is exactly one. Along E, rb=-x+O(x^3)
and h(b) has order at least six, so the leading coefficient is -a.
At the generic point of Gamma, g=b(a r+h/b) and a r is a unit.

To cancel the double pole of (6) at E, H must have a pole of exactly
order two at the value B(0). But then H(F) has a double pole along
Gamma as well, where G0 is regular. There is no cancellation there.
Consequently G0+H(F) cannot be polynomial on U0. This proves the claim
for the last possible case and completes all of (1).

The failure is not absence of a rational primitive. It is incompatible
principal parts on two reduced components of one actual fibre. Other
components of that fibre, when present, are unnecessary for the proof.

## 5. A sharp quartic control and the remaining question

Take a=lambda=1 and h=b^4. Then

    F=rb+b^4=b^4-v,
    H=b^2,              F-H^2=-v in L_1,
    G0=-r^2/2-2rb^3-(4/3)b^6,
    J_(r,b)(F,G0)=r^2,  J_(x,t)(F,G0)=1.

Both F and H are global; F has source t-degree four. On Uinf the only
critical point of F is (r,b)=(0,0), which is on D. On E, F_x=-1.
Hence F has no critical point on the original affine source. Its
boundary restriction F|D=b^4 is nonconstant. Nevertheless (6) has a
double pole on E, and the special-fibre proof rules out every polynomial
mate. The boundary point is not falsely declared noncritical on W.

More generally F=a rb+B(b) is source-critical-free exactly when B_1=0:
on Uinf the possible critical point is b=0,r=-B_1/a, and the source
derivative on E is -a. The B_2!=0, B_1=0 subcase is already excluded by
the rational residue, while B_1=B_2=0 needs the two-component obstruction.

This connects directly to [the global quartic root
theorem](planar_jc48_sep08_quartic_boundary.md), without using it as
a dependency. Even when its proposed conclusion H in L_2 and
F-H^2 in L_1 is granted, this example shows why a terminal regularity
argument is still necessary. The present theorem covers the entire
class linear in r, including arbitrarily high source t-degrees; it
does not cover a general H^2+L, or a coordinate with degree at least
two in r. The next question is whether an analogous relative-primitive
description can preserve the distinct component valuations for those
larger classes.

## 6. Reproduction and scope

Run `python3 04-computation/planar_jc48_sep08_boundary_linear.py` and
the same command with `python3 -O`. The finite-exact controls check the
literal chart, the full symbolic residue and primitive formulas, bounded
families for globality and criticality, and the quartic hostile by two
Jacobian calculations. They do not infer the unbounded theorem from a
finite search. The pole-degree argument and the two component valuations
are analytic proofs above. The [independent analytic/source audit](planar_jc48_sep08_boundary_linear_audit.md)
passes, with 149 always-active gates in each mode and byte-identical output.

Source SHA256: `aa57669f5b20202adb7e9e39e10b2c7d600a32d188e1303963dbf68e0d068d53`.

Output SHA256: `30a014d4d7850aae83ce1b54c13b026236916e79edb0dd7d667168adb690b193`.
