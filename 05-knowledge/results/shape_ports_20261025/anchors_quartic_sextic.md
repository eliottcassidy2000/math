# The nonzero-defect quartic/sextic source lane is excluded

25 October 2026. **PROVED elementary polynomial obstruction, independently audited.** See [the independent audit](anchors_audit.md). Exact symbolic and finite controls are supplemental. This closes one lane of the [genus-six source frontier](../genus_six_source_frontier_20261020/FRONTIER.md), not the planar Jacobian conjecture and not the whole generic conductor-twelve boundary.

## 1. Exact claim and inherited frontier

Let k be a characteristic-zero field. Let

    P = y^4 + a(x)y^2 + b(x)y + d(x),
    Q = y^6 + c(x)y^4 + B(x)y^3 + e(x)y^2 + h(x)y + f(x)

belong to k[x,y]. There is no such pair with

    J = P_x Q_y - P_y Q_x in k*,
    beta = 2B - 3b != 0.

The constant-Jacobian equations themselves force beta to be constant. Write it as kappa below, reserving k for the field. In fact the proof gives a sharper statement: if J lies in k[x] and beta!=0, then either J=0 or deg_x(J) is congruent to 2 modulo 3. Every positive degree in that congruence class occurs in the explicit constant-a controls below. The missing y^5 coefficient in Q is a hypothesis. The inherited genus-six frontier proves it in its normalized monic quartic/sextic source lane; this note does not assume it for every degree-six coordinate.

The closest proved mechanism is the audited genus-five quartic/sextic first-integral obstruction in [the conductor-ten packet](../conductor_ten_source_exclusion_20261020/CONDUCTOR_TEN_SOURCE_EXCLUSION.md). The corrected near miss is setting beta=0 at generic codimension six. The canonical hostile is P=y^4, Q=y^6+y^3+y: its generic one-variable algebra is immersed and has semigroup <4,6,9>, but its planar Jacobian is zero. The least-used sidecar is the degree at infinity of the full constant coefficient of J, after retaining both scalar terms.

The current boundary before this note is explicit in genus-six FRONTIER sections 6--8. This proof handles its beta!=0 branch. Combined with the separately inherited beta=0 exclusion, the terminal (4,6) pair is removed in that setting. The pairs (3,7), (4,5), (4,10), and (6,9) remain obligations. No conductor-fourteen theorem follows.

## 2. Integrate the top rows, retaining actual scalar terms

The coefficients of y^7 down through y^3 give constants C,kappa,E,H,F in k such that

    c = 3a/2 + C,
    B = (3b+kappa)/2,
    e = 3d/2 + 3a^2/8 + C a + E,
    h = 3ab/4 + C b + 3kappa a/8 + H,
    f = 3ad/4 + C d + 3b^2/8 + 3kappa b/8 - a^3/16 + E a/2 + F.

The legal constant target shear Q -> Q-C P-F removes C and F and preserves J exactly. From now on C=F=0. No x-dependent scalar is subtracted.

Put

    I = 3bd + (3kappa/2)d - (3/4)a^2b - (3kappa/16)a^2 + 2Eb + Ha,
    K = (3/2)d^2 - (3/4)a^2d + 2Ed - (3/4)ab^2
        - (3kappa/8)ab + Hb + (3/32)a^4 - (E/2)a^2.

Direct coefficient expansion gives the exact identity

    J = I' y^2 + K' y + h d' - b f'.

Thus J independent of y forces I=L and K=M for constants L,M in k. The actual scalar Jacobian is still

    J = h d' - b f'.                                      (2.1)

## 3. The exact plane relation

Write

    r=2b+kappa, U=L+E kappa.

Eliminating d gives N(a,b)=0, where

    N = 9a[kappa^2 a^3 - 16b r^3]
        + 96kappa a^2(Ha-U) + 256(Ha-U)^2
        + r^2(384Hb - 256E^2 - 384M).                      (3.1)

For example the exact Sylvester resultant is

    Res_d(I-L,K-M) = 3N/512.

This is a necessary relation only. Its zeros are not automatically Keller pairs. Conversely, in the r!=0 fraction-field chart, I=L gives exactly

    d = a^2(4b+kappa)/(8r) - 2E/3 - 2Ha/(3r) + 2U/(3r).
                                                                  (3.2)

## 4. Constant-coordinate branches

If b is constant and a were nonconstant, (3.1) has unique highest term 9kappa^2 a^4, impossible. So a is constant. If r!=0, I=L then makes d constant. If r=0, K=M is a quadratic polynomial equation in d with nonzero leading coefficient 3/2, again making d constant. Therefore J=0.

Suppose a is constant and b is nonconstant. Using 4b+kappa=2r-kappa, equation I=L makes d=a^2/4-2E/3+(-kappa a^2/8-2Ha/3+2U/3)/r, a constant plus a constant divided by r. Since d is polynomial and r is nonconstant, that numerator vanishes and d is constant. The y coefficient K'=0 becomes

    b' (H - (3/2)ab - (3kappa/8)a) = 0.

Characteristic zero and nonconstant b imply a=H=0. Then d=-2E/3, h=0, and

    J = -(3/8)b(2b+kappa)b'.

Its degree is 3 deg(b)-1>=2, so it cannot be a nonzero constant.

This disposes of zero polynomials as well as nonzero constants; the remaining case has a and b both nonconstant.

## 5. Forced degrees and the rational comparison curve

Let A=deg(a)>0 and B0=deg(b)>0. In N the two possible top degrees are 4A and A+4B0. Every other term has strictly smaller degree than their maximum: in particular a^3 is below a^4, and b^3 is below ab^4. Their leading coefficients are nonzero since kappa!=0. Hence they must tie:

    3A=4B0,
    deg(a)=4m, deg(b)=3m, m>=1.                          (5.1)

Use deg_infinity(R)=deg numerator-deg denominator for a nonzero rational function, and deg_infinity(0)=-infinity. Define

    z=a/r; deg_infinity(z)=m.

Divide (3.1) by 9ar^3. The correction numerator has degree at most 12m (the a^2 Ha term), while the denominator has degree 13m. Consequently

    deg_infinity(b-kappa^2 z^3/16) <= -m.                (5.2)

Define a rational comparison triple

    b0 = kappa^2 z^3/16,
    a0 = z(2b0+kappa),
    d0 = z^2(2b0+kappa)(4b0+kappa)/8.

With E=H=L=M=0 this triple satisfies I=K=0 identically. Its degrees are 3m,4m,8m. Equations (3.2) and (5.2) give

    deg_infinity(b-b0) <= -m,
    deg_infinity(a-a0) <= 0,
    deg_infinity(d-d0) <= 4m.                           (5.3)

The last inequality follows by expanding the quadratic expression in b in (3.2): changing b by degree at most -m changes that degree-8m expression by at most 4m; the remaining terms have degrees at most m, 0, and -3m.

This is a rational comparison at infinity, not a polynomial coordinate change or a claim that z is polynomial. Finite poles of z are permitted.

## 6. The scalar Jacobian has unavoidable positive degree

For the comparison triple put

    h0 = 3a0(2b0+kappa)/8,
    f0 = 3a0d0/4 + 3b0^2/8 + 3kappa b0/8 - a0^3/16.

From (5.3), the actual h,f obey

    deg_infinity(h-h0) <= 3m,
    deg_infinity(f-f0) <= 8m.

Rational differentiation lowers degree at infinity by at least one. Using (2.1), expand the difference from J0=h0 d0'-b0 f0'. Each term has degree at most 11m-1:

    (h-h0)d0': 3m+(8m-1),
    h0(d-d0)': 7m+(4m-1),
    (h-h0)(d-d0)': <=7m-1,
    (b-b0)f0': -m+(12m-1),
    b0(f-f0)': 3m+(8m-1),
    (b-b0)(f-f0)': <=7m-1.

Meanwhile direct differentiation of the comparison triple gives

    J0 = (3kappa^4/16384) z^2(kappa z^3+8)
         (10kappa^2 z^6+61kappa z^3+64) z'.               (6.1)

Because deg_infinity(z)=m>0 in characteristic zero, z' has degree m-1 and nonzero leading coefficient. Thus

    deg_infinity(J0)=12m-1,
    leading coefficient = (15/8192) kappa^7 m lc(z)^12 !=0.

The error degree 11m-1 is strictly smaller. Therefore the actual polynomial J has degree exactly 12m-1>0, contradicting J in k*. This proves the claim.

## 7. Equality control and the precise remaining boundary

The comparison triple is itself a useful non-Keller control. Set z=x and kappa=1, with E=H=L=M=0, and form P,Q from section 2. All positive y coefficients of J vanish, beta=1, and J is exactly the degree-eleven polynomial (6.1). Thus the scalar degree obstruction is attained; the first-integral curve is nonempty and rational. Its rationality or I=K=0 cannot replace the missing nonzero-constant-Jacobian condition.

More generally any nonconstant polynomial z of degree m gives a polynomial control with deg J=12m-1. The argument only uses J independent of y until its final contradiction. Thus its cases give the stated congruence law: 12m-1 when a,b are both nonconstant, 3 deg(b)-1 when a is constant and b is not, and J=0 when b is constant. Every degree 3s-1, s>=1, is attained by a=H=E=d=0, b=x^s and any nonzero kappa, with J=-(3/8)b(2b+kappa)b'. The proof rules out constancy, not existence of solutions to the two conserved quantities. It does not classify every polynomial solution of I=L,K=M, establish an automorphism theorem without the displayed shape, or address any of the four other genus-six terminal degree pairs.

## 8. Research board and anchor-language transition

Anchor: the actual polynomial Jacobian, including its scalar coefficient. Niche: conserved-quantity curves and exact degree-at-infinity repair. Wildcard: compare with guarded Collatz affine blocks and the shape atlas.

| Live concept | Object/operation | Preserved data | Lost data and cheapest test |
|---|---|---|---|
| Source presentation | monic (4,6) coefficients | actual y-degrees, beta | a semigroup alone loses the presentation; keep alternate (6,9,4) hostile |
| First integrals | (a,b,d) -> (I,K) | vanishing y²,y rows | scalar J; test the degree-eleven comparison family |
| Rational chart | (a,b) -> z=a/(2b+beta) | exact fraction-field identities | finite poles; do not claim a polynomial coordinate change |
| Degree completion | comparison curve plus error | leading scalar Jacobian | lower jets; verify error strictly below 12m-1 |
| Collatz fuel | fixed 27-block -> affine residue map | a guarded endpoint congruence | magnitude/future legality; retain binary fuel from Session53 |
| Shape charts | marked ray in simplex/sphere/box | invertible normalized coordinates | scale/integrality if marker dropped; retain arithmetic sidecar |

The precise common move is to retain the consumer after a quotient. The planar source map to a first-integral curve loses scalar J; Session53's 27-block map to the C2 affine quotient modulo seven loses height and the dyadic guard. The comparison family above and the inherited two-copy 27-block source are separate hostile controls. There is no mathematical reduction between the open Jacobian and Collatz conjectures here.

The exact Collatz next test remains native transition after the repeated-block fuel expires, or a content-sensitive extension beyond the already closed endpoint families. No new Collatz convergence claim is made by this planar proof.

## 9. Reproduction and scope

From the repository root:

    python3 04-computation/shape_anchors_20261025/verify_quartic_sextic.py
    python3 -O 04-computation/shape_anchors_20261025/verify_quartic_sextic.py

The checker verifies the symbolic identities directly, independently forms the Sylvester resultant, tests rational comparison functions with finite poles, and replays polynomial sharpness controls of several degrees and nonzero beta values. Explicit exceptions are raised instead of assertions. No finite census is used to justify the all-degree exclusion; sections 2--6 provide the proof.
