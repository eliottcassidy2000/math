# The all-finite (5,3) quartic boundary class

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This concerns the complete global square-prefix entry on the fixed DG
surface. Neither placement at infinity is included in this theorem.

## 1. Statement, inherited mechanism, and coefficient universe

Keep `W=(P1_x x P1_z)\{z=x^2}`, `t=1/(z-x^2)`, and
`omega=dx wedge dt`. Its second chart is
`x=1/r`, `t=-r^2-r^4b_D`, so `omega=r^2 dr wedge db_D`.
Let `H in L_2`, `L in L_1`, `F=H^2+L`, and assume the
leading binary octic of `H` has two distinct finite zeros of
multiplicities five and three.

**Theorem.** A rational mate of `F` forces `L`
constant. Consequently no polynomial mate of any degree exists.
Section 9 supplies actual rational mates with constant `L` in
this exact partition for every finite fivefold-point position.

Actual source/surface scalings normalize the separation and
leading scalar. Write `u=x-p`, keeping arbitrary complex `p`.
This ordinary polynomial translation is not a surface
automorphism. The complete original global section constraints
give

    N=u^5(u-1)^3,
    P=u^2[2u^4+(-4p-6)u^3+A u^2+B u+b],
    Q=u^4+(-4p-3)u^3+(A+4p^2-3)u^2
        +(-2Ap+B+12p^2+12p+1)u+c,
    M=u^2(mu4 u^2+mu3 u+mu2),
    R=mu4 u^2+(mu3-2p mu4)u+r0,
    H=N t^2+P t+Q,  L=M t+R.                         (1)

The active-point jets responsible for the factors `u^2` are
proved in Section 2. Before those jets, the constant and linear
coefficients of `P,M` are free. All remaining lower rows in
(1) are induced by the full global filtration; none is a frozen
prefix or selected subfamily. In local boundary coordinates
`s=z-x^2`, use `calN=N+sP+s^2Q`, `calM=M+sR` to distinguish
sections from their leading polynomials.

Dependencies are [the full section filtration](planar_jc48_sep08_dg_quadratic.md),
[inverse leading exactness](planar_jc48_sep08_leading_exactness.md),
[M-unit Newton/Morse analysis](planar_jc48_sep08_quartic_common_root.md),
[finite first-jet obstructions](planar_jc48_sep08_shared_roots.md),
and the compact primitive-space mechanism in
[the all-finite (4,3,1) proof](planar_jc48_sep08_four_three_one.md).
The proof derives the additional inverse coefficients it uses.

Live concepts are the full normal jets, normalized pole orders,
inverse-series trace, zeros of the relative form, and complete
Riemann--Roch spaces. The source-to-target maps retain the
original fibre constant and volume. Local principal parts alone
lose ordinary affine points on the same vertical line. Section 7
records a concrete failed correction and the replacement that
preserves those points. The sharp constant-`L` family prevents
turning the polynomial conclusion into blanket rational exclusion.

## 2. Active jets and the moving-root residue

If `L` is constant, the nonconstant factor `2H` in `J(F,G)`
already excludes polynomial mates. Suppose `L` is nonconstant
and a rational mate exists. The fivefold point must be active,
meaning shared with `calM` and with vanishing normal derivative
`calN_s`. Otherwise its relative forms are regular: an M-unit
has `5!=3j`, and a normal unit is regular. At the triple point
an inactive branch is regular or has a forbidden logarithm;
an active triple has exactly two normalized differential poles
of order two, hence primitive capacity two. Since `N` has
degree eight, a generic point on `D` has differential order two
and primitive local degree three. This exceeds that total
capacity even on the same component. If no active point exists,
the primitive is constant on every compact component instead.

Thus the finite first-jet obstruction requires
`P(0)=P'(0)=M(0)=M'(0)=0`, giving precisely (1). The complete
inverse coefficient `T2=-M/(4N)` makes `M du/N` rational-exact
by normalized trace. Its residue at zero gives

    mu4+3mu3+6mu2=0.                                (2)

Do not yet assert that the triple point is M-unit: a shared
normal-unit triple is permitted at this stage.

Use the actual chart `w=u^2t`, with volume `u^-2 du wedge dw`.
The first two rows of the complete function are

    f(w)=(bw+c)^2+mu2 w+r0,
    g(w)=2(bw+c)(-w^2+B w+Q1)+mu3 w+mu3-2p mu4,
    F=f(w)+u g(w)+O(u^2).                            (3)

For every simple root of the original equation `f(w)=zeta`,
the moving-root residue is `(g'f'-g f'')/(f')^3`. If `b!=0`,
there are two such actual smooth outer branches generically;
the residue identity would force `g=C f'`. But the cubic
coefficient of `g` is `-2b`, while `f'` is linear. Thus

    b=0.                                            (4)

If `mu2!=0`, `f` is now linear and exactness forces `g` to be
constant. Its quadratic and linear coefficients give
`c=0`, `mu3=0`. By (2) this is **case I**:

    M=lambda u^2(1-6u^2),  R=-6lambda u^2+12p lambda u,
    b=c=0,  lambda!=0.                              (5)

If `mu2=0`, (2) gives **case II**:

    M=lambda u^3(1-3u),  R=-3lambda u^2+(1+6p)lambda u,
    b=0,  lambda!=0.                                (6)

The alternative `mu2=mu3=0` would make `M=0`, hence `L`
constant. An additive constant `r0` is absorbed in the generic
fibre value. Now, and only now, the triple point is M-unit:
its values are `-5lambda` in (5) and `-2lambda` in (6).

## 3. Complete local divisors and the nonsquare-leading sidecar

Under (4), the normal order at the fivefold point is at least
three. In case I its local equation has one smooth outer branch
`s~u^2` with differential order `-2`, and three high
determinations `s~u^(8/3)`, forming one actual branch
`u=tau^3`, `s~tau^8`, with differential order `-4`.
The leading high face is `a^2+lambda Z^3`, with simple roots;
the generic outer face is linear after removing its zero factor.
One plus three exhausts the local Weierstrass degree four.
The primitive divisor therefore has degrees one and three,
total four.

If the triple point has `P(1)!=0`, its normal order is zero.
The M-unit cancellation pair forms one branch on which the
differential vanishes to order four, forcing primitive local
degree five. This is incompatible with total capacity four.
Thus in case I

    P(1)=0,  B=4p+4-A.                              (7)

In case II the complete fivefold face is obtained with
`u=tau^2`, `s=tau^5 Z`. It has four simple nonzero
determinations, forming two normalized branches, each with
differential order `-4`. Their primitive divisor has degrees
three and three, total six. This capacity does not by itself
force `P(1)=0`.

At an M-unit triple with `P(1)=0`, normal order at least two
gives unit differentials. Normal order one has the balanced
cubic with simple roots or one double root. The double-root
Morse contact is at most two generically: contact one gives
a unit on the ramified branch, whereas contact two gives
forbidden nonzero logarithms. Under the mate all such triple
forms are therefore units. If `P(1)!=0`, there is instead the
single zero of order four described above.

All other differential zeros are the four generic points on
`D`, each of order two. There is no further boundary point
at infinity because `N8=1`. On the original smooth generic
fibre in `W\D`, the form is a unit. These inventories are
complete and refer to the original volume and original fibre.

The geometric generic fibre is integral. A nontrivial polynomial
composition of `F` must have outer degree two or four because
`deg_t F=4`. Outer degree four would make `N` square up to
scalar, contradicted by multiplicities five and three. Outer
degree two would give `(K-H)(K+H)=L-constant`; the second
factor has `t` degree two, forcing `L` constant. Thus `F` is
noncomposite. The classical noncomposite/closed-polynomial
criterion and relative algebraic closure give geometric generic
integrality in characteristic zero; see
[Arzhantsev--Petravchuk, Theorem 1 and Lemma 3](https://arxiv.org/pdf/math/0608157v2)
and [THM-3827, generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases](../../01-canon/theorems/THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases.md).

The resulting genera and complete primitive-space dimensions are

| Case | Fivefold poles of eta | Triple | Genus | Primitive divisor E | dim L(E) |
|---|---|---|---:|---|---:|
| I | 2,4 | units after (7) | 2 | 1,3 | 3 |
| II, P(1)=0 | 4,4 | units | 1 | 3,3 | 6 |
| II, P(1)!=0 | 4,4 | zero of order 4 | 3 | 3,3 | 4 |

Each dimension follows from Riemann--Roch with `deg E>2g-2`.
In particular the genera cannot be interchanged or inferred
from a pole table that omitted the normal-unit triple zero.

## 4. Case I: a complete three-dimensional primitive space

Put

    U=1/u,
    J=u^2(u-1)^2t+u^2-(2+2p)u.                       (8)

The global `L_1` function `J` is regular at `D`. At the
active outer branch `U` has a simple pole and `J` is regular;
at the high branch they have pole orders three and two.
At the triple, the factor `(u-1)^2` cancels the pole of `t`.
The only affine denominator is `u`, whose line has constant
`F` and is absent from the generic fibre. Thus `1,U,J`
belong to `L(E)`. Independence follows by their `t` degrees
and then independence of `1,1/u`; they form a complete basis.

Reduce the actual brackets `u^2(F_u G_w-F_w G_u)` modulo
the original `F-zeta`, with (5),(7). For a primitive
`a(zeta)U+b(zeta)J+constant`, the `w^3` coefficient is

    4a u^2(u-1)^6.

Hence `a=0`. The remaining `w` coefficient is
`-2b lambda u(u-1)^2`, so `b=0`, contradicting bracket one.

## 5. Two further exact inverse coefficients in case II

For completeness derive the needed coefficients, rather than
assume every inverse row is exact. For `F(u,T)=v^-4`, a
rational mate expanded in the actual radical field satisfies

    partial_u G(u,T)=(v^5/4)partial_v T,
    g_(k+4)'=(k/4)T_k.

Thus every `T_k du` with `k!=0` is exact in that field; when
`T_k` is rational, normalized trace descends exactness to
`C(u)`. Center the quadratic by `t=y-P/(2N)` and put

    D0=Q-P^2/(4N),  E0=R-MP/(2N).

The centered quartic is `(Ny^2+D0)^2+My+E0`.
Factor it formally into two monic quadratics with linear
coefficients `sigma,-sigma`. For the pair whose roots begin
with `+N^(-1/2)v^(-1)` and `-N^(-1/2)v^(-1)`, write
`sigma=M v^2 V(v^4)/(2N)`. The coefficient identities give

    (1-E0 z)V^2 +D0 M^2 z^2 V^4/(4N)
        +M^4 z^3 V^6/(64N^2)=1,
    V=1+(E0/2)z+(3E0^2-D0 M^2/N)z^2/8+... .

The two inverse branches are `y(v),y(-v)` and their sum is
`-sigma`. Consequently

    T6=-ME0/(8N),
    T10=-M(3E0^2-D0 M^2/N)/(32N).                    (9)

These rational coefficients retain all complete lower rows.
In case II, their residues at zero are, successively,

    Res T6 du=lambda^2(B+12p+2)/16,
    Res T10 du=lambda^3 c/32
         after B=-12p-2.                            (10)

Therefore

    B=-12p-2,  c=0,
    d=P(1)=A-16p-6.                                 (11)

No low coefficient of the original `Q` was independently
altered. The remaining cases are distinguished by `d`.

## 6. Case II with d=0: all six primitive coefficients are controlled

Use `U,J` from (8) and put `V=(u-1)H/u`. A full basis is

    1, U, J, UJ, V, VJ.                              (12)

At each active branch the pole bounds are respectively
`0,2,1,3,2,3`, since `H` is finite there. At the triple,
`J` is regular and the factor `u-1` cancels the pole of `H`.
They are regular at `D`, and have no other affine denominator
than `u`, which is absent generically. Independence follows
from `t` degrees three for `VJ`, two for `V`, then the distinct
coefficients of `J,UJ`, and finally `1,U`. Hence (12) is the
complete space from the genus-one row of the table.

For the five nonconstant basis functions in the order
`U,J,UJ,V,VJ`, form their actual brackets and reduce modulo
`F-zeta`. The five homogeneous coefficient rows at

    u^8 w^3, u^6 w^2, u^5 w^2, u^4 w, u^3 w

have determinant

    192 lambda zeta(-lambda^2+32p^2 zeta).            (13)

It is nonzero in `C(zeta)` for every `lambda!=0` and every
`p`, including zero. All five primitive coefficients vanish,
so the target bracket one is impossible. This is a determinant
of the complete primitive space, not a degree cutoff for mates.

## 7. Case II with d!=0: an affine-safe genus-three correction

A tempting fourth function is
`(J+d/(u-1))/u`. It cancels the boundary triple's rational
principal part but creates poles at the two ordinary affine
points over `u=1` on the generic fibre. It is therefore **not**
in the primitive space. The first failed implication was that
repairing a boundary principal part paid the complete vertical
line. This rejected function is not used in any matrix below.

Instead put

    J0=u(u-1)^3t+u^2-(3+2p)u,
    K(u)=-u^2+(2p+3)u-A+(12p+2)/u,
    k0=K(1)=-A+14p+4,  k1=K'(1)=-10p-1,
    Adj=J0-k0+k1(U-1),
    W0=H Adj.                                       (14)

These functions have only the denominator `u`. Near the
normal-unit triple, the actual pole branch satisfies

    t=-P/N-(H-Q)/P+O(N(H-Q)^2/P^3),
    H=(unit)(u-1)^(-3/2).

Therefore `J0=K(u)+O((u-1)^(3/2))`. The explicit subtraction
in (14) kills both the constant and linear analytic jets of
`K`. Thus `Adj` has order at least `3/2`, and `W0` is regular
at that boundary branch. Unlike the rejected correction, it is
also regular at every ordinary affine point over `u=1`.

At each active point `J0` has pole order three, `U` has order
two, and `H` is finite, so `W0` has at most order three.
All functions in (14) are regular at `D`, since `H,J0` are
global and `U` tends to zero. No other affine poles are
introduced. Hence

    1,U,J0,W0                                       (15)

is a basis of the four-dimensional space. Independence follows
from the unique `t` degree three of `W0`, then degree one of
`J0`, and the remaining `1,U`.

For the three nonconstant functions in (15), take the homogeneous
bracket rows `u^8w^3,u^6w^2,u^5w^2`. Their determinant is

    16p{lambda[A^2-(56p+12)A+640p^2+336p+36]
              -12zeta(A-8p-6)}.                    (16)

For `p!=0`, this is nonzero: a zero fibre slope forces
`A=8p+6`, leaving `256lambda p^2` inside braces. For `p=0`,
use rows `u^8w^3,u^6w^2,u^4w`; their determinant is

    -20lambda[(6-A)lambda+12zeta],                  (17)

also nonzero. Thus all three coefficients vanish in either
case. This finishes the genus-three branch without introducing
any new affine pole, and completes the rational conclusion.

## 8. Exact controls and audit boundary

The producer verifies complete shifted global rows, the moving
residue reductions, all inverse coefficient identities, both
residues in (10), the local leading faces and normalized
orders, the corrected genus-three analytic jets, and the full
matrices (13),(16),(17). The proof supplies branch exhaustion,
geometric integrality and completeness of the primitive spaces.
Every matrix row is extracted directly from a polynomial
remainder; no denominator clearing changes the displayed
monomial addresses. The source checks this polynomial property
separately for every basis function.

Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_five_three.py
python3 -B -O 04-computation/planar_jc48_sep08_five_three.py
```

Both runs pass **117 gates** and reproduce the same 437 bytes
as [the frozen output](planar_jc48_sep08_five_three.out).
The [source](../../04-computation/planar_jc48_sep08_five_three.py)
has 13,167 bytes, SHA256
`0ceac3b9198d66241883179a4fc0b7fc598f2e838af87561f42d1d966951dba4`.
Output SHA256:
`3013be961d8e1c884adb298f858ef01e4d721d8b9521e08ffa4faf2d92e40c9b`.
Semantic gate digest:
`a1fc2ea5ca623ebc16b2c697ee095b8c019938f119491a01e669a98e99a95bcf`.

The [independent complete audit](planar_jc48_sep08_five_three_audit.md)
accepts every branch, genus and full primitive space. It reconstructs
both inverse residues and every determinant from original-coordinate
Jacobians and separate polynomial long division, without importing
the producer. All117 normal and optimized gates match the frozen output.

## 9. Same-partition rational sharpness and remaining placements

For every finite `p`, put

    A0=u(u-1),  Z0=u^2(u-1)t+u-2p-1,
    H=A0 Z0^2+c,
    G_H=[1/(3u^2)+5/(3u)+1/(u-1)]/Z0.

The original global section constraints hold for every `p`,
the leading polynomial is precisely `u^5(u-1)^3`, and
`J(H,G_H)=1`. Hence `H^2+c0` has rational mate `G_H/(2H)`.
This is an actual all-parameter hostile to removing the
constant-`L` boundary; it is not a polynomial mate.

If the fivefold or triple point is infinity, the leading
polynomial is a pure cubic or pure quintic and passes the
leading differential gate. Those placements remain separate
OPEN tasks here. This note makes only the all-finite claim.
