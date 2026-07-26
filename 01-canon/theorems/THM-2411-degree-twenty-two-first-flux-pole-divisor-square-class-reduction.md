---
id: THM-2411
title: "Degree-twenty-two first-flux pole-divisor emptiness and projected square-class atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the polynomial exact-square-prefix, genuine nonsplit degree-twenty-two
  terminal branch, legal target translation leaves five weighted parameters
  B,C,D,E,W. On the divisor where the first flux loses its Z coefficient,
  the residual first flux forces the fixed constant coefficient E to equal
  a degree-five polynomial in the centered coordinate y. Hence y is constant;
  the second flux then makes T and q constant, contradicting the genuine
  deck. Thus the entire first-flux pole divisor is empty. The exact quadratic
  F_2, sextic R_6, perfect-square classification, and H_2 S_2^2 example remain
  valid only as projected coefficient-cone controls after forgetting the
  fixed-E sidecar. This closes only the A=0 chart, not degree twenty two,
  JC(2), or DC(2).
source: codex-2026-07-26-degree-twenty-two-pole-divisor
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2217-square-prefix-pole-alternative-and-odd-leading-degree-terminal-wall
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2406-degree-eighteen-H4-weighted-pole-deep-wall-collapse
script: 04-computation/jc2_degree22_first_flux_pole_divisor_thm2411.py
output: 05-knowledge/results/jc2_degree22_first_flux_pole_divisor_thm2411.out
script_sha256: 7d14a16aab791db2da9dc2749117db6bdfe539096fe8d7670a9c852d44d956e3
output_sha256: 916848b045d78ada4ecd6e911014668dc511c8d9c0b019bf82ca8ca2f2467000
hash_basis: working-tree bytes (LF)
---

# THM-2411 -- the degree-twenty-two first-flux pole divisor is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The degree-eighteen program uses a trigonal spectral curve because the first
constant Faber flux eliminates `Z=T^2` generically. At degree twenty two, one
important exceptional divisor behaves differently: the coefficient of `Z`
in the first flux vanishes. The first version of this theorem projected the
remaining equations to a quadratic cover of the centered coefficient line.
That projection forgot a decisive sidecar: `E` is a fixed coefficient in the
constant field, not a function allowed to vary with `y`.

Retaining `E` closes the divisor immediately. The quadratic, sextic, and
square-class calculations remain exact and useful as an atlas of the
projected universal coefficient family, but they do not describe live Keller
trajectory branches.

## 1. Inherited nonsplit terminal coordinates

Work over

```text
R=C[x],                         K=C(x),
```

where `C` is algebraically closed of characteristic zero. In the polynomial
exact-square-prefix and genuine nonsplit terminal branch of THM-2214 and
THM-2217, depress the quartic over the usual quadratic function-field
extension:

```text
P=w^4+2d w^2+q w+(d^2-s),       q^2=T!=0.             (1)
```

The deck involution fixes `d,s,T` and sends

```text
(w,q) -> (-w,-q).                                     (2)
```

It is genuine, so the constant field of the extension is still `C` and
`q` is not a deck-fixed constant.

After the exact target-shear quotient of THM-2230 and the full Faber gauge,
write a reduced degree-twenty-two mate as

```text
Q=E_22+alpha E_18+beta E_14+gamma E_10
      +delta E_6+epsilon E_2.                         (3)
```

Here `E_j` is the polynomial part at fibre infinity of `P^(j/4)`, and all
five displayed coefficients are in `C`. THM-2129's Hamiltonian identity
and nonsplit parity give the two constant flux equations

```text
Phi_Q=0,                     Psi_Q=Psi0 in C.          (4)
```

The third flux and the Keller one-form remain available. If

```text
F=R_Q/q,       q=A_src/U,       U^2=V,
```

then, in the notation of the exact-square prefix,

```text
A_src(2T F'+F T')=2kappa T,       kappa in C*.         (5)
```

Nothing below spends (5). The first two fluxes already close `mathcal A=0`;
the one-form remains available for the complementary chart.

## 2. The degree-twenty-two target-translation normal form

For the legal target translation

```text
P_c=P+c,                       c=2alpha/11,            (6)
```

the finite Faber translation law is

```text
E_(4j-2)(P_c)
 =sum_(i=0)^(j-1) binom(j-1/2,i)c^i E_(4(j-i)-2)(P).
                                                            (7)
```

Rewriting the fixed polynomial `Q` in the basis attached to `P_c` kills
the `E_18` coefficient. Define

```text
B =beta-9alpha^2/22,

C =gamma-7alpha beta/11+21alpha^3/121,

D =delta-5alpha gamma/11+35alpha^2 beta/242
        -315alpha^4/10648,

E =epsilon-3alpha delta/11+1815alpha^2 gamma/29282
        -385alpha^3 beta/29282+63alpha^5/29282,

W =Psi0+(105alpha^6-770alpha^4 beta+4840alpha^3 gamma
        -31944alpha^2 delta+234256alpha epsilon)/644204.
                                                            (8)
```

Then the exact whole-polynomial identity is

```text
Q=E_22(P_c)+B E_14(P_c)+C E_10(P_c)+D E_6(P_c)
       +E E_2(P_c),                                     (9)
```

and the second flux in these coordinates has constant value `W`. The first
and third fluxes are unchanged. Since `P_c` replaces `s` by `s-c`, put

```text
y=11s-2alpha,          u=dT,          Z=T^2.           (10)
```

From now on `s` means the translated coordinate, so equivalently
`alpha=0` and `y=11s`.

The natural weights are

```text
wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).               (11)
```

Every equation below is homogeneous for these weights.

## 3. The first two normalized fluxes

Define

```text
mathcal A
 =616B-1089u+63y^2,                                  (12)

mathcal K
 =-745360Buy+6160By^3+2342560Cu-58080Cy^2
   +511104Dy-3748096E
   +922383u^2y-25410uy^3+63y^5.                      (13)
```

The first flux is

```text
Phi_Q/q=-N_1/7496192,

N_1=1331 mathcal A Z+4 mathcal K.                    (14)
```

The second flux is

```text
Psi_Q-W=N_2/1319329792,                              (15)
```

where

```text
N_2
 =15944049Z^2
  +65591680BZy-206145280CZ-162339408Zuy+2236080Zy^3
  +1443016960Bu^2-71554560Buy^2+98560By^4
  +449771520Cuy-1239040Cy^3
  -1978994688Du+16355328Dy^2
  -239878144Ey-1319329792W
  -1190488992u^3+147581280u^2y^2
  -1219680uy^4+672y^6.                               (16)
```

The exact companion reconstructs (14)--(16) twice: once by the quartic
Laurent recurrence and once by a direct finite multinomial expansion. It
checks all six Faber rows of degrees `2,6,10,14,18,22`.

## 4. The first-flux pole divisor

Consider the divisor

```text
mathcal A=0.                                         (17)
```

Equation (17) is precisely the locus where the first flux loses its `Z`
coefficient. It gives

```text
u=(616B+63y^2)/1089.                                 (18)
```

The remaining part of the first equation is `mathcal K=0`, which gives

```text
E=[-71148B^2y+745360BC+5082By^3+43560Cy^2
    +287496Dy+945y^5]/2108304.                       (19)
```

Solving (17) and the residual first flux for `u` and `E` is algebraically
reversible in the universal coefficient variety. Substitution in (16) gives

```text
N_2=F_2/243,                                         (20)
```

where

```text
F_2(Z,y)
 =3874403907Z^2
  -6375511296BZy-50093303040CZ-1738775808Zy^3
  +59838693376B^3+10491199488B^2y^2+41215426560BCy
  -272021815296BD+587575296By^4+4817387520Cy^3
  -31794757632Dy^2-320597139456W+20901888y^6.
                                                            (21)
```

Thus every degree-twenty-two Keller trajectory on (17) lies on the explicit
quadratic curve

```text
F_2(Z,y)=0.                                          (22)
```

No root choice or squaring has entered this identity. Projecting away the
solved value of `E`, however, loses its fixed-constant condition; Section 5
restores exactly that sidecar.

## 5. The fixed-`E` sidecar empties the divisor

The normalized coefficients `B,C,D,E,W` all lie in the fixed constant field
`C`. This fact must be retained when using (19). Define

```text
P_5(Y)
 =945Y^5+5082B Y^3+43560C Y^2
   +(-71148B^2+287496D)Y+745360BC.                    (22a)
```

Direct substitution of (18) in the residual first flux gives the scalar
associate

```text
mathcal K
 =(16/9)[P_5(y)-2108304E].                            (22b)
```

Thus `mathcal K=0` makes `y` algebraic over `C`. Since `C` is algebraically
closed and `y in C(x)`, it follows that `y in C`; the fixed leading
coefficient `945!=0` prevents (22b) from becoming the zero polynomial.

Equation (18) now makes `u` constant. Equation (21), whose leading
`Z` coefficient is the nonzero constant `3874403907`, makes `Z` algebraic
over `C`, hence constant. Therefore `T^2=Z`, then `T`, and finally `q` are
constants in the algebraic function field. The genuine deck fixes constants
but sends `q` to `-q`, contradicting (2) and `q!=0`.

Consequently:

> **First-flux pole-divisor closure.** There is no genuine nonsplit
> degree-twenty-two Keller trajectory on `mathcal A=0`.

```text
{genuine degree-22 trajectories on mathcal A=0}=empty. (23)
```

The first failed inference in the earlier square-class route was the
projection

```text
(B,C,D,E,W;Z,y) -> (B,C,D,W;Z,y).
```

It preserves the formal equation `F_2=0` but destroys the condition that `E`
is one fixed constant for the whole trajectory.

## 6. The projected sextic square-class atlas

The remaining calculations in Sections 6--8 concern that projected universal
coefficient family. They are exact algebraic controls, not surviving Keller
branches.

The quadratic discriminant is

```text
Disc_Z(F_2)=63478233612288 R_6(y),                   (24)
```

with

```text
R_6
 =42525y^6+205821By^4+1568160Cy^3
  +(-1920996B^2+7762392D)y^2
  -14609056B^3+66411576BD+39530700C^2+78270786W.
                                                            (25)
```

In particular

```text
deg_y R_6=6,                  lc_y R_6=42525!=0.      (26)
```

Completing the square in `Z` identifies the normalization of (22), up to
removing square factors, with

```text
v^2=R_6(y).                                          (27)
```

### 6.1 What projected rationality would force

If one deliberately forgets the fixed-`E` sidecar and supplies a nonconstant
rational parametrization `(Z,y)` of an irreducible component of (22), it
extends by properness to a nonconstant morphism from `P^1`. Riemann--Hurwitz
then forces that projected component to have genus zero.

If `R_6` is a square, (22) splits and this is the degree-zero squarefree
case. Otherwise (27) is a connected double cover. Write

```text
R_6=H S^2,                 H squarefree.              (28)
```

Because `deg R_6=6`, `deg H` is even. The smooth projective model of
`v^2=H(y)` has genus

```text
g=(deg H-2)/2.                                      (29)
```

Genus zero and (26) therefore force exactly

```text
deg H=0                  or                  deg H=2. (30)
```

Equivalently, the only projected square classes that can admit such a
rational parametrization are

```text
R_6=H_0 S_3^2,             H_0 in C*,

or

R_6=H_2 S_2^2,             H_2 squarefree quadratic.
                                                            (31)
```

This statement is retained only as an atlas of the quotient that forgot `E`;
Section 5 proves that neither class lifts to a genuine trajectory on
`mathcal A=0`.

## 7. Exact classification of the projected perfect-square locus

The first alternative in (31) can be written without an existential cubic:

> **Perfect-square criterion.** The sextic (25) is a square in `C[y]` if and
> only if
>
> ```text
> BC=0,
>
> D=22141B^2/79200,
>
> W=-(2080981B^3+13186800C^2)/41164200.              (32)
> ```
>
> Under (32),
>
> ```text
> R_6
>  =42525[y^3+(121/50)By+(1936/105)C]^2.             (33)
> ```

To prove necessity, normalize the leading cubic in a putative square to be
monic. The missing `y^5` term in (25) kills its quadratic coefficient.
The `y^4` and `y^3` coefficients then force the two lower coefficients in
(33). The `y` coefficient forces `BC=0`; the `y^2` coefficient forces the
displayed value of `D`; and the constant coefficient forces the displayed
value of `W`. Substitution proves sufficiency. All pivots are nonzero in
characteristic zero.

The axes are included, not divided away:

```text
C=0:  42525[y^3+(121/50)By]^2,

B=0:  42525[y^3+(1936/105)C]^2,

B=C=0: 42525y^6.
```

Criterion (32) classifies the projected perfect-square coefficient locus.
None of these points lifts to a `mathcal A=0` Keller trajectory because the
fixed-`E` obstruction of Section 5 applies first.

## 8. Projected quadratic and squarefree controls

The second alternative in (31) is not a formal artifact. For example,

```text
B=C=0,              D=-175/10648,       W=-175/161051
                                                            (34)
```

gives

```text
R_6=42525(y^2-2)(y^2+1)^2.                            (35)
```

Thus the projected `(B,C,D,W)` coefficient cone actually meets the
`H_2S_2^2` locus. It does not supply a constant `E` or a Keller trajectory.
Conversely, at

```text
(B,C,D,W)=(1,1,1,1),
```

the exact Euclidean algorithm gives

```text
gcd(R_6,R_6')=1,                                    (36)
```

so the projected universal family also contains generic genus-two curves.
Equations (34)--(36) are coefficient-cone controls only.

## 9. The retained whole-polynomial sidecar

The flux quotient forgets that the reduced mate is a polynomial in the
original exact-square-prefix coordinates. Preserve that information as
follows. Write

```text
P=H^2+L
```

and, for `j>=1`, define the finite binomial truncation

```text
mathcal R_j(P,H)
 =sum_(i=0)^(j-1) binom(j-1/2,i)
      H^(2j-1-2i)L^i.                                (37)
```

Put `S_j=E_(4j-2)-mathcal R_j`. Independent Faber expansion gives the top
sidecar

```text
E_22-mathcal R_6(P,H)
 =33T/2048 [
      14L^4-28L^3s-14L^2Td+42L^2s^2
      -LT^2+56LTds-56Ls^3
      +14T^2d^2+6T^2s-140Tds^2+70s^4
   ].                                                (38)
```

and the lower sidecars

```text
S_4=35T(L^2-2Ls-Td+3s^2)/128,

S_3=-5T(-L+2s)/16,

S_2=3T/8,                         S_1=0.              (39)
```

The individual `E_(4j-2)` belong a priori to `K[z]`; polynomiality does not
apply to `E_22` alone. What is known to lie in `R[z]` is the full-mate
difference

```text
Q-[mathcal R_6+B mathcal R_4+C mathcal R_3
       +D mathcal R_2+E mathcal R_1]

 =S_6+B S_4+C S_3+D S_2
 in R[z].                                             (40)
```

Here `Q` and every `mathcal R_j(P,H)` are polynomials. Equation (40), not
the top-seed identity (38) by itself, is the exact whole-polynomial
pole-cancellation sidecar for a later parametrization.

## 10. Exact verification and scope

Run

```text
python 04-computation/jc2_degree22_first_flux_pole_divisor_thm2411.py
python -O 04-computation/jc2_degree22_first_flux_pole_divisor_thm2411.py
```

The companion:

1. reconstructs the Laurent rows of degrees `2,6,10,14,18,22` both by
   recurrence and direct multinomial expansion;
2. checks the whole-polynomial translation identity and all five invariants
   in (8);
3. verifies (14)--(25), including the exact scalar in (24), and checks that
   the fixed-`E` obstruction (22b) has degree five with leading coefficient
   `945`;
4. checks both directions of the projected criterion (32)--(33), including
   the `B`-axis, `C`-axis, and their intersection;
5. checks the projected squarefree and `H_2S_2^2` controls (34)--(36); and
6. derives (38)--(40), including every lower Faber sidecar, independently of
   the flux elimination.

The theorem is scoped to the divisor `mathcal A=0` inside the genuine
nonsplit, polynomial exact-square-prefix degree-twenty-two terminal branch.
It proves that divisor empty. It does not treat `mathcal A!=0`, does not turn
the projected controls (31) into trajectories, and does not prove `JC(2)` or
`DC(2)`.

## 11. Independent hostile audit

An independent audit reconstructed the degree-twenty-two Laurent bank from
the recurrence and checked it against the separate finite multinomial
formula. It rechecked the legal target translation, all five normalized
constant coefficients, both flux numerators, and the exact associate

```text
mathcal K|_(mathcal A=0)
 =(16/9)[P_5(y)-2108304E].
```

The audit attacked the repaired step directly. Because `E` is one fixed
element of the algebraically closed constant field and the leading
coefficient of `P_5` is the nonzero constant `945`, the equation cannot
become an identity in a nonconstant `y`. Once `y` is constant, (18) makes
`u` constant, and the nonzero `Z^2` coefficient in (21) makes `Z`, then
`T` and `q`, constant. This contradicts the genuine deck. No use of the
projected square-class atlas is needed for the emptiness conclusion.

The auditor also checked the quotient boundary: (24)--(36) remain exact
after forgetting `E`, but the displayed `H_2S_2^2` control does not lift to
a trajectory and is never used as one. The full-mate polynomial sidecar,
including all lower Faber rows, was rederived; no top-seed polynomiality
claim survives by itself.

Normal, optimized, and stored transcripts byte-match after LF normalization,
and the declared hashes agree with the working-tree bytes. No mathematical,
typing, dependency, or scope defect remains. **QED.**
