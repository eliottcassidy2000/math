---
id: THM-2347
title: "Degree-eighteen double-zero wall saturation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On the
  degree-eighteen wall 20BC+21W=0, every point whose structured Mordell
  polynomial has square-class degree at most four either lies in an
  already-closed two-sparse plane or satisfies 126D=25B^2. In the
  B!=0, C!=0 chart, three coefficients of the terminal subresultants,
  after exact factor stripping and the substitution X=C^2, have a
  grevlex basis containing (126D-25)^2. THM-2345 closes the latter
  common-root wall, so the complete double-zero wall contains no Keller
  trajectory. Other parameter walls, the residual off-wall H_2/H_4
  strata, and JC(2) remain open.
source: codex-2026-07-25-double-zero-wall-saturation
depends_on:
  - THM-2314-degree-eighteen-bd-linear-ratio-closure
  - THM-2316-degree-eighteen-cd-linear-ratio-closure
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2345-degree-eighteen-common-root-wall-saturation
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
  - THM-2342-degree-eighteen-deep-wall-first-flux-cover-obstruction
script: 04-computation/jc2_degree18_double_zero_wall_saturation_thm2347.py
output: 05-knowledge/results/jc2_degree18_double_zero_wall_saturation_thm2347.out
script_sha256: 597041ac116c1ec9d0cca4d2c5d3815320922dc10c9ce8f48f8726c0fb009f4a
output_sha256: 7937687edb626dff625dd360bdab937068e105fc13cdcdc3c0a7807981cde4ae
hash_basis: working-tree bytes (LF)
---

# THM-2347 -- the double-zero wall saturates to the common-root wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2332 reduces every remaining degree-eighteen Keller trajectory to

```text
F=4P^3+49Q^2=H S^2,                  deg(H) in {0,2,4}.       (1)
```

THM-2345 closed the codimension-one wall on which `P(0)=0`.  This theorem
closes its natural dual: the wall on which `Q` gains a second zero at the
origin,

```text
20BC+21W=0.                                             (2)
```

The mechanism is symmetric only at first sight.  Condition (2) does not
give `F` a visible square factor.  Instead, the square-class multiplicity
forces three terminal subresultant coefficients to vanish, and their
stripped ideal remembers the missing common-root coordinate.

## 1. The double-zero wall

Retain THM-2332's exact covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                            (3)
```

The last coefficient in (3) is

```text
2361960BC+2480058W=118098(20BC+21W).                 (4)
```

Thus (2) is exactly the wall

```text
Q=y^2 R_4(y).                                        (5)
```

It is weighted homogeneous for

```text
wt(B,C,D,W)=(2,3,4,5).                               (6)
```

The cases omitted by the main normalization are already closed.  If
`B=0`, equation (2) gives `W=0`, leaving only the `C`--`D` plane; its
axes are closed by THM-2297 and its complete singular ratio bank by
THM-2316.  If `B!=0` but `C=0`, equation (2) again gives `W=0`, and the
corresponding `B`--`D` branch is closed by THM-2314.

It remains to consider

```text
B C !=0.                                             (7)
```

Over `C`, weighted scaling normalizes `B=1`.  Equation (2) becomes

```text
W=-20C/21.                                           (8)
```

In this chart,

```text
P=245y^4+1890y^2-24300+122472D,

Q=y^2(539y^4+11340y^2+183708Cy+72900-367416D).       (9)
```

The degree and leading coefficient of `F` remain

```text
deg_y(F)=12,                 lc_y(F)=73060029!=0.    (10)
```

## 2. Low square class forces four common derivative roots

Equation (1) and (10) give

```text
deg(S)=(12-deg(H))/2 >=4.                            (11)
```

Differentiating (1) shows

```text
S divides gcd(F,F'),             deg gcd(F,F')>=4.   (12)
```

This implication does not assume a generic root pattern.  It remains valid
if `H` or `S` meets a higher-multiplicity collision.

Because the leading coefficients of `F` and `F'` are nonzero constants,
the principal-subresultant specialization theorem applies without a
degree-drop exception.  The complete subresultant degree profile over
`Q(C,D)` is

```text
12,11,10,9,8,7,6,5,4,3,2,1,0.                      (13)
```

Write the degree-three and degree-two members as

```text
Sres_3=a_3y^3+a_2y^2+a_1y+a_0,

Sres_2=b_2y^2+b_1y+b_0.                              (14)
```

Condition (12) makes both polynomials in (14) vanish identically.  In
particular,

```text
a_3=a_0=b_2=0.                                      (15)
```

## 3. Three stripped coefficients recover the missing wall

Put

```text
J=126D-25.                                          (16)
```

Exact division gives, up to nonzero rational contents,

```text
a_3=J^3 f_0(C^2,J),

a_0=C J^5 f_1(C^2,J),

b_2=J^6 f_2(C^2,J).                                 (17)
```

The powers in (17) are exact: after the displayed divisions the three
quotients are not divisible by `J`, and the middle quotient is not
divisible by another `C`.  The parity in `C` follows from the involution

```text
(y,C) -> (-y,-C).                                   (18)
```

Set

```text
X=C^2.                                              (19)
```

After substituting `D=(J+25)/126`, clearing only nonzero rational contents,
and using grevlex order in `(X,J)`, the exact Groebner basis of

```text
(f_0,f_1,f_2)                                       (20)
```

is

```text
[
 -(30618X+361)^2 L(X,J),

 J(30618X+361)^3
   (199644669X^2+7654500X+62500),

 J^2
],                                                  (21)
```

where

```text
L(X,J)
 =1233383522515969671 JX^2
  +27991525644715500 JX
  +158573241187500 J
  -659000910244935988096920 X^5
  -54559736202362496221340 X^4
  -1758935316808435500000 X^3
  -27710180945527500000 X^2
  -213905002500000000 X
  -648671875000000.                                 (22)
```

The last basis element is the decisive certificate:

```text
J^2 belongs to (f_0,f_1,f_2).                       (23)
```

Suppose first that `J!=0`.  Since `C!=0`, equations (15)--(17) force all
three generators in (20) to vanish.  Equation (23) then forces `J^2=0`, a
contradiction.  Therefore every low-square-class point in the remaining
chart satisfies

```text
J=0,

126D=25                         when B=1.            (24)
```

Undoing the weighted normalization turns (24) into the invariant equation

```text
126D=25B^2.                                         (25)
```

## 4. Close the complete wall

Let a Keller trajectory satisfy (2).  The `B=0` and `C=0` cases were
closed in Section 1.  In the remaining chart THM-2332 gives (1), hence
(11)--(15), and Section 3 gives (25).  The trajectory therefore lies on
the complete common-root wall closed by THM-2345.  This is again
impossible.

Consequently

```text
no degree-eighteen Keller trajectory satisfies 20BC+21W=0.            (26)
```

This closes a second codimension-one weighted wall in the higher-support
degree-eighteen cone.  It does not close the residual off-wall `H_2` or
`H_4` strata, any other parameter wall, `JC(2)`, or `DC(2)`.

## 5. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_double_zero_wall_saturation_thm2347.py
python3 -O 04-computation/jc2_degree18_double_zero_wall_saturation_thm2347.py
```

Both transcripts are byte-identical to the stored output.  The companion
checks (4)--(5), (9)--(10), the complete profile (13), every exact factor
and maximal factor order in (17), the `C`-parity reduction, the three
generator signatures, every element of (21)--(22), and the zero remainder
in (23).  A generic deep point has gcd degree five, while two off-deep
hostile controls have gcd degree zero.  No executable check uses Python
`assert`.

## 6. Independent hostile audit

An independent pass checked the weighted `B=1` normalization and the
`B=0`/`C=0` handoff to the proved two-sparse closures. It rederived that
`F=HS^2` with `deg H<=4` forces `deg gcd(F,F')>=4`, including higher
root collisions, and verified that the constant nonzero leading
coefficients of `F,F'` remove every subresultant specialization
degree-drop exception. It then checked that the three selected
subresultant coefficients are used only after the valid `C*J!=0` factor
stripping, so `J^2` in the stripped ideal forces the contradiction exactly
in the stated chart.

The full exact referee was replayed independently under ordinary and
optimized Python. Both several-minute computations are byte-identical to
the stored transcript after LF normalization and reproduce the complete
subresultant profile, exact factor orders, three-element grevlex basis,
terminal `J^2`, and deep/off-deep gcd controls. The frontmatter LF hashes
were also reproduced.

QED.
