---
id: THM-2314
title: "Degree-eighteen B--D linear-ratio closure"
status: >
  PROVED + VERIFIED-EXACT. In the genuine nonsplit polynomial
  exact-square-prefix degree-eighteen branch of THM-2262/2297, the two
  rational linear-factor points D/B^2=4075/85176 and D/B^2=25/126 in
  THM-2311's B--D bank are empty. After the weighted normalization B=1,
  both trigonal spectral curves are absolutely irreducible. At the first
  ratio there are eight simple branch points and two smooth totally
  ramified cubic fibres, giving normalization genus four. At the second,
  y=0 is an ordinary triple point with three unramified normalization
  branches, while six other simple branch points give normalization genus
  one. A rational Keller trajectory is therefore constant and yields the
  inherited nonsplit-deck contradiction. Four algebraic B--D ratios and 29
  two-sparse ratios overall remain; this does not prove JC(2).
source: codex-2026-07-25-degree18-bd-linear-ratios
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
script: 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
output: 05-knowledge/results/jc2_degree18_bd_linear_ratio_closure_thm2314.out
script_sha256: 99e7a5a92be4cb75446d08b973f19fe1dfba7ad182b5ad095e8324613ce61714
output_sha256: 414d2cb756f8e3b8144582d84762f9430128c469b2a378932b8e7ed238a638a5
hash_basis: working-tree bytes (LF)
---

# THM-2314 -- the two rational B--D ratios have positive-genus spectra

THM-2311 reduces every exactly two-sparse degree-eighteen survivor to one of
`31` weighted-projective ratio points. On the `B`--`D` line, two of its six
points are rational:

```text
D/B^2=4075/85176,                 D/B^2=25/126.      (1)
```

The repeated branch values behind these factors have quite different local
meanings. The first ratio replaces two ordinary branch values by smooth
total ramification and does not lower the genus at all. The second creates
an ordinary triple point and lowers genus four to genus one. Both remain
too curved to carry a rational Keller trajectory.

## 1. Exact normalized curves

Use the scope and normalized invariants of THM-2262/2297:

```text
(B,C,D,W),                       weights (2,3,4,5),

G_0(u,y;B,C,D,W)=0,             wt(u,y)=(2,1).       (2)
```

On the off-axis `B`--`D` line, choose a constant weighted scaling with
`B=1`. This gives an isomorphic spectral curve and sends `D` to the
invariant ratio

```text
t=D/B^2.                                             (3)
```

With `C=W=0`, write

```text
G_t
 =-26040609u^3+(49601160+1607445y^2)u^2

  +(-20995200-2857680y^2-52907904t-138915y^4)u

  +777600y^2+33592320t+78120y^4+1959552ty^2+1127y^6.
                                                            (4)
```

The scaling in (3) is used only to identify isomorphic algebraic curves.
It is not treated as a target-preserving quotient of the retained Keller
one-form.

## 2. Both spectral curves are absolutely irreducible

The leading `u` coefficient in (4) is a nonzero constant. If `G_t` were
reducible over `C(y)`, its cubic degree would give a root `r(y)`. After
dividing by the leading coefficient the root is integral over `C[y]`;
because `C[y]` is integrally closed,

```text
r(y) in C[y].                                        (5)
```

If `deg r=d>2`, the term `u^3` after substitution has `y`-degree `3d`,
strictly larger than the degrees from `y^2u^2`, `y^4u`, and `y^6`.
Its leading term cannot cancel. Hence

```text
deg r<=2.                                            (6)
```

The curve is even in `y`. Thus `r(-y)` is also a root. If these roots are
distinct, the third root is their difference from the even sum of all
three roots, and is itself even. If they coincide, `r` is already even.
Consequently any polynomial root would force one of the form

```text
u=ay^2+b.                                            (7)
```

Substitute (7) in (4) and set its four coefficients in `y^2` equal to
zero. Exact Buchberger reduction over `Q` gives

```text
t=4075/85176:       reduced Groebner basis {1},

t=25/126:           reduced Groebner basis {1}.      (8)
```

The coefficient ideals remain the unit ideal after extending constants to
`C`; (7) is impossible. Therefore both `G_t` are absolutely irreducible.
Their projective normalizations are connected degree-three covers of the
`y`-line.

## 3. The common infinity fibre is unramified

At infinity use

```text
r=1/y,                         v=u/y^2.              (9)
```

For either ratio,

```text
r^6G_t(v/r^2,1/r)=L_infinity(v)+O(r^2),

L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (10)
```

Its discriminant is

```text
Disc_v L_infinity
 =-153384762202971019112448 !=0.                     (11)
```

The implicit-function theorem therefore gives three distinct smooth
branches with local parameter `r`. They are the three points over infinity,
and all are unramified for projection to the `y`-line.

## 4. The ratio 4075/85176 retains genus four

Put

```text
t_1=4075/85176.
```

The exact branch discriminant is

```text
Delta_1(y)
 =-56684737689882624/4826809

  *(91y^2+1215)^2

  *(1577224103y^8+203464170y^6-147517112925y^4
      +1389005982000y^2-3037628182500).              (12)
```

Call the degree-eight factor `h_8`. Exact Euclidean reduction gives

```text
gcd(h_8,h_8')=1,             gcd(h_8,91y^2+1215)=1. (13)
```

Thus the eight roots of `h_8` are distinct simple branch values. At each,
the normalization has one ramification point of index two.

The squared quadratic in (12) is not a pair of singular collisions. Put
`x=y^2`. At its fibre

```text
x_0=-1215/91
```

the whole cubic becomes

```text
G_t(u,x_0)=-729/15379*(819u-295)^3.                 (14)
```

At the triple root `u_0=295/819`,

```text
partial_x G_t(u_0,x_0)=-16329600/169 !=0.           (15)
```

Both roots `y_0` of `91y^2+1215` are nonzero, so
`partial_yG_t=2y_0 partial_xG_t` is nonzero. The curve is smooth there.
Solving locally for `y-y_0` shows that it has order three in the parameter
`u-u_0`; each fibre is one totally ramified point of index three.

There are no other finite branch values by (12), and Section 3 removes
infinity. The total ramification is therefore

```text
8*(2-1)+2*(3-1)=12.                                 (16)
```

Riemann--Hurwitz gives

```text
2g-2=3*(-2)+12=6,                  g=4.              (17)
```

## 5. The ratio 25/126 has an elliptic normalization

Put

```text
t_2=25/126.
```

This time

```text
Delta_2(y)
 =-19442865027629740032 y^6

  *(7889y^6+211680y^4+1814400y^2+2916000).           (18)
```

Call the sextic factor `h_6`. Exact reduction gives

```text
gcd(h_6,h_6')=1,                    h_6(0)!=0.       (19)
```

Its six roots are distinct simple branch values, each with ramification
index two.

The factor `y^6` records a singular crossing rather than ramification.
Shift

```text
U=u-40/63.
```

The exact local equation at `y=0` is

```text
G_t
 =-26040609U^3+1607445U^2y^2-816480Uy^2
   -138915Uy^4-10080y^4+1127y^6.                   (20)
```

Its tangent cone is

```text
-5103U(5103U^2+160y^2).                             (21)
```

Over `C`, (21) is a product of three distinct lines, none equal to
`y=0`. Thus `(40/63,0)` is an ordinary triple point with delta invariant
`3`. Its normalization has three local branches, and `y` is a uniformizer
on each. All three are unramified over `y=0`.

Consequently the only ramification consists of the six simple roots in
(19):

```text
sum_P(e_P-1)=6.                                     (22)
```

Riemann--Hurwitz now gives

```text
2g-2=3*(-2)+6=0,                   g=1.              (23)
```

The ordinary-triple delta `3` independently explains the drop from the
genus-four generic spectrum to genus one.

## 6. Positive genus closes both Keller ratios

A putative Keller trajectory supplies

```text
(u,y) in C(x)^2,                   G_t(u,y)=0.        (24)
```

If nonconstant, (24) extends to a nonconstant morphism from `P^1` to the
connected projective normalization. Riemann--Hurwitz forbids such a map
when the target genus is `4` or `1`. Hence `u` and `y` are constant.

The wall `y=0` is already empty by THM-2262, so the inherited first-flux
identity

```text
Z=T^2=-2N_2/(5103y)                                 (25)
```

makes `Z`, then nonzero `T`, and then `q` constant. The genuine nonsplit
deck fixes the algebraically closed constant field but sends `q` to `-q`.
This contradicts `q!=0`.

Therefore neither ratio in (1) can carry a degree-eighteen Keller
trajectory in the stated branch.

## 7. Consequence, loss ledger, and scope

The `B`--`D` bank of THM-2311 shrinks from six points to the four roots of

```text
22656250-772734375t+7600635000t^2
 -30805790400t^3+46376717184t^4.                    (26)
```

Across all exactly two-sparse planes, the unresolved bank shrinks from
`31` to `29`.

```text
source:
  THM-2311's two rational B--D ratio points;

map:
  weighted normalization B=1, followed by normalization of G_t;

preserved:
  the labelled B--D orbit, spectral-curve isomorphism class,
  normalization components, ramification, and genus;

temporarily forgotten:
  the scaling of the third flux, Keller one-form, and whole-Faber sidecar;

why no restoration is needed here:
  positive genus already forces every rational spectral trajectory to be
  constant, before any differential scaling can matter;

hostile control:
  the nearby representative D=1 has squarefree branch discriminant,
  detecting an accidentally generic or identically repeated computation.
                                                            (27)
```

This theorem closes only the two ratios in (1), inside the genuine
nonsplit polynomial exact-square-prefix degree-eighteen branch. The four
quartic `B`--`D` ratios, the other `25` two-sparse ratios, every
three-/four-sparse singular stratum, split/even-leading descent, other
Newton edges, `JC(2)`, and `DC(2)` remain open.

## 8. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
python3 -O 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
```

Both runs are byte-identical to the stored output. The companion verifies
the two complete discriminant factorizations; squarefreeness and
coprimality; the common separable infinity cubic; both exceptional local
charts; the ordinary-triple tangent cone; absolute irreducibility through
unit Groebner bases; and a squarefree `D=1` hostile control. The local
normalization, Riemann--Hurwitz, and deck arguments are the mathematical
proof above rather than delegated computer conclusions.
