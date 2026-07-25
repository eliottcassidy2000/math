---
id: THM-2314
title: "Degree-eighteen B--D ratio-bank closure"
status: >
  PROVED + VERIFIED-EXACT. In the genuine nonsplit polynomial
  exact-square-prefix degree-eighteen branch of THM-2262/2297, all six
  weighted ratios in THM-2311's B--D bank are empty. After B=1
  normalization, every spectral curve is absolutely irreducible. The two
  rational ratios have normalization genera four and one: their repeated
  branch factors encode respectively smooth total cubic ramification and
  one ordinary triple point. The four roots of the irreducible quartic
  factor uniformly have two ordinary nodes and eight simple branch points,
  hence normalization genus two. Positive genus makes every rational
  Keller trajectory constant and yields the inherited nonsplit-deck
  contradiction. Exactly 25 two-sparse ratios on the other four planes
  remain; this does not prove JC(2).
source: codex-2026-07-25-degree18-bd-linear-ratios
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
script: 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
output: 05-knowledge/results/jc2_degree18_bd_linear_ratio_closure_thm2314.out
script_sha256: 9d81ea499269344d23aded1dd1bacc6faad729d37dc2ef5e3d63e00e68e48a29
output_sha256: fe597f0aedf5a1fb8a4ec800cded5e7471fb45ebc04e0effefd62603adc0dcd6
hash_basis: working-tree bytes (LF)
---

# THM-2314 -- the full B--D ratio bank has positive-genus spectra

THM-2311 reduces every exactly two-sparse degree-eighteen survivor to one of
`31` weighted-projective ratio points. On the `B`--`D` line, two of its six
points are rational:

```text
D/B^2=4075/85176,                 D/B^2=25/126.      (1)
```

The other four are the roots of one irreducible quartic. The repeated branch
values behind the three factors have different local meanings. The first
rational ratio replaces two ordinary branch values by smooth total
ramification and does not lower the genus. The second creates an ordinary
triple point and lowers genus four to genus one. Each quartic ratio creates
two ordinary nodes and lowers genus four to genus two. All six remain too
curved to carry a rational Keller trajectory.

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

## 2. All six spectral curves are absolutely irreducible

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
zero. At the rational points, exact Buchberger reduction over `Q` gives

```text
t=4075/85176:       reduced Groebner basis {1},

t=25/126:           reduced Groebner basis {1}.      (8)
```

For the remaining points, let the primitive polynomial

```text
p(t)
 =46376717184t^4-30805790400t^3+7600635000t^2
   -772734375t+22656250.                             (8a)
```

Modulo `11`, its monic reduction is

```text
t^4+4t^2-t+2.
```

Writing `bar p` for this monic reduction, the Rabin certificate

```text
gcd(bar p,t^(11^2)-t)=1,   t^(11^4)-t=0 mod bar p   (8b)
```

proves that this reduction, and hence (8a), is irreducible. Put
`K=Q(alpha)` for one root. The same four coefficient equations from (7)
have reduced Groebner basis `{1}` in `K[a,b]`. A unit ideal stays a unit
after passing to an algebraic closure and under every complex embedding of
`K`. Thus (7) is impossible at all four conjugates as well.

All six `G_t` are therefore absolutely irreducible. Their projective
normalizations are connected degree-three covers of the `y`-line.

## 3. The common infinity fibre is unramified

At infinity use

```text
r=1/y,                         v=u/y^2.              (9)
```

For every ratio in the bank,

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

## 6. Every quartic ratio has a genus-two normalization

Let `alpha` be a root of the irreducible polynomial (8a). Its discriminant
is the nonzero integer

```text
-499128191381233551206575897907721679687500000000000000.
                                                            (23a)
```

The four complex embeddings of `K=Q(alpha)` therefore give exactly the four
remaining `B`--`D` ratios. In the power basis of `K`, put

```text
X
 =36/212384375

  *(40119541056alpha^3-21007917000alpha^2
      +3603521250alpha-206640625).                  (23b)
```

This is nonzero: its displayed representative is a nonzero polynomial of
degree three, while the minimal polynomial of `alpha` has degree four.
Exact division in `K[y]` gives

```text
Delta_alpha(y)
 =-153384762202971019112448

  *(y^2-X)^2 h_8(y),                                (23c)

deg h_8=8,

gcd(h_8,h_8')=1,                 gcd(h_8,y^2-X)=1.  (23d)
```

Thus `h_8` gives eight distinct simple branch values. It remains to
interpret the squared quadratic. Write `H(u,x)=G_alpha(u,y)` with `x=y^2`
and put

```text
R
 =4/1911459375

  *(139498220736alpha^3-60035887800alpha^2
      +8517363750alpha-89421875).                   (23e)
```

At `x=X`, the cubic has one simple root and the double root `R`. Exact
reduction gives

```text
H(R,X)=H_u(R,X)=H_x(R,X)=0.                         (23f)
```

Hence each of the two points

```text
(u,y)=(R,+sqrt(X)),              (R,-sqrt(X))
```

is singular. In local coordinates `U=u-R`, `Y=y-y_0`, the coefficient of
`U^2` in the tangent cone is

```text
A
 =-472392/4334375

  *(108948478464alpha^3-37534039200alpha^2
      +3767242500alpha-72640625),                   (23g)
```

and its tangent discriminant is

```text
Theta
 =78364164096/173375

  *(13250490624alpha^3-29352342600alpha^2
      +9881156250alpha-904703125).                  (23h)
```

Both are nonzero for the same degree reason as `X`. Therefore the tangent
cone has two distinct lines, and neither is `Y=0`. Each singularity is an
ordinary node of delta invariant one; its two normalization branches are
unramified for the `y`-projection.

The node factor in (23c) contributes no normalization ramification.
The eight simple roots in (23d) each contribute one, and Section 3 removes
infinity. Therefore, uniformly under all four embeddings,

```text
sum_P(e_P-1)=8,

2g-2=3*(-2)+8=2,                   g=2.              (23i)
```

Equivalently, the two ordinary nodes account for the genus drop from four
to two.

## 7. Positive genus closes all six Keller ratios

A putative Keller trajectory supplies

```text
(u,y) in C(x)^2,                   G_t(u,y)=0.        (24)
```

If nonconstant, (24) extends to a nonconstant morphism from `P^1` to the
connected projective normalization. Riemann--Hurwitz forbids such a map
when the target genus is `4`, `1`, or `2`. Hence `u` and `y` are constant.

Undoing the constant weighted curve isomorphism of Section 1 shows that the
original retained `u` and `y` are constant as well. The wall `y=0` is
already empty by THM-2262, so its inherited first-flux identity

```text
Z=T^2=-2N_2/(5103y)                                 (25)
```

makes `Z`, then nonzero `T`, and then `q` constant. The genuine nonsplit
deck fixes the algebraically closed constant field but sends `q` to `-q`.
This contradicts `q!=0`.

Therefore none of the six `B`--`D` ratios can carry a degree-eighteen
Keller trajectory in the stated branch.

## 8. Consequence, loss ledger, and scope

The entire six-point `B`--`D` bank of THM-2311 is empty. Across all exactly
two-sparse planes, the unresolved bank shrinks from `31` to `25`.

```text
source:
  THM-2311's full six-point B--D ratio bank;

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
                                                            (26)
```

This theorem closes only the `B`--`D` bank, inside the genuine nonsplit
polynomial exact-square-prefix degree-eighteen branch. The other `25`
two-sparse ratios, every three-/four-sparse singular stratum,
split/even-leading descent, other Newton edges, `JC(2)`, and `DC(2)` remain
open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
python3 -O 04-computation/jc2_degree18_bd_linear_ratio_closure_thm2314.py
```

Both runs are byte-identical to the stored output. The companion verifies
the two rational discriminant factorizations; the mod-`11` irreducibility
certificate and uniform quartic-field factorization; every squarefree and
coprime residual; the common separable infinity cubic; the smooth total
ramification, ordinary triple, and ordinary-node local charts; absolute
irreducibility through unit Groebner bases over `Q` and `K`; and a
squarefree `D=1` hostile control. The local normalization,
Riemann--Hurwitz, and deck arguments are the mathematical proof above
rather than delegated computer conclusions.
