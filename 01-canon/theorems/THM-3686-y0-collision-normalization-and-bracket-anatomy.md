---
id: THM-3686
title: "The y=0 collision surface: normalization, boundary two-form, and bracket anatomy"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT.  On y=0, the explicit weighted
  quartic Keller map of THM-3438 has image ring
  R=C[A,B,C]/(16A^3C^2+27B^4-36AB^2).  Its normalization is the smooth
  surface obtained by adjoining g=1-x^2z, with three explicit relations;
  the source plane is the complement of two boundary lines.  The retained
  collision consists of the g=+1 and g=-1 normalization points over
  (A,B,C)=(0,0,1), while g=0 is a third, missing point of that fibre.  The
  source parametrization is unramified with everywhere-injective differential,
  1 is a sum of two Poisson
  brackets of elements of R, and the remaining one-bracket equation is
  exactly a planar Jacobian-counterexample construction problem.  No
  homogeneous pair works in any degree, and an exact Groebner gate excludes
  every pair through total target degree three.  The transported source
  two-form has divisor the missing g=0 boundary line, giving two exact
  boundary-jet laws for any survivor.  Arbitrary mixed-weight degree-at-
  least-four one-bracket pairs, JC(2), and arbitrary quartic C3 data remain
  OPEN; no counterexample is claimed.
source: codex-jc-quartic-c3-construction / 2026-08-22
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
related:
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3276-derivative-normalized-inverse-different-line-and-all-finite-jet-hostile
  - THM-3441-weighted-quartic-jelonek-components-and-inertia-parity
  - THM-3685-weighted-quartic-native-polynomial-graph-descent-no-go
scripts:
  - 04-computation/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.py
  - 04-computation/jacobian_y0_target_degree3_one_bracket_no_go_thm3686.py
outputs:
  - 05-knowledge/results/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.out
  - 05-knowledge/results/jacobian_y0_target_degree3_one_bracket_no_go_thm3686.out
script_sha256:
  - f1b92216adc03101227861caccb5e02b19f7d8e52e9be7aa14cfa5b50b84f1cb
  - 04b31009f7808d434e89a8ea99ac7d40250ac16f9920d60e17236bbef12612f8
output_sha256:
  - e989173b9afabbbcbda2ba06a596ed422fc3cc9094ed4c618874c48b20714b3e
  - aed11227c2637b890632f8c54f55ec596ff93550fcdbbf1ff8757c23ffd1bf90
hash_basis: LF-normalized bytes
---

# THM-3686 -- the retained collision is a normalization-and-one-bracket problem

**PROVED + VERIFIED-EXACT + FINITE-EXACT.**

This theorem follows the simplest source surface containing the known
collision of the THM-3438 three-dimensional Keller map.  The restriction is
not itself a planar Keller map.  Instead, it produces a singular target
surface whose normalization remembers the lost sheet coordinate.  That
turns the counterexample search into one precise algebraic question:

```text
does this collision-retaining target subalgebra contain P,Q with {P,Q}=1? (1)
```

We determine the normalization, its missing boundary, the canonical-form
debt in `(1)`, the additive bracket width, the complete homogeneous no-go,
and the first finite mixed-weight gate.  Work over `C`, and use

```text
{P,Q}=P_x Q_z-P_z Q_x.                                  (2)
```

## 1. Exact image ring on `y=0`

Restrict the weighted quartic map of THM-3438 to the source plane `y=0`.
With

```text
t=x^2z,
A=3z(2-t),
B=2xz(2-t),
C=x(1-t),                                                (3)
```

direct expansion gives

```text
H(A,B,C)=16A^3C^2+27B^4-36AB^2=0.                       (4)
```

Let `phi:C[A,B,C]->C[x,z]` be the substitution `(3)`.  The bracket

```text
{A,B}=-12z(t-2)(t-1)                                    (5)
```

is nonzero, so `A,B` are algebraically independent and `ker(phi)` has
height one.  The polynomial `(4)` is irreducible.  Indeed, as a primitive
quadratic in `C` over `C[A,B]`, its discriminant is

```text
576 A^3 B^2(4A-3B^2).                                   (6)
```

The last factor has odd valuation in `C(A,B)`, so `(6)` is not a square;
Gauss's lemma applies.  Therefore

```text
R:=im(phi)=C[A,B,C]/(H).                                 (7)
```

The singular locus of `Spec(R)` is exactly the union of the two axes

```text
L_C: A=B=0,                    L_A: B=C=0.               (8)
```

This follows immediately from

```text
H_A=48A^2C^2-36B^2,
H_B=36B(3B^2-2A),
H_C=32A^3C.                                             (9)
```

Thus the retained collision lands on a singular target axis; it is not an
unexplained failure of a smooth parametrization.

## 2. The source differential is nevertheless injective everywhere

The three source minors are

```text
{A,B}=-12z(t-2)(t-1),
{A,C}=-6(2t^2-4t+1),
{B,C}=-2x(3t^2-6t+2).                                  (10)
```

They generate the unit ideal in `C[x,z]`.  More precisely, two already
suffice:

```text
1=(-1/6+2t-t^2){A,C}+2xz(t-2){B,C}.                    (11)
```

So the map `A^2_(x,z)->Spec(R)` is unramified, equivalently its differential
has rank two everywhere.  It is **not** an immersion in the
algebraic-geometric sense, because `(24)` makes it noninjective.  The
collision must be explained by normalization branches and nonproper
boundary, not by a critical point on the source plane.

## 3. Exact smooth normalization

Adjoin the source sheet coordinate

```text
g=1-x^2z=1-t.                                          (12)
```

In the common function field it satisfies

```text
g=2AC/(3B),                     g^3-g+BC/2=0.           (13)
```

It is therefore integral and birational over `R`.  Put

```text
S=C[A,B,C,g]/I,                                        (14)

I=(BC-2g(1-g^2),
   3Bg-2AC,
   4A(1-g^2)-3B^2).                                   (15)
```

Equations `(15)` give the exact presentation of `R[g]`.  One quick domain
proof is to project to the `g`-line.  Off `g(1-g^2)=0`, choose `B!=0`; then

```text
C=2g(1-g^2)/B,
A=3B^2/[4(1-g^2)],                                    (16)
```

and the third equation follows.  This open set is irreducible of dimension
two.  Every special fibre at `g=0,+1,-1` has dimension one, and each of its
axis branches is a limit of `(16)`.  At `g=0`, holding `B` nonzero reaches
the parabola `C=0, A=3B^2/4`; holding `C` nonzero and using
`B=2g(1-g^2)/C` reaches the other line.  The same formula reaches the
`C`-axis branches at `g=+1,-1`.  For an `A`-axis branch at `g=epsilon` with
`epsilon=+1` or `-1`, hold `A=a!=0`, take `B=s`, and choose the formal branch
`g(s)` through `epsilon` satisfying `1-g(s)^2=3s^2/(4a)`; then
`C=3sg(s)/(2a)`.  The remaining origins follow by closure.  Hence `(14)` is
an irreducible two-dimensional affine variety set-theoretically.

The surface `(14)` is smooth.  The Jacobian of `(15)` is

```text
[ 0,          C,       B,  6g^2-2 ]
[-2C,        3g,     -2A,       3B ]
[4(1-g^2),  -6B,       0,     -8Ag ].                  (17)
```

All `3x3` minors vanish modulo `I`.  If its rank were below two, the `2x2`
minors `2C^2` and `6B^2` would give `B=C=0`.  Equations `(15)` then give
`g=0,+1,-1`; at each value the first and third, or first and second, rows of
`(17)` are visibly independent.  Thus the rank is exactly two everywhere.
The Jacobian criterion also makes the coordinate ring in `(14)` reduced;
together with irreducibility, it is a domain.  The natural surjection from
that ring to the literal subring `R[g]` is between two-dimensional domains
and is an isomorphism on `B!=0`, so its kernel has height zero and vanishes.
Thus `(15)` is the promised exact presentation.  Smoothness makes `S`
normal, and

```text
S is the normalization of R.                            (18)
```

## 4. The source plane is a normalization open, and the collision separates

In `S`, define the disjoint boundary lines

```text
D_0: g=0,  A=B=0              (the C-line),
D_-: g=-1, B=C=0              (the A-line).             (19)
```

The source formulas refine `(3)` to

```text
A=3(1+g)z,
B=2(1+g)xz,
C=gx.                                                    (20)
```

They give an isomorphism

```text
A^2_(x,z)  ~=  U:=S minus (D_0 union D_-).              (21)
```

For completeness, the inverse `x` glues from

```text
C/g,              2(1-g^2)/B,              3B/(2A),    (22)
```

on the corresponding principal opens.  The inverse `z` glues from

```text
A/[3(1+g)],                    (1-g)g^2/C^2.             (23)
```

The complements of those opens are exactly the lines in `(19)`, and
relations `(15)` make the formulas agree.

The known source collision is now completely transparent:

```text
(x,z)=( 1,0)  -> (A,B,C,g)=(0,0,1, 1),
(x,z)=(-1,2)  -> (A,B,C,g)=(0,0,1,-1).                 (24)
```

Over the common point `(A,B,C)=(0,0,1)`, the normalization equation is

```text
g(g-1)(g+1)=0.                                         (25)
```

The two source points are the `g=+1` and `g=-1` branches.  The third branch
`g=0` is the missing boundary `D_0`.  Thus the cubic sheet coordinate
algebraizes the collision data exactly; forgetting it is precisely what
retains the collision in `R`.

## 5. Canonical-form divisor and two boundary-jet laws

Transport the source form

```text
omega=dx wedge dz                                             (26)
```

to the normalization.  Near the generic point of `D_0`, `(C,g)` are smooth
coordinates and

```text
x=C/g,                    z=g^2(1-g)/C^2,
omega=-(g/C^2)dC wedge dg.                              (27)
```

Near the generic point of `D_-`, use `(A,B)`.  Differentiating the third
relation in `(15)` and simplifying gives

```text
omega=-(1/(4Ag))dA wedge dB.                            (28)
```

It follows that

```text
div_S(omega)=D_0.                                      (29)
```

There is a zero of order one on the missing cubic branch and no pole; the
form is nonvanishing on `U` and at the generic point of `D_-`.  This rules
out the hoped-for immediate pole contradiction, but replaces it with an
exact global incidence condition.

Indeed, if `P,Q in R` solve `(1)`, then on `S`

```text
dP wedge dQ=omega.                                     (30)
```

Hence the morphism `(P,Q):S->A^2` has ramification divisor exactly `D_0`
and is etale at the generic point of `D_-`.

The corresponding target-axis jet laws are explicit.  On `L_C`, write

```text
p(C)=P(0,0,C),       q(C)=Q(0,0,C),                    (31)
```

and evaluate all displayed partial derivatives at `(A,B)=(0,0)`.  Expansion
of

```text
B=2g(1-g^2)/C,
A=3g^2(1-g^2)/C^2                                    (32)
```

in `(30)` gives

```text
p' Q_B-P_B q'=0,                                      (33)

2[p'(3Q_A+2Q_BB)-(3P_A+2P_BB)q']
 +4[(P_B)'Q_B-P_B(Q_B)']=-1.                          (34)
```

On `L_A`, put `p(A)=P(A,0,0)`, `q(A)=Q(A,0,0)` and evaluate at `B=C=0`.
Since `C=-3B/(2A)+O(B^3)` near `D_-`, equation `(30)` gives

```text
4A[p'Q_B-P_Bq']-6[p'Q_C-P_Cq']=1.                     (35)
```

Equations `(33)`--`(35)` are the exact inverse-different/cofactor sidecar
for this collision surface: first-order dependence on the missing branch,
second-order unit response there, and a transverse unit response on the
other boundary.  They are necessary, not yet sufficient.

## 6. Grading blocks every homogeneous Darboux pair in all degrees

Use `u=g`, `h=1-u^2`, and weights

```text
wt(x)=1,               wt(z)=-2.                       (36)
```

Then

```text
A=3x^-2h,              B=2x^-1h,              C=xu.   (37)
```

Every weight-`r` element of `R` has the form `x^r p(u)`.  Direct
differentiation gives the weight-Wronskian law

```text
{x^r p(u),x^s q(u)}
 =x^(r+s+1)[s p'(u)q(u)-r p(u)q'(u)].                 (38)
```

For a homogeneous bracket to be a nonzero constant, `r+s=-1`.  After
swapping the two entries, write

```text
r=-R,                    s=R-1,             R>=1.      (39)
```

A monomial `A^iB^jC^k` of weight `-R` contains
`h^(i+j)` with `2i+j>=R`, so

```text
h^ceil(R/2) divides p.                                 (40)
```

A monomial of weight `R-1` has `k=R-1+2i+j`, so

```text
u^(R-1) divides q.                                     (41)
```

Equation `(38)` becomes

```text
R p q'+(R-1)p'q.                                      (42)
```

If `R>=2`, its leading coefficient is

```text
[R deg(q)+(R-1)deg(p)] lc(p)lc(q),                    (43)
```

which is nonzero, while `(40)`--`(41)` make its degree at least two.  If
`R=1`, `(42)` is `p q'`; either `q` is constant and the bracket is zero, or
`h|p` makes its degree at least two.  Therefore

```text
no two homogeneous elements of R have bracket 1, in any degree.          (44)
```

This is an elementary subalgebra-level manifestation of the general fact
that planar graded Keller maps cannot be counterexamples; see T. Shaska,
*Graded Keller maps and the Jacobian Conjecture*, arXiv:2607.20210v2.  That
paper is context, not a dependency of `(44)`.

## 7. Exact finite gate through target degree three

Let `V_N` be the span of all nonconstant target monomials

```text
A^iB^jC^k,                 1<=i+j+k<=N.                (45)
```

The companion places independent coefficients on arbitrary `P,Q in V_N`,
expands `{P,Q}-1` in `C[x,z]`, and computes the exact coefficient ideal over
`Q`.  Its reduced Groebner bases are

```text
N=2:  9 basis monomials, 18 variables,  55 equations, basis [1],
N=3: 19 basis monomials, 38 variables, 131 equations, basis [1].          (46)
```

Thus

```text
P,Q of total target degree at most three cannot solve {P,Q}=1.           (47)
```

This is a decomposability obstruction, not merely a failure of additive
span.  For `N=3`, the linear bracket map

```text
Lambda^2(V_3) -> C[x,z]                                (48)
```

has 171 columns, rank 75, nullity 96, and adjoining the constant `1` leaves
the rank 75.  So `1` is already a linear combination of degree-three
brackets, but is not one decomposable wedge.  For `N=2`, the rank rises from
26 to 27 after adjoining `1`, so even additive span fails there.

One explicit degree-three positive control is

```text
1={C,A/6-B^2/4}
  +(1/12){C^2,A^2C}
  -(13/6){C^3,AB^2}
  +5{BC^2,ABC}.                                        (49)
```

## 8. One is already a sum of two target-subalgebra brackets

The additive gap collapses further in target degree four.  A short
three-bracket identity is

```text
Q_0=-B^2/4+A/6+A^2C^2/6,

1={C,Q_0}+(7/4){BC^2,ABC}-(13/8){AC^3,B^2}.            (50)
```

More sharply, choose `a,b in C` satisfying

```text
108a^2+2040a-403=0,
23324a+920b^2-4459=0.                                  (51)
```

These equations have solutions and force `b!=0`.  Put

```text
alpha=(6a+104)/17,             beta=(84a-91)/(102b),

P_1=C+alpha AC^3,
Q_1=A/6-B^2/4+(1/6-3a/2)A^2C^2-(12b/7)ABC,

P_2=BC^2+beta AC^3,
Q_2=(8b/7)A+(3b/28)B^2-(15b/14)A^2C^2+(7/4)ABC.       (52)
```

Then direct exact reduction modulo `(51)` gives

```text
{P_1,Q_1}+{P_2,Q_2}=1.                                (53)
```

For provenance, use row support

```text
U=(C,BC^2,AC^3),                  V=(A,B^2,A^2C^2,ABC) (54)
```

and coefficient matrix

```text
[1/6,   -1/4,  1/6-3a/2, -12b/7]
[8b/7, 3b/28,    -15b/14,    7/4]
[a,    -13/8,           0,      b].                    (55)
```

Under `(51)`, row three is `alpha` times row one plus `beta` times row two.
Factoring this rank-two tensor gives exactly `(52)`--`(53)`.

There is no cosmetic one-bracket compression inside the four functions just
displayed.  Scale

```text
(F_0,F_1,F_2,F_3)=(P_1,Q_1,102bP_2,Q_2).               (56)
```

Over the coefficient field in `(51)`, the exact linear system

```text
sum_(i<j) w_ij {F_i,F_j}=1                             (57)
```

has the unique solution

```text
w_01=1,            w_23=1/(102b),            all cross w_ij=0.          (58)
```

The companion represents that field in the `Q`-basis `(1,a,b,ab)` and
checks uniqueness by a rank-24 rational system.  Its Pluecker scalar is

```text
w_01w_23-w_02w_13+w_03w_12=1/(102b)!=0.                (59)
```

Thus the alternating form in `(58)` has rank four and is not decomposable.
No constant linear combinations of this particular four-function span can
turn `(53)` into one bracket.  This is a minimal hostile for compression,
not a no-go after adjoining other functions of `R`.

Define the additive bracket length of `1` in `R` to be the least `ell` for
which it is a sum of `ell` brackets of elements of `R`.  Equation `(53)`
proves only

```text
bracket_length_R(1)<=2.                                (60)
```

It does **not** prove that the length is two, because an arbitrary
mixed-weight one-bracket identity remains open.

## 9. Exact counterexample frontier

If `P,Q in R` satisfy `{P,Q}=1`, then `(P,Q):A^2_(x,z)->A^2` is a planar
Keller map.  Both entries take the same value at the two points in `(24)`,
because every element of `R` factors through `(A,B,C)`.  It would therefore
be a planar Jacobian counterexample.

Conversely, this normalization lane asks for exactly such a pair.  The
current entry conditions are now rigid:

```text
* both entries must mix weights;
* at least one entry has target degree >=4;
* their two-form must vanish simply on D_0 and nowhere on D_-;
* their C-axis and A-axis jets must satisfy (33)--(35).                   (61)
```

This is the first genuinely admissible support region.  The quartic
three-dimensional collision has not disappeared; it has been concentrated
into a decomposable-wedge problem with an exact normalization boundary and
two exact jet sidecars.

The result does not exclude arbitrary mixed-weight degree-at-least-four
pairs, prove `JC(2)`, or exclude general quartic `C3` inverse-different/
cofactor data.  Native polynomial graph descents are separately closed by
THM-3685.

## 10. Exact reproduction

Run

```bash
python3 -B 04-computation/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.py
python3 -B -O 04-computation/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.py

python3 -B 04-computation/jacobian_y0_target_degree3_one_bracket_no_go_thm3686.py
python3 -B -O 04-computation/jacobian_y0_target_degree3_one_bracket_no_go_thm3686.py
```

Both modes byte-match their stored transcripts.  The second computation is
the expensive gate (about forty seconds on the recorded SymPy 1.14 setup).
All truth checks use explicit exceptions and survive optimized execution.

**QED.**
