---
id: THM-2071
title: "Quadratic-fiber rigidity and tame normal form for planar Keller pencils"
status: >
  PROVED. If one member P of a complex planar Keller pair has degree two
  along a linear source fiber, then its leading fiber coefficient is constant,
  its complementary component has reduced fiber degree one, and the pair has
  an explicit tame normal form. The intermediate square/parity descent is
  strengthened by a centered parity decomposition and a noncancellation law
  with central-binomial coefficient. Consequently no hypothetical JC(2)
  counterexample has a quadratic member in any output-pencil/source-foliation
  direction. This is a pencil-direction theorem, not a statement about generic
  cover degree or a Jelonek classification.
source: codex-2026-07-21-JC2-quadratic-fiber
depends_on: []
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - MISTAKE-237
script: 04-computation/jc2_quadratic_fiber_square_gate_codex_20260721.py
output: 05-knowledge/results/jc2_quadratic_fiber_square_gate_codex_20260721.out
script_sha256: 3d5ce81db8601a3035db28ae63bf3d003d4f72372d304c3250098a64a8efb267
output_sha256: 2b692ddf3606a0226bb22b88b0a2060a5ecff951b41ba02a8dee3af111821030
hash_basis: repository blobs with LF line endings
---

# THM-2071 -- quadratic-fiber rigidity and tame normal form

Let

```text
P=A(x)y^2+B(x)y+C(x),       A!=0,
Q in C[x,y],                J(P,Q)=P_x Q_y-P_y Q_x=kappa in C*.
```

For this fixed source-fiber direction define the reduced complementary degree

```text
mu_y(P,Q) = min_H deg_y(Q-H(P)),       H in C[T].        (1)
```

Every polynomial in the minimum is nonzero, because `Q=H(P)` would have zero
Jacobian with `P`. Thus (1) is an ordinary minimum of a nonempty subset of
the natural numbers.

The conclusions are:

1. `A` is constant.
2. `mu_y(P,Q)=1`.
3. The pair is tame and has the explicit normal form in Section 3 below.

The proof first gives two useful intermediate constraints: the reduced degree
is odd, and `A=aU^2` is a square up to a nonzero scalar. A centered parity
argument then rules out every nonconstant `U` simultaneously. Thus not only a
nonsquare `A`, but every nonconstant top fiber coefficient, admits no
polynomial Jacobian mate in this quadratic-fiber setting.

## 1. Top-coefficient law

Choose `H` attaining (1), put `Q_0=Q-H(P)`, and write

```text
Q_0=sum_(j=0)^n q_j(x)y^j,       q_n!=0,       n=mu_y(P,Q).
```

Target shear does not change the Jacobian, so `J(P,Q_0)=kappa`. For `n>=1`,
the coefficient of `y^(n+1)` is

```text
n A'(x)q_n(x)-2A(x)q_n'(x)=0.                         (2)
```

Equivalently,

```text
(q_n^2/A^n)'=0,
```

and characteristic zero gives

```text
q_n^2=c A^n,                     c in C*.             (3)
```

The case `n=0` is impossible: then `Q_0=S(x)` and

```text
J(P,S)=-(2Ay+B)S'.                                    (4)
```

The coefficient of `y` in (4) forces `S'=0`, contradicting `kappa!=0`.

## 2. Even descent, odd stopping degree, and the square gate

Suppose `n=2k` is positive and even. From (3),

```text
q_n=d A^k                                  (d in C*),  (5)
```

because `C(x)` is a field and `C` contains a square root of `c`. But `P^k`
has top fiber term `A^k y^(2k)`. Therefore

```text
Q_0-dP^k
```

has smaller `y`-degree and the same Jacobian with `P`, contradicting the
minimality of `n`. Hence `n` is odd.

Now factor (3) in the UFD `C[x]`. For every irreducible `pi`,

```text
2 ord_pi(q_n)=n ord_pi(A).
```

Since `n` is odd, every `ord_pi(A)` is even. Thus `A=aU^2`. This proves the
square gate.

The same argument is an explicit terminating descent, not just a minimality
argument: whenever the current complementary degree is even, (5) identifies
the unique top target shear `dP^k`; subtract it and continue. The descent
cannot stop at degree zero, so it stops at an odd degree. This is one rigorous
quadratic-face instance of the leading-form propagation sought in THM-1345;
it says nothing about arbitrary Newton faces.

## 3. Reduced degree one: exact tame normal form

Suppose `n=1`, so after a target shear

```text
Q_0=L(x)y+M(x).
```

The coefficients of `y^2`, `y`, and `1` in `J(P,Q_0)=kappa` are respectively

```text
A'L-2AL'                         =0,                   (6)
B'L-2AM'-BL'                    =0,                   (7)
C'L-BM'                         =kappa.               (8)
```

Here `L!=0`, since `L=0` makes (7)--(8) incompatible with `kappa!=0`. Equation
(6) gives `L^2/A` constant, so write

```text
A=aL^2,                         a in C*.              (9)
```

After division by `L^2`, equation (7) becomes

```text
(B/L)'=2aM'.                                           (10)
```

Thus `B/L-2aM` is constant in `C(x)`. In particular

```text
V=B/L in C[x],                 M=V/(2a)+m_0.          (11)
```

Substituting (9)--(11) into (8) yields

```text
kappa=L (C-V^2/(4a))'.                                (12)
```

Both factors on the right are polynomials and their product is a nonzero
constant. Therefore

```text
L=ell in C*,       A=a ell^2 in C*,
C-V^2/(4a)=(kappa/ell)x+d_0.                          (13)
```

Set

```text
Y=y+B(x)/(2A).
```

Equations (11)--(13) give

```text
P=A Y^2+(kappa/ell)x+d_0,
Q-H(P)=ell Y+m_0.                                     (14)
```

For target coordinates `(u,v)=(P,Q)`, the inverse is

```text
Y=(v-H(u)-m_0)/ell,
x=(ell/kappa)(u-A Y^2-d_0),
y=Y-B(x)/(2A).                                        (15)
```

It is polynomial. The construction consists of a triangular target shear, a
triangular source shear, and affine/elementary triangular maps, so the map is
tame. This degree-one endpoint can also be fed to THM-2063 after the nonlinear
target shear; (6)--(15) record the stronger normal form native to the
quadratic descent.

## 4. Centered parity rigidity closes every odd residual degree

Absorb the scalar in the square gate into `U`, so `A=U^2`. Work temporarily
over `K=C(x)` and set

```text
h=B/(2U),             z=Uy+h,             D=C-h^2.
```

Then `P=z^2+D`. If a polynomial `Q_0(x,y)` is written as an element of
`K[z]`, the chain rule gives the exact identity

```text
J_xy(P,Q_0)=U L(Q_0),
L=D' partial_z-2z partial_x,                             (16)
```

where `partial_x` in `L` holds `z` fixed. Notice that `L(P)=0`, `L` sends the
even part in `z` to an odd polynomial, and it sends the odd part to an even
polynomial.

Split `Q_0=Q_ev+Q_odd` by `z`-parity. Since `kappa/U` is even, (16) gives

```text
L(Q_ev)=0,                 L(Q_odd)=kappa/U.            (17)
```

The first equation loses no polynomial information. Indeed, write the even
part uniquely as

```text
Q_ev=sum_i c_i(x)P^i,                 c_i in K.
```

Because `L(P)=0`,

```text
L(Q_ev)=-2z sum_i c_i'(x)P^i.
```

The powers of `P=z^2+D` are independent over `K`, so every `c_i'=0` and
every `c_i` lies in `C`. Thus `Q_ev=H(P)` for an honest `H in C[T]`. It can be
removed by a target shear even though the parity split was performed over
`K`. After this shear, rename the nonzero odd part as `Q_0` and write

```text
Q_0=sum_(j=0)^r a_j(x)z^(2j+1),          a_r!=0.        (18)
```

The highest coefficient in `L(Q_0)` says `a_r=c in C*`. The coefficients of
`z^(2j)`, followed by the constant coefficient, give the full triangular
system

```text
2a_(j-1)'=(2j+1)D'a_j,       j=r,r-1,...,1,
D'a_0=kappa/U.                                        (19)
```

In particular `D'` is not identically zero. Starting from `a_r=c`, (19)
shows recursively that every `a_j` is a polynomial in `D`; more precisely,

```text
a_(r-k)(D)=c binom(r+1/2,k)D^k + lower powers of D,
                                      0<=k<=r.          (20)
```

Now use the information that was discarded by merely working over `K[z]`:
the original constant fiber coefficient

```text
q_0(x)=Q_0(x,0)=sum_(j=0)^r a_j(D)h^(2j+1)             (21)
```

is a polynomial in `x`. If `h` had a finite pole, then `D=C-h^2` would have
leading polar part `-h^2`. By (20), every top-degree contribution in (21)
has the same polar order, and their total coefficient is

```text
c S_r,
S_r=sum_(k=0)^r (-1)^k binom(r+1/2,k)
   =(-1)^r binom(2r,r)/4^r !=0.                        (22)
```

The last equality is the finite alternating-binomial identity. Hence the
leading pole cannot cancel. This contradicts `q_0 in C[x]`, so `h` has no
finite pole. Therefore

```text
h in C[x],                    equivalently U divides B. (23)
```

It follows that `D=C-h^2` and `a_0(D)` are polynomials. The last equation in
(19) now says that the polynomial `D'a_0(D)` equals `kappa/U`. Consequently
`1/U` is a polynomial, so `U` is constant. Let `d=deg D`; the same equation
forces `d>=1`. Formula (20) says `deg_T a_0=r`, and therefore

```text
0=deg(D'a_0(D))=(r+1)d-1.
```

Thus `r=0` and `d=1`. The reduced complementary degree is `2r+1=1`, and
Section 3 supplies the tame normal form and polynomial inverse. This also
shows directly that the apparent nonconstant-square residual was empty at
every odd degree, not only in the cubic case.

## 5. The cubic system and the rational-integrability audit

For `r=1`, retain an inessential constant target term and write

```text
Q_0=c z^3+E(x)z+F(x).
```

Equation (16) is equivalent coefficient by coefficient to

```text
E'=(3c/2)D',             F'=0,             D'E=kappa/U.
```

Hence, for constants `e,f in C`,

```text
E=(3c/2)D+e,             F=f,
((3c/4)D^2+eD)'=kappa/U.                           (24)
```

On returning to `y`, the linear coefficient is

```text
q_1=U((3c/2)C+e)+(3c/8)B^2/U.                     (25)
```

Thus polynomiality of `q_1` already gives `U | B^2`. The constant coefficient
is stronger:

```text
q_0=-(c/2)h^3+((3c/2)C+e)h+f.                     (26)
```

At a finite pole of `h`, the first term in (26) cannot cancel, so `h` is a
polynomial and `U | B`, exactly the `r=1` instance of (22).

Equation (24) also exposes an independent rational-integrability obstruction.
For nonconstant `U in C[x]`, a rational function `R in C(x)` satisfying

```text
R'=1/U                                                   (27)
```

exists if and only if `U=u(x-alpha)^m` with `u in C*` and `m>=2`. To prove
necessity, write the distinct-root multiplicities of `U` as `m_i`. A rational
derivative has no simple pole, so every `m_i>=2`, and `R` has pole order
`m_i-1` at that root. If `d=deg U` and there are `s` distinct roots, the
rational-map degree of `R` is therefore `d-s`. At infinity, integrating
`R'=1/U` shows that `R-R(infinity)` has a zero of order `d-1`. A fiber
multiplicity cannot exceed the map degree, so `d-1<=d-s`, whence `s=1`.
The converse is immediate by integrating a negative power. Applied to (24),
this lemma reduces the cubic branch to a single-root power even before (26);
the pole noncancellation then eliminates that last family.

## 6. Pencil and coordinate invariance

The gate is intrinsic to a pair consisting of an output direction and a
linear source foliation.

Fix `P` and replace a linear target complement by

```text
Q'=alpha Q+beta P,                  alpha!=0.
```

As `H` ranges over `C[T]`, so does

```text
K(T)=(H(T)-beta T)/alpha,
```

and

```text
Q'-H(P)=alpha (Q-K(P)).
```

Thus `mu_y(P,Q')=mu_y(P,Q)`. Scaling `P` also changes neither the algebra
`C[P]` nor the square class of `A`, because every nonzero complex scalar is a
square. Affine changes preserving the chosen source foliation send

```text
x -> ax+b,                 y -> cy+linear(x),
```

so they preserve fiber degree, reduced degree, and the square property of the
top coefficient. For a different source direction, one simply applies the
theorem in a new affine coordinate system.

Consequently, for an arbitrary planar Keller map and any nonzero output-pencil
member `R`, any linear complement defines the same reduced degree. Combining
this theorem with THM-2063 gives the following all-directions gate for a
hypothetical counterexample:

```text
deg_fiber R <= 1    is impossible by THM-2063;

deg_fiber R = 2     is impossible by THM-2071.                         (28)
```

## 7. Newton/Jelonek frontier after closing the quadratic face

In the Newton picture, `A(x)` is the coefficient polynomial on the top
`y=2` row. Sections 2 and 4 show successively that its zero divisor would
have to be even and then that it must be empty. Thus the entire quadratic
fiber face is closed; there is no surviving nonconstant-square/cubic cell.
The next unresolved source-fiber degree for a hypothetical counterexample is
degree at least three in the chosen direction.

This remains a rank-one-at-infinity coefficient theorem, not VC(4). Neither
`mu_y` nor the degree of `P` along one linear fiber has been identified with a
generic function-field cover degree. THM-1330's Jelonek/monodromy constraints
remain separate necessary data, and MISTAKE-237 forbids calling this descent
equivalent to the Jelonek, Vanishing-Conjecture, or Lame-for-polygons programs.

## 8. Exact symbolic referee

The companion script verifies:

- (2) for complementary degrees `1` through `7`;
- the even target-shear cancellation for powers `P^1` through `P^4`;
- the three equations (6)--(8);
- the centered cubic system (24), including (25);
- the triangular odd recurrence and central-binomial noncancellation (20)--(22)
  for `r=0` through `6`;
- a nontrivial rational instance of (14), its constant Jacobian, and both
  compositions of (15);
- bounded hostile linear-system searches against a nonsquare top coefficient
  and a square-but-nonconstant control.

The bounded searches are controls only; the centered proof above is
degree-uniform. The stored output ends in `RESULT: PASS`. QED.
