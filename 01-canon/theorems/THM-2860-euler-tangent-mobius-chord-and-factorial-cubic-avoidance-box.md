---
id: THM-2860
title: "Euler-tangent Möbius chord and factorial cubic-avoidance box"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
  REQUESTED.  For a real two-plane in the mean-zero hyperplane of a
  four-slot space, failure of Euler transversality is exactly a Segre
  coplanarity determinant.  Its generic normal is a fractional-linear
  function of the four Euler eigenvalues; the sole denominator-degenerate
  branch is a three-slot coordinate face.  On every normalized factorial
  support 0<=a_0<a_1<a_2<a_3<=30, the resulting two homogeneous septic
  cubic-divisibility remainders are coprime modulo 1000033.  Thus no cubic
  Maxwell line in these 31465 cells is Euler-tangent.  This does not
  exclude a quartic-harmonic zero or a shared Euler-transverse line, and
  it does not make the Euler multiplier observable from scalar moments.
source: root/gmc-euler-tangent-mobius-chord-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
related:
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2849-four-slot-first-window-macaulay-box
script: 04-computation/gmc_euler_tangent_mobius_chord_thm2860.py
output: 05-knowledge/results/gmc_euler_tangent_mobius_chord_thm2860.out
script_sha256: 379a441ff3e2d64b4047f22112cb28e14a577f0bc65bbbad7d12af8bda0c72db
output_sha256: 098bb831662b6b6d5ee4d14c56bbbd2bb4365a61c64cf0992fef3be94f4d846b
hash_basis: LF-normalized bytes
---

# THM-2860 -- Euler-tangent Möbius chord and factorial cubic avoidance

**PROVED + FINITE-EXACT + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
REQUESTED.**

THM-2848 reduces a factorial four-slot common point to a real chord line
of the quadratic null conic.  This theorem measures whether that chord is
tangent to its own Euler drift.

The answer is exact.  Euler-tangent chords form one dual conic, whose
generic points are fractional-linear functions of the four support
exponents.  More importantly, an exact census proves that none of the
factorial **cubic** multipole lines through exponent thirty lies on that
conic.  The result supplies a transverse coordinate, not its scalar
observation: the quartic and multiplier-access boundaries remain open.

## 1. The intrinsic Euler determinant

Let

```text
V=R^4,                 epsilon(x)=sum_(j=0)^3 x_j,
h=ker(epsilon),
A=diag(a_0,a_1,a_2,a_3),              a_0<a_1<a_2<a_3. (1)
```

For a real two-plane `E subset h`, choose an ordered basis `u,v` and put

```text
c=u+iv,
Delta_A(E;c)=det_C[c,cbar,Ac,A cbar].                  (2)
```

Direct column conversion gives

```text
Delta_A(E;c)=-4 det_R[u,v,Au,Av] in R.                 (3)
```

Thus `(2)` is real.  Its gauge laws are

```text
Delta_A(E;lambda c)=|lambda|^4 Delta_A(E;c),           (4)

det[u',v',Au',Av']
 =det(B)^2 det[u,v,Au,Av]                              (5)
```

when `lambda in C*` and `(u',v')=(u,v)B`,
`B in GL_2(R)`.  Hence the zero predicate, and after fixing the ordered
slot orientation also its sign, depends only on the plane.

Let `pi:V -> V/E`.  Equation `(3)` says exactly that

```text
Delta_A(E)=0
 iff det(pi o A|_E)=0
 iff there is 0!=T in E with AT in E.                  (6)
```

The vectors `T,AT` are independent.  Otherwise `T` would be an
`A`-eigenvector.  Since the eigenvalues are distinct, such a vector is a
coordinate vector, and no nonzero coordinate vector lies in `h`.
Consequently the last condition in `(6)` has the sharper form

```text
E=span_R{T,AT},
epsilon(T)=epsilon(AT)=0.                              (7)
```

This is the Euler-tangent chord normal form.

## 2. Segre coplanarity and the Möbius normal

Write a plane in `h` as

```text
E=ker(epsilon) intersect ker(b),                       (8)
```

where the real covector `b=(b_0,b_1,b_2,b_3)` is taken modulo scaling
and addition of a constant covector.  A vector `T` as in `(7)` exists
exactly when

```text
det
 [ 1        1        1        1       ]
 [ a_0      a_1      a_2      a_3     ]
 [ b_0      b_1      b_2      b_3     ]
 [ a_0 b_0  a_1 b_1  a_2 b_2  a_3 b_3 ] =0.           (9)
```

Indeed, the four rows in `(9)` impose respectively

```text
epsilon(T)=0,       epsilon(AT)=0,
b(T)=0,             b(AT)=0.                           (10)
```

This proves both directions of `(9)`.

Equivalently, the four Segre points

```text
[1:a_j:b_j:a_jb_j] in P^3                            (11)
```

are coplanar.  Thus there are real
`alpha,beta,gamma,delta`, not all zero, with

```text
alpha+beta a_j+(gamma+delta a_j)b_j=0
                                                    for j=0,1,2,3. (12)
```

The pair `(gamma,delta)` cannot vanish: otherwise one affine-linear
polynomial in `a` would vanish at four distinct points.

There are now exactly two branches.

1. If `gamma+delta a_j!=0` for all four slots, then

   ```text
   b_j=-(alpha+beta a_j)/(gamma+delta a_j).            (13)
   ```

   Thus the normal values are an honest fractional-linear image of the
   Euler eigenvalues.  When the four `b_j` are distinct, `(13)` is
   equivalently equality of the corresponding cross ratios.

2. A linear denominator can vanish at at most one exponent, say `a_k`.
   Equation `(12)` forces the numerator to vanish there too.  The
   numerator and denominator are therefore proportional, so the other
   three `b_j` are equal.  Since `b` is not constant modulo `epsilon`,
   `b_k` is different and `(8)` reduces to

   ```text
   E={x in h:x_k=0}.                                   (14)
   ```

   This is precisely a three-slot coordinate face.

Conversely every normal in `(13)`, and every degenerate normal `(14)`,
satisfies `(9)`.  This completes the global iff classification.  In the
factorial cubic problem THM-2824 already excludes `(14)` as a multipole
line; the finite computation below also retains it as a hostile boundary
point of the parameter line.

## 3. One projective line of tangent chords

Put `[r:t] in P^1` and define

```text
T(r,t)=(
 (a_2-a_1)r+(a_3-a_1)t,
-(a_2-a_0)r-(a_3-a_0)t,
 (a_1-a_0)r,
 (a_1-a_0)t).                                         (15)
```

Exact cancellation gives

```text
epsilon(T)=0,                 epsilon(AT)=0.            (16)
```

Conversely `(15)` parametrizes the two-dimensional kernel of the two
covectors `epsilon` and `epsilon o A`.  Hence every Euler-tangent chord
is

```text
E_[r:t]=span_R{T(r,t),AT(r,t)}.                        (17)
```

The points `[1:0]` and `[0:1]` are the two displayed three-slot faces
missing slots three and two.  The other two coordinate-face points occur
at the corresponding finite roots of the first two entries in `(15)`.

## 4. Factorial cubic divisibility gives two septics

Fix normalized factorial slots

```text
f_a(s)=s^a/a!,                    L(s^n)=n!,            (18)
```

and use the coefficient vectors

```text
U=T(r,t),                         V=AT(r,t).            (19)
```

Define

```text
g11=L(U^2),       g12=L(UV),       g22=L(V^2),

t111=L(U^3),      t112=L(U^2V),
t122=L(UV^2),     t222=L(V^3).                         (20)
```

The real quadratic on `(19)` is positive definite.  Comparing
coefficients in a possible factorization of its binary cubic gives the
exact iff

```text
Q|C on E_[r:t]
 iff I1(r,t)=I2(r,t)=0,                                (21)
```

where

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22,

I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2.             (22)
```

Every entry of `U,V` is linear in `(r,t)`.  The `g` entries have degree
two and the `t` entries degree three, so `I1,I2` are homogeneous binary
forms of degree seven.

A real projective line in the mean-zero three-space is a Maxwell
multipole line of the cubic harmonic exactly when its two conjugate
intersections with `Q=0` are cubic roots.  Equivalently, its restricted
positive quadratic divides its restricted cubic.  Thus `(21)--(22)`
identify Euler-tangent cubic multipole lines with common projective roots
of these two septics.

## 5. Exact avoidance through exponent thirty

The exact universe is

```text
0<=a_0<a_1<a_2<a_3<=30,
number of supports=binom(31,4)=31,465.                  (23)
```

Reduce the normalized factorial tensors modulo

```text
p=1,000,033.                                           (24)
```

Every factorial degree used in `(20)` is at most `3a_3<=90<p`, so all
normalizing factorials are units modulo `p`.  The companion constructs
`I1,I2` directly from ordered factorial tensors.  On three independent
control supports, a second constructor expands the rational polynomials
in `s`, applies `L(s^n)=n!`, and reduces the result modulo `p`; the two
constructors agree.

For every support in `(23)`, exact finite-field Euclidean division gives

```text
deg I1(1,z)=deg I2(1,z)=7,
gcd_Fp[z](I1(1,z),I2(1,z))=1.                          (25)
```

The two degree statements certify nonzero leading terms

```text
[z^7]I_nu(1,z)=I_nu(0,1)!=0,             nu=1,2,       (26)
```

so no common point is lost at `[0:1]`.  Because reduction preserves both
degrees and all denominators are units, `(25)` says the rational
resultant reduces to a nonzero element of `F_p`.  It is therefore nonzero
over `Q`.  Hence

```text
I1(r,t)=I2(r,t)=0 has no point in P^1_C                (27)
```

in any of the `31,465` cells.

Combining `(17)`, `(21)`, and `(27)` proves the finite exact conclusion:

```text
Every factorial cubic Maxwell line on a support with a_3<=30
is Euler-transverse:

Delta_A(E)!=0.                                         (28)
```

The modular gcd engine has two hostile controls.  A planted pair with
common factor `1+z+z^2` returns gcd degree two.  The same
`has_full_septic_degrees` predicate used on every census cell rejects a
planted degree-six/degree-seven pair and accepts a full-degree positive
control.  Thus neither `(25)` nor the projective-infinity gate is an
unchecked default of the implementation.

Finally, a SHA-256 stream binds, in lexicographic support order, every
support, every coefficient of both septics, and every computed gcd
degree:

```text
341b01f418f373633762c2289da11ad02c0678ef39cb18097ffcb73dec2104d7.
                                                               (28a)
```

This is a coefficient-level census binding, not a digest of only the
support labels or final ranks.

## 6. What the Euler eigenvalues still do not reveal

At a hypothetical four-slot common point

```text
Z=sum_j c_j f_(a_j),                                   (29)
```

every `c_j` is nonzero.  Otherwise the first three scalar moments would
give a forbidden three-slot common point by THM-2824.  Consequently

```text
det[Z,AZ,A^2Z,A^3Z]
 =(product_j c_j) product_(i<j)(a_j-a_i)
 !=0.                                                  (30)
```

Thus `Z` is cyclic for the four-eigenvalue Euler operator.

For fixed `m`, put

```text
x_r^(m)=L(Z^(m-1) A^r Z).                              (31)
```

The characteristic polynomial of `A` gives the fourth-order recurrence

```text
x_(r+4)-e1 x_(r+3)+e2 x_(r+2)-e3 x_(r+1)+e4 x_r=0,    (32)
```

where `e_i` are the elementary symmetric functions of the four support
exponents.  But its initial values are `x_0,x_1,x_2,x_3`; the recurrence
does not recover `x_1` from the scalar value `x_0`.

Integration by parts identifies the missing first seed:

```text
x_0^(m)=L(Z^m),

m x_1^(m)
 =L(Theta(Z^m))
 =L((s-1)Z^m).                                        (33)
```

Equation `(30)` shows that no lower-degree Euler operator identity can
collapse this four-seed observability boundary.  Therefore the support
eigenvalues organize the missing multiplier tower, but they do not turn
it into scalar moment data.

## 7. Exact companion and scope

Run

```text
python3 04-computation/gmc_euler_tangent_mobius_chord_thm2860.py
python3 -O 04-computation/gmc_euler_tangent_mobius_chord_thm2860.py
```

Both modes byte-match the stored transcript.  The companion verifies:

1. the real determinant, complex gauge law, Segre determinant, generic
   Möbius branch, degenerate coordinate-face branch, and tangent pencil;
2. independent rational and modular remainder constructors on three
   controls;
3. the planted common-factor hostile and the census's own full-degree
   predicate on positive and degree-drop controls; and
4. all `31,465` supports in `(23)`, including both septic leading terms
   and the exact modular gcd, with the full coefficient/gcd stream bound
   by `(28a)`.

The theorem is finite in `a_3`.  A separate random modular probe above
this box is not part of the theorem or its evidence.  Even inside the
box, `(28)` does **not** exclude THM-2848's branch `F^o=0` or a shared
cubic--quartic line: it only proves that every such cubic line would be
Euler-transverse.  It proves neither four-slot SFC/GMC nor access to the
multiplier in `(33)`.

**QED, pending independent hostile audit.**
