---
id: THM-3465
title: "Nonreal cyclic-character Keller rigidity and all-degree HFC separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For any finite
  orientation-preserving Euclidean rotation of the real plane, a real
  polynomial map whose complex coordinate
  has a nonreal rotation character and whose real Jacobian is a nonzero
  constant is linear; its character must already occur in the linear source
  representation.  Consequently, in the Fourier coordinates of THM-3310,
  the entire cyclic-
  eigenvector HFC(3) lane, in every degree, is disjoint from THM-3303's
  constant-Jacobian sector: the only Keller survivors are scalar multiples
  of z or zbar, and their third normalized simplex moment is nonzero.  The
  anti-linear real-structure mate, equal component degrees, and
  eta!=eta^(-1) is load-bearing.  This proves neither HFC(3), FC(3), nor
  JC(2); mixed-character HFC cells and the cyclic lane with nonconstant
  Jacobian remain open.
source: root/factorial-jacobian-alternation/2026-08-15
audit: >
  Two independent agents rederived the top-bracket argument, equal-degree
  binary-form lemma, character contradiction, coordinate-Jacobian factor,
  and third-moment exit.  The exact companion checks the binary-form formula,
  finite-rotation character/star bookkeeping, the full cyclic-quartic bracket
  table, triangle normalization, and sharp independent-mate and order-two
  hostiles.
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
related:
  - THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow
  - THM-3016-jacobian-pair-cross-term-rigidity-at-subleading-order
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3328-boundary-cone-overlap-and-anti-tangent-keller-passport
script: 04-computation/factorial_cyclic_character_keller_rigidity_thm3465.py
output: 05-knowledge/results/factorial_cyclic_character_keller_rigidity_thm3465.out
script_sha256: d072afb9ee508ae874fdf253084f17b206f55a030bd009382b1005c39fabe8ca
output_sha256: 3ffba457df22503a986f20d188c3e6b68491dd643fe448a55aee13357c184a96
semantic_sha256: 91d4afcb83b4fa6e0ca54da2946b6afd6db335888ffc9199c1351f2080811db6
hash_basis: raw bytes
---

# THM-3465 -- nonreal cyclic characters force a Keller map to be linear

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. The common carrier

Let `rho` be a nontrivial finite orientation-preserving Euclidean rotation of
a real affine plane.  Recenter at its fixed point and choose conjugate complex
linear coordinates `z,w` so that

```text
rho(z,w)=(xi z,xi^(-1)w),                         |xi|=1.           (0)
```

For `g in C[z,w]`, let `g^dagger` conjugate every coefficient and exchange
`z,w`.  On the real plane it is the ordinary complex conjugate of `g`.  Write

```text
p=(g+g^dagger)/2,              q=(g-g^dagger)/(2i).               (0a)
```

The source-to-target map in this theorem is

```text
g  |->  (g,g^dagger)  |->  {g,g^dagger}_{z,w}
   |->  F_g=(p,q).                                                  (0b)
```

It retains a finite rotation character and the real Keller predicate on the
same polynomial.  The anti-linear dagger in `(0b)` is not optional.

For the HFC specialization, on the standard triangle put

```text
z=s_1+omega s_2+omega^2 s_3,
w=s_1+omega^2 s_2+omega s_3,             omega^3=1, omega!=1.       (1)
```

Thus `w=conjugate(z)` on the real triangle.  Its cyclic rotation acts in the
THM-3310 convention by

```text
rho(z,w)=(omega^2 z,omega w).                                    (2)
```

## 2. Rigidity theorem

**Theorem.**  In the finite-rotation setup `(0)`, suppose

```text
g o rho=eta g,                 eta!=conjugate(eta).                 (3)
```

If

```text
det D_(x,y)(p,q)=c in R\{0},                                     (4)
```

then `deg g=1`, `eta` is one of `xi,xi^(-1)`, and `g` is a nonzero
multiple of the corresponding linear coordinate.  In the triangle convention
`(2)`, more precisely,

```text
eta=omega    => g=A w,
eta=omega^2  => g=A z,                  A!=0.                       (5)
```

*Proof.*  The dagger partner has character `eta^(-1)=conjugate(eta)` and the
same ordinary degree `d`.  Put `K=Jac_(x,y)(z,w)`, a fixed nonzero purely
imaginary constant.  Directly from `(0a)`,

```text
Jac_(x,y)(g,g^dagger)=-2i Jac_(x,y)(p,q),
Jac_(x,y)(g,g^dagger)=K {g,g^dagger}_{z,w},
```

so

```text
Jac_(x,y)(p,q)=(iK/2){g,g^dagger}_{z,w}.                           (6)
```

Hence `(4)` makes the bracket in `(6)` a nonzero constant.  On the triangle,
`K=-3 sqrt(3)i`, recovering the factor `3 sqrt(3)/2`.

Let `g_d` be the top homogeneous part.  If `d>=2`, the component of degree
`2d-2` in the bracket is

```text
{g_d,g_d^dagger}_{z,w}=0.                                         (7)
```

We use the elementary equal-degree binary-form lemma.  If nonzero homogeneous
forms `f,h` both have degree `d`, write

```text
f=w^d F(z/w),                 h=w^d H(z/w).
```

Direct differentiation gives

```text
{f,h}=d w^(2d-2)(F'H-FH').                                        (8)
```

In characteristic zero, `(8)` vanishes exactly when `(H/F)'=0`, hence when
`h=lambda f`.  Applying this to `(7)` makes `g_d^dagger` proportional to
`g_d`.  But these two nonzero forms have the distinct characters `eta^(-1)`
and `eta`, impossible by `(3)`.  Thus `d<=1`.

The nonreal character excludes a constant term.  The only degree-one source
characters are `xi,xi^(-1)`, proving the general assertion.  In `(2)` their
spaces are `Cz,Cw`, giving `(5)`.  The bracket is `-|A|^2` in the first
displayed case and `|A|^2` in the second, so `(4)` is indeed nonzero.  QED.

The mechanism is the leading-form dependence already native to planar
Jacobian pairs, but the finite character plus real-structure sidecar turns
dependence into a contradiction.  A finite `C3` eigenspace is not a
`C^*`-graded space: it mixes weights such as `w` and `z^2`.  Thus this result
is not an application of positive weighted-homogeneous Keller rigidity.

## 3. All-degree HFC consequence

THM-3300 reduces the cyclic-eigenvector HFC(3) obligation to moments whose
orders are multiples of three.  With normalized uniform triangle measure,
THM-3310 gives

```text
<z^3>=<w^3>=1/10.                                                 (11)
```

Every triangle survivor in `(5)` therefore has

```text
<g^3>=A^3/10!=0.                                                  (12)
```

Consequently:

> No cyclic eigenvector on the triangle, in any polynomial degree, can both
> have all positive simplex moments zero and have constant nonzero real
> Jacobian.

This closes the intersection of two large search lanes, not either lane
itself.  In particular:

- THM-3310/3321/3323's full-support-five cyclic quartic remains open as an
  HFC moment cell, but every point in it has zero or nonconstant real
  Jacobian;
- THM-3303's boundary-collision implication remains useful for
  mixed-character candidates, but its hypothesis is empty on every pure
  cyclic character cell; and
- no statement is made about nonhomogeneous `FC(3)` or arbitrary planar
  Keller pairs.

## 4. Quartic control and the shorter proof it reveals

For the current degree-four chart

```text
g=A w+B z^2+C z w^2+D z^3w+E w^4,                                  (13)
```

the exact companion expands `{g,g^dagger}` into ten monomials.  Its top
coefficients include

```text
[w^6]=-4E conjugate(D),
[z^3w^3]=8(|D|^2-2|E|^2),
[z^6]=-4D conjugate(E).                                            (14)
```

Constancy forces `D=E=0`; then `[z^2w^2]=-3|C|^2` forces `C=0`, and
`[zw]=4|B|^2` forces `B=0`.  This directly closes the quartic slice and
  independently checks `(5)`.  The top-form proof in Section 2 explains why
the coefficient cascade persists in every degree.

## 5. Sharp failure boundaries

**Reality/dagger is load-bearing.**  If the opposite-character mate is
independent, then

```text
g=w+Bz^2,                    h=z,                    {g,h}=-1        (15)
```

for arbitrary `B`.  Unequal component degrees permit a nonlinear triangular
Keller pair.  Equation `(15)` is not a real map of the form `(0a)`.

**The nonreal character is load-bearing.**  For order two,
`chi=chi^(-1)`.  The real odd shear

```text
(x,y) |-> (x,y+x^3)                                             (16)
```

is nonlinear, equivariant under `(x,y)|->(-x,-y)`, and has Jacobian one.

**A two-step operation kernel is too short.**  With the univariate factorial
functional, the companion verifies

```text
f=x^2+(-4+2i)x+(2-2i),
L(f)=L(f^2)=0,                 L(f^3)=32+80i,
L(f f^dagger)=8.                                                  (17)
```

Thus `ker L intersect M_f^(-1)(ker L)` is only a depth-two address.  Neither
the full FC orbit nor the dagger energy follows from it.

## 6. Information contract and next live cell

```text
source:      a real-plane polynomial in a nonreal finite-rotation character
target:      its real planar map (Re g,Im g); then the triangle HFC cell
operation:   Hermitian Poisson bracket with g^dagger
preserved:   degree, character, simplex moments, constant-real-J predicate
destroyed:   mixed-character top layers, boundary ancestry, inverse effectivity
sidecars:    dagger, source area form/orientation, character label
hostiles:    independent mate (15), order-two shear (16)
```

The next joint experiment should not revisit the now-empty pure character
cells.  Decompose a general top layer into

```text
V_1 direct_sum V_omega direct_sum V_(omega^2).                       (18)
```

Keller forces the two equal-degree top coordinates to be proportional, so
after a target phase the live top form is dagger-real and necessarily mixes
paired characters.  The cheapest bounded probe is degree at most four:
alternate the next bracket layer with simplex moments through six and the
multiplier bank `{1,z,w,z^2,zw,w^2}`.  Any modular moment reduction must obey
THM-3310's denominator guard.

## 7. Exact verification

Run

```text
python3 04-computation/factorial_cyclic_character_keller_rigidity_thm3465.py
python3 -O 04-computation/factorial_cyclic_character_keller_rigidity_thm3465.py
```

and compare raw bytes with the declared output.  The companion uses exact
symbolic and rational arithmetic, explicit failure gates, and no assertion-
  sensitive checks.  It verifies `(8)`, the triangle controls `(11)--(17)`,
  and finite-rotation character/star bookkeeping through degree 24 for orders
  three through sixteen.  The universal result is proved in Section 2, not
  inferred from those finite controls.

**QED.**
