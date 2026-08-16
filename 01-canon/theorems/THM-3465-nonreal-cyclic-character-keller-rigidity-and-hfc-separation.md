---
id: THM-3465
title: "Nonreal cyclic-character Keller rigidity and all-degree HFC separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the Fourier
  coordinates of THM-3310, every real planar polynomial map whose complex
  coordinate lies in either nontrivial C3 character and whose real Jacobian
  is a nonzero constant is linear.  Consequently the entire cyclic-
  eigenvector HFC(3) lane, in every degree, is disjoint from THM-3303's
  constant-Jacobian sector: the only Keller survivors are scalar multiples
  of z or zbar, and their third normalized simplex moment is nonzero.  The
  anti-linear real-structure mate, equal component degrees, and
  chi!=chi^(-1) are load-bearing.  This proves neither HFC(3), FC(3), nor
  JC(2); mixed-character HFC cells and the cyclic lane with nonconstant
  Jacobian remain open.
source: root/factorial-jacobian-alternation/2026-08-15
audit: >
  Two independent agents rederived the top-bracket argument, equal-degree
  binary-form lemma, character contradiction, coordinate-Jacobian factor,
  and third-moment exit.  The exact companion checks the binary-form formula,
  character/star bookkeeping, the full cyclic-quartic bracket table, triangle
  normalization, and sharp independent-mate and order-two hostiles.
depends_on:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3310-degree-four-cyclic-eigenspace-on-the-triangle
related:
  - THM-3016-jacobian-pair-cross-term-rigidity-at-subleading-order
  - THM-3303-keller-simplex-null-moments-force-a-boundary-collision
  - THM-3328-boundary-cone-overlap-and-anti-tangent-keller-passport
script: 04-computation/factorial_cyclic_character_keller_rigidity_thm3465.py
output: 05-knowledge/results/factorial_cyclic_character_keller_rigidity_thm3465.out
script_sha256: e74f1ffdfe240e4a848fce72f5e997faedd9b10a9f95ac7b0c159b60b2bc2671
output_sha256: edc1dc8ec628d5ebe3af15f1ba791137a6b39e8f8a632b7441d7fd516cf981b0
semantic_sha256: 77ce350716ae0245cbb99ba4d4fdb79faf86011fab6eb2324b139e2ef85bfeee
hash_basis: raw bytes
---

# THM-3465 -- nonreal cyclic characters force a Keller map to be linear

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. The common carrier

On the standard triangle put

```text
z=s_1+omega s_2+omega^2 s_3,
w=s_1+omega^2 s_2+omega s_3,             omega^3=1, omega!=1.       (1)
```

Thus `w=conjugate(z)` on the real triangle.  Its cyclic rotation acts in the
THM-3310 convention by

```text
rho(z,w)=(omega^2 z,omega w).                                    (2)
```

For `g in C[z,w]`, let `g^dagger` conjugate every coefficient and exchange
`z,w`.  On the real plane it is the ordinary complex conjugate of `g`.  Write

```text
p=(g+g^dagger)/2,              q=(g-g^dagger)/(2i).                (3)
```

The source-to-target map in this theorem is

```text
g  |->  (g,g^dagger)  |->  {g,g^dagger}_{z,w}
   |->  F_g=(p,q).                                                   (4)
```

It keeps the HFC simplex moments and the real Keller predicate on the same
polynomial.  The anti-linear dagger in `(4)` is not optional.

## 2. Rigidity theorem

**Theorem.**  Suppose `g!=0` has one of the two nontrivial characters

```text
g o rho=chi g,                 chi in {omega,omega^2}.              (5)
```

If

```text
det D_(x,y)(p,q)=c in R\{0},                                    (6)
```

then `deg g=1`.  More precisely,

```text
chi=omega    => g=A w,
chi=omega^2  => g=A z,                  A!=0.                       (7)
```

*Proof.*  The dagger partner has character `chi^(-1)` and the same ordinary
degree `d`.  In the affine triangle coordinates used by THM-3310,

```text
Jac_(x,y)(z,w)=-3 sqrt(3)i,
Jac_(x,y)(g,g^dagger)=-2i Jac_(x,y)(p,q),
```

so

```text
Jac_(x,y)(p,q)=(3 sqrt(3)/2){g,g^dagger}_{z,w}.                     (8)
```

Hence `(6)` makes the bracket in `(8)` a nonzero constant.

Let `g_d` be the top homogeneous part.  If `d>=2`, the component of degree
`2d-2` in the bracket is

```text
{g_d,g_d^dagger}_{z,w}=0.                                         (9)
```

We use the elementary equal-degree binary-form lemma.  If nonzero homogeneous
forms `f,h` both have degree `d`, write

```text
f=w^d F(z/w),                 h=w^d H(z/w).
```

Direct differentiation gives

```text
{f,h}=d w^(2d-2)(F'H-FH').                                       (10)
```

In characteristic zero, `(10)` vanishes exactly when `(H/F)'=0`, hence when
`h=lambda f`.  Applying this to `(9)` makes `g_d^dagger` proportional to
`g_d`.  But these two nonzero forms have the distinct characters `chi^(-1)`
and `chi`, impossible because `omega!=omega^2`.  Thus `d<=1`.

The nontrivial character excludes a constant term.  The corresponding
degree-one character space is `Cw` or `Cz`, giving `(7)`.  The bracket is
`-|A|^2` in the first case and `|A|^2` in the second, so `(6)` is indeed
nonzero.  QED.

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

Every survivor in `(7)` therefore has

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
independently checks `(7)`.  The top-form proof in Section 2 explains why
the coefficient cascade persists in every degree.

## 5. Sharp failure boundaries

**Reality/dagger is load-bearing.**  If the opposite-character mate is
independent, then

```text
g=w+Bz^2,                    h=z,                    {g,h}=-1        (15)
```

for arbitrary `B`.  Unequal component degrees permit a nonlinear triangular
Keller pair.  Equation `(15)` is not a real map of the form `(3)`.

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
source:      a triangle polynomial in a pure nontrivial C3 character
target:      its real planar map (Re g,Im g)
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
sensitive checks.  It verifies `(8)--(17)` and the character/star bookkeeping
through degree 24.  The universal result is proved in Section 2, not inferred
from those finite controls.

**QED.**
