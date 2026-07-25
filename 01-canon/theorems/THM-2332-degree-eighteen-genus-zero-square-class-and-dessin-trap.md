---
id: THM-2332
title: "Degree-eighteen genus-zero square-class and Hurwitz-signature trap"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. In every
  remaining three- or four-sparse stratum of the genuine nonsplit
  polynomial exact-square-prefix degree-eighteen branch, the trigonal
  spectral curve is absolutely irreducible. Its branch discriminant is
  exactly -2^12*3^21*7^2 times 4P(y)^3+49Q(y)^2 for an explicit quartic
  P and sextic Q. Any putative rational Keller trajectory forces the
  normalization to have genus zero. The order/maximal-order discriminant
  formula and Riemann--Hurwitz then force
  4P^3+49Q^2=H S^2 over C, where H is squarefree of degree 0, 2, or 4
  and S has degree 6, 5, or 4 respectively. These are exactly the
  connected degree-three branch signatures (3,3), (3,2,2), and
  (2,2,2,2). Thus the higher-support singular cone reduces to three
  low-branch square-class/Hurwitz strata. None is proved empty here;
  JC(2) remains open.
source: codex-2026-07-25-degree18-square-class-dessin
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2328-degree-eighteen-bw-ratio-bank-closure
related:
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
  - THM-2324-degree-eighteen-bc-rational-ratio-closure
script: 04-computation/jc2_degree18_genus_zero_square_class_thm2332.py
output: 05-knowledge/results/jc2_degree18_genus_zero_square_class_thm2332.out
script_sha256: 4612b4f81ff35194a9ec4bde26e117afbc138a7011778fcb63980646e0fa5743
output_sha256: f75376c3d9ea8a1e91426da158d771d69831a80aa165e36323f5462cf07914a9
hash_basis: working-tree bytes (LF)
---

# THM-2332 -- the residual degree-eighteen locus is a low-branch square class

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2297 puts the degree-eighteen spectral equation in the weighted form

```text
G_0(u,y;B,C,D,W)=0,

wt(y,u,B,C,D,W)=(1,2,2,3,4,5).                     (1)
```

THM-2311 and THM-2314/2316/2320/2324/2328 close every support of size at
most two in `(B,C,D,W)`.  The remaining repeated-branch cone therefore
has support size three or four.

The raw repeated-branch resultant is not the right residual object.  A
putative rational trajectory would make the normalization rational, and
a rational degree-three cover has only four units of ramification.  The
right object is therefore the square class of the branch discriminant,
which forgets singular-model index squares but retains the odd inertia.

## 1. Exact depressed-cubic coordinates

Write

```text
G_0=a_3u^3+a_2u^2+a_1u+a_0
```

and translate

```text
u=v-a_2/(3a_3).
```

After division by `a_3`, the equation is

```text
v^3+p(y)v+q(y)=0,                                  (2)
```

where

```text
p=(16/964467)P,
q=(64/703096443)Q,                                 (3)

P
 =245y^4+1890By^2-24300B^2+122472D,                (4)

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                          (5)
```

Thus `P` is a structured quartic and `Q` is a sextic with zero constant
term.  Exact expansion of the cubic discriminant gives the unexpectedly
small Mordell form

```text
Delta_0(y):=Disc_u G_0

 =-2099434729254912(4P(y)^3+49Q(y)^2)

 =-2^12*3^21*7^2(4P(y)^3+49Q(y)^2).                (6)
```

The polynomial in parentheses has degree twelve and leading coefficient

```text
4*245^3+49*539^2=73060029 !=0.                     (7)
```

Equation (6), rather than `Disc_y Delta_0` by itself, is the exact object
carried forward.

## 2. Uniform irreducibility on the higher-support locus

Every remaining curve is absolutely irreducible.  This follows from the
same coefficient mechanism that independently audited THM-2324 and
THM-2328.

If the cubic in `u` were reducible over an algebraic closure of the
constant field, it would have a root in the rational function field in
`y`.  Dividing by the constant leading coefficient makes this root
integral over the integrally closed polynomial ring, so it is a
polynomial.  Weighted degree comparison bounds it by

```text
u=ay^2+by+c.                                        (8)
```

Let

```text
L(a)=1127-138915a+1607445a^2-26040609a^3.          (9)
```

The `y^6` and `y^5` coefficients after (8) are

```text
L(a),                       bL'(a).                 (10)
```

Since

```text
Disc L=-153384762202971019112448 !=0,               (11)
```

equations (10) force `b=0`.  The `y^3` coefficient is then

```text
-435456C.                                          (12)
```

If `C!=0`, this is already impossible.  If `C=0`, the `y` coefficient
becomes

```text
-5878656W.                                         (13)
```

Every support of size at least three among four coordinates contains
`C` or `W`.  Hence (12)--(13) exclude a polynomial root on all five
higher-support strata.  The normalization is connected and the map

```text
pi:C_tilde -> P^1_y
```

has degree three.

## 3. Infinity contributes no ramification

In the weighted infinity chart

```text
z=1/y,                         v=u/y^2,
```

all parameter terms vanish at `z=0`, leaving (9) in the coordinate `v`.
By (11), there are three distinct smooth points over infinity and `z` is
a local parameter at each.  Thus `pi` is unramified at infinity.

This also explains why the leading coefficient of `Delta_0` is exactly
the constant in (11): all twelve polynomial-discriminant degrees occur
at finite `y`.

## 4. A surviving rational trajectory forces genus zero

A putative retained Keller branch supplies

```text
(u,y) in C(x)^2,                  G_0(u,y)=0.         (14)
```

If `(u,y)` were constant, THM-2262's wall exclusion and first-flux formula

```text
Z=T^2=-2N_2/(5103y)
```

would make nonzero `T` and `q` constant, contradicting the genuine
nonsplit deck `q -> -q`.  Hence a survivor gives a nonconstant rational
map

```text
P^1_x -> C_tilde.                                  (15)
```

Riemann--Hurwitz applied to (15) forces `C_tilde` to have genus zero.
Applying it again to the degree-three map `pi` gives

```text
-2=3*(-2)+R_pi,

R_pi=4.                                             (16)
```

All ramification is tame.  A finite branch value has one of two inertia
types:

```text
transposition:    (2,1),     contribution 1;
three-cycle:      (3),       contribution 2.        (17)
```

Let `n_2` and `n_3` count these values.  Equations (16)--(17) give

```text
n_2+2n_3=4.                                        (18)
```

The complete possibilities are

```text
(n_2,n_3)=(0,2), (2,1), (4,0),                    (19)

branch signatures:
  (3,3),       (3,2,2),       (2,2,2,2).           (20)
```

Each occurs for a connected product-one tuple in `S_3`; the companion
enumerates the transitivity and product-one conditions rather than only
the integer partition.  The first two are two- and three-branch
Belyi-type signatures.  The last is a genuine four-point Hurwitz
signature, so “dessin” is only a mnemonic for the finite monodromy
carrier, not a claim that every case is a classical three-value Belyi
map.

## 5. The order discriminant forces the square-class trap

Divide `G_0` by its constant leading coefficient.  Then

```text
A=C[y][u]/(G_0)
```

is a free cubic `C[y]`-order in the function field of `C_tilde`.  Let
`O` be the integral closure.  Since `C[y]` is a PID, the standard
order/maximal-order discriminant formula is

```text
disc(A)=disc(O) I(y)^2                             (21)
```

up to a nonzero constant, where `I(y)` generates the index ideal.

Because infinity is unramified, the finite field discriminant is, again
up to a nonzero constant,

```text
disc(O)
 =product_(transposition values a)(y-a)
  *product_(three-cycle values b)(y-b)^2.           (22)
```

Put

```text
H(y)=product_(transposition values a)(y-a).         (23)
```

Then `H` is squarefree and

```text
deg H=n_2 in {0,2,4}.                              (24)
```

Combining (6), (21), and (22), and absorbing nonzero complex constants
into the square factor, yields the exact necessary condition

```text
4P(y)^3+49Q(y)^2=H(y)S(y)^2,                       (25)

H squarefree,

(deg H,deg S) in {(0,6),(2,5),(4,4)}.              (26)
```

Conversely, (25) alone is not sufficient for a Keller trajectory: it
records only the discriminant square class of the polynomial order.  It
does not impose the Faber equations, flux, one-form, deck, or the exact
local normality data.

## 6. Frontier effect

The degree-eighteen residual problem is no longer

```text
classify an unspecified singular weighted trigonal cone.
```

It is the three-branch algebraic classification

```text
cyclic branch:
  4P^3+49Q^2=S_6^2;                                (27a)

mixed branch:
  4P^3+49Q^2=H_2 S_5^2,       H_2 squarefree;      (27b)

four-simple branch:
  4P^3+49Q^2=H_4 S_4^2,       H_4 squarefree.      (27c)
```

Here `P,Q` are not arbitrary: they have exactly the coefficient pattern
(4)--(5), and at least three of `(B,C,D,W)` are nonzero.  Coefficient
comparison in (27a)--(27c), with the Keller one-form and Faber sidecars
restored only after this square-class sieve, is the next exact target.

The intrinsic finite combinatorial carrier is a product-one transitive
tuple in `S_3`, with vertices the three sheets and labels the local
inertia cycles.  It is not a tournament: transpositions are involutive
ties, three-cycles are ternary rotations, and imposing a head-to-head
orientation would destroy precisely the parity information retained by
the discriminant square class.

No three- or four-sparse stratum is proved empty here.  Split/even
descent, other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_genus_zero_square_class_thm2332.py
python3 -O 04-computation/jc2_degree18_genus_zero_square_class_thm2332.py
```

Both transcripts are byte-identical to the stored output.  The companion
verifies the depressed coefficients (3)--(5), the full identity (6), the
degree and leading term (7), the irreducibility coefficients
(9)--(13), all five higher supports, separable infinity, and every
transitive product-one `S_3` signature in (19)--(20).  The
order-discriminant formula, the two Riemann--Hurwitz uses, and the Keller
implication are the mathematical proof above, not computer assumptions.
QED.
