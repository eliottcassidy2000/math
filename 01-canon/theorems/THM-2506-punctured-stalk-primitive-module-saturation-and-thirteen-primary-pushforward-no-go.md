---
id: THM-2506
title: "Punctured-stalk primitive-module saturation and thirteen-primary pushforward no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every essential
  THM-2436 punctured 91-stalk parent, all 72 primitive C_13 x C_7
  Fourier coefficients of the integral defect are nonzero, not merely
  one. Thus every preassigned primitive mode has nonvanishing locus
  exactly the essential locus, improving the fixed-mode parent-mass
  floors from 1/252 and 1/168 to 2/7 and 3/7. The primitive group-algebra
  summand is the field Q(zeta_91); convolution by any THM-2474
  squarefree collision kernel is an automorphism there, so abstract
  primitive spectral cancellation is impossible. Conversely, every
  affine homomorphic pushforward from C_13 x C_7 to a finite 13-group
  kills the defect by the row-sum law. Hence complete colour saturation
  still does not produce a THM-2471 owner/deep-root current. The packet
  belongs to the already empty high-septimal branch, no live row is
  excluded, and LRC(14) remains open.
source: codex-2026-07-27-punctured-stalk-primitive-module
depends_on:
  - THM-2435-top-blocker-essential-parent-and-punctured-stalk-carrier
  - THM-2436-punctured-ninety-one-stalk-repeated-step-spectrum
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
related:
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2474-squarefree-first-collision-primitive-character-saturation
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
script: 04-computation/lrc14_punctured_stalk_primitive_module_thm2506.py
output: 05-knowledge/results/lrc14_punctured_stalk_primitive_module_thm2506.out
script_sha256: 98324f21d835dca9154d5c5ab9090229282627f743bc64b054d1dd11e494afc2
output_sha256: a82eb92c151654e6b4680158ea1597e88249ee1e61e4a15e5f9bd62698065e01
hash_basis: working-tree bytes (LF)
---

# THM-2506 -- the punctured stalk fires every primitive colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2436 proves that every nonflat punctured `91`-stalk defect has *some*
mixed `F_13 x F_7` Fourier coefficient.  It then uses a union bound over the
`72` possible coefficients to select one fixed mode.  Rationality makes that
loss unnecessary.  The entire mixed block is one irreducible rational module:

```text
one nonzero mixed coefficient
  -> one nonzero element of Q(zeta_91)
  -> all 72 Galois embeddings are nonzero.                    (1)
```

This is the same rational-module mechanism that drives THM-2449 and
THM-2474.  It removes the remaining abstract colour-selection and spectral
cancellation questions on this stalk.  It also reveals a sharp categorical
boundary: the `C_7` factor carrying the defect is annihilated by every affine
homomorphic pushforward to a `13`-primary root.

## 1. Pointwise primitive saturation

Retain all notation and hypotheses of THM-2436.  Thus

```text
P=T_(7^(M+1))(A)                                             (2)
```

is the quotient-parent image, `mathcal E` is its flat six-unit tiling locus,
and for almost every `Y in P`,

```text
d_Y:F_13 x F_7 -> Z,

sum_(r in F_7)d_Y(h,r)=0                 for every h.         (3)
```

For primitive roots `zeta_13,zeta_7`, use THM-2436's unnormalized transform

```text
dtilde_Y(alpha,beta)
 =sum_(h,r)d_Y(h,r)
   zeta_13^(-alpha h)zeta_7^(-beta r).                       (4)
```

Then, almost everywhere, exactly one alternative holds:

```text
Y in mathcal E:
  d_Y=0 and every mixed transform vanishes;

Y in P minus mathcal E:
  dtilde_Y(alpha,beta)!=0
  for every alpha in F_13^*, beta in F_7^*.                  (5)
```

### Proof by the rational group algebra

CRT and the cyclotomic decomposition give

```text
Q[C_13 x C_7]
  isomorphic to
Q direct_sum Q(zeta_13) direct_sum Q(zeta_7)
  direct_sum Q(zeta_91).                                    (6)
```

The row-sum law (3) kills the trivial and pure-`13` summands.  Therefore

```text
d_Y lies in Q(zeta_7) direct_sum Q(zeta_91).                 (7)
```

The pure-`7` summand consists exactly of functions independent of `h`.
THM-2436 Section 3 proves, uniformly for one source, two distinct sources,
and coincident labelled sources, that

```text
d_Y is h-independent             iff d_Y=0.                  (8)
```

It also identifies `d_Y=0` with the flat locus and proves `d_Y!=0` on
`P minus mathcal E`.  Thus every essential defect has a nonzero projection to
the last field in (6).

That projection is one nonzero element of `Q(zeta_91)`.  Its `72` embeddings
are precisely the values in (4) with

```text
alpha!=0,                         beta!=0.                    (9)
```

A field embedding cannot send a nonzero element to zero, proving (5). QED.

Equivalently, THM-2436 Section 4 supplies one nonzero mixed coefficient and
THM-2449's coprime Galois all-or-all law supplies the other `71`.  The module
proof explains why these are the same statement.

## 2. Fixed-mode locus and quantitative invoices

Fix *before seeing `Y`* any one primitive pair `(alpha,beta)`.  Equation (5)
gives the exact nonvanishing locus

```text
{Y in P:dtilde_Y(alpha,beta)!=0}
   =P minus mathcal E                         almost everywhere. (10)
```

THM-2435's essential-parent invoice is

```text
mu(P minus mathcal E)>=(1+k)/7.                          (11)
```

Therefore every preassigned primitive mode has the sharp inherited floors

```text
k=1:       nonvanishing parent mass >=2/7,

k=2:       nonvanishing parent mass >=3/7,                    (12)
```

with uniform floor `2/7`.  This replaces the union-bound floors

```text
(2/7)/72=1/252,

(3/7)/72=1/168.                                               (13)
```

There is also a tiny but exact coefficient floor.  THM-2436's complete atlas
has at most nine holes.  Positive holes and negative duplicates balance, so

```text
||d_Y||_1<=18.                                                (14)
```

For one primitive value `theta=dtilde_Y(alpha,beta)`, (5) makes `theta` a
nonzero algebraic integer in the degree-`72` field `Q(zeta_91)`.  Its field
norm is a nonzero integer, while every conjugate is bounded by (14).  Hence

```text
|dtilde_Y(alpha,beta)|>=18^(-71)             on P minus E,    (15)
```

and the fixed-mode unnormalized `L^2` energy obeys

```text
integral_P |dtilde_Y(alpha,beta)|^2 dY
 >=((1+k)/7)18^(-142).                                       (16)
```

For the normalized Fourier convention, divide (15) by `91` and (16) by
`91^2`.  These floors are structural rather than numerically competitive.

## 3. Squarefree collision convolution is an automorphism

Let

```text
G=C_13 x C_7,

K_prim=e_prim Q[G] isomorphic to Q(zeta_91).                  (17)
```

If `c in Q[G]` is any THM-2474 squarefree collision correlation on `G`, then
every primitive transform of `c` is nonzero.  Thus its primitive projection

```text
c_prim=e_prim c                                               (18)
```

is a nonzero element of the field `K_prim`.  Convolution by `c` restricts on
`K_prim` to multiplication by `c_prim`, and is therefore an automorphism:

```text
x_prim!=0                  implies (c*x)_prim!=0.             (19)
```

More strongly, each of the `72` primitive transforms of `c*x` is the product
of two nonzero transforms.  Applying this to an essential THM-2436 defect
shows that squarefree collision convolution cannot erase any of its primitive
colours.  THM-2474's normalized currents insert the common factor `91^(-2)`;
that scalar does not affect invertibility.

Equation (19) closes the **abstract spectral-cancellation** question.  It is
not yet a physical composition theorem: no current result realizes this
group-algebra convolution on one owner-supported source/arrival/deep fibre.

## 4. Every affine thirteen-primary pushforward vanishes

Let `V` be a finite abelian `13`-group and let

```text
L:C_13 x C_7 -> V                                             (20)
```

be affine with homomorphic linear part.  The restriction of that linear part
to `C_7` is zero, because `V` has no nonidentity element of order dividing
`7`.  Hence `L(h,r)` is independent of `r`.

Define the signed pushforward

```text
(L_*d_Y)(v)=sum_(L(h,r)=v)d_Y(h,r).                            (21)
```

For each fixed `h`, all seven terms lie in the same fibre of `L`, and (3)
gives their sum as zero.  Therefore

```text
L_*d_Y=0                         for every such L.             (22)
```

This is an exact no-go, not a dimension heuristic.  The punctured-stalk
defect can have every one of its `72` primitive colours nonzero while every
affine pushforward to a `13`-primary root vanishes.  A bridge to THM-2471 must
retain a non-affine physical sidecar or a common `C_7` coordinate; matching
cardinalities or discarding the septimal factor cannot work.

## 5. A two-row hostile with all colours and zero pushforward

The boundary already occurs in one explicit THM-2436 cover.  Normalize the
guard to

```text
G={0,...,25},                                                  (23)
```

take one blocker source column `B_5`, and take five step-one APs of length
thirteen with centres

```text
84, 71, 33, 47, 58.                                           (24)
```

The guard, APs, and blocker cover all of `Z/91Z`.  Before adding the blocker,
the only nonzero defect rows are

```text
d(0,-)=e_5-e_3,

d(1,-)=e_5-e_4.                                               (25)
```

Put

```text
omega=zeta_13^(-alpha),             xi=zeta_7^(-beta).
```

Direct substitution in (4) factors every primitive value as

```text
dtilde(alpha,beta)
 =xi^3(xi-1)(xi+1+omega xi).                                  (26)
```

The first two factors are nonzero.  If the last vanished, then
`|1+xi|=1`; equivalently `xi+conjugate(xi)=-1`, which would make `xi` a
third root of unity.  This is impossible for primitive seventh-root `xi`.
Thus all `72` values in (26) are nonzero.

Nevertheless (22) kills every affine `13`-primary pushforward.  This is the
minimal mechanism behind the boundary:

```text
primitive 91-colour saturation
  does not imply
physical 13-root current.                                    (27)
```

## 6. Anti-Kakeya and temporal scope

The rational module is the rigorous content of the anti-Kakeya reframe.
The static coordinates are

```text
h = quotient-stalk row,              r = blocker-source column. (28)
```

They are not THM-2461/2471's prescribed-time target phase and source atom.
Three further coordinates remain absent:

1. the affine normalization carry `kappa(Y)`, needed for coherent absolute
   Fourier phase;
2. the ancestry sheet `b in Z/13^K Z`, distinguishing source refinement from
   arrival refinement; and
3. the residue `a mod d_0` needed for uniform descent to the deep root.

The carry cancels from centre differences, which is why THM-2436's rounding
amplifier works, but it does not cancel from an integrated endpoint current.
Likewise, an `F_7`-only Dvir argument would merge distinct `91`-directions:
steps `1` and `8` already have the same slope modulo seven.  Any polynomial
replacement for the atlas must retain the full CRT/carry object.

Finally, THM-2436 proves the high-septimal branch containing this packet is
empty.  Reusing (5) on the live `165` low-septimal rows requires a genuine
transplant theorem.  The module upgrade does not silently move the packet
between branches.

## 7. Exact gain and stopping boundary

THM-2506 proves

```text
essential punctured-stalk defect
  -> all 72 primitive colours pointwise
  -> any fixed colour on the whole essential locus
  -> squarefree collision convolution cannot cancel it.       (29)
```

It removes the factor-`72` fixed-mode loss in THM-2436 and closes abstract
primitive spectral cancellation with THM-2474.  It does not physicalize that
convolution, transport absolute phase, recover the discarded `C_7` coordinate
after a `13`-primary pushforward, couple source and arrival atoms, descend to
the deep sheet, exclude a live scalar row, or prove LRC(14).

The next proof object is therefore a same-parent ancestry cospan carrying the
primitive module together with `kappa(Y),b,a mod d_0`, or an independent proof
that the owner loop has nonzero THM-2365 drift.  More colour or a larger scalar
energy estimate is no longer the missing step.

## 8. Exact companion

Run

```text
python3 04-computation/lrc14_punctured_stalk_primitive_module_thm2506.py
python3 -O 04-computation/lrc14_punctured_stalk_primitive_module_thm2506.py
```

The dependency-free exact companion:

- reconstructs the guard, five centred APs, blocker cover, and two-row defect
  in (23)--(25);
- constructs `Phi_91` and reduces all `72` primitive transforms exactly;
- verifies the factorization (26) coefficientwise in `Q(zeta_91)`;
- convolves with the `91`-unit mask and verifies its primitive multiplier
  `mu(91)=1` on all colours;
- checks the old and new fixed-mode floors exactly; and
- exhausts all `2,366` affine maps to `Z/13Z` and `Z/169Z`, finding zero
  nonzero pushforwards.

Normal and optimized runs must reproduce

```text
05-knowledge/results/lrc14_punctured_stalk_primitive_module_thm2506.out
```

byte-for-byte.

An independent audit rederived the rational-module decomposition, pointwise
all-colour law, locus floors, algebraic-norm invoice, collision-convolution
automorphism, affine-pushforward no-go, hostile factorization, and temporal
scope before promotion.

QED.
