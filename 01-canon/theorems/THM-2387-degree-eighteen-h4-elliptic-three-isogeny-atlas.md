---
id: THM-2387
title: "Degree-eighteen H4 elliptic three-isogeny atlas"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Every
  genuine coprime H_4 survivor has a smooth elliptic quadratic
  resolvent E:t^2=H_4 and a connected unramified Cardano three-cover
  X->E. The divisor of A=7Q+S_4t is three times a degree-zero divisor
  selecting one lift over each P-root; its Abel-Jacobi class is a
  nonzero point of Pic^0(E)[3], and complementary selections negate the
  class. The compatible selections form an explicit quotient atlas,
  not a canonical 16/2 bijection. Tate normal form, the Velu/Lattes
  quotient, both j-maps, and the binary-quartic degree-four sector
  polynomial are exact. Over C the four unoriented/eight oriented
  torsion sectors always exist, so this is structure rather than an H_4
  obstruction. The structured P,Q,S_4 trace, H_4 emptiness, degree
  eighteen, JC(2), and DC(2) remain open.
source: codex-2026-07-26-h4-elliptic-three-isogeny-atlas
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2386-degree-eighteen-h4-common-root-elimination
related:
  - THM-2341-degree-eighteen-deep-wall-local-genus-split
  - THM-2373-degree-eighteen-rational-charged-section-atlas
script: 04-computation/jc2_degree18_h4_elliptic_three_isogeny_atlas_thm2387.py
output: 05-knowledge/results/jc2_degree18_h4_elliptic_three_isogeny_atlas_thm2387.out
script_sha256: bbee4041003db38d24580eb3d139d52ef585a883d329f9ab3da134030052b631
output_sha256: 62b844a08e4245d486a30d477934a60636be90fd917e5b563213658ce6d79b04
hash_basis: working-tree bytes (LF)
---

# THM-2387 -- the H4 survivor is an elliptic three-isogeny sector

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

Retain THM-2332's structured covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y,                          (1)

F=4P^3+49Q^2.
```

Let a **genuine coprime H4 survivor** mean a genuine genus-zero
four-transposition survivor in the sense of THM-2332 and THM-2386:

```text
F=H4*S4^2,

H4 squarefree,                    deg(H4,S4)=(4,4),

Res_y(P,Q)!=0.                                    (2)
```

THM-2386 proves the last condition for every genuine H4 survivor. This
theorem identifies the elliptic and modular structure forced by (2). It
does not eliminate that structure.

## 1. The quadratic resolvent is elliptic

Over `C`, let

```text
E: t^2=H4(y).                                      (3)
```

The squarefree quartic gives a smooth genus-one curve with two points
`infinity_+`, `infinity_-` over infinity. Define the conjugate Cardano
factors

```text
A=7Q+S4*t,                  Abar=7Q-S4*t.           (4)
```

Equations (1)--(3) give the exact norm identity

```text
A*Abar
 =49Q^2-H4*S4^2
 =-4P^3.                                           (5)
```

Adjoin a cube root

```text
z^3=-4A                                             (6)
```

to `C(E)`. Put

```text
w=-4P/z.
```

Then (5)--(6) give

```text
w^3=-4Abar,               zw=-4P.
```

For `u=z+w`,

```text
u^3
 =z^3+w^3+3zwu
 =-56Q-12Pu.                                      (7)
```

THM-2332's depressed cubic is

```text
v^3+p(y)v+q(y)=0,

p=(16/964467)P,            q=(64/703096443)Q.       (8)
```

The exact reconstruction

```text
v=(2/1701)(z-4P/z)                                 (9)
```

turns (7) into (8), since

```text
(16/964467)(2/1701)=12(2/1701)^3,

64/703096443=56(2/1701)^3.                         (10)
```

Multiplying `z` by the three cube roots of unity gives the three roots
of (8). Thus the normalization `X` of (6) is the cyclic degree-three
base change of the connected trigonal cover to its quadratic
discriminant resolvent.

## 2. The Cardano divisor is nonzero three-torsion

Let `alpha` be a root of `P` of multiplicity `m`. Coprimality in (2)
gives `Q(alpha)!=0`. Hence the two points of `E` above `alpha` are
distinct, and

```text
A+Abar=14Q
```

is a unit at both. Equation (5) therefore puts a zero of order `3m` in
exactly one of `A,Abar`; the other factor is a unit. Let
`R_alpha,epsilon` denote the lift selected by the zero of `A`, and set

```text
Z_epsilon=sum_(P(alpha)=0) m_alpha R_alpha,epsilon.
```

This is a degree-four divisor, with multiplicity retained if `P` is not
squarefree.

At each point over infinity, `P^3` has pole order twelve. Each summand
of `A` and `Abar` has pole order at most six. Their nonzero product in
(5) forces each factor to have pole order exactly six at both infinity
points. Consequently

```text
div_E(A)=3D_epsilon,

D_epsilon
 =Z_epsilon-2*infinity_+-2*infinity_-.             (11)
```

In particular `deg(D_epsilon)=0`, and

```text
T_epsilon:=AJ(D_epsilon)
 belongs to Pic^0(E)[3].                           (12)
```

This class is nonzero. Indeed, if `D_epsilon=div(g)`, then (11) gives
`A=cg^3` for some `c!=0`. Constants have cube roots over `C`, so the
Kummer algebra (6) would split. But THM-2332 gives an
absolutely irreducible cubic with nonsquare discriminant square class
`H4`; after adjoining the discriminant square root, its cyclic cubic
Galois closure remains connected. Therefore

```text
T_epsilon in Pic^0(E)[3] minus {0}.                (13)
```

Since every valuation of `A` is divisible by three, (6) is unramified.
Connectedness makes it an etale isogeny

```text
pi:X -> E,                         deg(pi)=3.        (14)
```

Riemann--Hurwitz gives `g(X)=1`.

The complementary lift selection `bar(epsilon)` is selected by
`Abar`. The divisor of `P` on `E` is

```text
div_E(P)
 =Z_epsilon+Z_bar(epsilon)
  -4*infinity_+-4*infinity_-.
```

Hence

```text
D_epsilon+D_bar(epsilon)=div_E(P),

T_bar(epsilon)=-T_epsilon.                         (15)
```

## 3. The correct selection atlas

For a fixed pair `(E,P)`, let `Lift(P/E)` be the set of choices of one
of the two points over each distinct root of `P`; its multiplicity is
then retained in `Z_epsilon`. Form `D_epsilon` as in (11). The
compatible sector atlas is

```text
A(P,E)
 ={
    (epsilon,T):
    epsilon in Lift(P/E),
    AJ(D_epsilon)=T in Pic^0(E)[3] minus {0}
   }
   /
   ((epsilon,T) ~ (bar(epsilon),-T)).               (16)
```

Equation (16) is the exact object carried forward. It is not the false
claim

```text
"the 16 lift selections modulo complement
 are canonically the eight nonzero 3-torsion points."                 (17)
```

Even when `P` has four distinct roots, the Abel--Jacobi map from its
sixteen selections can collide, and a general selection need not be
three-torsion. When `P` has fewer than four distinct roots, there are
fewer than sixteen rootwise selections. The actual norm factor `A`
supplies one complementary compatible pair; no cardinality statement
about (16) is needed.

Over `C`,

```text
Pic^0(E)[3] is isomorphic to (Z/3Z)^2.
```

It has eight nonzero **oriented** points and four cyclic subgroups

```text
(Pic^0(E)[3] minus {0})/{+1,-1}
 is P^1(F3),                         cardinality 4. (18)
```

Thus nonzero three-torsion is automatic on every complex elliptic
curve. The content of (11)--(16) is how the structured Cardano divisor
lands in that torsion, not the existence of torsion itself.

## 4. Tate normal form and the incoming isogeny

The subgroup `<T_epsilon>` is the kernel of the dual isogeny `E->X`.
Equivalently, (14) is an incoming cyclic isogeny to `E`. To parameterize
it, choose a generator `T_X` of `ker(pi)` on the source. After choosing
origins and coordinates, the pointed source has Tate normal form

```text
X:
  Y^2+aXY+bY=X^3,

T_X=(0,0),

b(a^3-27b)!=0.                                    (19)
```

The line `Y=0` is tangent at `T_X` and has intersection multiplicity
three, so `T_X` has exact order three. The discriminant is

```text
Delta_X=b^3(a^3-27b).                              (20)
```

The distinction between `T_epsilon` on the target and `T_X` on the
source is duality data. Identifying their orientations requires a deck
character, equivalently a choice of primitive cube root of unity; it is
not a canonical equality of points.

Velu's quotient by `<T_X>` is

```text
E'=X/<T_X>:

Y^2+aXY+bY
 =X^3-5abX-b(a^3+7b),                             (21)
```

and `E'` is isomorphic to the resolvent `E` in (3). If `(x,y)` is a
point of the source, translation by `T_X` and `-T_X` gives

```text
x((x,y)+T_X)    =-by/x^2,

x((x,y)-T_X)    =b(y+b+ax)/x^2.
```

Therefore the quotient's `x`-coordinate is the degree-three
Velu/Lattes map

```text
f(x)
 =x+x(P+T_X)+x(P-T_X)
 =(x^3+abx+b^2)/x^2.                              (22)
```

Put

```text
lambda=a^3/b.                                      (23)
```

This is invariant under Weierstrass scaling and forgets the sign of the
kernel generator. Direct invariant calculation gives

```text
j_X=lambda(lambda-24)^3/(lambda-27),

j_E=lambda(lambda+216)^3/(lambda-27)^3.            (24)
```

The second formula is the one relevant to the quartic resolvent (3);
the first belongs to the incoming Cardano cover.

## 5. The binary-quartic sector polynomial

Homogenize the resolvent quartic as

```text
H4(U,V)
 =h0 U^4+h1 U^3V+h2 U^2V^2+h3 UV^3+h4 V^4.
```

Use the binary-quartic invariants

```text
I=12h0h4-3h1h3+h2^2,

J=72h0h2h4+9h1h2h3
  -27h0h3^2-27h1^2h4-2h2^3.                       (25)
```

With this convention,

```text
27 Disc(H4)=4I^3-J^2,

j(E)=6912I^3/(4I^3-J^2).                           (26)
```

Equating (26) with the target formula in (24) gives the exact sector
polynomial

```text
M_H4(lambda)
 =6912I^3(lambda-27)^3
  -(4I^3-J^2)lambda(lambda+216)^3.                 (27)
```

Smoothness says `4I^3-J^2!=0`. Thus

```text
deg_lambda(M_H4)=4,

M_H4(27)
 =-(4I^3-J^2)*27*243^3 !=0.                       (28)
```

The four roots counted with modular multiplicity are the four
unoriented incoming three-isogeny sectors. At elliptic curves with extra
automorphisms, distinct subgroups may give the same `lambda` value, so
(27) is not a claim that all four roots are simple or canonically
labeled. Choosing a generator over each projective line restores the
eight oriented sectors in (18).

Equation (27) is automatic rather than obstructive over `C`: every
elliptic curve has the four lines in (18), and every degree-four
polynomial has four roots with multiplicity. Its value is that it turns
the Cardano divisor into an exact modular coordinate which can be
combined with the special coefficient patterns in (1).

What remains open is the structured trace:

```text
lambda-sector
 + selected P-lift divisor
 + A=7Q+S4*t
 + the sparse degree-eighteen P,Q patterns
 + the Keller/Faber/one-form sidecars.              (29)
```

Matching only `j(E)` in (27) forgets every item after the first.

## 6. Weil-pairing tournament and its loss boundary

There is an intrinsic binary relation here, but only after retaining its
gauge. The vertices are the four lines

```text
L in P^1(F3)=projective lines of E[3].              (30)
```

Choose a generator `t_L` of each line and a primitive cube root
`zeta`. For distinct lines, nondegeneracy of the Weil pairing gives

```text
e_3(t_L,t_M) in {zeta,zeta^(-1)}
```

with no ties. Orient

```text
L -> M  iff  e_3(t_L,t_M)=zeta.                    (31)
```

Negating `t_L` reverses exactly the three arcs incident to `L`: it is a
vertex switch. Replacing `zeta` by `zeta^(-1)` reverses every arc. In
the representatives

```text
infinity=(0,1),  0=(1,0),  1=(1,1),  -1=(1,-1),
```

with determinant as symplectic form, one tournament is

```text
0 -> infinity,       1 -> infinity,      -1 -> infinity,

0 -> 1,              1 -> -1,            -1 -> 0. (32)
```

It is a directed triangle feeding a sink, with score multiset
`{0,2,2,2}`. Only its switching class, up to global converse and
relabeling, is intrinsic.

The exact preservation/loss sidecar is:

```text
vertices:
  four cyclic order-three subgroups;

pairwise observable:
  the exponent +/-1 of the Weil pairing;

orientation gauge:
  generator sign at each vertex;

ties:
  none;

preserved:
  projective sector incidence and alternating-pairing switching class;

lost:
  oriented point T versus -T,
  the primitive root zeta,
  lambda and j,
  the lift selection epsilon and divisor positions,
  P,Q,H4,S4,A and their scales,
  the isogeny coordinates,
  every Keller/Faber/one-form sidecar.              (33)
```

Accordingly, the tournament is a compact carrier for the four-sector
symplectic gauge. It is neither equivalent to (16) nor an H4
certificate.

## 7. Scope and exact reproduction

Every genuine coprime H4 survivor therefore supplies:

```text
a smooth elliptic resolvent E;

a connected etale three-isogeny X->E;

a nonzero Cardano class in Pic^0(E)[3];

a compatible lift-selection sector in (16);

a root lambda of M_H4;

a four-vertex Weil switching class with the loss sidecar (33).       (34)
```

The first five objects are exact geometric structure, not evidence that
the locus is empty. The next mathematical target is to reinsert the
structured trace (29), preferably in one of THM-2373's rational charged
sections, without discarding the scale or divisor sidecars.

Run

```bash
python3 04-computation/jc2_degree18_h4_elliptic_three_isogeny_atlas_thm2387.py
python3 -O 04-computation/jc2_degree18_h4_elliptic_three_isogeny_atlas_thm2387.py
```

Both transcripts are byte-identical to the stored output. The companion
checks the norm and Cardano constants, general Weierstrass invariant
identities, the flex at `(0,0)`, both exact discriminants and j-maps,
the translation formulas and (22), the full symbolic binary-quartic
discriminant identity, (27)--(28), smooth and singular controls, all
eight nonzero vectors and four projective lines over `F3`, every edge in
(32), and the switching action of all four generator negations. The
divisor, connectedness, etaleness, dual-isogeny interpretation, and
atlas scope are the mathematical proof above, not computer assumptions.
QED.
