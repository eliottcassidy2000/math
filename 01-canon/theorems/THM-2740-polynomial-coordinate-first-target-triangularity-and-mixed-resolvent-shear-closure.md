---
id: THM-2740
title: "Polynomial-coordinate first-target triangularity and mixed resolvent shear closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over a
  characteristic-zero field, if U is a polynomial coordinate
  of k[x,y] and Jac(U,V) is a nonzero constant, then relative to a conjugate
  coordinate S one has V=(kappa/Jac(U,S))S+G(U); hence (U,V) is a polynomial
  automorphism.  On every polynomial graph in the THM-2696 quotient,
  A=x^2-2y=-2t is a coordinate.  Therefore `(A,V)` is triangular for every
  restricted polynomial second target V, including every quotient-compatible
  V=Phi(A,B,d) and every mixed shear B+H(A,d).  The mixed primitive and
  explicit inverse are exact.  Neither the noncoordinate-first-target lane,
  nongraph resolvent surfaces, JC(2), nor DC(2) is closed.
source: a4-resolvent-next-gate-scout-2026-07-28
audit: thm2705-2709-audit-2026-07-28-coordinate-first (independent all-degree, mixed-shear, characteristic-boundary, and exact replay audit)
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
related:
  - THM-2699-affine-plane-linear-projection-keller-slice-classification
  - THM-2702-polynomial-graph-coordinate-projection-keller-classification
  - THM-2705-linear-target-planes-containing-d-polynomial-graph-keller-classification
  - THM-2709-complete-linear-target-plane-polynomial-graph-keller-classification
  - THM-2715-nonlinear-d-target-shear-polynomial-graph-keller-classification
script: 04-computation/jacobian_coordinate_first_mixed_shear_thm2740.py
output: 05-knowledge/results/jacobian_coordinate_first_mixed_shear_thm2740.out
script_sha256: c152534cbfc0356b98de84ee09da31a17143c40b694d25b923dea75dcef1a3dd
output_sha256: d5a10c83a85f19a9c155d806c35e5a8ef44ec5b1866bbb1c85d4bb3fda4049c4
hash_basis: LF-normalized bytes
---

# THM-2740 -- polynomial-coordinate first-target triangularity and mixed resolvent shear closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The coordinate-pair calculations in THM-2702 and the nonlinear shear
classification in THM-2715 share a simpler invariant.  On their graph
plane, the first target `A` is already a polynomial coordinate.  Once one
component of a planar Keller pair is a coordinate, the second cannot hide a
counterexample: constant Jacobian makes it affine-linear in the conjugate
coordinate.

This closes a much larger target family than one-variable shears.  The
second target may be any polynomial after restriction to the graph,
including arbitrary joint dependence on all three quotient coordinates
`A,B,d`.

## 1. Coordinate-first triangularity

Let `k` be a field of characteristic zero.  Suppose

```text
U,S in k[x,y],                 k[U,S]=k[x,y].             (1)
```

Then

```text
delta=Jac(U,S) in k^*.                                    (2)
```

Let `V in k[x,y]` satisfy

```text
Jac(U,V)=kappa in k^*.                                    (3)
```

Express `V=P(U,S)` with `P in k[U,S]`.  The chain rule gives

```text
Jac(U,V)=delta partial_S P.                               (4)
```

Hence

```text
partial_S P=kappa/delta.                                  (5)
```

Characteristic zero is load-bearing: integrating `(5)` in the polynomial
ring yields

```text
V=(kappa/delta)S+G(U),             G in k[U].             (6)
```

Conversely every polynomial in `(6)` has Jacobian `kappa`.  The inverse is
explicit.  From target values `(u,v)`, recover

```text
s=(delta/kappa)(v-G(u)),                                 (7)
```

and then apply the polynomial inverse of the coordinate automorphism
`(x,y)->(U,S)`.  Therefore

```text
U a polynomial coordinate and Jac(U,V) in k^*
   ==> (U,V) is a polynomial automorphism.                (8)
```

In particular, every hypothetical counterexample to the planar Jacobian
conjecture must have **neither** component a polynomial coordinate.  This is
a reduction of the search space, not a proof that one component of every
Keller pair is a coordinate.

## 2. The graph coordinate in the reflection quotient

Work over `C` with the THM-2696 coordinates

```text
A=x^2-2y,             B=y^2-2xz,             d=z,         (9)
```

and restrict to any polynomial graph

```text
Gamma_f={z=f(x,y)} ~= A2_(x,y).                          (10)
```

Put

```text
t=y-x^2/2.                                                (11)
```

The change `(x,y)<->(x,t)` is a polynomial automorphism and

```text
A=-2t,
Jac_(x,y)(A,x)=2.                                         (12)
```

Thus `(A,x)` is a coordinate pair on every graph, independently of `f`.
For **any** polynomial `V in C[x,y]`, equation

```text
Jac(A,V)=kappa!=0                                         (13)
```

is, in the coordinates `(x,t)`, simply

```text
2 partial_x V=kappa.                                      (14)
```

Sections 1 and 2 give the complete response form

```text
V=(kappa/2)x+G(t),                 G in C[t].             (15)
```

The inverse is

```text
t=-A/2,
x=2(V-G(t))/kappa,
y=t+x^2/2.                                                (16)
```

Consequently every such pair is a triangular polynomial automorphism.

The word "any" in `(13)` is source-ring exact.  In particular it includes
the restriction of every ambient polynomial

```text
Phi(A,B,d) in C[A,B,d]                                   (17)
```

to `Gamma_f`.  Hence all quotient-compatible target pairs

```text
(A,Phi(A,B,d))                                            (18)
```

are closed on arbitrary polynomial graphs: if their Jacobian is a nonzero
constant, they are automorphisms of the graph plane.

## 3. The arbitrary mixed-shear primitive

The family named as the next target in THM-2715 is

```text
V=B+H(A,d),                  H in C[A,d].                 (19)
```

It is already covered by `(18)`, but retaining its primitive explains why
the earlier one-variable classification generalized without a new degree
argument.

Write

```text
F(x,t)=f(x,t+x^2/2),
c=-kappa/4.                                               (20)
```

Since `A=-2t` is fixed under `partial_x|_t`, direct differentiation of
`(19)` and `(13)` gives

```text
xF_x+F-xt-x^3/2
  -(1/2)H_d(-2t,F)F_x=c.                                  (21)
```

The left side is an exact derivative in `x`.  Therefore there is a unique
`K in C[t]` for the chosen `F` such that

```text
xF-H(-2t,F)/2=x^2t/2+x^4/8+cx+K(t).                      (22)
```

Conversely, differentiating `(22)` recovers `(21)`.  Using
`y=t+x^2/2` in `B=y^2-2xF`, equation `(22)` cancels every occurrence of
`F` and gives

```text
B+H(A,d)=t^2-2cx-2K(t)
          =t^2+(kappa/2)x-2K(t).                         (23)
```

This is `(15)` with

```text
G(t)=t^2-2K(t).                                           (24)
```

Its explicit inverse is

```text
t=-A/2,
x=(2/kappa)(V-t^2+2K(t)),
y=t+x^2/2.                                                (25)
```

No classification of the possible graph polynomials `F` is required for
invertibility.  THM-2715 remains stronger in a different direction: for
`H=H(d)` it classifies every graph that can occur, including the sharp
quadratic coefficient wall and its two sections.

## 4. Sharp controls and failure boundaries

### 4.1 Genuine mixed dependence and two graph sections

Take

```text
H(A,d)=-d^2+2d+A^3+2A.                                   (26)
```

Both graphs

```text
f_1=-x-y,
f_2=y-x+2                                                (27)
```

satisfy

```text
Jac(A,B+H(A,d))=-4.                                      (28)
```

Despite being distinct sections, they induce the same target map.  In
`(x,t)` coordinates it is

```text
(A,V)=(-2t,-2x-2t+(-2t)^3+2(-2t)),                       (29)
```

which has the inverse `(16)`.  Thus joint `A,d` dependence is genuinely
present but supplies no new escape.

### 4.2 Nonzero response is essential

If `kappa=0`, then every `V=G(U)` satisfies `Jac(U,V)=0`; the pair has
one-dimensional image.  Division by `kappa` in `(7)` and `(16)` is therefore
not cosmetic.

### 4.3 Characteristic zero is essential

In characteristic `p>0`,

```text
(x,y) -> (x,y+y^p)                                        (30)
```

has Jacobian one, but is not a polynomial automorphism: an automorphism of
`k[x][y]` fixing `x` must have degree one in `y`, whereas the second
coordinate in `(30)` has degree `p`.  Equation `(5)` fails to imply `(6)`
because `partial_S(S^p)=0`.

### 4.4 The live complex boundary

Over `C`, no counterexample is supplied by dropping the coordinate
hypothesis.  Rather, that is exactly where the planar Jacobian conjecture
remains.  Proving from the Keller condition that one component is a
polynomial coordinate would itself settle the conjecture.

Likewise an arbitrary resolvent surface need not be a global graph, and the
restriction of `A` to it need not be a polynomial coordinate.  The theorem
does not transport `(12)` through a birational chart or a normalization.

## 5. Exact companion

Run

```text
python 04-computation/jacobian_coordinate_first_mixed_shear_thm2740.py
python -O 04-computation/jacobian_coordinate_first_mixed_shear_thm2740.py
```

and compare both transcripts with

```text
05-knowledge/results/jacobian_coordinate_first_mixed_shear_thm2740.out.
```

The companion uses exact symbolic arithmetic and explicit failure raises.
It checks:

- `(14)` on the generic fifteen-coefficient polynomial of total degree at
  most four and recovers exactly five free coefficients of `G(t)`;
- the nonlinear coordinate pair `U=x+y^2,S=y`, a degree-four `G(U)`, and
  the exact inverse at `kappa=7`;
- the mixed polynomial `(26)`, both sections `(27)`, their common target,
  Jacobian `-4`, and the exact inverse;
- the zero-response boundary and the characteristic-three hostile `(30)`.

Normal and optimized transcripts byte-match the stored twelve-line output.
The script contains no Python `assert` nodes.  The all-degree proof is
Sections 1--3, not the bounded generic coefficient check.

## 6. Consequence and next target

The closed graph hierarchy is now

```text
literal coordinate pairs                 THM-2702;
all linear target planes                 THM-2705/2709;
one-variable nonlinear d-shears          THM-2715;
every second target paired with A        THM-2740.        (31)
```

The next meaningful graph target is not another polynomial added to the
second coordinate.  It is a pair whose first component is a nonlinear
combination `U(A,B,d)` that is not already known to be a polynomial
coordinate after graph restriction.  A useful finite decision problem is:

```text
classify which low-degree U(A,B,d)|_(Gamma_f) are coordinates,
then apply (8) before solving any second-target PDE.       (32)
```

This is a coordinate-detection problem, not a new Jacobian calculation.

The theorem proves no statement about a Weyl-algebra endomorphism or the
unsettled `JC(2)`/`DC(2)` provenance.  It closes exactly the coordinate-first
graph lane.

An independent hostile audit rederived the all-degree coordinate-ring proof,
the graph-coordinate chain rule, and the mixed-shear primitive and inverse.
It checked the zero-response and positive-characteristic boundaries, matched
both declared hashes, and replayed normal, optimized, and stored transcripts.

QED.
