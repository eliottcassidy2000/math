---
id: THM-3622
title: "Compiler one-observable graph closure, normalization, and arm debt"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the THM-3561 compiler with the single separating observable z=xq,
  the affine graph closure has an exact saturated four-relation presentation.
  The source map is injective, birational, and quasi-finite with image the
  complement of one smooth line, but the closure is nonnormal along a second
  covered line. Its smooth normalization has three boundary arms: the source
  covers d=1, omits d=-2 over the covered singular line, and omits d=0 over
  the absent smooth line. The missing coordinate x has a simple pole on both
  omitted arms.
source: root / graph_closure_geometry normalization wildcard, 2026-08-21
audit: >
  PASS -- an independent hostile reconstruction recovered both saturated
  presentations, the exact image and singular tangent debt, the smooth
  normalization and its three arms, the glued inverse charts, and every
  divisorial valuation; normal, optimized, and stored 115-gate transcripts
  are byte-identical.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3618-compiler-one-graph-observable-fibre-separator-no-embedding
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3614-russell-cylinder-collision-free-full-linear-projection-rigidity
script: 04-computation/jc2_compiler_one_observable_graph_closure_normalization_thm3622.py
output: 05-knowledge/results/jc2_compiler_one_observable_graph_closure_normalization_thm3622.out
script_sha256: 8f1f338f1cfc5dc132564593c69a83ac82adaee023739b4a59a241c33fcbcd7e
output_sha256: 7eb9821770b8831d35b81ef3fa503df4fdb7a94cd7dadcec6c89842936ffc39c
hash_basis: raw LF bytes
---

# THM-3622 -- compiler one-observable graph closure, normalization, and arm debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

All rings and closed points are over `C`. Retain the THM-3561 compiler and
the first separator from THM-3618:

```text
d=1+x^2q,
b=(d-1)(d+2)^2,
c=xd(d+2),
e=q(d+3),
z=xq,
A=C[x,q],                       R=C[b,c,e].             (1)
```

Since

```text
e=4q+z^2,                       Q:=(e-z^2)/4=q,         (2)
```

the one-observable graph ring is

```text
S=R[z]=C[b,c,Q,z]  subset  A,
Z_graph=Spec(S).                                             (3)
```

The symbol `Z_graph` denotes the surface; uppercase `Z` below denotes the
coordinate corresponding to `z`.

## 0. Statement

The closure `Z_graph` has the exact finite presentation in section 1. The
map induced by `S subset A`,

```text
psi:A2_(x,q) -> Z_graph,                                  (4)
```

is injective, birational, and quasi-finite. Its exact set-theoretic image is

```text
psi(A2)=Z_graph minus L_infinity,

L_infinity={B=-4,Q=Z=0},             C free.             (5)
```

The complement is a smooth affine line. Nevertheless `(4)` is not an open
immersion: the entire singular locus is the covered affine line

```text
L_0={B=0,Q=Z=0},                      C free.             (6)
```

At every point of `L_0`, the target tangent space has dimension three while
the source differential spans only a two-plane.

Adjoining the integral element `d` gives the normalization

```text
T=S[d],                         Z_tilde=Spec(T).          (7)
```

The surface `Z_tilde` is smooth. Above `Q=Z=0` it has three disjoint lines

```text
N_r={d=r,Q=Z=0},        C free,       r in {1,-2,0}.      (8)
```

The finite normalization map satisfies

```text
N_1 and N_-2 -> L_0,            N_0 -> L_infinity.       (9)
```

Moreover,

```text
A2 is isomorphic to Z_tilde minus (N_0 union N_-2).       (10)
```

Thus `N_1` is the arm covered by the source, `N_-2` is a missing arm hidden
pointwise behind the same line `L_0`, and `N_0` is the missing arm whose image
is visibly absent. In the common function field,

```text
polar_divisor_(Z_tilde)(x)=N_0+N_-2.                    (11)
```

Both poles are simple.

## 1. Exact saturated presentation of the graph closure

Let `P=C[B,C,Q,Z]`. Define

```text
j1=B^2+4B-4C^2Q-C^2Z^2,

j2=BQ+BZ^2-3CQZ-CZ^3,

j3=BZ^3+2CQ^2-2CQZ^2-CZ^4-6QZ-2Z^3,

j4=CQ^3-3Q^2Z-4QZ^3-Z^5.                              (12)
```

Then the exact kernel of

```text
P -> A,              (B,C,Q,Z) |-> (b,c,q,z)           (13)
```

is

```text
J=(j1,j2,j3,j4),              S=P/J.                   (14)
```

The four displayed polynomials are the reduced lexicographic Groebner basis
for the order `B>C>Q>Z`. Their mechanism is more transparent from the two
generic denominator equations

```text
BQ^3=Z^2(3Q+Z^2)^2,
CQ^3=Z(Q+Z^2)(3Q+Z^2).                                 (15)
```

Precisely,

```text
J=(the two polynomials in (15)):Q^infinity.             (16)
```

This saturation is load-bearing. The unsaturated pair `(15)` retains
boundary torsion and is not an exact presentation.

For proof of exactness, after inverting `Q`, equations `(15)` give

```text
B=Z^2(3Q+Z^2)^2/Q^3,
C=Z(Q+Z^2)(3Q+Z^2)/Q^3,                                (17)
```

so

```text
(P/J)_Q = C[Q,Q^(-1),Z].                               (18)
```

Equation `(16)` says that `P/J` has no `Q`-torsion. The map `(13)` becomes
an isomorphism after inverting `Q`; therefore its kernel is `Q`-torsion and
must vanish. In particular `J` is prime, `Z_graph` is an irreducible affine
surface, and `(14)` is not merely a set-theoretic equation list.

The first equation is the original smooth compiler-surface relation in the
new coordinates:

```text
j1=B(B+4)-C^2E,                    E=4Q+Z^2.            (19)
```

## 2. Birationality, exact image, and the missing line

On `Q!=0`, `(17)` and

```text
x=Z/Q,
d=1+Z^2/Q                                                     (20)
```

give an exact inverse to `(4)`. On the boundary `Q=0`, the reduced
Groebner basis becomes

```text
B^2+4B-C^2Z^2,          BZ^2,          Q,          Z^3. (21)
```

Hence, set-theoretically,

```text
Q=0  implies  Z=0 and B(B+4)=0,                         (22)
```

so the boundary is exactly `L_0 union L_infinity`.

The source boundary is explicit:

```text
q=0  implies  (b,c,Q,z)=(0,3x,0,0).                    (23)
```

It covers all of `L_0` and none of `L_infinity`. Together with the inverse
on `Q!=0`, this proves `(5)`. THM-3618 proves that `z` separates every
closed compiler fibre, so `(4)` is injective. Formula `(20)` proves equality
of function fields, hence birationality. Its closed fibres contain at most
one point, so this finite-type morphism is quasi-finite.

The image is dense and open but not closed. Consequently `(4)` is not
finite. Section 6 gives the sharper algebraic reason: `A=S[x]`, while `x`
has divisorial poles on the normalization and therefore is not integral over
`S`.

## 3. Singular locus and the tangent direction missed by separation

Let `M_J` be the `4x4` Jacobian matrix of `(j1,j2,j3,j4)` with respect to
`(B,C,Q,Z)`. Exact reduction gives

```text
J+(all 2x2 minors of M_J)=(B,Q,Z).                      (24)
```

Thus

```text
Sing(Z_graph)=L_0                                                  (25)
```

scheme-theoretically. In particular the singular locus has codimension one,
so `Z_graph` is not normal. At

```text
p_kappa=(0,kappa,0,0) in L_0,
```

the Jacobian has only one nonzero row,

```text
(4,0,-4kappa^2,0).                                    (26)
```

Therefore

```text
T_(p_kappa)Z_graph={delta B=kappa^2 delta Q},
dim T_(p_kappa)Z_graph=3.                               (27)
```

At the source point `(x,q)=(kappa/3,0)`, the two columns of `d psi` are

```text
(0,3,0,0),
(kappa^2,4kappa^3/27,1,kappa/3).                       (28)
```

They span the two-plane in `(27)` satisfying the additional equation

```text
delta Z=(kappa/3)delta Q.                              (29)
```

The pure `Z` tangent direction belongs to `(27)` but not to `(29)`. At the
generic point of `L_0`, a transverse slice has tangent cone

```text
Q(kappa Q-3Z)=0.                                      (30)
```

Indeed, fixing `C=kappa`, equation `j1` gives
`B=kappa^2Q+O((Q,Z)^2)`; the degree-two initial form of `j2` is then
`kappa Q(kappa Q-3Z)`, and `kappa` is a unit at the generic point.

The branch `3Z=kappa Q` is the covered arm `N_1`; the branch `Q=0` is the
omitted arm `N_-2`. Pointwise separation detects neither this second branch
nor the failure of the inverse to be regular there.

By contrast, at `(-4,kappa,0,0) in L_infinity`, two independent Jacobian rows
are

```text
(-4,0,-4kappa^2,0),             (0,0,-4,0),            (31)
```

so `L_infinity` lies in the smooth locus.

## 4. The exact smooth normalization

In `C[d,C,Q,Z]`, put

```text
k1=d^3+d^2-2d-CZ,
k2=d^2Z+2dZ-CQ,
k3=dQ-Q-Z^2,
k4=dZ^3+3QZ+3Z^3-CQ^2,
k5=CQ^3-3Q^2Z-4QZ^3-Z^5.                              (32)
```

Then

```text
T=S[d]
 =C[d,C,Q,Z]/(k1,k2,k3,k4,k5).                         (33)
```

These five polynomials are the reduced lexicographic Groebner basis for
`d>C>Q>Z`. The first three attractive relations are

```text
CZ=d(d-1)(d+2),
CQ=d(d+2)Z,
Q(d-1)=Z^2.                                           (34)
```

They are only a generic presentation. The exact ideal in `(33)` is their
`Q`-saturation; `k4,k5` remove the boundary torsion.

As in section 1, this also proves exactness of `(33)`: after inverting `Q`,
the relations give `d=1+Z^2/Q` and determine `C`, while the saturation says
that the quotient has no `Q`-torsion. The quotient therefore embeds in its
domain localization and maps injectively onto the subring `S[d]`.

The map to `(14)` is

```text
B=(d-1)(d+2)^2=d^3+3d^2-4,
(C,Q,Z)=(C,Q,Z).                                      (35)
```

It is finite because `d` satisfies the monic equation

```text
d^3+3d^2-(B+4)=0.                                     (36)
```

It is birational because `d=1+Z^2/Q` when `Q!=0`.

The ring in `(33)` is smooth. On `Q!=0`, `(34)` gives the coordinate chart
`C[Q,Q^(-1),Z]`. When `Q=0`, equations `(32)` force `Z=0` and

```text
d(d-1)(d+2)=0.                                        (37)
```

At `d=0,1,-2`, respectively, explicit `2x2` Jacobian minors have values

```text
2,                9,                -18,              (38)
```

independently of `C`. These charts cover the whole surface. Equivalently,
the ideal generated by `(k1,...,k5)` and all Jacobian `2x2` minors is the
unit ideal. Thus `T` is normal. Being finite, birational, and normal over
`S`, it is exactly the integral closure of `S` in `C(x,q)`.

## 5. The three arms and the normalization-open source

The three lines in `(8)` have images

```text
d=1:       B=0,        N_1 -> L_0;
d=-2:      B=0,        N_-2 -> L_0;
d=0:       B=-4,       N_0 -> L_infinity.              (39)
```

The source has `d=1` whenever `q=0`, so it covers `N_1`. It omits `N_0`
and `N_-2`.

This is an exact scheme statement, not only a point inventory. The open
complement of `N_0 union N_-2` is covered by

```text
Q!=0,                         d(d+2)!=0.               (40)
```

On these charts define, respectively,

```text
x=Z/Q,                        x=C/[d(d+2)].             (41)
```

The formulas agree on the overlap because `CQ=d(d+2)Z`. The other two
relations in `(34)` give

```text
xQ=Z,                         xZ=d-1,                  (42)
```

so `d=1+x^2Q` and `C=xd(d+2)`. This constructs the inverse to the source
map and proves `(10)`.

The normalization explains two logically different debts:

| debt | normalization arm | target effect | pointwise visibility |
|---|---|---|---|
| visible closure debt | `N_0` | its smooth image `L_infinity` is absent | visible as missing points |
| invisible tangent debt | `N_-2` | its image `L_0` is already covered by `N_1` | invisible to a fibre separator |

A curve approaching `N_-2` has target limit on the covered line `L_0`, while
its source coordinate satisfies `x -> infinity`. Hence the set-bijection

```text
A2 -> Z_graph minus L_infinity                              (43)
```

does not have a continuous or regular inverse along `L_0`.

One closed hostile witness in the source is

```text
d=-2,                 x^2q=-3.                         (44)
```

Writing `x=u^(-1)` gives

```text
(B,C,Q,Z)=(0,0,-3u^2,-3u),              u!=0.          (45)
```

Its target image has closure containing the covered origin of `L_0`, whereas
the lift has `x=u^(-1)` and escapes to infinity. Thus the image of this closed
source curve is not closed in `Z_graph minus L_infinity`, proving directly
that the inverse set-map is not continuous.

## 6. Exact valuations and the pole divisor

At the generic point of each `N_r`, the residue of `C` is a unit. The
equation

```text
CZ=d(d-1)(d+2)                                         (46)
```

and the nonzero derivatives at `r=0,1,-2` show that `Z` is a uniformizer.
The leading terms are:

| arm | leading terms in uniformizer `Z` | source status |
|---|---|---|
| `N_1` | `d-1~CZ/3`, `Q~3Z/C`, `B~3CZ`, `E~12Z/C`, `x~C/3` | covered |
| `N_-2` | `d+2~CZ/6`, `Q~-Z^2/3`, `B~-C^2Z^2/12`, `E~-Z^2/3`, `x~-3/Z` | omitted over `L_0` |
| `N_0` | `d~-CZ/2`, `Q~-Z^2`, `B+4~3C^2Z^2/4`, `E~-3Z^2`, `x~-1/Z` | omitted over `L_infinity` |

Thus `x` is regular on the normalization-open source, has a simple pole on
each omitted arm, and has no other pole. This proves `(11)` and also proves
that `x` is not integral over `S`. Hence `A=S[x]` is not finite over `S`.

The divisorial valuation on the visibly missing curve can be read directly
on `Z_graph`. At the generic point of `L_infinity`, the denominators below
are units and the relations give

```text
Q =Z^2(CZ-B)/(B-3CZ),
E =Z^2(CZ-3B)/(B-3CZ),
x =(B-3CZ)/[Z(CZ-B)],
d =2CZ/(B-CZ),
B+4=C^2E/B.                                            (47)
```

Therefore `Z` is a uniformizer and

```text
v_Linfinity(Z)=1,
v_Linfinity(Q)=v_Linfinity(E)=v_Linfinity(B+4)=2,
v_Linfinity(x)=-1.                                     (48)
```

The corresponding leading residues are

```text
Q/Z^2 -> -1,
E/Z^2 -> -3,
(B+4)/Z^2 -> 3C^2/4,
xZ -> -1,
d/Z -> -C/2.                                           (49)
```

THM-3618's punctured `d=0` path is the specialization `C=0` of this whole
missing divisorial arm. The closure debt is an affine line, not merely the
single limiting point visible on that one path.

## 7. Relation to THM-3618 and scope ledger

THM-3618 proves that `z` separates every closed fibre and that `R[z]` is
strictly smaller than `A`. The present theorem identifies the exact geometry
behind both facts:

* separation removes source collisions but does not supply the second
  normalization branch over `L_0`;
* `N_-2` is lost without changing the point-set image;
* `N_0` is lost and exposes the missing line `L_infinity`;
* adjoining `x` pays both simple-pole debts, consistently with
  `R[z,x]=A`.

What is proved:

* the exact saturated presentations `(14)` and `(33)`;
* injectivity, birationality, quasi-finiteness, and exact image `(5)`;
* the scheme-theoretic singular locus and tangent-plane loss;
* the smooth finite normalization, its three arms, and the open
  identification `(10)`;
* the exact valuations and polar divisor `(11)`.

What is not claimed:

* that pointwise injectivity implies a closed embedding, properness, or
  finiteness;
* a classification of closures for every separator `e^(m-1)z` with `m>1`;
* a classification of arbitrary multi-observable graph rings;
* any counterexample to, proof of, or reduction of `JC(2)`.

## 8. Exact companion contract

The companion verifies, without truth-bypassing assertion statements:

* compiler identities and `Q=q`;
* the four-relation reduced lex basis, the `Q`-saturation `(16)`, and the
  generic inverse;
* the exact `Q=0` boundary, both boundary lines, and the source line;
* the singular Jacobian ideal, target tangent space, and source differential;
* the five-relation normalization basis, saturation, monic integrality, and
  unit singular ideal;
* all three boundary arms, the two inverse charts, and their gluing;
* the three leading-term valuation packets and the exact local identities
  `(47)--(49)`;
* nonfiniteness and syntax hygiene.

The saturation-to-kernel implication, finite-birational normalization
argument, and divisorial interpretation are proof-driven from those exact
identities. The companion contains no finite degree box or probabilistic
substitute.
