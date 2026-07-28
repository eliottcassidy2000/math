---
id: THM-2768
title: "Modular C2-C3 quotients to A4/S4 and Bass-Serre cycle ranks"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  PSL2(Z)=C2*C3 has explicit tetrahedral and octahedral quotients
  A4/S4 obtained by adding the (2,3,3)/(2,3,4) face relation.  Their kernels
  are torsion-free free groups F3/F5.  The finite (2,3)-biregular Bass--Serre
  quotients have cycle ranks 3/5 and suppress to K4/the cube.  Under the
  matching quotient, the binary generator collapses in A4/V4=C3 but survives
  in S4/V4=S3.  These are finite quotients of the modular free product, not
  identifications with the A4/S4 monodromy of THM-2766.
source: root/modular-binary-ternary-a4-s4-2026-07-28
depends_on:
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
related:
  - THM-2646-braid-modular-exponent-fibre-product-and-conjugacy-classification
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
script: 04-computation/modular_c2_c3_a4_s4_bass_serre_thm2768.py
output: 05-knowledge/results/modular_c2_c3_a4_s4_bass_serre_thm2768.out
script_sha256: bdfedadcca55c215a1218b86ab91857f0e950f005602b483a266bcb01cee9595
output_sha256: 27cd3383ab3fa3a293101ce837ab8d9b37e13506d3bd845358112a967673b91a
hash_basis: LF-normalized bytes
---

# THM-2768 -- modular binary/ternary quotients and their lost relations

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

The dyadic and ternary grammars do coexist on one exact object, but the
object is the modular free product and the live finite monodromy groups are
its quotients.  The extra face relation is load-bearing.  It turns the
infinite `(2,3)` Bass--Serre tree into the tetrahedral or octahedral finite
frame and changes the fate of the binary generator under the `V4` matching
quotient.

## 1. The common modular carrier

Write

```text
Gamma=PSL_2(Z)=<x,y | x^2=y^3=1>=C2*C3.                (1)
```

Its Bass--Serre tree has vertices

```text
Gamma/<x>   disjoint_union   Gamma/<y>                  (2)
```

and one edge for every element of `Gamma`.  The two vertex types have
degrees two and three.  This is the canonical object on which the binary and
ternary free factors coexist without an additional relation.

The Farey tessellation and the familiar fraction trees are related views,
not literal synonyms for `(2)`: a rooted Stern--Brocot or Calkin--Wilf tree
chooses a positive sector and an orientation, while `(2)` is the full
unoriented Bass--Serre carrier.  Nothing below identifies an LRC dyadic tower
or a `C3` chronology with that sector without an explicit action map.

## 2. The tetrahedral and octahedral quotients

On `{1,2,3,4}`, define

```text
x_A=(12)(34),       y_A=(123),
x_S=(12),           y_S=(234).                         (3)
```

With composition read from right to left, exact permutation multiplication
gives

```text
ord(x_A),ord(y_A),ord(x_A y_A)=(2,3,3),
ord(x_S),ord(y_S),ord(x_S y_S)=(2,3,4).                (4)
```

The first pair generates `A4`; the second generates `S4`.  The classical
spherical triangle presentations therefore give

```text
Gamma/<<(xy)^3>> = Delta(2,3,3) = A4,
Gamma/<<(xy)^4>> = Delta(2,3,4) = S4.                  (5)
```

The spherical area formula gives orders `12` and `24`, respectively, and
the explicit images in `(3)` attain those orders.  Thus `(5)` is equality,
not merely a surjection.

Equation `(5)` is the first stopping rule for the proposed modular analogy:
`A4` and `S4` are not `C2*C3`.  They are obtained by adding a face relation
of length three or four.  Any transfer from modular normal forms to finite
monodromy must retain which relation was imposed.

## 3. Free kernels and finite Bass--Serre frames

Let `K_A,K_S` be the kernels in `(5)`.  Every finite-order element of a free
product is conjugate into a factor.  Since the two factor generators retain
their exact orders in both quotients, neither kernel meets a conjugate of a
factor.  Hence each kernel acts freely on the Bass--Serre tree and is a free
group.

The Euler characteristic is

```text
chi(Gamma)=1/2+1/3-1=-1/6.                              (6)
```

A torsion-free subgroup of index `d` is free of rank `1+d/6`.  Therefore

```text
K_A=F_3,                 K_S=F_5.                       (7)
```

The quotient graph `K\T` has vertices `G/<x>` and `G/<y>` and edges `G`.
Consequently its exact censuses are

```text
G=A4:   6 degree-two vertices, 4 degree-three vertices,
        12 edges, beta_1=12-6-4+1=3;

G=S4:  12 degree-two vertices, 8 degree-three vertices,
        24 edges, beta_1=24-12-8+1=5.                  (8)
```

Suppressing every degree-two vertex turns the first graph into `K4` and the
second into the cube graph.  Thus the tetrahedral and octahedral names are
literal incidence geometry.  The cube is a partial cube; `K4` is not.  The
partial-cube feature is therefore a property of the `S4` quotient frame,
not a consequence of the bare free product.

The cycle ranks in `(8)` agree with the kernel ranks in `(7)` because the
finite graphs are precisely the quotient graphs whose fundamental groups
are `K_A,K_S`.

## 4. The `V4` matching quotient treats the binary factor differently

Let `V4` act as the kernel of the permutation action of `S4` on the three
perfect matchings of four letters.  Then

```text
A4/V4=C3,                    S4/V4=S3.                  (9)
```

For the generators `(3)`, one has

```text
x_A in V4,      y_A mod V4 has order 3;
x_S mod V4 has order 2,      y_S mod V4 has order 3.   (10)
```

Thus the `A4` matching quotient collapses the binary generator completely,
while the `S4` matching quotient retains both free-factor images.  This is
the exact group-theoretic reason the cyclic cubic branch sees only a ternary
quotient, whereas the full cubic branch sees `S3`.

The preimages of `V4` inside `Gamma` sharpen the distinction.  In the `A4`
branch the preimage is the kernel of `Gamma -> C3`, so Reidemeister--Schreier
gives

```text
<x> * <yxy^-1> * <y^2xy^-2> = C2*C2*C2,                (11)
```

and

```text
1 -> F3 -> C2*C2*C2 -> V4 -> 1.                        (12)
```

In the `S4` branch the preimage is the kernel of the surjection
`Gamma -> S3`.  Both factor maps are injective, so this index-six subgroup
is torsion-free and

```text
1 -> F5 -> F2 -> V4 -> 1.                              (13)
```

Equations `(12)--(13)` are the honest `V4`-torsor statements.  The middle
groups differ; there is no branch-independent `V4` tree move.

## 5. Relation to the quadratic-cubic pullback

THM-2766 proves that a rank-two, product-oriented Kummer plane over a cyclic
cubic base gives `A4`, while the same plane over a full `S3` cubic base gives
`S4`.  Equations `(3)--(13)` provide abstract modular quotient presentations
of those two groups and explain the behavior of their matching quotients.

They do not produce a canonical modular action on the cubic splitting field,
the Kummer torsor, a Keller source, or a Jelonek boundary.  Choosing the
generators `(3)` uses a labeling and, in the `S4` branch, a reflection class.
The field-level monodromy of THM-2766 does not supply that choice as an
affine polynomial operation.  In particular, the theorem gives neither the
unit/class-group carrier of THM-2655 nor a `JC(2)` or `DC(2)` conclusion.

## 6. Exact verification and boundary ledger

Run

```bash
python 04-computation/modular_c2_c3_a4_s4_bass_serre_thm2768.py
python -O 04-computation/modular_c2_c3_a4_s4_bass_serre_thm2768.py
```

The companion uses explicit exceptions, no floating point, and no
truth-bearing Python assertions.  It generates the two permutation groups,
checks all triangle orders and matching actions, constructs every coset and
edge of both finite Bass--Serre quotients, and exhausts all `4!`/`8!`
relabellings to identify the suppressed graphs with `K4` and the cube.  Both
modes byte-match the stored transcript.

```text
PROVED HERE (candidate):  explicit modular quotients to A4 and S4;
                          triangle relations (2,3,3) and (2,3,4);
                          torsion-free kernels F3 and F5;
                          finite Bass--Serre censuses and cycle ranks;
                          K4/cube suppression;
                          exact V4 matching-quotient behavior;
                          preimage extensions (12)--(13).

NOT PROVED:               A4 or S4 equals PSL2(Z);
                          a canonical modular action on THM-2766 data;
                          a literal identification of rooted fraction trees
                          with every dyadic/C3 tower in the repository;
                          a Keller, JC(2), DC(2), graceful-tree, or LRC(14)
                          consequence.                                  (14)
```

QED (candidate).
