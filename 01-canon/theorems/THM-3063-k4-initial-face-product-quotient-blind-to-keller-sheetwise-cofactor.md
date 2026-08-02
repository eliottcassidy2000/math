---
id: THM-3063
title: "K4 initial-face product quotient is blind to the sheetwise Keller cofactor"
status: >
  PROVED + VERIFIED-EXACT (independent audit pending).  Insert four
  branchwise primitive-element cofactors into the quartic root-difference K4
  by the natural vertex-unit gauge.  Every matching monomial then sees only
  the product of the four cofactors, and the oriented contraction remains the
  identically zero Pluecker sum.  Keller constancy instead requires four
  separate cofactor-times-derivative values to agree.  An exact product-one
  cofactor twist preserves the full matching triple, hence every THM-3058
  initial face and filtered contraction jet, while destroying Keller
  constancy.  The hostile can be chosen C3-inertia-equivariant; a residue-7
  packet gives the exact minimum face and residues.  This refutes only the
  vertex-gauge instantiation, not every possible branch-labelled
  augmentation, and excludes no Keller map or JC branch.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch
  - THM-3058-k4-hafnian-initial-face-augmentation-and-unbounded-cancellation-jet
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
related:
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion
  - THM-3049-k4-matching-monomial-tropical-root-extraction-clutch
script: 04-computation/quartic_c3_hafnian_cofactor_blindness_thm3063.py
output: 05-knowledge/results/quartic_c3_hafnian_cofactor_blindness_thm3063.out
script_sha256: 4bca80fdca42fa928fc17e56aacd75635fb169076d022080acf4238bd524b85f
output_sha256: 76680f699ab6b49ac9a0d8fbe9da0478c6d4cc81d06260b4e9312363698f9bf1
hash_basis: LF-normalized bytes
---

# THM-3063 -- the hafnian face forgets the sheetwise Keller equation

**PROVED + VERIFIED-EXACT (independent audit pending).**

## 1. The proposed bridge and its exact contraction

[THM-3058](THM-3058-k4-hafnian-initial-face-augmentation-and-unbounded-cancellation-jet.md)
repairs valuation-only matching data by retaining the minimum matching face
and its initial-residue sum.  [THM-3059](THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample.md)
shows that the first missing coordinate on a quartic `C3` Jelonek face is
instead the branchwise primitive-element Jacobian cofactor.  This theorem
tests the most direct way to combine them and proves that it loses exactly
the needed branch identity.

Let `K` be a field of characteristic different from two, let
`z_0,z_1,z_2,z_3` be four distinct elements of a splitting field, and put

```text
e_ij=z_i-z_j.                                             (1)
```

Let `c_i` be four nonzero branchwise cofactor values.  Insert them by
THM-3058's vertex-unit gauge:

```text
B_ij=c_i c_j e_ij.                                       (2)
```

For the three perfect matchings, use the oriented amplitudes

```text
A_1=B_01B_23,
A_2=-B_02B_13,
A_3=B_03B_12.                                            (3)
```

Write

```text
C=c_0c_1c_2c_3,
P_1=e_01e_23,       P_2=e_02e_13,       P_3=e_03e_12.    (4)
```

Then, exactly rather than only in the associated graded,

```text
(A_1,A_2,A_3)=C(P_1,-P_2,P_3),                           (5)
A_1+A_2+A_3=C(P_1-P_2+P_3)=0.                            (6)
```

Equation (6) is the oriented Pluecker identity already isolated in
THM-3046 and THM-3058.  It has two immediate consequences:

1. every matching valuation is shifted by the common scalar `v(C)`, the
   minimum face is unchanged, and every initial residue is multiplied by
   the same nonzero residue of `C`;
2. the THM-3058 initial-face sum is zero, and the full filtered contraction
   jet vanishes at every depth because the contraction itself is identically
   zero.

Thus the proposed cofactor insertion does not move the quartic packet off
THM-3058's cancellation wall.

## 2. The product quotient

More is true: even before contracting, the complete matching triple (5)
depends on `(c_0,c_1,c_2,c_3)` only through `C`.  The product-one torus

```text
T={(d_0,d_1,d_2,d_3) in (K*)^4: d_0d_1d_2d_3=1}         (7)
```

acts by

```text
c_i -> d_i c_i                                           (8)
```

and fixes every `A_j` exactly.  Hence it fixes not merely the first
THM-3058 augmentation but all matching valuations, residues, units, and
arbitrarily deep jets.

This is the exact contraction loss:

```text
four labelled cofactors -> their product C.               (9)
```

The lost fibre has dimension three before inertia constraints.

## 3. What the Keller equation actually asks for

Let `a in K*` and

```text
f(T)=a product_i(T-z_i),
D_i=f'(z_i)=a product_(j!=i)(z_i-z_j).                   (10)
```

In a primitive-element graph chart, absorb the common leading coefficient,
the sign, and the coordinate-change determinant into `c_i`.  The four
physical Jacobian values are then

```text
J_i=c_iD_i.                                              (11)
```

This is the same chain-rule cofactor identified in THM-3059; the common
leading scalar `a` may equally be absorbed into every `c_i`.  A Keller chart
requires the sheetwise predicate

```text
J_0=J_1=J_2=J_3=kappa in K*.                             (12)
```

The product quotient sees only

```text
product_i J_i=C product_i D_i,                            (13)
```

which is a discriminant scalar.  It does not see the three ratios needed
for (12).

The failure is exact.  Start with the abstract Keller-compatible packet

```text
c_i=kappa/D_i.                                           (14)
```

For any `d in T`, set `c'_i=d_i c_i`.  Equations (5)--(9) give the identical
full matching triple, whereas

```text
J'_i=d_i kappa.                                          (15)
```

Unless all four `d_i` agree, (15) is not constant.  Therefore the Keller
predicate is not constant on a fibre of the THM-3058 vertex-gauge quotient.
No statistic of that quotient, including its unbounded contraction jet, can
decide (12).

## 4. Exact `C3` minimum face from THM-3059

The no-go survives the inertia symmetry relevant to the live frontier.
Use THM-3059's planar `m=1` reciprocal packet

```text
M(X)=-X^4+X^3+u,
u=s^3.                                                    (16)
```

At `s=0`, the three small roots have leading coefficients solving
`a^3=-1`, while the fixed unit root has residue one.  Over residue
characteristic seven choose

```text
(a_0,a_1,a_2)=(3,5,6).                                  (17)
```

In the edge order

```text
01,02,03,12,13,23,                                       (18)
```

the root-difference valuations are

```text
(1,1,0,1,0,0).                                          (19)
```

Every perfect matching therefore has valuation one.  The signed initial
residues of `(P_1,-P_2,P_3)` are

```text
(2,4,1) in F_7,                     2+4+1=0.             (20)
```

Thus the exact THM-3058 face is

```text
F={1,2,3},                         sigma=0.               (21)
```

The derivative valuations at the three small roots and the fixed root are

```text
v_s(D_i)=(2,2,2,0),                                      (22)
```

and their normalized residues are

```text
(6,5,3,6).                                               (23)
```

The Keller-compatible cofactors (14), with `kappa=1`, consequently have

```text
valuations (-2,-2,-2,0),
initial residues (6,3,5,6).                              (24)
```

Their product has valuation `-6` and residue one.  Hence the cofactor-gauged
matching face has valuation vector

```text
(-5,-5,-5),                                              (25)
```

with the same signed residue word `(2,4,1)` and zero initial sum.

## 5. Minimal inertia-equivariant hostile

In the `C3` local quartic algebra, an inertia-equivariant product-one twist
may be constant on the cubic orbit and take a separate value on the fixed
factor:

```text
d(t)=(t,t,t,t^(-3)).                                     (26)
```

Its product is one.  It preserves all three matching amplitudes exactly,
but it changes the Jacobian packet from

```text
(kappa,kappa,kappa,kappa)
```

to

```text
(t kappa,t kappa,t kappa,t^(-3)kappa).                   (27)
```

The latter is constant exactly when `t^4=1`.  The first prime-residue hostile
outside characteristics two and three is `F_7` with `t=2`:

```text
d=(2,2,2,1/8),                    product d_i=1,          (28)
residue(d)=(2,2,2,1).                                    (29)
```

The original Jacobian residues are `(1,1,1,1)`; after (28) they are
`(2,2,2,1)`.  Nevertheless the full matching triple is byte-for-byte the
same algebraic triple, not merely the same face or clutch.

The twist (26) is legitimate in the normalized local quartic algebra: the
linear fixed factor and cubic factor are comaximal, and (26) is constant on
each factor.  This does **not** assert that both cofactor packets arise from
polynomial maps with the same target data.  It proves the logically required
quotient statement: the proposed matching observable forgets a coordinate
on which Keller constancy depends.

## 6. Exact stopping boundary

The natural THM-3058 bridge therefore has the connection ledger

```text
source:       labelled quartic roots plus four branch cofactors;
map:          vertex-gauged K4 matching monomials;
preserved:    root matching packet and product of cofactors;
destroyed:    three sheetwise cofactor ratios;
target test:  THM-3058 minimum face and contraction jet;
hostile:      d=(2,2,2,1/8) on the C3/fixed decomposition. (30)
```

After `C3` symmetry, the cheapest missing scalar is the pointed ratio

```text
R_J=J_C3/J_fixed.                                       (31)
```

Keller requires `R_J=1`.  THM-3058's `V4`-blind matching quotient has no
pointed fixed sheet and cannot form (31).  A valid next augmentation must
retain:

1. the actual inertia-fixed sheet;
2. the branchwise primitive-element cofactor in the full graph order;
3. its normalized unit residue on the fixed and cubic factors; and
4. affine-source regularity of that cofactor/section.

Items 1--4 are exactly the pointed order/owner sidecars of THM-3038,
THM-3042, THM-3057, and THM-3059.  Another symmetric matching contraction is
not a substitute.

```text
PROVED HERE:       exact product-only cofactor dependence;
                   product-one torus blindness;
                   incompatibility with sheetwise Keller constancy;
                   exact C3 face, valuations, and F7 residues;
                   C3-equivariant minimal hostile.

REFUTED HERE:      the natural vertex-gauge use of THM-3058 as a Keller-C3
                   discriminator.

NOT PROVED:        impossibility of every alternative branch-labelled
                   augmentation; physical realization of the twisted packet;
                   exclusion of a Keller C3 component, A4, S4, G1, JC(2),
                   or DC(2).                              (32)
```

## 7. Exact companion

Run

```text
python3 04-computation/quartic_c3_hafnian_cofactor_blindness_thm3063.py
python3 -O 04-computation/quartic_c3_hafnian_cofactor_blindness_thm3063.py
```

Both modes must LF-byte-match the stored transcript.  The companion checks
the symbolic Pluecker and vertex-gauge factorizations, product-one torus,
abstract Keller packet and its twist, exact rational `C3` twist, residue-7
minimum face, derivative/cofactor valuations and residues, and the unchanged
gauged matching packet.  Every truth-bearing check uses explicit runtime
exceptions rather than Python assertions.
