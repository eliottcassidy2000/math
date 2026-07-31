---
id: THM-2726
title: "a21-transverse integral split response three-pole closure"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  Every
  geometrically integral member of the full chosen-sheet split degree-22
  Faber response family with a21 nonzero has three distinct normalization
  poles of the physical function q.  A nonconstant physical trajectory would
  pull all three pole fibres back to the source P1, contradicting THM-2723's
  two-pole capacity.  No nonzero first-flux or generic-genus hypothesis is
  used.  The a21=0 and geometrically reducible/nonreduced loci remain open.
source: root/a21-transverse-three-pole-closure-2026-07-28
audit: coordinate-first-audit-2026-07-28 (independent local-branch derivation; quotient, pole-pullback, and final-text line audit)
depends_on:
  - THM-2719-full-split-odd-faber-generic-normalization-genus-four-hundred-nineteen
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2718-split-prime23-five-pole-rational-primitive-closure
---

# THM-2726 -- transverse integral split responses have three forbidden poles

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2719 used the three branches at infinity to compute generic genus.  The
same branches carry a uniform divisor obstruction: the physical affine
coefficient `q` has a pole on each one.  This removes every geometrically
integral `a_21!=0` response member, including exceptional genus-drop members.

## 1. Statement

Use THM-2719's full chosen-sheet split degree-twenty-two response family

```text
F_23=Phi_22+sum_j a_j h^(22-j)Phi_j-lambda h^23,
G_24=Psi_22+sum_j a_j h^(22-j)Psi_j-W h^24            (1)
```

in `P(1,2,3,4)_[h,d,q,s]`.  Let `C_a` be any geometrically integral member
and assume

```text
a_21!=0.                                               (2)
```

Then `C_a` supports no physical polynomial Keller trajectory arising from a
split polynomial exact-square prefix and the displayed constant-coefficient
degree-22 Faber gauge.

No condition is imposed on `lambda` or `W`, and no generic smoothness or
genus assumption is used.  The conclusion is uniform over every
geometrically integral specialization satisfying `(2)`.

## 2. Three coarse infinity branches

The only possible infinity singularity is

```text
P_infty=[h:d:q:s]=[0:1:0:0].                          (3)
```

On the `d=1` index cover the local involution is

```text
(h,q,s) |-> (-h,-q,s).                                (4)
```

The transverse coefficient is the nonzero constant

```text
Phi_21(1,0,0)=88179/131072.                            (5)
```

Hence `(2)` lets `F_23=0` solve analytically for `h`.  The ordinary
order-six faces of the two equations are

```text
Phi_(22,6)=-(231/128)q s(q^2-3s^2)(3q^2-s^2),

Psi_(22,6)=-(231/256)(q^2-s^2)
                         (q^4-14q^2s^2+s^4).          (6)
```

After solving for `h`, every `h`-dependent contribution to `G_24` has
ordinary order at least seven, so its order-six face remains exactly the
second polynomial in `(6)`.

The six tangent lines of the second face have

```text
q=alpha s,                 alpha^2 in
{1, 7+4sqrt(3), 7-4sqrt(3)}.                          (7)
```

None of the three values in `(7)` belongs to `{0,3,1/3}`, so the first face
in `(6)` is nonzero on every one of those lines.  Equivalently the two faces
are coprime; their resultant is `-2^24` after removing their displayed
nonzero scalars.  It follows on every lifted branch that

```text
ord_s(q)=1,                       ord_s(h)=6.          (8)
```

The involution `(4)` freely pairs the six lines.  Thus the coarse
normalization has three distinct smooth branches, with invariant equations

```text
Q=q^2=c s^2+O(s^3),
c in {1, 7+4sqrt(3), 7-4sqrt(3)}.                     (9)
```

The physical affine coefficient is the invariant function

```text
q_aff=q/h^3.                                          (10)
```

Indeed `(4)` changes both numerator and denominator by `-1`.  Equations
`(8)`--`(10)` give on each of the three coarse branches

```text
ord(q_aff)=1-3*6=-17.                                 (11)
```

Equivalently the invariant square `Q/h^6` has order `-34`.  Therefore the
normalization of `C_a` has at least three distinct pole points of `q_aff`.

## 3. Every physical map is nonconstant and surjective

A physical trajectory gives at the generic source point a rational map

```text
gamma_0:P1_x ---> C_a.                                (12)
```

Its generic image lies in the chart `h=1`, where `q_aff` is exactly the
source coefficient

```text
q_source=A_src/U.                                     (13)
```

The map cannot be constant.  If it were, the affine response coordinates
`d,q,s` would all be constant.  Every Faber coefficient `a_j` is already in
`C`, so the exact third observable `R_Q` would be constant.  THM-2723 instead
gives

```text
R_Q'=kappa/U!=0.                                      (14)
```

Because `C_a` is geometrically integral, `(12)` factors through its
normalization `nu:X->C_a`.  Properness extends the resulting rational map
across the omitted source points.  A nonconstant morphism between projective
integral curves is finite onto, hence surjective:

```text
gamma:P1_x -> X.                                      (15)
```

## 4. Pullback exceeds the source capacity

Let `O_1,O_2,O_3` be the three normalization points from `(9)`.  Each fibre
`gamma^(-1)(O_i)` is nonempty, and the three fibres are pairwise disjoint.
Pullback of `(11)` makes `q_source=gamma^*(q_aff)` have a pole at every point
in their union.  Thus `q_source` has at least three distinct pole points on
`P1_x`.

THM-2723 proves that the same function `A_src/U` has at most two pole points
on the source projective line.  This contradiction proves the theorem.

## 5. Scope and surviving boundary

The proof uses neither the generic genus `419` nor a nonzero value of the
first flux.  It removes the entire geometrically integral `a_21!=0` locus,
including singular and genus-drop specializations.

It does not remove `a_21=0`.  It also does not automatically remove a
reducible or nonreduced response member: its three infinity branches may lie
on different image components, and the component containing a physical
trajectory must be audited separately.  Non-exact-prefix source charts,
Faber gauges outside `(1)`, the remaining split branch, `JC(2)`, and `DC(2)`
remain open.

An independent hostile audit recomputed both order-six faces and their
stripped resultant, checked the index-two quotient and pole order `-17`, and
rederived nonconstancy, normalization lifting, surjectivity, and divisor
pullback.  It specifically verified that no condition on `lambda` or `W` is
used and certified this final text line by line.
