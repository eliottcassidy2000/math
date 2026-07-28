---
title: "LRC14 origin-return, endpoint deconvolution, and rainbow-thinning pipeline"
date: 2026-07-28
status: PROOF INTERFACE / OPEN PHYSICAL SERVICE; all theorem inputs are separately scoped below
source: root-2026-07-28-origin-return-rainbow-synthesis
---

# The algebraic endgame is assembled; the physical table is not

The incoming THM-2642/2644/2645/2647/2648 chain changes the last-lane LRC
question.  The remaining problem is not whether a saturated eleven-sheet
pair contains charged information or thin witnesses.  It does.  The missing
object is one lawful nonnegative same-base transition table on which origin,
composition, and edge restriction coexist.

## 1. Inheritance pass

- **Closest proved mechanism:** THM-2644.  A nonnegative weight on an odd
  torsor is decoded by two different quadratic rungs: Gram purity and an
  oriented common-middle return.
- **Canonical hostile:** THM-2642.  Every affine clutch has positive support,
  so support alone cannot select a carry or origin.
- **Closest physical selector:** THM-2640.  The predecessor carry gives a
  private singleton root in every refined cell, but it is target-neutral;
  covariance would require the unsupported clutch `c->c+7 delta`.
- **Corrected near miss:** THM-2645.  Exact representation multiplicity has
  all twelve charged `C13` characters, but retains a free common-origin
  gauge and is not a physical transition table.
- **Least-used sidecar:** THM-2647.  Once one absolute labelled endpoint
  two-set is supplied, the other endpoint is uniquely reconstructed by the
  signed inverse of `I+T_d`; its unavoidable total variation is `13/2`.
- **New thin atlas:** THM-2648.  Two rainbow matchings cover all thirteen
  carries, minimally, but the matching restriction is not yet a lawful LRC
  event.

No item in this list lowers the `165`-row ledger by itself.

## 2. The exact proposed composition

The desired service has three typed rungs.

```text
physical same-base table
  -> nonnegative common-origin weight mu on C13
  -> absolute endpoint pair (A,B)
  -> two rainbow edge charts covering every carry.             (1)
```

The theorem stack would act as follows.

1. Compute `M=sum mu`, `E=sum mu^2`, `delta=M^2-E`, and the **oriented**
   return `R=(mu*mu)(0)` through one common physical middle fibre.
   `E=M^2>0` makes `mu` a singleton; `R>delta` then forces that singleton
   to be the identity by THM-2644 and oddness of `C13`.
2. The selected origin turns one relative two-hole endpoint into an absolute
   labelled set `A`.  From the full multiplicity `r`, THM-2647 gives

   ```text
   1_B=(I+T_d)^(-1) T_(-a0)[r-9*1],
   ||(I+T_d)^(-1)||_1=13/2.                                (2)
   ```

3. With `(A,B)` fixed, THM-2648 supplies two eleven-edge rainbow matchings
   whose two-point carry holes are disjoint.  Their union has multiplicity
   one on four carries and two on nine, hence retains all charged colours.

This is an implication, not a construction.  The first arrow must preserve
nonnegativity, a common middle, and the physical endpoint gauge.  The last
arrow must make restriction to the selected matching edges measurable on
the same physical object.  Current coefficient fibres supply neither.

## 3. Why the modular `2/3` sidecar cannot choose the `13`-origin

THM-2646 proves that the missing coordinate over modular `C6` is an integral
braid height.  That is the correct model of a sidecar, but it cannot be
imported abstractly into (1):

```text
Hom(C13,C6)=0,                 Hom(C13,Z)=0.                (3)
```

Therefore no translation-equivariant homomorphic `C6` colour or integer
height separates a free `C13` origin torsor.  A mixed `2/3/13` coordinate
must come from a literal common ancestor, carry, or endpoint action; CRT
cardinality alone does nothing.  This is the sharp interaction between the
user's two special free factors and the live prime thirteen: `2` and `3`
organize the modular quotient, while `13` remains an orthogonal origin
gauge until a physical map couples them.

## 4. Matching, not tournament

THM-2648's two charts are matchings with symmetric carry multiplicity.  The
right external analogy is inherited vertex colouring of perfect matchings:
endpoint labels survive a matching contraction only when they are carried
with the edges.  An uncoloured union retains support but loses allocation.
There is no intrinsic antisymmetric pair observable here, so orienting the
two charts as a tournament would discard rather than recover information.

The terminology is from Krenn--Gu--Soltesz,
[*Questions on the Structure of Perfect Matchings inspired by Quantum
Physics*](https://arxiv.org/abs/1902.06023); only the label-retention analogy
is used here.

The role of colouring is only explanatory.  It does not make the abstract
matching graph a measurable LRC subset.

## 5. Cheapest decisive test

For one exact surviving row, refine a THM-2549 same-base chronology cell by
the physical predecessor carry of THM-2640, and then:

1. enumerate the induced common-origin histogram `mu(u)` on `C13`, retaining
   the same carry observer at both adjacent clocks;
2. print exact `M,E,delta,R`, retaining the common middle and orientation;
3. test `E=M^2` and `R>delta` before any Fourier scalarization;
4. verify that the selected root commutes with the lawful target action,
   rather than imposing the missing `c->c+7 delta` clutch;
5. if the branch lands, reconstruct `(A,B)` by (2) and test whether both
   THM-2648 edge indicators are events in the original sigma-algebra;
6. if it fails, preserve the first failed predicate: impurity, absent return,
   gauge mismatch, or nonmeasurable edge restriction.

A positive result would compose existing theorems into a genuine carry
cover.  A negative result would identify the exact missing physical
coordinate.  Recomputing support, charged energy, or an unlabelled rainbow
atlas cannot decide this test.
