---
id: THM-3794
title: "Constant-unit affine surfaces have no quadratic etale plane map"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  An integral affine surface
  over an algebraically closed field of characteristic zero whose only units
  are constants admits no dominant etale morphism to A2 of generic degree
  two.  The proof passes to the finite normalization, uses normal-surface
  Cohen--Macaulayness to obtain a flat double cover, and turns its branch
  polynomial into a forbidden nonconstant unit on the original surface.
  The constant-unit hypothesis is sharp.  Consequently every Darboux map on
  the THM-3783 and THM-3785 factor surfaces has surface degree at least three;
  their pulled-back planar Keller degrees are at least six and nine,
  respectively.
source: root / quadratic-normalization branch-unit lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-23).  The audit
  independently rederived the normalization inclusion and normalization-form
  Zariski Main step; height and Cohen--Macaulay hypotheses for miracle
  flatness; the trace-zero rank-one summand and Cayley--Hamilton quadratic
  law; the entire non-etale fibre over h=0; and the branch-unit contradiction.
  It also verified both degree-floor applications and the punctured squaring
  hostile.  The normalization/open-immersion direction,
  local miracle-flatness hypotheses, trace-zero line summand, quadratic
  algebra law, entire ramification-fibre deletion, and unit pullback were
  checked separately.  The punctured squaring map is an exact hostile control
  showing that nonconstant units are the missing coordinate.
depends_on: []
related:
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
external:
  - "Stacks Project, Tags 05K0 (Zariski Main) and 00R4 (miracle flatness)."
---

# THM-3794 -- constant units forbid quadratic etale maps to the plane

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be algebraically
closed of characteristic zero, let `X=Spec B` be an integral affine surface
with

```text
B*=k*,                                                               (1)
```

and let

```text
f=(A,C):X -> A2=Spec R,                 R=k[A,C],                     (2)
```

be a dominant etale morphism.  Then

```text
[Frac(B):Frac(R)] != 2.                                               (3)
```

This is a statement about the intermediate affine surface, not only about a
polynomial source plane.  Its load-bearing object is the branch equation of
the finite normalization that Zariski Main places around a quasi-finite map.

## 1. Put the missing sheets back

Assume for contradiction that the degree in `(3)` is two.  Let `K=Frac(B)`
and let `S` be the integral closure of `R` in `K`.  Since `R` is a finitely
generated algebra over a field, the normalization

```text
pi:Xbar=Spec S -> A2                                             (4)
```

is finite, integral, and generically of degree two.

Every element of `S` is integral over `R`, hence integral over `B` after the
inclusions `R subset B subset K`.  Etaleness makes `B` normal, so `S subset
B`.  The induced map

```text
j:X -> Xbar,                    f=pi o j                           (5)
```

is birational.  It is quasi-finite because `pi` is finite and `f` is
quasi-finite.  Zariski's Main Theorem, in its normalization form, therefore
identifies `j` with an open immersion.  Thus `Xbar` is the finite two-sheeted
completion of the sheets that may escape from `X`.

## 2. The completion is a flat quadratic algebra

A normal Noetherian surface is Cohen--Macaulay.  At every point of `Xbar`,
the finite local map in `(4)` has regular target, Cohen--Macaulay source,
equal local dimensions, and zero-dimensional fibre.  Miracle flatness makes
`pi` flat.  It is consequently finite locally free of constant rank two.

The trace splitting, valid because `2` is invertible, is

```text
S=R direct_sum M,                    M=ker(Tr_(S/R)),              (6)
```

where `M` is a rank-one projective `R`-module.  Since `Pic(A2)=0`, choose a
generator `u` of `M`.  Cayley--Hamilton for the trace-zero element `u` gives

```text
S=R[u]/(u^2-h),                      h in R.                       (7)
```

The ring `S` is a domain of field degree two, so `h` is nonzero and is not a
square in `Frac(R)`.

## 3. Ramification becomes a forbidden unit

At every point of `Xbar` above `h=0`, equation `(7)` has `u=0` and

```text
Omega_(S/R)=S du/(2u du)                                             (8)
```

is nonzero.  Hence `pi` is not etale at any point of

```text
pi^(-1)V(h)=V(u).                                                     (9)
```

But `j` is an open immersion and `f=pi o j` is etale.  Therefore

```text
j(X) intersection V(u)=empty.                                      (10)
```

Equivalently, the regular function

```text
h(A,C)=u^2 in B                                                     (11)
```

has no zero on the affine variety `X`.  It is therefore a unit.  By `(1)`
it lies in `k*`.  Dominance makes `R -> B` injective, so already `h in k*`
inside `R`.  Since `k` is algebraically closed, `u^2-h` then splits and
`S` is not a domain, contradicting its construction.  This proves `(3)`.

The mechanism can be summarized as

```text
quadratic missing sheet
 -> finite flat double cover
 -> one branch divisor and no unramified companion sheet above it
 -> the whole branch fibre must be deleted
 -> branch polynomial becomes a unit
 -> constant-unit contradiction.                                  (12)
```

The third arrow is special to degree two.  In higher degree a branch fibre
may also contain unramified sheets, so its discriminant can still vanish on
the etale open and the same one-line unit argument does not apply.

## 4. Sharp hostile boundary

The constant-unit hypothesis cannot be dropped.  On

```text
X_*=Spec k[u,u^(-1),v]=G_m x A1
```

the map

```text
(u,v) |-> (s,t)=(u^2,v) : X_* -> A2_(s,t)                         (13)
```

is dominant, etale, and generically of degree two.  Its finite completion is
`A2_(u,v) -> A2_(s,t)`; the branch line `u=0` has been removed, and the branch
polynomial `s` pulls back to the nonconstant unit `u^2`.  This is exactly the
failure forbidden by `(1)`.

## 5. Counterexample-design consequences

Any affine factor surface whose coordinate ring embeds in a polynomial ring
`k[x,y]` automatically has only constant units: a unit and its inverse would
both be source polynomials.  Therefore a factorization of a planar Keller map

```text
A2 -> X -> A2                                                     (14)
```

can never have second arrow of generic degree two.

For the THM-3783 quadratic factor surface, the already proved birational
exit and `(3)` force every Darboux map to have surface degree at least three.
Its degree-two source atlas then yields a planar Keller map of degree at least
`2*3=6`.  For the THM-3785 cubic pseudo-plane, the same argument raises the
surface degree to at least three and the pulled-back planar degree to at least
`3*3=9`.  In particular the residual-divisor lane of THM-3790 cannot close at
surface degree two, regardless of the topology of its residual curve.

The scheme-theoretic inputs are the normalization form of
[Zariski's Main Theorem](https://stacks.math.columbia.edu/tag/05K0) and
[miracle flatness](https://stacks.math.columbia.edu/tag/00R4).  **QED.**
