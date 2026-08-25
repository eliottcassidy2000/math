---
id: THM-3675
title: "Russell-cylinder critical-fold formal-conjugacy closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For either
  THM-3642 zero-debt collision polynomial Q_6 or Q_*, and every nonzero
  critical vertical displacement H in t^2 C[t], no arbitrary target two-form
  has nonzero constant pullback on q=Q(x)+H(t), w=t.  A formal target/source
  coordinate conjugates H to its leading monomial.  The transformed constant
  Jacobian is an x-independent series; THM-3673's retained cokernel kills all
  of that series except the nonzero J_0 debt, giving a contradiction.  Hence
  no target pair exists in these families.  H'(0)!=0, other Q, and
  nonordinary tangent collisions remain open.
source: kps-s194 / THM-3673 formal-reparameterization continuation, 2026-08-21
audit: >
  PASS -- kps-s195 independently checked formal monomialization and inverse
  series, the simultaneous source/target coordinate change, completed-ring
  two-form typing, transformation of a constant source density to the
  x-independent series kappa*psi'(u), and annihilation of every positive
  coefficient row except the nonzero retained J_0 debt.  Normal and
  optimized companions returned all 50 gates with the pinned hashes.  No
  correction was required.
depends_on:
  - THM-3673-russell-cylinder-monomial-ramification-debt-dilation
related:
  - THM-3623-russell-cylinder-even-general-vertical-fold-all-order-closure
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3641-russell-cylinder-principal-noneven-curvature-debt-boundary
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
script: 04-computation/jc2_russell_cylinder_critical_fold_formal_conjugacy_thm3675.py
output: 05-knowledge/results/jc2_russell_cylinder_critical_fold_formal_conjugacy_thm3675.out
script_sha256: 543dd998a41924e16c1f6d1cf1439432fe050cef84e7c9301ecd25fb977bfd00
output_sha256: 699d1494ad63078014daec6e9efd03244a9f5afbd28435a0d6738728af5115f8
hash_basis: raw LF bytes
---

# THM-3675 -- every critical vertical fold over Q_6 and Q_* is closed

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The
monomial obstruction in THM-3673 survives arbitrary higher terms of the
vertical displacement.  The reason is not that those terms disappear: they
are absorbed by a formal coordinate change, after which a hypothetical
constant Jacobian becomes a stable-only series.  The retained debt row is
blind to stable-only variation but not to its nonzero constant term.

All rings, germs, and two-forms are over `C`.

## 0. Statement

Fix either collision polynomial

```text
Q_6=Q_1-(259/36)x^2(x^2-1)^2,

Q_*= -x^7-(27/4)x^6+3x^5+18x^4-3x^3
       -(27/2)x^2+x-3/4,                                (1)
```

with `Q_1` as in THM-3642.  Let

```text
0!=H(t) in t^2 C[t],
E_(Q,H)(x,t)=(x,Q(x)+H(t),w=t).                         (2)
```

For every target two-form `Omega`, including every decomposable form
`dF^# wedge dG^#` from an actual target pair,

```text
E_(Q,H)^*Omega != kappa dx^dt             for kappa in C*. (3)
```

Consequently no target pair on either fold has nonzero constant source
Jacobian.  Any such pair would have identified the retained triple and hence
would have been a planar Jacobian-conjecture counterexample; `(3)` excludes
these two critical-fold families.

## 1. Formal monomialization is target-compatible

Write

```text
k=ord_0(H)>=2,                alpha=[t^k]H!=0.           (4)
```

The series

```text
phi(t)=t (H(t)/(alpha t^k))^(1/k)=t+O(t^2)              (5)
```

is a formal coordinate.  Let `t=psi(u)` be its inverse.  Then

```text
psi(0)=0,              psi'(0)=1,
H(psi(u))=alpha u^k.                                  (6)
```

This is simultaneously a source and target change.  On the completed target
ring replace the stable coordinate `w` by

```text
W=phi(w).                                               (7)
```

Transporting an arbitrary target two-form through `(7)` gives another
arbitrary formal target two-form.  Along the source, `w=t=psi(u)` and hence
`W=u`; the transformed fold is exactly

```text
q=Q(x)+alpha u^k,             W=u.                      (8)
```

THM-3673 applies to formal two-form jets in the completed regular target
coordinates, so no polynomiality of `phi` or `psi` is needed.  An actual
global polynomial form maps into this formal universe, which is enough for
an obstruction.

## 2. What a constant Jacobian becomes

Assume for contradiction that the pullback in `(3)` equals
`kappa dx^dt`.  Under `t=psi(u)`, it becomes

```text
kappa psi'(u) dx^du=sum_(n>=0) u^n K_n(x) dx^du,

K_n(x)=kappa [u^n]psi'(u).                              (9)
```

Every `K_n` is constant in `x`, and `K_0=kappa` by `(6)`.  Therefore

```text
K_n'=K_n''=0,
Lambda(K_n)=0,

Lambda(P)=(5P(-1)-18P(0)+13P(1))/18.                  (10)
```

Apply the THM-3673 dilated identities to `(8)`.  In both identities the
`K_k` value block also kills constants:

```text
Q_6:   2012/2187-2012/2187=0,
Q_*:   4/9-4/9=0.                                      (11)
```

Thus the left side `Lambda(K_(2k))`, every derivative term, and the entire
`K_k` block vanish.  The remaining `K_0` block gives respectively

```text
0=alpha^2 kappa (365888/6561),
0=alpha^2 kappa (5440/81).                              (12)
```

Both are impossible because `alpha,kappa!=0`.  This proves `(3)`.

## 3. Relation to the even-fold theorem

THM-3623 already closes every critical vertical displacement when `Q` is
even, using an all-order moving-triple recurrence.  Neither `Q_6` nor `Q_*`
is even.  The present theorem reaches these two non-even folds by a different
mechanism: a finite retained debt identity, character dilation, and formal
conjugacy.  It therefore does not subsume THM-3623 and does not assume its
side parity.

## 4. Strict scope

The theorem closes exactly

```text
Q in {Q_6,Q_*},              0!=H in t^2 C[t].          (13)
```

It does not cover

- `H'(0)!=0`, including THM-3629's genuinely mixed nonlinear pairs;
- another ordinary zero-debt polynomial on a THM-3641 curvature plane;
- a nonordinary tangent-collision slope packet;
- implicit source planes not of the form `(2)`; or
- arbitrary polynomial maps of `A^2`.

No Keller pair is constructed.  `JC(2)` remains **OPEN**.  The sharp next
ordinary-fold target is to derive the fourth-debt functional over a Hermite
parameter family and test whether its zero set is empty, while the sharp
orthogonal target remains the `H'(0)!=0` mixed-pair/Darboux boundary.

## 5. Exact controls

The companion constructs the monomializing inverse `psi` recursively over
`Q` for three deliberately nonmonomial controls of orders `2,3,4`, verifies
`H(psi(u))=alpha u^k` through the order needed by THM-3673, and evaluates the
constant-series cancellation and surviving debt exactly.  These are controls
for `(5)--(12)`; the all-`H` proof is formal and not replaced by the three
examples.

```bash
python3 04-computation/jc2_russell_cylinder_critical_fold_formal_conjugacy_thm3675.py
python3 -O 04-computation/jc2_russell_cylinder_critical_fold_formal_conjugacy_thm3675.py
```

Normal and optimized transcripts are byte-identical to the stored output.
The companion reports `50` active gates and contains no Python assertion
statements.
