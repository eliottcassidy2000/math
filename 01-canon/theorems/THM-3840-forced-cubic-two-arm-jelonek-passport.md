---
id: THM-3840
title: "A nonlinear cubic plane atlas forces an origin/branch Jelonek pair"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Every polynomial Keller atlas of the THM-3811 surface has the
  target origin and at least one of three explicit smooth cubic-branch values
  in its nonproper-value set.  THM-3836 supplies both Laurent arms; the row
  unit k is nonconstant on each etale source component, and a pole on its
  smooth completion gives the two escaping valuations.  The branch
  normalization parameter and surviving companion root are exact.
source: jc_quartic_c3_construct / cubic two-arm nonproper-value lane, 2026-08-23
audit: >
  PROVISIONAL EXACT CANDIDATE.  The companion verifies all slope resultants,
  both arm limits, discriminant incidence and smoothness, exclusion of triple
  roots, branch-normalization values, the companion Vieta packet, and
  pairwise distinctness of the three candidates.  The curve-completion and
  valuative nonproperness argument is human proof.  Normal and optimized runs
  byte-match the frozen transcript; independent hostile audit remains.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
  - THM-3836-cubic-factor-cofactor-darboux-packet
related:
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
  - THM-3578-zariski-main-boundary-rank-and-sheet-debt
script: 04-computation/jc2_forced_cubic_two_arm_jelonek_passport_thm3840.py
output: 05-knowledge/results/jc2_forced_cubic_two_arm_jelonek_passport_thm3840.out
script_sha256: 82b821f1201e2ec2acdd6d77ce92a74b29c7eae5a3868603b25d0dea46baea39
output_sha256: 99a71f761cc2eee49b26d0c468974de40d33c2d331509876626acdb76d448310
semantic_sha256: fb59cb14b46a4edafbe622a82ab9c7d42b5af615d160a4b249ddbbd28f5017b9
hash_basis: raw LF bytes
---

# THM-3840 -- the forced two-arm fibre gives a Jelonek pair

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**  Work over an algebraically closed field `K` of characteristic zero.
Let

```text
psi:A2_(x,y) -> U                                                (1)
```

be a dominant etale plane atlas of the THM-3811 nonlinear cubic surface, and
write

```text
F=(A,C):A2 -> A2,                 Jac(A,C)=lambda in K*.         (2)
```

Let `J_F` denote the nonproper-value, or Jelonek, set of `F`.  There is a
root `a` of

```text
r(Z)=3Z^3+7Z^2+1                                                (3)
```

such that

```text
O=(0,0) in J_F,
T_a=(-2/[a(7a^2+3)], -1/a^2) in J_F.                            (4)
```

The point `T_a` is a smooth nontriple point of the cubic discriminant.  Its
ramified double root has normalization parameter `q=1/a`, while the surviving
companion sheet has marked-root ratio

```text
beta=2a/(7a^2-1).                                                (5)
```

The three possible points `T_a` are pairwise distinct and all differ from
`O`.  The theorem forces `O` and at least **one** of them; it does not force
all three or identify the whole Jelonek curve.

## 1. The two source curves carry a nonconstant unit

THM-3836 supplies a root `a` of `(3)` for which the pulled-back cubic pencil
member `h-ak=0` has nonempty comaximal minus and plus factors.  Choose an
irreducible source curve `E^-` on the minus side and `E^+` on the plus side.
THM-3831 identifies their intrinsic targets as

```text
U_a^-: C=0,
U_a^+: C=-1/a^2,                                                (6)
```

with each arm isomorphic to `G_m` using the intrinsic function `k`.

The restriction of the etale morphism `psi` to either source component is
again etale over the corresponding smooth arm.  A nonempty etale morphism of
curves has open one-dimensional image.  Consequently

```text
k|_(E^-), k|_(E^+) are nonconstant units.                        (7)
```

This proves the nonconstancy needed below rather than assuming that every
source component surjects set-theoretically onto the full arm.

## 2. Poles of the row unit are source-infinite valuations

Let `bar(E)^+` be the smooth projective completion of the normalization of
either source curve.  A nonconstant rational function on a projective curve
has a pole.  Since `k` is a unit on the affine curve, every zero and pole lies
in the boundary `bar(E)^+ minus E`.

Such a pole is genuinely source-infinite.  If the corresponding DVR map
extended to a point of `A2_(x,y)`, every source polynomial, in particular the
polynomial `k`, would be regular there, contradicting its negative valuation.
Thus it is a valid valuative witness to failure of properness of `F`.

Put

```text
Q(a)=7a^2+3.                                                     (8)
```

The exact arm formulas in THM-3831 are

```text
on E^-: (A,C)=(a/[Q(a)k], 0),
on E^+: (A,C)=(a/[Q(a)k]-2/[aQ(a)], -1/a^2).                    (9)
```

At a pole of `k`, the local function `1/k` tends to zero.  The two target
limits in `(9)` are therefore exactly `O` and `T_a`.  This proves `(4)` by
the valuative criterion for properness.  Notice that `aQ(a)!=0`, so every
displayed address is defined, and `T_a!=O` because its second coordinate is
nonzero.

## 3. The plus endpoint is a smooth simple branch value

The THM-3811 cubic discriminant is

```text
Delta=
 A(C+5A)(4C+19A)(3C-17A)
 +C^2(162A^3+126A^2C-4C^3)-27A^2C^4.                           (10)
```

Exact substitution of `(4)` and reduction modulo `r(a)` gives

```text
Delta(T_a)=0,                 Delta_A(T_a)=0,
Delta_C(T_a)!=0.                                                   (11)
```

For the last inequality, the numerator resultant with `r(a)` is

```text
27703695576000!=0.                                               (12)
```

Thus `T_a` is a smooth point of the branch curve.  The branch normalization
from THM-3811 is

```text
R(q)=(q-3)(q+1)(q+2),
A(q)=-2q^2R(q)/(3q^2+7)^2,
C(q)=-qR(q)/(3q^2+7).                                           (13)
```

Substituting `q=1/a` yields exactly `T_a`.  A triple root would require
`q^2=7/3`, but

```text
Res_a(r(a),7a^2-3)=5245!=0.                                    (14)
```

Hence the branch cubic has one double and one simple companion root.  Vieta's
relations give the companion `q`-coordinate and marked-root ratio

```text
q_comp=(7a^2-1)/(2a),              beta=1/q_comp,                (15)
```

which is `(5)`; its denominator is nonzero because the corresponding
resultant is `1363`.

Finally, equality `T_a=T_b` would imply `a^2=b^2` from the second coordinate.
Distinct roots would then satisfy `b=-a`, but

```text
Res_Z(r(Z),r(-Z))=216!=0.                                       (16)
```

Thus the three candidates are pairwise distinct.

The theorem ties the mandatory bichromatic fibre to the precise `S_3`
branch/companion-sheet transition.  It supplies two forced nonproper values,
not a polynomial atlas, not the full nonproper curve, and not a Jacobian
counterexample.  **QED, pending independent hostile audit.**
