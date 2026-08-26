---
id: THM-4209
title: "Generic mixed-B K-zero exact-weight-nine planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-4147/4176/4192 + VERIFIED-EXACT + INDEPENDENTLY
  SOURCE-AUDITED. On the live exact-weight-nine b=d=0 reduced (2,3) seam at
  K=0 and Delta=5696/105, every genuinely mixed row with
  eta*zeta*(eta+zeta)!=0 is excluded for arbitrary Phi,Theta, without a
  critical-discriminant or projection-squarefreeness hypothesis. The direct
  anti-diagonal calculation is an independently verified redundancy check:
  THM-4176 had already excluded the whole zeta=-eta, eta!=0 wall, including
  Theta=-Delta at K=0. The genuinely mixed source resultant has residual
  length 21 and total critical length 25; its complete packet and pure cubic
  carrier force strict finite/full monodromy deficits. THM-4205 subsequently
  closes the Y-only critical walls and hence the whole exact-weight-nine K=0
  coefficient face. Entry, M>=10, JC(2), and DC(2) remain OPEN.
source: codex-pair-entry-jc-mixed-20260826
depends_on:
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
  - THM-4192-complete-p-only-k-zero-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
  - THM-4205-complete-exact-weight-nine-k-zero-planar-jacobian-exclusion
script: 04-computation/jc23_k0_generic_mixed_b_exclusion_thm4209.py
output: 05-knowledge/results/jc23_k0_generic_mixed_b_exclusion_thm4209.out
independent_audit_script: 04-computation/jc23_k0_generic_mixed_b_exclusion_independent_audit_thm4209.py
independent_audit_output: 05-knowledge/results/jc23_k0_generic_mixed_b_exclusion_independent_audit_thm4209.out
script_sha256: 1af42b8de216c56196a716e5e45b4bb143eeaa1d35c4faea3f55d1c8c6c79a5d
output_sha256: 21584a5339909af0ce7d96eedddae1df364f4cb2497e30bd51de217cc4147a11
independent_audit_script_sha256: acc0696e768fbc2480398a14f7af3bdd56b09d7f92a61ad53e9f9edb62bf65d0
independent_audit_output_sha256: 64cb1cad1357158864a1ca68a0bf8624de20f6ebc08ea27fab5274586ea9b7fc
hash_basis: raw LF bytes
primary_audit: >
  PASS. A rational (s,p) source calculation reconstructs the polynomial
  critical pair, its Hessian bridge, the generic p^6 R21 resultant, the
  directly recomputed anti R19 resultant, all endpoint units, both universal
  fibres, the Newton polygon and packet, the pure cubic carrier, and all four
  strict response inequalities. It also records the R17/R15 degree drops on
  the next repeated-tangent wall without promoting them.
independent_audit: >
  ACCEPT. A separate source-coordinate referee re-expands the complete source,
  reconstructs both normalized projections, verifies the pt^2 gradient
  change and typed Hessian bridge, recomputes all resultants from the
  specialized source pairs, audits source and normalized infinity, restores
  the four universal Morse points, and independently recovers both packets,
  carriers, and response budgets. Projection collisions are not assumed away.
---

# THM-4209 -- generic mixed-B K-zero exact-weight-nine planar Jacobian exclusion

**PROVED RELATIVE TO THM-4147/4176/4192 + VERIFIED-EXACT + INDEPENDENTLY
SOURCE-AUDITED; JC(2) AND DC(2) REMAIN OPEN.**

## 1. Exact statement and inheritance pass

Work over `C` in the live exact-`M=9`, `b=d=0` reduced `(2,3)` seam. Retain

```text
K=2848/45-(7/6)Delta.                                  (1)
```

Thus the wall `K=0` is exactly

```text
Delta=5696/105.                                       (2)
```

> **Theorem.** No normalized exact-weight-nine source on the following
> coefficient stratum is a nonautomorphic planar Keller pair:
>
> ```text
> G: eta*zeta*(eta+zeta)!=0, Phi and Theta arbitrary.   (3)
> ```

No critical-resultant discriminant, squarefree projected-root condition, or
separation condition at `T=-1/6` is imposed. Multiple roots of a critical
projection are allowed.

THM-4147 already excludes the critical-open part of this row. The advance
here is to replace projection genericity by the complete reduced source
critical scheme, closing every projection-collision wall inside `(3)`.
THM-4192 supplies the closest rational-source precedent on the `zeta=0`
wall.

The anti-diagonal `zeta=-eta`, `eta!=0` is not a second new stratum.
THM-4176 already excludes that entire repeated-top wall for every
`Phi,Theta` when `Delta!=0`; condition `(2)` makes `Delta` nonzero. The
anti-diagonal source calculations retained below are an independent
redundancy check, not a scope extension. MISTAKE-519 records the correction.

The inheritance pass is:

- closest proved mechanism: THM-4147's complete packet and prime-carrier
  finite/full response obstruction;
- canonical hostile: `eta+zeta=0`, where the normalized leading row and top
  boundary branch collide;
- corrected near miss: a repeated root of a resultant projection does not
  imply a non-Morse source critical point;
- least-used sidecar: the source ideal `(A,B)` together with the Hessian
  bridge, before eliminating `s`.

The live concept board was

```text
source ideal | projection resultant | universal fibres | infinity rows
Newton packet | cubic carrier | response defect.       (4)
```

## 2. Complete source and polynomial critical pair

Put

```text
s=XT,                 p=T+s^2,                 t=p-s^2,
y=sp.                                                    (5)
```

On `(2)`, the complete source is

```text
G=-s^2/(2t)+H,

H=-3p+(8/3)p^2-(1376/135)p^3
  +Phi*s*p^3+(5696/105)*p^4+Theta*s^2*p^3
  +eta*s*p^4+zeta*s^3*p^3.                            (6)
```

There is no omitted exact-weight-nine monomial. Define

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.                                       (7)
```

Direct cancellation proves that `A,B` are polynomials and gives

```text
t^2 G_s=pA,
2t^2 G_p=t^2B-sA.                                    (8)
```

The gradient transformation has determinant `pt^2`. Differentiating `(8)`
on the source critical ideal gives the correctly typed bridge

```text
p det D_(s,p)(A,B)
 =2t^2 det Hess_(s,p)(G)                 mod (A,B).    (9)
```

There are no points hidden on the excluded source divisors:

```text
p=0: A=-s, B=-6;
t=0: A=-s, and the only candidate s=p=0 has B=-6.      (10)
```

In the generic mixed row, the `s`-degrees and leading rows are

```text
(deg_s A,deg_s B)=(6,3),
(LC_s A,LC_s B)=(3zeta p^2,9zeta p^2).                (11)
```

The corresponding anti rows are `(-3eta p^2,-9eta p^2)`. Thus no finite
`p!=0` source point is lost at `s=infinity` in either stratum.

## 3. Direct critical lengths without projection genericity

### 3.1 Generic mixed row

Eliminating `s` from the complete pair `(A,B)` gives identically

```text
Res_s(A,B)=p^6 R_21(p),                                (12)

R_21(0)=2^6*3^9*zeta^3,
[p^21]R_21=2^2*3^12*eta^3*zeta^2*(eta+zeta)^4.        (13)
```

Both endpoints are units under row `G` of `(3)`, for every `Phi,Theta`.
The residual source critical scheme therefore has length `21`, counted with
intersection multiplicity.

### 3.2 Anti-diagonal redundancy check

One must specialize the source pair first and recompute its resultant;
specializing only a frozen generic leading-degree determinant would be
invalid. Direct recomputation at `zeta=-eta` gives

```text
Res_s(A,B)=p^6 R_19(p),                                (14)

R_19(0)=-1,259,712 eta^3,
[p^19]R_19=1,327,104 eta^5(Delta+Theta)^4.             (15)
```

These are units when `eta*(Delta+Theta)!=0`. The residual source critical
length is `19`, with `Phi` arbitrary. This independently recovers the generic
row-A ledger of THM-4176 after imposing `K=0`; it is not needed for the new
theorem `(3)`.

### 3.3 Reducedness and the four universal points

For a hypothetical Keller realization, the inherited Keller-Hessian
congruence and `(9)` make every full source critical point Morse. Hence the
critical scheme is reduced even if several points have the same `p`- or
`T`-coordinate. No discriminant of `R_21`, `R_19`, or a normalized
projection is needed.

The rational source chart collapses two inherited coordinate fibres. In
normalized `(X,T)` coordinates they restore exactly

```text
T=0:    X^2=-6, G=0,   det Hess(G)=+6;
T=-1/6: X^2= 6, G=1/2, det Hess(G)=-6.                 (16)
```

Thus there are two Morse points on each fibre. The exact affine critical
lengths are

```text
L_G=21+2+2=25,
L_A=19+2+2=23.                                        (17)
```

The independent normalized audit finds degrees `(8,9)` and the factor
`T^56(6T+1)^2 Q_21` in row `G`, and degrees `(7,8)` with
`T^42(6T+1)^2 Q_19` in row `A`. Its infinity rows are respectively

```text
9T^8(eta+zeta),                 8T^7(Delta+Theta),     (18)
```

which independently recover the exact gates in `(3)`.

## 4. Generic mixed boundary packet

For `Q=q^-1`, clearing the denominator of `G=q` gives

```text
F_Q(s,p)=(s^2-p)(1-QH)-Q s^2/2.                       (19)
```

The lower Newton polygon has counterclockwise vertices

```text
(0,1),(2,0),(5,3),(1,5),(0,5),                        (20)
```

and Pick data

```text
(2Area,boundary,genus)=(31,11,11).                    (21)
```

Its complete nonvertical packet is

```text
e_G=(8,8,4,2,2,2,1),
n_G=sum e_G=27,
delta_G=sum(e-1)=20.                                  (22)
```

The four primitive faces are

```text
s^2(1-Q/2)-p,
s^2(1-Q/2-Q*zeta*s^3*p^3),
Q*p^3*s*(p-s^2)*(eta*p+zeta*s^2),
Q*p^5*(Delta+eta*s).                                  (23)
```

In particular, `Phi,Theta` are interior to this boundary calculation. The
finite critical scheme and THM-3827 as inherited through THM-4147 give
geometric connectedness. The Newton genus is at most `11`, while the Keller
residue identity and defect `20` force genus at least `11`. Equality makes
`(22)` complete: there is no hidden affine ramification or unrecorded genus
loss.

## 5. Pure cubic carrier and generic response contradiction

At `K=0`, the unique nonrational boundary equation is

```text
zeta W^3=q-1/2.                                       (24)
```

It is separable in characteristic zero. It is irreducible over `C(q)`
because the divisor `q-1/2` has valuation one at `q=1/2`, which cannot be the
valuation of a cube. Thus `(24)` is one prime cubic carrier. THM-4120/4122
and THM-4147's finite-separable-carrier transport apply.

If the carrier responds finitely, its three conjugate index-two branches are
removed from the origin response, so

```text
n=21, beta=3.                                         (25)
```

The maximum handle-plus-carrier merger capacity is

```text
2n-L_G-1+beta=42-25-1+3=19<n-1=20.                   (26)
```

If both handles are identities, the total index is only `beta=3<20`. If the
carrier responds fully at the origin, then `n=27` and the commutator bound is

```text
ind([X,Y])<=2(n-L_G)=2(27-25)=4<delta_G=20.            (27)
```

Both responses are impossible. This proves row `G` of `(3)` for every
`Phi,Theta`, including every critical-projection collision wall.

## 6. Anti-diagonal redundancy response

On `zeta=-eta`, the two raw top branches merge. THM-4147's local coordinate
calculation has tangent cone

```text
Q a(eta*a-(Delta+Theta)u).                             (28)
```

Under `eta(Delta+Theta)!=0`, this is one ordinary node. Its two normalized
branches have index seven, so the complete packet is

```text
e_A=(7,7,4,2,2,2,1),
n_A=25,
delta_A=18.                                           (29)
```

The same genus squeeze makes this packet complete. The carrier is now

```text
-eta W^3=q-1/2,                                       (30)
```

again one prime separable cubic. In the finite response,

```text
n=19, beta=3,
2n-L_A-1+beta=38-23-1+3=17<n-1=18.                   (31)
```

In the full response,

```text
ind([X,Y])<=2(n_A-L_A)=2(25-23)=4<delta_A=18.          (32)
```

The identity-handle case again has index only `3<18`. This reproduces the
`K=0`, `Theta!=-Delta` slice of THM-4176 by a separate source calculation.
It adds no coefficient scope to `(3)`.

## 7. Repeated-tangent audit and the true `K=0` frontier

The anti specialization

```text
zeta=-eta, eta!=0, Theta=-Delta                         (33)
```

is not a mere endpoint failure: the two tangent-cone factors in `(28)`
repeat. Direct source recomputation gives the exact degree staircase

```text
Phi!=0: Res_s(A,B)=p^6 R_17,
         [p^17]R_17=777,924 Phi^4 eta^5;

Phi=0:  residual degree 15,
         leading coefficient
         =229,431,851,352,064 eta^5/50,625.            (34)
```

These resultants are **VERIFIED-EXACT**. They are the `K=0` specializations
of THM-4176's complete row-B/row-C degree towers, whose normalized packets
are respectively

```text
(6,6,4,2,2,2,1),              (5,5,4,2,2,2,1).
```

Thus `(33)` is already excluded by THM-4176. A proof confined to this file
would still have to normalize the repeated tangent before using `(34)` in a
monodromy budget; the dependency supplies exactly that missing step.

Other hostile specializations have the expected type changes:

```text
eta=0:      residual degree 20;
zeta=0:     residual degree 18;
eta=zeta=0: exits exact weight nine.                   (35)
```

The `zeta=0`, `eta!=0` specialization is already closed by THM-4192, while
`eta=zeta=0` exits exact weight nine. Consequently the genuinely unresolved
`K=0` coefficient locus is confined to the Y-only row

```text
eta=0,                         zeta!=0.                 (36)
```

THM-4147 closes its critical-open part. THM-4205 subsequently uses the
source-ideal specialization of this audit to close every residual critical
wall in `(36)`. None of these specializations is silently inferred from the
generic resultant.

## 8. Audit and replay

The primary source audit builds `(6)--(15)` symbolically, verifies `(9)`,
reconstructs the universal fibres, Newton hull, packet, pure carrier, response
budgets, and both degree drops in `(34)`. The independent referee expands the
source separately, uses both rational and normalized projections, and
recomputes the anti resultant only after specialization. It verifies source
infinity and normalized infinity independently and never assumes projected
squarefreeness.

Replay from the repository root:

```bash
python3 -B 04-computation/jc23_k0_generic_mixed_b_exclusion_thm4209.py
python3 -B \
  04-computation/jc23_k0_generic_mixed_b_exclusion_independent_audit_thm4209.py
```

Compare stdout with the corresponding frozen outputs. Normal and `python3
-O -B` primary replays byte-match.

## 9. Strict scope

This theorem newly removes the critical-open hypothesis only on the genuinely
mixed `K=0` stratum `(3)`. Its anti-diagonal computations are independently
verified redundancy for THM-4176. Together, THM-4192 closes the P-only row,
THM-4176 closes the anti-diagonal row, this theorem closes the genuinely mixed
row, and THM-4147 closes the critical-open Y-only row. THM-4205 subsequently
closes all remaining Y-only critical walls, hence the whole exact-weight-nine
`K=0` coefficient face. Entry into the reduced seam, `M>=10`, JC(2), and
DC(2) remain open.
