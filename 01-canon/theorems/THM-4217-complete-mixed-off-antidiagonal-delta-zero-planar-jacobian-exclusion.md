---
id: THM-4217
title: "Complete mixed off-antidiagonal Delta-zero planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4147/4159
  + VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED. In the inherited b=d=0
  reduced (2,3) seam, the complete exact-M=9 coefficient locus Delta=0 and
  eta*zeta*(eta+zeta)!=0 contains no nonautomorphic planar Keller pair. Its
  five exhaustive critical lengths are 25,24,23,22,21, and every row has
  genus 11 and packet (8,8,4,2,2,2,1). JC(2), DC(2), entry, other cells,
  the remaining coefficient walls, and M>=10 remain OPEN.
source: codex-planar-jacobian-delta-zero-session-20260826
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
related:
  - THM-4180-repeated-top-delta-zero-planar-jacobian-exclusion
  - THM-4183-p-only-delta-zero-planar-jacobian-exclusion
  - THM-4209-generic-mixed-b-k-zero-exact-weight-nine-planar-jacobian-exclusion
script: 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217.py
output: 05-knowledge/results/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217.out
script_sha256: d4a6e8594b93842ede059017857e5f40d4a70da26d4a7de08084da71a6a2cd61
output_sha256: 9822ca6c47aede26ab2088bb6c60eba958971367d373351329edcbd98c74ec99
independent_audit_script: 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217_independent_audit.out
independent_audit_script_sha256: 0d8e64906574fd92dbbdde0891fbc55941359c7c423e85022dcf4a88e1a410ca
independent_audit_output_sha256: 93ca99a16230cb608366322e0df8ab79c8b057a217b4493c97db29c572115f07
---

# THM-4217 -- complete mixed off-antidiagonal Delta-zero planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4147/4159
+ VERIFIED-EXACT + INDEPENDENTLY SOURCE-AUDITED; JC(2) REMAINS OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam. Put

```text
K0=2848/45.
```

> **Theorem.** The complete exact-weight-nine coefficient locus
>
> ```text
> Delta=0,                 eta*zeta*(eta+zeta)!=0            (1)
> ```
>
> contains no nonautomorphic planar Keller pair.

The words *complete coefficient locus* refer only to `(1)` inside the named
reduced seam. They do not assert entry into that seam or the planar Jacobian
conjecture.

The inheritance pass is:

- the closest proved mechanism is THM-4147's exact critical-length, labelled
  packet, prime-carrier, and response obstruction;
- the canonical hostile is the `Delta=0` face contraction: the old
  `Delta`-bearing face becomes a monomial and cannot simply retain its old
  label;
- the corrected near miss is MISTAKE-519: repeated-top and `K=0` rows are
  already closed by THM-4176/4205, so neither was the open frontier here;
- the least-used sidecar is the replacement diagonal variable `W=sp`, paired
  with a source-first resultant before any projected-root test.

THM-3992/3997 and THM-4007 supply the complete source and force

```text
K=2848/45-(7/6)Delta=K0.                                  (2)
```

No critical-resultant discriminant, squarefree projected-root hypothesis, or
condition on `Phi` is assumed.

## 2. Complete source and the five critical rows

Use

```text
s=XT,             p=T+s^2,             t=p-s^2=T,
P=T+X^2T^2,                              Y=XTP,

G=-s^2/(2t)+H,

H=-3p+(8/3)p^2-(1376/135)p^3+K0s^2p^2
  +Phi sp^3+Theta s^2p^3+eta sp^4+zeta s^3p^3.          (3)
```

This is the complete exact-`M=9` source on `(1)`. Define

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.                                         (4)
```

They are polynomials and direct differentiation gives

```text
t^2G_s=pA,             2t^2G_p=t^2B-sA.                 (5)
```

Moreover,

```text
(deg_s A,deg_s B)=(6,3),
(LC_s A,LC_s B)=(3zeta p^2,9zeta p^2).                 (6)
```

Thus `zeta!=0` prevents a critical intersection at `s=infinity` for finite
nonzero `p`. At `p=0`, one has `A=-s,B=-6`; at `t=0`, `A=-s` and the only
candidate again has `B=-6`. These are exactly the two coordinate fibres
treated separately in Section 3.

The source Hessian bridge is the ideal identity

```text
p det D(A,B) = 2t^2 det Hess_(s,p)(G)       modulo (A,B). (7)
```

Hence a hypothetical Keller realization makes every residual source-critical
point Morse. Repeated roots of a projection resultant may merge `p`-values,
but cannot lower the critical-scheme length.

The general source-first elimination is

```text
Res_s(A,B)=p^6 R_21(p),

R_21(0)=-46656 zeta D,
[p^21]R_21=2^2 3^12 eta^3 zeta^2(eta+zeta)^4,           (8)

D=4K0^2 Theta-27zeta^2.
```

This closes the `D!=0` row. If `D=0`, condition `(1)` gives the unique chart

```text
zeta=(2K0/3)u=(5696/135)u,       Theta=3u^2,       u!=0. (9)
```

There is no sign quotient. Put

```text
J=8544Phi-1215u^3-22784u,

E=2460375u^4-204543360u^2+5580439552,
S=547499520eta+uE,

T0=27064125u^4-5739517440u^2+47239069696.             (10)
```

Each wall was substituted in `(A,B)` before elimination. The exact tower is

| coefficient row | source resultant | live bottom endpoint |
|---|---|---|
| `D!=0` | `p^6R_21` | `-46656zeta D` |
| `D=0,J!=0` | `p^7R_20` | `(8305770496/1125)u^2J` |
| `D=J=0,S!=0` | `p^8R_19` | `(2916352/16875)u^2S` |
| `D=J=S=0,T0!=0` | `p^9R_18` | `-(512/375)u^5T0` |
| `D=J=S=T0=0` | `p^10R_17` | first strict-transform row is a unit |

On `J=0` and `S=0`, respectively,

```text
Phi=u(1215u^2+22784)/8544,
eta=-uE/547499520.                                    (11)
```

The leading endpoint in the first three wall rows is

```text
(129777664/11390625)eta^3u^2(135eta+5696u)^4,          (12)
```

and on `S=0` it becomes

```text
-u^9 A0^4 E^3/
  11731773052904797284477266998693409587200000,

A0=18225u^4-1515136u^2-129777664.                     (13)
```

Here `E=0` is exactly the `eta=0` loss and `A0=0` is exactly the
`eta+zeta=0` repeated-top loss in the `S=0` chart. Both are outside `(1)`.

For the terminal row, `T0` is squarefree and exact Euclidean algorithms give

```text
gcd(T0,u)=gcd(T0,E)=gcd(T0,A0)=1,
gcd(T0,[p^1]R_18)=1.                                  (14)
```

Thus every root of `T0`, not merely a numerical control, preserves the scope,
the top endpoint, and the next bottom row. Equations `(8)--(14)` prove that
the five rows are exhaustive and disjoint and that their residual degrees are
exactly `21,20,19,18,17`.

## 3. Coordinate fibres and exact affine lengths

In normalized coordinates, `(3)` is

```text
G_N=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3+K0Y^2
    +Phi P^2Y+Theta PY^2+eta P^3Y+zeta Y^3.            (15)
```

Put `f=(G_N)_X/T` and `h=(G_N)_T`. Direct reduction gives two universal
pairs:

```text
T=0,      X^2=-6,       det Hess(G_N)=+6,
T=-1/6,   X^2= 6,       det Hess(G_N)=-6.              (16)
```

They are the two `t=0` and two `p=0` points omitted by the rational source
chart. Restoring all four points to the residual degrees in Section 2 gives
the exact affine critical lengths

```text
L=25,24,23,22,21.                                      (17)
```

The nonzero Hessians in `(16)` and the bridge `(7)` are why no projected-root
discriminant is needed.

## 4. Contracted boundary, replacement face, and packet

For the generic source fibre write

```text
F_Q=(s^2-p)(1-QH)-Qs^2/2.                              (18)
```

Its lower Newton polygon is

```text
(0,1),(2,0),(5,3),(1,5),(0,4),
(2Area,boundary,g_Pick)=(30,10,11).                    (19)
```

The four nonvertical primitive faces, in order, are

```text
s^2(1-Q/2)-p,

s^2[(1-Q/2)-K0Q(sp)^2-zeta Q(sp)^3],

Qsp^3(p-s^2)(eta p+zeta s^2),

Qp^4[-1376/135+eta sp].                               (20)
```

The vertical face

```text
p[-1+Q(-3p+(8/3)p^2-(1376/135)p^3)]
```

lies on the affine divisor `s=0` and is not an infinity puncture.

The fourth face in `(20)` is the replacement discarded by a naive
`Delta -> 0` specialization. In `W=sp` it has the simple nonzero root

```text
W=1376/(135eta),                                       (21)
```

so it supplies one rational index-four place.

For the length-two top edge use `s=z^-1,p=(1-a)z^-2`. The leading strict
transform is, up to a unit,

```text
a(1-a)^3[eta(1-a)+zeta].                              (22)
```

Its two boundary roots are `a=0` and `a=(eta+zeta)/eta`. They are distinct,
finite torus roots precisely under the three units in `(1)`. The derivative
in `a` is a unit at both roots, and the residue differential gives index
eight at each. Thus this face contributes two rational index-eight places,
with no hidden `Phi` or `Theta` wall.

The second face in `(20)` is the cubic carrier

```text
q-1/2=K0W^2+zeta W^3.                                 (23)
```

It defines a prime degree-three extension `C(W)/C(q)`: the polynomial map in
`W` has degree three, and characteristic zero makes it separable. Explicitly,
its discriminant is

```text
(q-1/2)[4K0^3-27zeta^2(q-1/2)] != 0 in C(q).          (24)
```

The three geometric conjugates have index two. The first face in `(20)` gives
one rational index-one place. THM-4103's primitive-normal formula therefore
gives the candidate labelled packet

```text
(8,8,4,2,2,2,1)
= rational (8,8,4,1) + one cubic orbit (2,2,2).        (25)
```

Its degree and defect are

```text
n=27,                 defect=20=2*11-2.               (26)
```

The inherited THM-3827/4147 connectedness argument applies: a nontrivial
closed-polynomial factor would put a complete curve inside the finite
critical scheme. Pick gives normalization genus at most `11`, while
Riemann--Hurwitz over the elliptic target and `(26)` gives genus at least
`11`. Hence the genus is exactly `11`, `(25)` is complete, and there is no
hidden affine ramification or normalization loss.

## 5. Full and finite response contradictions

THM-4120 and the finite-separable-carrier transport of THM-4147 leave exactly
two responses.

In the full response every place in `(25)` lies above the target origin. The
cover degree and origin index are

```text
n=27,                    ind(mu_O)=20.                 (27)
```

The two handle permutations generate transitively. If `L` is one of `(17)`,
the transported fixed sheets and the commutator-support lemma give

```text
ind(mu_O)<=2(n-L)=4,6,8,10,12.                         (28)
```

Every number in `(28)` is strictly below `20`.

In the finite response, the prime residue extension forces all three carrier
conjugates together. Removing their six sheets leaves

```text
n=21,       m=beta=3,
origin packet=(8,8,4,1),       ind(mu_O)=17.           (29)
```

After the splitting base change, the three carrier meridians are
transpositions on this same action. THM-4159's carrier-orbit lemma says that
the union of the two handle supports has size at least `n-m=18` and, writing
`k` for their intersection,

```text
k<=n+m-L=24-L,
ind(mu_O)<=2k+m.                                       (30)
```

For `L=25`, the fixed-sheet bound makes the total handle support at most
`2n-L=17`, already contradicting the union floor `18`. For
`L=24,23,22,21`, equation `(30)` gives

```text
ind(mu_O)<=3,5,7,9,                                    (31)
```

again all strictly below `17`. Thus neither response exists in any of the
five exhaustive rows. This proves the theorem.

## 6. Exact boundary and next frontier

The three excluded coefficient walls are structural and already typed:

```text
zeta=0:        P-only or lower-weight route;
eta=0:         Y-only route;
eta+zeta=0:    repeated-top route.                     (32)
```

They are handled by earlier theorems only in their stated coefficient
regions; `(32)` is not a global coverage claim. The result closes the mixed
off-antidiagonal `Delta=0` face inside the inherited exact-`M=9` seam. It does
not close generic nonzero-`Delta`, nonzero-`K` critical walls, another reduced
cell, entry, exact `M>=10`, `JC(2)`, or `DC(2)`.

The next dense exact-`M=10` scout has a candidate packet
`(9,9,6,2,2,2,1)` and one exact critical control of length `29`; those data
remain a proof candidate because uniform symbolic endpoints and an independent
audit have not yet been supplied.

## 7. Verification

The primary certificate checks the complete source, source-first resultant
tower, terminal all-roots Euclidean gates, Hessian bridge, universal fibres,
Newton polygon and faces, carrier discriminant, packet, and every response
inequality. The independent certificate reimplements every direct `(A,B)`
wall specialization, derives the distinct `(A,C_0)` projection from
resultant multiplicativity, proves both Hessian bridges, and independently
rebuilds the chart-loss ledger.

Replay with

```bash
python3 -B 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217.py
python3 -B -O 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217.py
PYTHONHASHSEED=4217 python3 -B 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217.py

python3 -B 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217_independent_audit.py
python3 -B -O 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217_independent_audit.py
PYTHONHASHSEED=4217 python3 -B 04-computation/jc23_mixed_offanti_delta_zero_complete_exclusion_thm4217_independent_audit.py
```

Compare stdout with the two frozen outputs.

**QED.**
