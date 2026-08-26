---
id: THM-4164
title: "Y-only triple-top-root planar Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On the Y-only eta=Delta=0
  exact-weight-nine wall, the triple-top-root locus with zeta*I_C!=0
  supports no nonautomorphic planar Keller pair. For J=15a^2+356 nonzero
  the genus is 9, affine critical length is 20, and packet is
  (8,7,2,2,2,1); on J=0 the genus is 8, length is 19, and packet is
  (8,4,2,2,2,2,1). THM-4165 subsequently closes the I_C=0 intersection.
  The zeta=0 wall, other cells, entry, M>=10, JC(2), and DC(2) remain OPEN.
source: codex-lrc14-jc-sharp-fronts-20260825b
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4007-live-two-three-third-normal-row-five-weight-floor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
related:
  - THM-4157-repeated-top-edge-wall-planar-jacobian-exclusion
  - THM-4159-inner-resultant-wall-planar-jacobian-exclusion
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4165-y-only-inner-top-triple-intersection-planar-jacobian-exclusion
script: 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164.py
output: 05-knowledge/results/jc23_y_only_triple_top_root_exclusion_thm4164.out
independent_audit_script: 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_y_only_triple_top_root_exclusion_thm4164_independent_audit.out
script_sha256: 2d64e99bd52789ec773f11eaf7ffe6f7bd88d67151e865b8af29c250c45d6c32
output_sha256: d8b3e5cc321e7a5309d909cb6376b85703e68241e2da8c5a991940ee19f14f1a
independent_audit_script_sha256: 2a906494eac680fa5a27bff5b30ad0fb75b8f3dc65373e0c38624381ddfc75cb
independent_audit_output_sha256: 2fc6d8d53a1b67f798d7f413e7eaf61fac8c45a9727c43ae14a7a943cde5908d
semantic_sha256: a9e1e53b621dabf0e5ebe569c228abfcbda27e65dfb78230415f0dafa763c1d9
independent_semantic_sha256: 65935648fd74d7d394a48bb4b13f9930259eb6be3f51b97634663a275cda731c
hash_basis: raw LF bytes
primary_audit: >
  PASS. The exact certificate parameterizes the complete triple-root chart,
  computes source and normalized critical projections with every endpoint
  and infinity gate, restores the four universal points, resolves both the
  smooth cubic tangency and its J=0 node, and checks all full and finite
  response inequalities. Normal, optimized, and hash-seeded outputs
  byte-match.
independent_audit: >
  ACCEPT. A clean-room source-chart referee uses the alternative pair
  (A,C_0), hence p^8 rather than p^6, a disjoint rational control, direct
  quadratic-field reduction, and an independent tangent-cone/Hensel branch
  audit. It reproduces both lengths, genera, packets, and response
  contradictions in normal, optimized, and hash-seeded runs.
weakest_link_audit: >
  The weakest inherited link is finite-separable cubic-carrier and fixed-sheet
  transport after the new boundary normalization. The carrier face is
  unchanged from THM-4147/4155, and both local normalizations are explicit.
  The newest local link, the J=0 Hensel split, has a squarefree tangent cone
  and is independently coefficient-audited; it remains the preferred hostile
  target for any further geometric referee.
---

# THM-4164 -- Y-only triple-top-root planar Jacobian exclusion

**PROVED RELATIVE TO THM-3827/3992/3997/4007/4103/4120/4122/4130/4147/4155
+ VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.** Work over `C`
in the live `b=d=0` reduced `(2,3)` seam.

## 1. Inheritance pass and statement

Put

```text
kappa=1376/135,                         K0=2848/45,
C(W)=zeta W^3+Theta W^2+Phi W-kappa,
I_C=4Theta K0^2-27zeta^2.
```

> **Theorem.** The exact-weight-nine locus
>
> ```text
> eta=Delta=0,             zeta!=0,
> C has one triple root,   I_C!=0
> ```
>
> contains no nonautomorphic planar Keller pair.

The closest proved mechanisms are THM-4155's Y-only critical/boundary
monodromy exclusion and THM-4161’s double-root strict transform.
The canonical hostile is exactly the triple-root locus omitted by THM-4161.
The corrected comparison is THM-4157: it treats a repeated slanted edge on
`zeta=-eta,Delta!=0`, whereas this is a horizontal-edge collision on
`eta=Delta=0`. The least-used sidecar is the paired critical/boundary endpoint
`J=15a^2+356`; retaining it reveals a second, nodal stratum rather than a
missing point.

The weakest proof link is geometric rather than computational: the final
finite response imports THM-4147/4155's finite-separable cubic-carrier and
fixed-sheet transport after normalization of the new boundary point. The
carrier face itself is unchanged, and both new local normalizations are
explicit, but this transport should remain an explicit audit target. On the
new algebraic side, the most delicate step is the `J=0` Hensel split into the
two branches with indices two and four; the tangent cone is squarefree and
both independent certificates verify the exact coefficients.

## 2. Complete triple-root chart

The constant and leading coefficients of `C` are nonzero, so its triple root
`r` is nonzero. With `a=r^-1`, every point in the statement has the unique
form

```text
C(W)=kappa(aW-1)^3,                  a!=0,
zeta=kappa a^3,
Theta=-3kappa a^2,
Phi=3kappa a.                                      (1)
```

Conversely `(1)` gives every triple-root point. Define

```text
J=15a^2+356,                 H_I=5805a^4+1013888.
```

Then

```text
I_C=-(44032/91125)a^2 H_I.                           (2)
```

Thus theorem scope is exactly `a*H_I!=0`, with no assumption on `J`.

## 3. Complete source and two critical projections on `J!=0`

Use the complete THM-4155 source

```text
s=XT, p=T+s^2, t=p-s^2, P=T+X^2T^2, Y=XTP,

H=-3p+(8/3)p^2-kappa p^3+K0s^2p^2
  +Phi sp^3+Theta s^2p^3+zeta s^3p^3,
G=-s^2/(2t)+H.
```

For

```text
A=(-sp+t^2H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2,
```

the exact identities `t^2G_s=pA` and `2t^2G_p=t^2B-sA` hold. Direct
elimination gives

```text
Res_s(A,B)=p^6R_16(p),

[p^16]R_16=
 (9563767153858210665857024/9341736328125)a^17J^2,
R_16(0)=
 (3877634048/16875)a^5(5805a^4+1013888)
 =-46656zeta I_C.                                     (3)
```

The leading `s`-rows of `(A,B)` are `3zeta p^2,9zeta p^2`, so no finite-`p`
intersection is lost at `s=infinity`.

Independently put `f=G_X/T`, `h=G_T`. Exact normalized elimination gives

```text
Res_X(f,h)=T^56(6T+1)^2Q_16(T),

[T^16]Q_16=
 (612081097846925482614849536/38306957530517578125)
 a^15J^2(5805a^4+1013888)^2,
Q_16(0)=-(3^15/2^7)zeta^7.                            (4)
```

The leading `X`-rows are both `9zeta T^8`, hence no residual root lies at
`X=infinity`. The `T^56` factor is a Sylvester degree-drop artifact. The
actual ideal restores

```text
T=0,    X^2=-6, G=0,   det Hess(G)=+6,
T=-1/6, X^2= 6, G=1/2, det Hess(G)=-6.                 (5)
```

Therefore

```text
L=16+2+2=20.                                           (6)
```

No residual discriminant or `Q_16(-1/6)` is a theorem hypothesis: a repeated
projected coordinate preserves scheme length, and a Keller realization makes
the actual critical scheme reduced by Hessian congruence. The exact rational
control `a=1` has `J=371`, `I_C!=0`, and both residuals are squarefree. The
independent `(A,C_0)` audit uses `a=2` and obtains `p^8R_16` with the same
endpoints.

## 4. Smooth cubic tangency, genus, and responses

Clear the generic fibre by

```text
F=q(s^2-p)-(s^2-p)H-s^2/2.
```

At the triple root put `u=1/p`, `v=s-1/a`, and
`K=u^4F(1/a+v,u^-1)`. Its first two rows are

```text
K=kappa a^3v^3+[8J/(45a^2)]u+higher.                  (7)
```

For `J!=0`, the point is smooth and

```text
u=-(172a^5/(3J))v^3+O(v^4).                           (8)
```

On the curve

```text
omega=ds/F_p=-u^2ds/K_u.
```

Thus `ord_v(omega)=6`, giving one rational index-seven place. It replaces
the three rational index-three places without lowering genus: both old and
new defects equal six.

The Newton polygon remains

```text
(0,1),(2,0),(5,3),(3,4),(0,4),       (2Area,B,I)=(27,11,9).
```

All other faces are unchanged from THM-4155, so

```text
packet=(8,7,2,2,2,1),
rational=(8,7,1) + cubic carrier=(2,2,2),
defect=16=2*9-2,                         genus=9.       (9)
```

The carrier equation remains

```text
q-1/2=K0W^2+zeta W^3.
```

Since `zeta!=0`, it is a prime separable degree-three residue extension.
THM-4120 and finite-separable transport leave exactly

```text
full n=22,                    finite (n,beta)=(16,3). (10)
```

In the full response `n-L=2`, hence the origin commutator has index at most
four, below packet defect sixteen. In the finite response the exact
capacities are `(13,14,3)`, all below `n-1=15`. Both responses are impossible.

## 5. Exact closure of `J=0`

Now work over

```text
Q[a]/(15a^2+356)=Q(sqrt(-1335)),
a=2sqrt(-1335)/15.                                    (11)
```

The inner gate stays live:

```text
I_C=335741565206528/6834375 !=0.                       (12)
```

Coefficientwise exact reduction of `(3),(4)` gives

```text
Res_s(A,B)=p^6R_15(p),
[p^15]R_15=
 -(3636055657731868916207553629845852858180471619584
   /1596123230438232421875)a,
R_15(0)=10525761290611774717952a/18984375,

Res_X(f,h)=T^56(6T+1)^2Q_15(T),
[T^15]Q_15=
 (299870607546709935491031806891677436336739204843700202793598976
  /10908504703026294708251953125)a,
Q_15(0)=
 -2385717670936132253981396577199885145595707392a
  /32842041778564453125,
Q_15(-1/6)=
 -71304646963571565929235399784047641547404476816773723491139584a
  /782625597463934606616497039794921875.               (13)
```

Both residuals are squarefree over the quadratic field, and the leading
coordinate rows remain units. Consequently

```text
L=15+2+2=19.                                           (14)
```

The boundary expansion now begins

```text
K=-3u^2-(16a/3)uv-(489856a/2025)v^3+higher.            (15)
```

Its tangent cone `u(-3u-(16a/3)v)` has two distinct lines, hence an ordinary
node of delta one. Hensel lifting gives the two branches

```text
u=-(16a/9)v+O(v^2),
u=-(30616/675)v^2+O(v^3).                              (16)
```

For the first, `u` and `K_u` have orders one, so `ord(omega)=1` and `e=2`.
For the second, `u` has order two while `K_u` has order one, so
`ord(omega)=3` and `e=4`. The node lowers genus to eight. Therefore

```text
packet=(8,4,2,2,2,2,1),
rational=(8,4,2,1) + cubic carrier=(2,2,2),
defect=14=2*8-2,                         genus=8,       (17)
full n=21,                    finite (n,beta)=(15,3).
```

The full commutator bound is again four, below fourteen. Finite capacities
are `(12,13,3)`, all below `n-1=14`. Thus `J=0` is also impossible, proving
the theorem in all stated scope.

## 6. Scope and replay

The inner intersection excluded from this theorem is

```text
I_C=0, equivalently 5805a^4+1013888=0.
```

THM-4165 subsequently closes that inner wall. No claim is made about the
zeta=0 wall, another reduced cell, entry, exact weight at least ten, `JC(2)`,
or `DC(2)`.

Primary replay:

```text
python3 -B 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164.py
python3 -B -O 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164.py
PYTHONHASHSEED=521 python3 -B 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164.py
```

Independent replay:

```text
python3 -B 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164_independent_audit.py
python3 -B -O 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164_independent_audit.py
PYTHONHASHSEED=613 python3 -B 04-computation/jc23_y_only_triple_top_root_exclusion_thm4164_independent_audit.py
```

**QED.**
