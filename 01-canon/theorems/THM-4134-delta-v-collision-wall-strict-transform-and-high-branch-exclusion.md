---
id: THM-4134
title: "Delta-V collision-wall strict transform and full-boundary branch exclusion"
status: >
  PROVED RELATIVE TO THM-3996/4053/4103/4120/4122/4130 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED. On the theta-only exact-M=8 collision wall
  Delta_V=0, the reciprocal source has two exact strata. Generically it has
  genus 8, packet (7,5,3,2,2,1), affine critical length 19 and degree 16 or
  20; on one secondary wall it has genus 7, packet (7,3,2,2,2,2,1), affine
  critical length 18 and degree 15 or 19. Critical strict transforms lose one
  or two roots respectively. Monodromy excludes the full-boundary degrees 20
  and 19. The lower degrees 16 and 15 with a finite horizontal BC carrier
  remain OPEN; the collision wall, JC(2), and DC(2) are not closed.
source: codex-frontier-synthesis-creative-20260825r
depends_on:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
script: 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134.py
output: 05-knowledge/results/jc23_delta_v_collision_wall_strict_transform_thm4134.out
independent_audit_script: 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_delta_v_collision_wall_strict_transform_thm4134_independent_audit.out
script_sha256: 4aad959ef19bd96015b219970368ebd241b5c57c534ea001238d796e4f322604
output_sha256: a0c3479c7ab4d56d210541d2992fd89fc3e0e5e0fc23ade07c0d5938bb6dffc0
independent_audit_script_sha256: 96f11131edfc834a533682dcce44ac8027cd27de31d2cebc5816475c82d73bcb
independent_audit_output_sha256: 07843acdc731cdb8a93b638a73d468d1c93900bfade8594af6ad1f2c8d47890b
semantic_sha256: a5d55d70a100e3f4c9413c9b4ab88fe5532b6d8c58740888042f31bc170ace2d
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy calculation rebuilds the complete (X,T) critical
  resultant, passes to reciprocal (S,u) coordinates, divides the exceptional
  multiplicity, computes both boundary normal forms and residue orders,
  freezes the degree and critical ledgers, checks every full-boundary support
  split, and retains a hostile residual and horizontal-section control.
  Normal, optimized and two hash-seeded replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A separate SymPy calculation imports no primary code and uses the
  rational (s,p) critical pair, p-resultant, z=1/p boundary chart, distinct
  dictionary permutations, a different mod-101 residual and elliptic
  addition reconstruction. It independently reproduces both strict-transform
  laws, packets, critical lengths, degree alternatives, exclusions and the
  shared semantic ledger. All replay modes byte-match its frozen output.
---

# THM-4134 -- Delta-V collision-wall strict transform

**PROVED RELATIVE TO THM-3996/4053/4103/4120/4122/4130 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED.**

THM-4130 empties the smooth nonresonant theta-only exact-`M=8` seam, but its
degree-sixteen critical eliminant loses its leading row on `Delta_V=0`.
Dividing that exceptional multiplicity reveals two wall strata. Their
full-boundary branches are impossible, but a horizontal quadratic carrier
prevents the smooth proof from crossing the wall completely.

## 1. Reciprocal wall atlas

Use THM-4130's normalized source polynomial `G(X,T)` and put

```text
u=T^(-1),             S=XT,             R=1+S^2u.          (1)
```

Then the exact reciprocal equation is

```text
u^3 G = K
      = A(S)R^3+uB(S)R^2-3u^2R-S^2u^4/2,                 (2)

A(S)=Theta S^2+Phi S-1376/135,
B(S)=8/3+(2848/45)S^2.                                   (3)
```

On the collision wall define

```text
W=135Phi^2+5504Theta=0,          Phi!=0,
S_0=2752/(135Phi),
D=273375Phi^2+2696167424.                                (4)
```

Exact substitution gives

```text
A(S)=Theta(S-S_0)^2,
B(S_0)=8D/(820125Phi^2).                                  (5)
```

### Generic stratum: `D!=0`

Writing `v=S-S_0`, the boundary equation begins

```text
Theta v^2+B(S_0)u.                                       (6)
```

It is therefore a smooth tangency with `u=-(Theta/B(S_0))v^2+...`.
The Keller residue differential has order four, so the two smooth-seam
index-three punctures collide into one index-five puncture. The source keeps
genus eight and has boundary packet

```text
(7,5,3,2,2,1).                                           (7)
```

### Secondary stratum: `D=0`

Here the linear `u` term vanishes. The first boundary form is

```text
Theta v^2+B'(S_0)uv-3u^2.                                (8)
```

Modulo `D=0`,

```text
Theta=489856/2025,              B'(S_0)^2=-91136/135,
disc(8)=501248/225 !=0.                                  (9)
```

Thus `(8)` is an ordinary node. Normalization lowers the genus to seven; the
two branches have residue order one and indices two. The packet is

```text
(7,3,2,2,2,2,1).                                         (10)
```

## 2. Critical strict transform

Let `Q_16(T)` be THM-4130's residual critical eliminant and set

```text
Qhat(u,W)=u^16 Q_16(1/u),
Theta=(W-135Phi^2)/5504.                                  (11)
```

On `W=0`,

```text
Qhat(u,0)=u Qhat_15(u),
Qhat_15(0)=513984438272 Phi^4 D/12075125625.              (12)
```

For `D!=0`, division by the proved multiplicity one gives

```text
u=-(11610/D)W+O(W^2),                                    (13)
S=S_0+
  2752(273375Phi^2-1348083712)W/(18225Phi^3D)+O(W^2),
G W^2 -> 2D^3/(82909778259375Phi^2).                     (14)
```

One critical point escapes both affine charts. Restoring the two inherited
critical pairs leaves exactly nineteen affine critical points under the
Keller-Morse hypothesis.

On `D=0`, the `u`-linear coefficient also vanishes but the quadratic row does
not. The divided normal law is

```text
u^2=-(5/501248)W+...,                                    (15)
```

so two critical roots escape and eighteen affine critical points remain.
The square-root transposition in `(15)` is critical-resultant inertia, not
cover inertia; no cover cycle is inferred from it.

## 3. Degree alternatives and high-branch exclusion

The target-pencil calculation `E_q(k(q))={O}` from THM-4120 is unchanged.
The rational punctures in `(7)` contribute `1+3+7+5=16` sheets over `O`,
while the BC quadratic point contributes either zero or four. Hence

```text
D!=0:             n in {16,20}.                          (16)
```

On the secondary wall, the rational contribution is `1+3+7+2+2=15`, giving

```text
D=0:              n in {15,19}.                          (17)
```

The full-boundary alternatives are impossible. At `n=20` or `n=19`, the
affine critical length is `n-1`. If the two nodal critical fibres have
defects `r_0+r_1=n-1`, THM-4130's fixed-sheet support argument gives

```text
|supp X|+|supp Y|<=n+1.                                  (18)
```

Transitivity forces the two supports to cover all sheets and meet in one
pivot. Thus `X,Y` are single cycles through that pivot and their commutator
is a three-cycle. Its cycle type contradicts both `(7)` and `(10)`. The
audits check all sixteen generic and fifteen secondary support splits.
Therefore

```text
n=20: CONTRADICTION,              n=19: CONTRADICTION.    (19)
```

## 4. Exact residual and hostile controls

Every hypothetical Keller pair on this wall is reduced to exactly one of

```text
D!=0: n=16, affine critical length 19, horizontal BC carrier;
D=0:  n=15, affine critical length 18, horizontal BC carrier. (20)
```

In either case the BC quadratic point maps to a finite target point over

```text
kappa Y^2=q-a^3/2.                                       (21)
```

This is precisely where the smooth finite-transport proof stops. The residual
is not vacuously removed by THM-4122's pole-ratio gate: after
`q=a^3/2+rho^2`, the quadratic base change has the polynomial section

```text
U=a/2+16rho^2/(9a^2),
V=-rho-64rho^3/(27a^3),                                  (22)
```

with intrinsic pole pair `(2,3)`.

At the exact hostile point

```text
(Phi,Theta)=(5504,-743040),                               (23)
```

the residual critical polynomial is squarefree of degree fifteen over `Q`.
Its primary and independent projections have squarefree reductions modulo
101 with distinct coefficient rows, providing a coordinate-independent
hostile control.

## 5. Scope and replay

This theorem excludes the two full-boundary branches and types the exact
survivors. It does not empty either survivor in `(20)`, prove the BC point
finite for an arbitrary map outside this chart, cross the other two collision
walls, treat `M>=9`, or prove `JC(2)` or `DC(2)`.

Run

```text
python3 -B 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134.py
python3 -B -O 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134.py
python3 -B 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134_independent_audit.py
python3 -B -O 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_delta_v_collision_wall_strict_transform_thm4134_independent_audit.py
```

All six streams match their corresponding frozen outputs, and both semantic
ledgers have digest
`a5d55d70a100e3f4c9413c9b4ab88fe5532b6d8c58740888042f31bc170ace2d`.
**QED.**
