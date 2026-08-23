---
id: THM-3815
title: "Quadratic r-profile constant-z2 nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the c=1
  cubic pseudo-plane, every canonical nodal carrier
  A=e^2-z/3+r(a_0+a_1e+a_2e^2)+h z^2 has a critical point.  A source-first
  slope chart gives two quartics and a degree-sixteen residual whose roots
  cannot all lie on the only forbidden linear divisor.  An independent
  Cr-first quartic-cubic elimination instead reduces a repeated-root-safe
  boundary hypothesis to five coefficients whose ideal contains a_2^2h^5.
  This closes the quadratic-r/constant-z2 mixed cell, not nonconstant z2
  profiles or cubic-and-higher mixed profiles.
source: root / all-frontier quadratic-mixed carrier lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (thm3798-hostile-audit, 2026-08-23).  The
  audit uses a different slope elimination: after removing the inherited
  a_2=0 and h=0 axes, it constructs a quartic p1 and cubic p2, factors their
  resultant as 11664 z^3 L B^2 H16, proves that boundary-only support would
  force a_2^2h^5=0, and reconstructs a finite source critical point without
  dividing by the harmless degree-drop factor.  It checks 2,164 finite-field
  hostiles as controls and has 39 active exact gates.  Normal and optimized
  executions LF-normalize exactly to their frozen transcripts; neither
  script contains Python assertions.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
related:
  - THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points
  - THM-3810-affine-r-profile-constant-z2-nodal-carriers-have-critical-points
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
script: 04-computation/jc2_quadratic_r_constant_z2_nodal_carrier_thm3815.py
output: 05-knowledge/results/jc2_quadratic_r_constant_z2_nodal_carrier_thm3815.out
script_sha256: 3d8537c90ba8346c65db66ea5d79c3c4d4407ecf7c53f76d9a70c3b6149b96d5
output_sha256: f7e75b7f1c7d492041141670bb4df4c8d8d24da84e8b3b8e6f414a83d5af2ed0
semantic_sha256: 2e90439d024488eefc6647b5a991d603276ed536a53e3fb5e3340fd27a565a2e
independent_script: 04-computation/jc2_quadratic_r_constant_z2_nodal_carrier_independent_audit_thm3815.py
independent_output: 05-knowledge/results/jc2_quadratic_r_constant_z2_nodal_carrier_independent_audit_thm3815.out
independent_script_sha256: f4fc2de02d98301edef10fa89c6b8ad061b666dd362920ccee5796ce50a7dc8a
independent_output_sha256: 62af0123d5e41c7f7f0b4d1bd2c1813ea333cb07f57ccf846ba6c4df144f1e6b
independent_semantic_sha256: 8bd378050e2497a2914cf202f9c0b51e0ca99835ec6f03826192fc56b77c464e
hash_basis: raw LF bytes
---

# THM-3815 -- every quadratic-r/constant-z2 nodal carrier is critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
algebraically closed of characteristic zero.  On

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re,       (1)
```

put, for arbitrary `q,b,c,h in k`,

```text
A=e^2-z/3+r(qe^2+be+c)+h z^2.                         (2)
```

Then `A` has a critical point on `Y`.  Since `Y` is a smooth symplectic
surface, `A` has no regular Darboux mate.

## 1. Hamiltonian packet and source-first slope

Write

```text
g=qe^2+be+c,              L=1-6hz,              K=1+2re.
```

The three Hamiltonian components are

```text
C_r=Lr^2-9z^2(2e+r(2qe+b)),
C_z=3gr^2-3K(2e+r(2qe+b)),
C_e=9gz^2-LK.                                      (3)
```

They obey the exact Casimir syzygy

```text
K C_r-3z^2 C_z+r^2 C_e=0.                          (4)
```

On the chart `z!=0`, set `u=r/z`.  The source equation itself, rather than
`C_r`, solves for

```text
w=z^2-u,                 e=w/(u^2z).               (5)
```

After substituting `(5)`, the exact denominator-free forms of `C_r,C_e` are

```text
p1=L u^4z-18w-18quzw-9bu^3z^2,
p2=9qw^2+9bzu^2w+9cz^2u^4-L(u^4+2u^3w),            (6)

C_r=z p1/u^2,                     C_e=p2/u^4.       (7)
```

Both polynomials in `(6)` are quartic in `u` before specialization.

## 2. The degree-sixteen residual cannot stay on the boundary

An exact fixed `8 by 8` Sylvester determinant gives

```text
Res_u(p1,p2)=-729z^8U(z),
deg_z U=16,             terms(U)=158,              U(0)=144.       (8)
```

The full residual is frozen structurally by SHA-256
`28b077c0775d3603163dddd01d078ba18feee257205a41e53d300684cbbbc449`.
Only the following small coefficient packet is needed for the proof.

Suppose first that `h!=0`.  The linear coefficient is

```text
[z]U=-3456h.                                         (9)
```

Thus `U` is nonconstant.  If every root of `U` lay on `L=0`, algebraic
closedness and `U(0)=144` would force `U=144L^d`.  Comparing `(9)` with
`[z](144L^d)=-864dh` forces `d=4`.  Put

```text
E_2=q^3+2qb+8c.
```

The successive exact coefficients of `U-144L^4` are

```text
[z^2]=324E_2,
c=-(q^3+2qb)/8       => [z^3]=144q,
q=c=0                => [z^5]=1296b^2,
q=b=c=0              => U-144L^4=128z^7L^5 !=0.    (10)
```

No parameter is divided in `(10)`, so every zero-axis seam is included.
This contradiction supplies a root `z_0` with `z_0L(z_0)!=0`.

If `h=0`, then `L=1`.  Were `U` constant, it would equal `144`.  The same
`z^2,z^3,z^5` cascade forces `q=c=b=0`, while the terminal specialization is

```text
U=144+128z^7,                                        (11)
```

again a contradiction.  Its nonzero root also has `z_0L(z_0)!=0`.

## 3. Projective gate and exact reconstruction

At either selected root,

```text
LC_u(p1)=Lz !=0,                 p1(0)=-18z^2 !=0.   (12)
```

Consequently the homogenized pair cannot meet at the point at infinity,
even if `p2` drops degree.  The vanishing resultant therefore supplies a
finite common root `u_0`, and `(12)` gives `u_0!=0`.  Define

```text
r_0=u_0z_0,              e_0=(z_0^2-u_0)/(u_0^2z_0). (13)
```

Direct identities, before imposing either quartic, give

```text
r_0^2e_0-z_0^3+r_0=0,
C_r(r_0,z_0,e_0)=z_0p1/u_0^2,
C_e(r_0,z_0,e_0)=p2/u_0^4.                           (14)
```

Thus `C_r=C_e=0`; equation `(4)` and `z_0!=0` give `C_z=0`.  The point in
`(13)` is the required source critical point.

## 4. Independent elimination and scope boundary

The hostile audit does not replay `(5)--(10)`.  On the genuine
`q h!=0` cell it solves `C_r` first, with

```text
W=1+quz,       M=Lu^2-9buz,       B=L+9bqz^2.
```

It obtains a quartic and a cubic whose resultant is

```text
11664 z^3 L B^2 H_16,             H_16(0)=144,
LC_z(H_16)=314928q^4h^4.          (15)
```

If every root of `H_16` lay on `LB=0`, repeated-root-safe multiplicity
counting would imply `H_16 | LB H_16'`.  After exact division, the ideal of
only the first five low primitive remainder coefficients contains

```text
q^2h^5,                                             (16)
```

contradicting the genuine-cell hypothesis.  Its reconstruction separately
checks the possible `W`, leading-coefficient, `u`, and `K` seams.  The axes
`q=0` and `h=0` agree with THM-3810 and THM-3805, while the primary proof
above already includes them uniformly.

This theorem closes exactly a quadratic polynomial `r` profile together
with a constant `z^2` coefficient.  It does **not** cover nonconstant
`z^2h(e)`, cubic-or-higher mixed `r` profiles, rational-pole mates, another
pseudo-plane arm, or arbitrary planar Keller maps.  THM-3812 independently
excludes a two-output bracket pair confined to arbitrary `r/z^2` arm
profiles, but does not replace this single-carrier critical-point theorem.
Normal and optimized executions of both exact companions LF-normalize to
the frozen outputs.  **QED.**
